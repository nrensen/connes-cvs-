"""Adversarial checkpoint, artifact, and input-contract tests."""

from __future__ import annotations

import json
from dataclasses import FrozenInstanceError, asdict

import mpmath as mp
import pytest

import connes_cvs.runner as runner
from connes_cvs.runner import CellConfig, GalerkinCell
from connes_cvs.sweep import run_sweep


def _fake_psi_point(index):
    x = mp.mpf(index)
    psi = x**3 + 2 * x
    psi_deriv = 3 * x**2 + 2
    return index, psi._mpf_, psi_deriv._mpf_, 0.0


def _cell(tmp_path, monkeypatch, *, checkpoint=True, artifact=False, tag="test"):
    monkeypatch.setattr(runner, "_worker_psi_point", _fake_psi_point)
    return GalerkinCell(
        CellConfig(c=2, N=2, T=10, dps=30, tag=tag),
        workers=1,
        progress_callback=False,
        ground_state="minimum",
        checkpoint_path=tmp_path / "cell.checkpoint.json" if checkpoint else None,
        artifact_path=tmp_path / "cell.json" if artifact else None,
    )


def _resign(payload):
    body = {key: value for key, value in payload.items() if key != "payload_sha256"}
    payload["payload_sha256"] = runner._payload_sha256(body)


def test_clean_and_resumed_runs_are_raw_identical(tmp_path, monkeypatch):
    clean = _cell(tmp_path, monkeypatch, artifact=True)
    config_before = asdict(clean.cfg)
    clean_artifact = clean.run()

    def unexpected_compute(_index):
        raise AssertionError("complete checkpoint attempted to recompute a point")

    monkeypatch.setattr(runner, "_worker_psi_point", unexpected_compute)
    resumed = GalerkinCell(
        clean.cfg,
        workers=1,
        progress_callback=False,
        ground_state="minimum",
        checkpoint_path=tmp_path / "cell.checkpoint.json",
    )
    resumed_artifact = resumed.run()

    assert asdict(clean.cfg) == config_before
    assert resumed.resumed_points == clean.cfg.N + 1
    assert clean.lambda_even._mpf_ == resumed.lambda_even._mpf_
    assert clean_artifact["lambda_even_raw_mpf"] == (
        resumed_artifact["lambda_even_raw_mpf"]
    )
    assert clean_artifact["matrix_sha256"] == resumed_artifact["matrix_sha256"]
    assert clean_artifact["eigenvector_sha256"] == (
        resumed_artifact["eigenvector_sha256"]
    )
    for i in range(clean.Q.rows):
        for j in range(clean.Q.cols):
            assert clean.Q[i, j]._mpf_ == resumed.Q[i, j]._mpf_


def test_checkpoint_and_artifact_hashes_are_valid(tmp_path, monkeypatch):
    cell = _cell(tmp_path, monkeypatch, artifact=True)
    artifact = cell.run()
    checkpoint = json.loads(
        (tmp_path / "cell.checkpoint.json").read_text(encoding="utf-8")
    )

    checkpoint_hash = checkpoint.pop("payload_sha256")
    assert checkpoint_hash == runner._payload_sha256(checkpoint)
    artifact_hash = artifact.pop("artifact_payload_sha256")
    assert artifact_hash == runner._payload_sha256(artifact)
    on_disk = json.loads((tmp_path / "cell.json").read_text(encoding="utf-8"))
    on_disk_hash = on_disk.pop("artifact_payload_sha256")
    assert on_disk_hash == runner._payload_sha256(on_disk)


def test_checkpoint_checksum_mismatch_is_rejected(tmp_path, monkeypatch):
    _cell(tmp_path, monkeypatch).run_phase_psi_cache()
    path = tmp_path / "cell.checkpoint.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["values"]["0"]["psi"]["mantissa_hex"] = "1"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        _cell(tmp_path, monkeypatch)._load_checkpoint()


def test_checkpoint_config_mismatch_is_rejected(tmp_path, monkeypatch):
    original = _cell(tmp_path, monkeypatch)
    original.run_phase_psi_cache()
    changed = GalerkinCell(
        CellConfig(c=2, N=2, T=10, dps=30, tag="different"),
        workers=1,
        progress_callback=False,
        checkpoint_path=tmp_path / "cell.checkpoint.json",
    )
    with pytest.raises(ValueError, match="configuration"):
        changed._load_checkpoint()


def test_checkpoint_backend_mismatch_is_rejected(tmp_path, monkeypatch):
    _cell(tmp_path, monkeypatch).run_phase_psi_cache()
    path = tmp_path / "cell.checkpoint.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["compute_fingerprint"]["mpmath_backend"] = "tampered"
    _resign(payload)
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="backend"):
        _cell(tmp_path, monkeypatch)._load_checkpoint()


def test_checkpoint_gmpy2_version_mismatch_is_rejected(tmp_path, monkeypatch):
    _cell(tmp_path, monkeypatch).run_phase_psi_cache()
    path = tmp_path / "cell.checkpoint.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["compute_fingerprint"]["gmpy2"] = "tampered"
    _resign(payload)
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="backend"):
        _cell(tmp_path, monkeypatch)._load_checkpoint()


def test_runner_rejects_arithmetic_fingerprint_change_mid_cell(tmp_path, monkeypatch):
    cell = _cell(tmp_path, monkeypatch, checkpoint=False)
    changed = dict(cell.compute_fingerprint)
    changed["mpmath_backend"] = "changed-mid-run"
    monkeypatch.setattr(runner, "_compute_fingerprint", lambda: changed)
    with pytest.raises(RuntimeError, match="changed during this cell"):
        cell.run()


def test_artifact_records_complete_mpmath_backend_provenance(tmp_path, monkeypatch):
    artifact = _cell(tmp_path, monkeypatch, artifact=True).run()
    environment = artifact["environment"]
    assert environment["mpmath_backend"] == mp.libmp.BACKEND
    for name, version in runner._gmpy2_versions().items():
        assert environment[name] == version
    assert artifact["compute_fingerprint"] == runner._compute_fingerprint()
    assert artifact["engine"] == (
        f"connes_cvs.runner v{artifact['compute_fingerprint']['connes_cvs']}"
    )
    timings = artifact["timings_seconds"]
    assert timings["total_wall_monotonic_s"] == timings["total_s"]
    assert timings["coordinator_process_cpu_s"] >= 0
    assert "excludes artifact hashing" in artifact["timing_scope"]


@pytest.mark.parametrize("content", ["not JSON", "[]", '{"x": 1, "x": 2}'])
def test_corrupt_or_duplicate_key_checkpoint_is_rejected(
    tmp_path, monkeypatch, content
):
    path = tmp_path / "cell.checkpoint.json"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError):
        _cell(tmp_path, monkeypatch)._load_checkpoint()


def test_noncanonical_checkpoint_records_and_fields_are_rejected(
    tmp_path, monkeypatch
):
    _cell(tmp_path, monkeypatch).run_phase_psi_cache()
    path = tmp_path / "cell.checkpoint.json"
    original = json.loads(path.read_text(encoding="utf-8"))

    noncanonical = json.loads(json.dumps(original))
    noncanonical["values"]["0"]["psi"]["mantissa_hex"] = "00"
    _resign(noncanonical)
    path.write_text(json.dumps(noncanonical), encoding="utf-8")
    with pytest.raises(ValueError, match="canonical|leading zero"):
        _cell(tmp_path, monkeypatch)._load_checkpoint()

    extra = json.loads(json.dumps(original))
    extra["unexpected"] = True
    _resign(extra)
    path.write_text(json.dumps(extra), encoding="utf-8")
    with pytest.raises(ValueError, match="field set"):
        _cell(tmp_path, monkeypatch)._load_checkpoint()


def test_path_collision_and_artifact_overwrite_policy(tmp_path, monkeypatch):
    same = tmp_path / "same.json"
    with pytest.raises(ValueError, match="different files"):
        GalerkinCell(
            CellConfig(c=2, N=1, T=10, dps=30),
            checkpoint_path=same,
            artifact_path=same,
        )

    artifact_path = tmp_path / "cell.json"
    artifact_path.write_text("do not overwrite", encoding="utf-8")
    cell = _cell(tmp_path, monkeypatch, checkpoint=False, artifact=True)
    with pytest.raises(FileExistsError, match="overwrite_artifact=True"):
        cell.run()
    assert artifact_path.read_text(encoding="utf-8") == "do not overwrite"

    replacement = GalerkinCell(
        cell.cfg,
        workers=1,
        progress_callback=False,
        ground_state="minimum",
        artifact_path=artifact_path,
        overwrite_artifact=True,
    )
    monkeypatch.setattr(runner, "_worker_psi_point", _fake_psi_point)
    replacement.run()
    assert json.loads(artifact_path.read_text(encoding="utf-8"))["schema"] == (
        "connes-cvs.cell-result.v2"
    )


def test_artifact_existence_is_rechecked_after_lock_acquisition(tmp_path, monkeypatch):
    artifact_path = tmp_path / "cell.json"
    cell = _cell(tmp_path, monkeypatch, checkpoint=False, artifact=True)
    original_lock = runner._exclusive_path_lock

    @runner.contextmanager
    def racing_lock(target):
        with original_lock(target):
            artifact_path.write_text("created by winner", encoding="utf-8")
            yield

    monkeypatch.setattr(runner, "_exclusive_path_lock", racing_lock)
    with pytest.raises(FileExistsError, match="overwrite_artifact=True"):
        cell.run()
    assert artifact_path.read_text(encoding="utf-8") == "created by winner"


def test_cell_config_is_frozen_and_validates_inputs():
    cfg = CellConfig(c=2, N=1, T=1, dps=15)
    with pytest.raises(FrozenInstanceError):
        cfg.N = 2

    invalid = [
        (dict(c=True, N=1, T=1, dps=15), TypeError),
        (dict(c=1, N=1, T=1, dps=15), ValueError),
        (dict(c=2, N=0, T=1, dps=15), ValueError),
        (dict(c=2, N=1, T=0, dps=15), ValueError),
        (dict(c=2, N=1, T=1, dps=14), ValueError),
        (dict(c=2, N=1, T=1, dps=15, tag=1), TypeError),
        (dict(c=2, N=1, T=1, dps=15, tag="x" * 257), ValueError),
        (dict(c=2, N=1, T=1, dps=15, flint_bits=1), ValueError),
    ]
    for kwargs, error in invalid:
        with pytest.raises(error):
            CellConfig(**kwargs)


def test_runner_and_sweep_reject_invalid_public_inputs(tmp_path):
    cfg = CellConfig(c=2, N=1, T=1, dps=15)
    assert GalerkinCell(cfg).ground_state == "minimum"
    with pytest.raises(TypeError, match="cfg"):
        GalerkinCell(object())
    with pytest.raises(TypeError, match="workers"):
        GalerkinCell(cfg, workers=True)
    with pytest.raises(TypeError, match="progress_callback"):
        GalerkinCell(cfg, progress_callback="yes")
    with pytest.raises(ValueError, match="ground_state"):
        GalerkinCell(cfg, ground_state="unknown")

    invalid_sweeps = [
        (dict(cutoffs=[]), ValueError),
        (dict(cutoffs="13"), TypeError),
        (dict(cutoffs=[2, 2]), ValueError),
        (dict(cutoffs=[1]), ValueError),
        (dict(cutoffs=[2], N=True), TypeError),
        (dict(cutoffs=[2], T=0), ValueError),
        (dict(cutoffs=[2], dps=14), ValueError),
        (dict(cutoffs=[2], workers=0), ValueError),
        (dict(cutoffs=[2], n_zeros=0), ValueError),
        (dict(cutoffs=[2], strict_zeros=1), TypeError),
        (dict(cutoffs=[2], overwrite_artifacts=1), TypeError),
    ]
    for kwargs, error in invalid_sweeps:
        with pytest.raises(error):
            run_sweep(**kwargs)

    not_directory = tmp_path / "file"
    not_directory.write_text("x", encoding="utf-8")
    with pytest.raises(OSError):
        run_sweep([2], checkpoint_dir=not_directory)


@pytest.mark.parametrize(("converged", "expected"), [(False, None), (True, mp.mpf("0.125"))])
def test_sweep_reports_gamma_error_only_for_accepted_root(
    tmp_path, monkeypatch, converged, expected
):
    import connes_cvs.sweep as sweep

    class FakeCell:
        def __init__(self, *args, **kwargs):
            self.eigvec_full = mp.matrix([[1], [1], [1]])
            self.lambda_even = mp.mpf(1)

        def run(self):
            return {
                "timings_seconds": {
                    "psi_cache_s": 1.0,
                    "matrix_assembly_s": 2.0,
                    "diagonalize_s": 3.0,
                    "total_s": 6.0,
                }
            }

    monkeypatch.setattr(sweep, "GalerkinCell", FakeCell)
    monkeypatch.setattr(
        sweep,
        "extract_zeros",
        lambda *args, **kwargs: [
            {"converged": converged, "error": mp.mpf("0.125")}
        ],
    )
    row = sweep.run_sweep([2], N=1, T=1, dps=15, workers=1)[2]
    assert row["gamma1_error"] == expected
    timing = row["timing"]
    assert timing["cache_sec"] == timing["psi_cache_s"] == 1.0
    assert timing["matrix_sec"] == timing["matrix_assembly_s"] == 2.0
    assert timing["diag_sec"] == timing["diagonalize_s"] == 3.0
    assert timing["zeros_sec"] == timing["zeros_s"]
    assert timing["total_sec"] == timing["total_with_zeros_s"]
    assert "excludes runner artifact hashing" in row["timing_scope"]
