"""Parallel, resumable compute-cell runner for the CvS Galerkin matrix.

The runner uses the same arithmetic kernels and matrix conventions as
``build_galerkin_matrix``. It parallelizes only the independent basis-index
quadratures, transports raw mpmath tuples losslessly, emits a live heartbeat,
and can checkpoint the expensive psi cache atomically for later resumption.

Examples
--------
Multiprocessing entry points must be protected on spawn platforms::

    from connes_cvs.runner import CellConfig, run_cell

    if __name__ == "__main__":
        result = run_cell(CellConfig(c=13, N=100, T=400, dps=80))
        print(result["lambda_even"])

Working precision is always explicit. Choose it by cross-precision
recomputation and reference comparison; this module does not auto-calibrate or
claim that a heuristic precision budget is sufficient.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import platform
import tempfile
import time
from contextlib import ExitStack, contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from multiprocessing import cpu_count, get_context
from numbers import Integral
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import mpmath as mp

from connes_cvs.operator import (
    HAS_FLINT,
    _compute_psi_pair,
    prime_powers_up_to,
)

if HAS_FLINT:
    from flint import ctx as flint_ctx


_C13_REFERENCE = "2.86545361493028029516151514986747977533e-59"
_CHECKPOINT_SCHEMA = "connes-cvs.psi-cache.v2"
_ARTIFACT_SCHEMA = "connes-cvs.cell-result.v2"


def _require_int(name: str, value: int, minimum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer, got {type(value).__name__}")
    value_int = int(value)
    if value_int < minimum:
        raise ValueError(f"{name} must be >= {minimum}, got {value_int}")
    return value_int


def _require_cutoff(c: int) -> int:
    return _require_int("c", c, 2)


def _mpf_record(value: mp.mpf) -> Dict[str, Union[int, str]]:
    """Encode a finite mpf without Python's decimal-integer digit limit."""
    value_mp = mp.mpf(value)
    if not mp.isfinite(value_mp):
        raise ValueError("only finite mpf values can be serialized")
    sign, mantissa, exponent, bitcount = value_mp._mpf_
    return {
        "sign": int(sign),
        "mantissa_hex": format(int(mantissa), "x"),
        "exponent_base_2": int(exponent),
        "bitcount": int(bitcount),
    }


def _mpf_from_record(record: object, dps: int) -> mp.mpf:
    """Decode and validate a canonical finite mpf checkpoint record."""
    if not isinstance(record, dict) or set(record) != {
        "sign", "mantissa_hex", "exponent_base_2", "bitcount"
    }:
        raise ValueError("invalid mpf record fields")
    sign = record["sign"]
    exponent = record["exponent_base_2"]
    bitcount = record["bitcount"]
    mantissa_hex = record["mantissa_hex"]
    for label, value in (("sign", sign), ("exponent", exponent), ("bitcount", bitcount)):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"mpf {label} must be an integer")
    if sign not in (0, 1):
        raise ValueError("mpf sign must be 0 or 1")
    if not isinstance(mantissa_hex, str) or not mantissa_hex:
        raise ValueError("mpf mantissa_hex must be a nonempty string")
    if mantissa_hex != mantissa_hex.lower() or any(
        ch not in "0123456789abcdef" for ch in mantissa_hex
    ):
        raise ValueError("mpf mantissa_hex must be canonical lowercase hexadecimal")
    if len(mantissa_hex) > max(64, int(dps * 0.84) + 64):
        raise ValueError("mpf mantissa exceeds the configured precision bound")
    if len(mantissa_hex) > 1 and mantissa_hex.startswith("0"):
        raise ValueError("mpf mantissa_hex must not contain leading zeroes")
    if abs(exponent) > 10_000_000:
        raise ValueError("mpf exponent is outside the supported checkpoint range")
    mantissa = int(mantissa_hex, 16)
    if mantissa == 0:
        if (sign, exponent, bitcount) != (0, 0, 0):
            raise ValueError("zero mpf records must use the canonical 0,0,0 encoding")
    else:
        if mantissa & 1 == 0 or bitcount != mantissa.bit_length():
            raise ValueError("mpf mantissa is not in canonical normalized form")
    value = mp.mpf((sign, mantissa, exponent, bitcount))
    if not mp.isfinite(value):
        raise ValueError("decoded mpf value is non-finite")
    return value


def _canonical_json_bytes(payload: object) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _payload_sha256(payload: object) -> str:
    return hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()


def _reject_duplicate_json_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _kernel_sha256() -> str:
    digest = hashlib.sha256()
    package_dir = Path(__file__).resolve().parent
    for name in ("kernels.py", "operator.py", "runner.py"):
        raw = (package_dir / name).read_bytes()
        digest.update(name.encode("utf-8") + b"\0" + len(raw).to_bytes(8, "big") + raw)
    return digest.hexdigest()


def _gmpy2_versions() -> dict:
    """Return wrapper and native-library versions for the gmpy2 stack."""
    try:
        import gmpy2
    except ImportError:
        return {"gmpy2": None, "gmp": None, "mpfr": None, "mpc": None}

    def version(name: str) -> Optional[str]:
        getter = getattr(gmpy2, name, None)
        return str(getter()) if callable(getter) else None

    return {
        "gmpy2": str(getattr(gmpy2, "__version__", "unknown")),
        "gmp": version("mp_version"),
        "mpfr": version("mpfr_version"),
        "mpc": version("mpc_version"),
    }


def _compute_fingerprint() -> dict:
    import connes_cvs

    flint_version = None
    native_flint_version = None
    if HAS_FLINT:
        import flint

        flint_version = getattr(flint, "__version__", "unknown")
        native_flint_version = getattr(flint, "__FLINT_VERSION__", None)
    gmpy_versions = _gmpy2_versions()
    return {
        "connes_cvs": connes_cvs.__version__,
        "kernel_sha256": _kernel_sha256(),
        "mpmath": getattr(mp, "__version__", "unknown"),
        "mpmath_backend": getattr(mp.libmp, "BACKEND", "unknown"),
        **gmpy_versions,
        "python_flint": flint_version,
        "native_flint": native_flint_version,
        "has_flint": HAS_FLINT,
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "system": platform.system(),
        "machine": platform.machine(),
    }


def _hash_mpf_values(values) -> str:
    digest = hashlib.sha256()
    for value in values:
        raw = _canonical_json_bytes(_mpf_record(value))
        digest.update(len(raw).to_bytes(8, "big"))
        digest.update(raw)
    return digest.hexdigest()


@dataclass(frozen=True)
class CellConfig:
    """One explicit CvS Galerkin compute cell."""

    c: int
    N: int
    T: int
    dps: int
    tag: str = ""
    flint_bits: Optional[int] = None

    def __post_init__(self) -> None:
        _require_cutoff(self.c)
        _require_int("N", self.N, 1)
        _require_int("T", self.T, 1)
        _require_int("dps", self.dps, 15)
        if not isinstance(self.tag, str):
            raise TypeError(f"tag must be a string, got {type(self.tag).__name__}")
        try:
            tag_size = len(self.tag.encode("utf-8"))
        except UnicodeEncodeError as exc:
            raise ValueError("tag must be valid UTF-8 text") from exc
        if tag_size > 256:
            raise ValueError("tag must encode to at most 256 UTF-8 bytes")
        if self.flint_bits is not None:
            minimum_bits = (3322 * self.dps + 999) // 1000
            bits = _require_int("flint_bits", self.flint_bits, minimum_bits)
            object.__setattr__(self, "flint_bits", bits)

    @property
    def flint_prec(self) -> int:
        """Effective python-flint precision in bits.

        The default preserves the historical package convention. Published
        c=100 production artifacts used ``flint_bits=4*dps`` and reproduction
        configurations must state that explicitly.
        """
        return self.flint_bits if self.flint_bits is not None else int(self.dps * 3.5)


_worker_globals: Dict[str, Any] = {}


def _init_worker(c: int, dps: int, flint_prec: int, T: int) -> None:
    mp.mp.dps = dps
    if HAS_FLINT:
        flint_ctx.prec = flint_prec
        flint_ctx.threads = 1
    L = mp.log(c)
    prime_data, _ = prime_powers_up_to(c)
    _worker_globals.clear()
    _worker_globals.update(
        c=c, dps=dps, flint_prec=flint_prec, T=T, L=L, prime_data=prime_data
    )


def _worker_psi_point(n_idx: int) -> Tuple[int, tuple, tuple, float]:
    started = time.perf_counter()
    psi, psi_d = _compute_psi_pair(
        n_idx,
        _worker_globals["L"],
        _worker_globals["T"],
        _worker_globals["dps"],
        _worker_globals["prime_data"],
    )
    return n_idx, psi._mpf_, psi_d._mpf_, time.perf_counter() - started


def _default_progress(done: int, total: int, elapsed_s: float) -> None:
    rate = done / elapsed_s if elapsed_s > 0 else 0.0
    eta = (total - done) / rate if rate > 0 else float("inf")
    print(
        f"psi-cache {done:4d}/{total} ({100 * done / total:5.1f}%) "
        f"elapsed={elapsed_s:8.1f}s ETA={eta:8.1f}s",
        flush=True,
    )


def _atomic_json_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
    except BaseException:
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        raise


@contextmanager
def _exclusive_path_lock(target: Path):
    """Portable single-writer lease for a checkpoint or artifact target."""
    lock_dir = Path(str(target) + ".lock")
    try:
        lock_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise RuntimeError(
            f"another run may be using {target}; lock exists at {lock_dir}. "
            "Remove a stale lock only after confirming no writer is active."
        ) from exc
    try:
        _atomic_json_write(
            lock_dir / "owner.json",
            {
                "pid": os.getpid(),
                "host": platform.node(),
                "created_utc": datetime.now(timezone.utc).isoformat(),
            },
        )
        yield
    finally:
        try:
            (lock_dir / "owner.json").unlink()
        except FileNotFoundError:
            pass
        try:
            lock_dir.rmdir()
        except FileNotFoundError:
            pass


def _decompose_decimal(value: str) -> tuple[int, int, str]:
    text = value.strip()
    sign = -1 if text.startswith("-") else 1
    text = text.lstrip("+-")
    exponent_suffix = 0
    if "e" in text.lower():
        mantissa_text, exponent_text = text.lower().split("e", 1)
        exponent_suffix = int(exponent_text)
    else:
        mantissa_text = text
    integer, _, fraction = mantissa_text.partition(".")
    digits = integer + fraction
    stripped = digits.lstrip("0")
    if not stripped:
        return 0, 0, ""
    leading_zeros = len(digits) - len(stripped)
    exponent = len(integer) - leading_zeros - 1 + exponent_suffix
    return sign, exponent, stripped


def _matching_significant_digits(reference: str, value: str) -> int:
    sign_a, exponent_a, digits_a = _decompose_decimal(reference)
    sign_b, exponent_b, digits_b = _decompose_decimal(value)
    if not digits_a or not digits_b or sign_a != sign_b or exponent_a != exponent_b:
        return 0
    for index, (digit_a, digit_b) in enumerate(zip(digits_a, digits_b)):
        if digit_a != digit_b:
            return index
    return min(len(digits_a), len(digits_b))


class GalerkinCell:
    """Run, resume, and emit one explicit ``(c, N, T, dps)`` cell."""

    def __init__(
        self,
        cfg: CellConfig,
        workers: Optional[int] = None,
        progress_callback: Optional[Union[Callable[[int, int, float], None], bool]] = None,
        ground_state: str = "minimum",
        checkpoint_path: Optional[Union[str, os.PathLike]] = None,
        artifact_path: Optional[Union[str, os.PathLike]] = None,
        overwrite_artifact: bool = False,
    ) -> None:
        if not isinstance(cfg, CellConfig):
            raise TypeError(f"cfg must be CellConfig, got {type(cfg).__name__}")
        if workers is None:
            workers = min(cpu_count(), 8)
        self.workers = _require_int("workers", workers, 1)
        if progress_callback is False:
            self.progress_callback: Optional[Callable[[int, int, float], None]] = None
        elif progress_callback is None:
            self.progress_callback = _default_progress
        elif callable(progress_callback):
            self.progress_callback = progress_callback
        else:
            raise TypeError("progress_callback must be callable, None, or False")
        if ground_state not in {"smallest_positive", "minimum"}:
            raise ValueError(
                "ground_state must be 'smallest_positive' or 'minimum', "
                f"got {ground_state!r}"
            )
        self.cfg = cfg
        self.ground_state = ground_state
        self.checkpoint_path = Path(checkpoint_path).expanduser().resolve() if checkpoint_path else None
        self.artifact_path = Path(artifact_path).expanduser().resolve() if artifact_path else None
        if self.checkpoint_path is not None and self.checkpoint_path == self.artifact_path:
            raise ValueError("checkpoint_path and artifact_path must be different files")
        if not isinstance(overwrite_artifact, bool):
            raise TypeError("overwrite_artifact must be a bool")
        self.overwrite_artifact = overwrite_artifact
        self.compute_fingerprint = _compute_fingerprint()
        self.timings: Dict[str, float] = {}
        self.psi_vals: Optional[Dict[int, mp.mpf]] = None
        self.psi_deriv_vals: Optional[Dict[int, mp.mpf]] = None
        self.Q: Optional[mp.matrix] = None
        self.lambda_even: Optional[mp.mpf] = None
        self.eigvec_full: Optional[mp.matrix] = None
        self.spurious_negatives: List[dict] = []
        self.resumed_points = 0
        self.effective_workers = 1

    def _checkpoint_payload(self) -> dict:
        if self.psi_vals is None or self.psi_deriv_vals is None:
            raise RuntimeError("psi cache has not been initialized")
        body = {
            "schema": _CHECKPOINT_SCHEMA,
            "config": asdict(self.cfg),
            "flint_prec": self.cfg.flint_prec,
            "compute_fingerprint": self.compute_fingerprint,
            "values": {
                str(index): {
                    "psi": _mpf_record(self.psi_vals[index]),
                    "psi_deriv": _mpf_record(self.psi_deriv_vals[index]),
                }
                for index in sorted(self.psi_vals)
                if index >= 0
            },
        }
        return {**body, "payload_sha256": _payload_sha256(body)}

    def _write_checkpoint(self) -> None:
        if self.checkpoint_path is not None:
            self._verify_compute_fingerprint()
            _atomic_json_write(self.checkpoint_path, self._checkpoint_payload())

    def _verify_compute_fingerprint(self) -> None:
        if _compute_fingerprint() != self.compute_fingerprint:
            raise RuntimeError(
                "compute backend, versions, platform, or kernel source changed "
                "during this cell"
            )

    def _load_checkpoint(self) -> None:
        self.psi_vals = {}
        self.psi_deriv_vals = {}
        if self.checkpoint_path is None or not self.checkpoint_path.exists():
            return
        max_hex = max(64, int(self.cfg.dps * 0.84) + 64)
        max_bytes = 65_536 + (self.cfg.N + 1) * (2 * max_hex + 2_048)
        if self.checkpoint_path.stat().st_size > max_bytes:
            raise ValueError(
                f"checkpoint {self.checkpoint_path} exceeds its size bound"
            )
        try:
            payload = json.loads(
                self.checkpoint_path.read_text(encoding="utf-8"),
                object_pairs_hook=_reject_duplicate_json_keys,
            )
        except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
            raise ValueError(f"invalid checkpoint {self.checkpoint_path}: {exc}") from exc
        expected_fields = {
            "schema",
            "config",
            "flint_prec",
            "compute_fingerprint",
            "values",
            "payload_sha256",
        }
        if not isinstance(payload, dict) or set(payload) != expected_fields:
            raise ValueError("checkpoint must use the exact canonical field set")
        if payload.get("schema") != _CHECKPOINT_SCHEMA:
            raise ValueError(f"unsupported checkpoint schema in {self.checkpoint_path}")
        expected_checksum = payload.get("payload_sha256")
        body = {key: value for key, value in payload.items() if key != "payload_sha256"}
        if not isinstance(expected_checksum, str) or not hmac.compare_digest(
            expected_checksum, _payload_sha256(body)
        ):
            raise ValueError("checkpoint payload checksum mismatch")
        if payload.get("config") != asdict(self.cfg):
            raise ValueError("checkpoint configuration does not match this cell")
        if payload.get("flint_prec") != self.cfg.flint_prec:
            raise ValueError("checkpoint flint precision does not match this cell")
        if payload.get("compute_fingerprint") != self.compute_fingerprint:
            raise ValueError(
                "checkpoint compute backend, versions, platform, or kernel source "
                "does not match this run"
            )
        values = payload.get("values")
        if not isinstance(values, dict):
            raise ValueError("checkpoint values must be an object")
        for key, row in values.items():
            try:
                index = int(key)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid checkpoint index {key!r}") from exc
            if (
                str(index) != key
                or not 0 <= index <= self.cfg.N
                or not isinstance(row, dict)
                or set(row) != {"psi", "psi_deriv"}
            ):
                raise ValueError(f"invalid checkpoint row {key!r}")
            try:
                psi = _mpf_from_record(row.get("psi"), self.cfg.dps)
                psi_d = _mpf_from_record(row.get("psi_deriv"), self.cfg.dps)
            except ValueError as exc:
                raise ValueError(f"invalid checkpoint row {key!r}: {exc}") from exc
            self.psi_vals[index] = psi
            self.psi_deriv_vals[index] = psi_d
        self.resumed_points = len(self.psi_vals)

    def _store_positive_point(self, output: Tuple[int, tuple, tuple, float]) -> None:
        if self.psi_vals is None or self.psi_deriv_vals is None:
            raise RuntimeError("psi cache has not been initialized")
        index, psi_raw, psi_d_raw, _wall_s = output
        psi = mp.mpf(psi_raw)
        psi_d = mp.mpf(psi_d_raw)
        self.psi_vals[index] = psi
        self.psi_deriv_vals[index] = psi_d
        if index:
            self.psi_vals[-index] = -psi
            self.psi_deriv_vals[-index] = psi_d

    def run_phase_psi_cache(self) -> None:
        started = time.perf_counter()
        mp.mp.dps = self.cfg.dps
        if HAS_FLINT:
            flint_ctx.prec = self.cfg.flint_prec
            flint_ctx.threads = 1
        self._load_checkpoint()
        psi_vals = self.psi_vals
        psi_deriv_vals = self.psi_deriv_vals
        if psi_vals is None or psi_deriv_vals is None:
            raise RuntimeError("checkpoint loader did not initialize the psi cache")
        pending = [index for index in range(self.cfg.N + 1) if index not in psi_vals]
        total = self.cfg.N + 1
        done = total - len(pending)
        for index in sorted(i for i in psi_vals if i > 0):
            psi_vals[-index] = -psi_vals[index]
            psi_deriv_vals[-index] = psi_deriv_vals[index]

        if done and self.progress_callback is not None:
            self.progress_callback(done, total, time.perf_counter() - started)
        if not pending:
            self.timings["psi_cache_s"] = time.perf_counter() - started
            return

        _init_worker(self.cfg.c, self.cfg.dps, self.cfg.flint_prec, self.cfg.T)
        self.effective_workers = min(self.workers, len(pending), cpu_count())
        if self.effective_workers == 1:
            outputs = map(_worker_psi_point, pending)
            for output in outputs:
                self._store_positive_point(output)
                done += 1
                self._write_checkpoint()
                if self.progress_callback is not None:
                    self.progress_callback(done, total, time.perf_counter() - started)
        else:
            context = get_context("spawn")
            with context.Pool(
                self.effective_workers,
                initializer=_init_worker,
                initargs=(self.cfg.c, self.cfg.dps, self.cfg.flint_prec, self.cfg.T),
            ) as pool:
                for output in pool.imap_unordered(_worker_psi_point, pending, chunksize=1):
                    self._store_positive_point(output)
                    done += 1
                    self._write_checkpoint()
                    if self.progress_callback is not None:
                        self.progress_callback(done, total, time.perf_counter() - started)
        self.timings["psi_cache_s"] = time.perf_counter() - started

    def run_phase_matrix_assembly(self) -> mp.matrix:
        if self.psi_vals is None or self.psi_deriv_vals is None:
            raise RuntimeError("run_phase_psi_cache must complete first")
        started = time.perf_counter()
        mp.mp.dps = self.cfg.dps
        if HAS_FLINT:
            flint_ctx.prec = self.cfg.flint_prec
            flint_ctx.threads = 1
        N = self.cfg.N
        dimension = 2 * N + 1
        Q = mp.matrix(dimension, dimension)
        for i in range(dimension):
            m_index = i - N
            for j in range(i, dimension):
                n_index = j - N
                if m_index == n_index:
                    value = self.psi_deriv_vals[n_index]
                else:
                    value = (
                        self.psi_vals[m_index] - self.psi_vals[n_index]
                    ) / (m_index - n_index)
                Q[i, j] = value
                if i != j:
                    Q[j, i] = value
        self.Q = Q
        self.timings["matrix_assembly_s"] = time.perf_counter() - started
        return Q

    def run_phase_diagonalize(self, Q: mp.matrix) -> None:
        started = time.perf_counter()
        mp.mp.dps = self.cfg.dps
        N = self.cfg.N
        dimension = 2 * N + 1
        V_even = mp.matrix(dimension, N + 1)
        V_even[N, 0] = 1
        inv_sqrt2 = 1 / mp.sqrt(2)
        for k in range(1, N + 1):
            V_even[N + k, k] = inv_sqrt2
            V_even[N - k, k] = inv_sqrt2
        Q_even = V_even.T * Q * V_even
        eigenvalues, eigenvectors = mp.eigsy(Q_even)
        negative_indices = [index for index in range(N + 1) if eigenvalues[index] < 0]
        self.spurious_negatives = [
            {"index": index, "value": mp.nstr(eigenvalues[index], 20)}
            for index in negative_indices
        ]
        if self.ground_state == "smallest_positive":
            candidates = [index for index in range(N + 1) if eigenvalues[index] > 0]
            if not candidates:
                raise RuntimeError("even-sector spectrum has no positive eigenvalue")
            selected = min(candidates, key=lambda index: eigenvalues[index])
        else:
            selected = min(range(N + 1), key=lambda index: eigenvalues[index])
        self.lambda_even = eigenvalues[selected]
        vector_even = mp.matrix(N + 1, 1)
        for index in range(N + 1):
            vector_even[index, 0] = eigenvectors[index, selected]
        vector_even /= mp.sqrt(sum(vector_even[index, 0] ** 2 for index in range(N + 1)))
        vector_full = V_even * vector_even
        vector_full /= mp.sqrt(sum(vector_full[index, 0] ** 2 for index in range(dimension)))
        self.eigvec_full = vector_full
        self.timings["diagonalize_s"] = time.perf_counter() - started

    def _c13_regression(self) -> dict:
        if self.lambda_even is None:
            raise RuntimeError("diagonalization must complete first")
        computed = mp.nstr(self.lambda_even, 39)
        matches = _matching_significant_digits(_C13_REFERENCE, computed)
        result = {
            "reference": _C13_REFERENCE,
            "computed": computed,
            "matching_digits": matches,
            "threshold": 22,
            "passed": matches >= 22,
        }
        if not result["passed"]:
            raise RuntimeError(
                "c=13 reference-cell regression failed: "
                f"{matches} matching significant digits, 22 required"
            )
        return result

    def _artifact(self) -> dict:
        if self.lambda_even is None or self.eigvec_full is None or self.Q is None:
            raise RuntimeError("cell has not completed")
        self._verify_compute_fingerprint()
        compute_fingerprint = self.compute_fingerprint
        artifact = {
            "schema": _ARTIFACT_SCHEMA,
            "config": {**asdict(self.cfg), "flint_prec": self.cfg.flint_prec},
            "lambda_even": mp.nstr(self.lambda_even, self.cfg.dps),
            "lambda_even_raw_mpf": _mpf_record(self.lambda_even),
            "matrix_shape": [self.Q.rows, self.Q.cols],
            "matrix_sha256": _hash_mpf_values(
                self.Q[i, j] for i in range(self.Q.rows) for j in range(self.Q.cols)
            ),
            "eigenvector_shape": [self.eigvec_full.rows, self.eigvec_full.cols],
            "eigenvector_sha256": _hash_mpf_values(
                self.eigvec_full[i, 0] for i in range(self.eigvec_full.rows)
            ),
            "timings_seconds": dict(self.timings),
            "timing_scope": (
                "compute pipeline through diagonalization; excludes artifact "
                "hashing, JSON serialization, and optional disk-write overhead"
            ),
            "workers": self.workers,
            "effective_workers": self.effective_workers,
            "resumed_points": self.resumed_points,
            "engine": f"connes_cvs.runner v{compute_fingerprint['connes_cvs']}",
            "compute_fingerprint": compute_fingerprint,
            "ground_state_selection": self.ground_state,
            "negative_eigenvalues": list(self.spurious_negatives),
            "environment": _environment(),
            "validation_scope": (
                "Finite-cell computation only. Working precision must be "
                "validated by independent reference or cross-precision runs."
            ),
        }
        if (self.cfg.c, self.cfg.N, self.cfg.T, self.cfg.dps) == (13, 100, 800, 150):
            artifact["c13_regression"] = self._c13_regression()
        artifact["artifact_payload_sha256"] = _payload_sha256(artifact)
        return artifact

    def run(self) -> dict[str, Any]:
        if (
            self.artifact_path is not None
            and self.artifact_path.exists()
            and not self.overwrite_artifact
        ):
            raise FileExistsError(
                f"artifact already exists: {self.artifact_path}; pass "
                "overwrite_artifact=True to replace it atomically"
            )
        lock_targets = sorted(
            {path for path in (self.checkpoint_path, self.artifact_path) if path},
            key=str,
        )
        with ExitStack() as stack:
            for target in lock_targets:
                stack.enter_context(_exclusive_path_lock(target))
            # Recheck after acquiring the writer lease: a competing run may
            # have created the artifact between the optimistic check above and
            # this lock acquisition.
            if (
                self.artifact_path is not None
                and self.artifact_path.exists()
                and not self.overwrite_artifact
            ):
                raise FileExistsError(
                    f"artifact already exists: {self.artifact_path}; pass "
                    "overwrite_artifact=True to replace it atomically"
                )
            started = time.perf_counter()
            process_started = time.process_time()
            self.run_phase_psi_cache()
            self.run_phase_diagonalize(self.run_phase_matrix_assembly())
            total_wall = time.perf_counter() - started
            self.timings["total_s"] = total_wall
            self.timings["total_wall_monotonic_s"] = total_wall
            self.timings["coordinator_process_cpu_s"] = (
                time.process_time() - process_started
            )
            artifact = self._artifact()
            if self.artifact_path is not None:
                _atomic_json_write(self.artifact_path, artifact)
            return artifact


def _environment() -> dict:
    import connes_cvs

    flint_version = None
    native_flint_version = None
    if HAS_FLINT:
        import flint

        flint_version = getattr(flint, "__version__", "unknown")
        native_flint_version = getattr(flint, "__FLINT_VERSION__", None)
    gmpy_versions = _gmpy2_versions()
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": cpu_count(),
        "mpmath": getattr(mp, "__version__", "unknown"),
        "mpmath_backend": getattr(mp.libmp, "BACKEND", "unknown"),
        **gmpy_versions,
        "python_flint": flint_version,
        "native_flint": native_flint_version,
        "connes_cvs": connes_cvs.__version__,
        "has_flint": HAS_FLINT,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }


def run_cell(
    cfg: CellConfig,
    workers: Optional[int] = None,
    progress_callback: Optional[Union[Callable[[int, int, float], None], bool]] = None,
    ground_state: str = "minimum",
    checkpoint_path: Optional[Union[str, os.PathLike]] = None,
    artifact_path: Optional[Union[str, os.PathLike]] = None,
    overwrite_artifact: bool = False,
) -> dict[str, Any]:
    """Run a cell and return its JSON-serializable artifact."""
    return GalerkinCell(
        cfg,
        workers=workers,
        progress_callback=progress_callback,
        ground_state=ground_state,
        checkpoint_path=checkpoint_path,
        artifact_path=artifact_path,
        overwrite_artifact=overwrite_artifact,
    ).run()
