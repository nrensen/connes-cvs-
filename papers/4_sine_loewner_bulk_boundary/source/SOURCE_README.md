# Bulk trace laws and algebraic boundary interaction for the sine Loewner operator

Author: Akiva Groskin. September 2026.

## Build the complete paper

Run from this directory with a standard TeX Live installation:

    pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
    pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex

The bibliography is included in references.tex. All mathematical statements and proofs are compiled from sections/. The figures in figures/ are TikZ inputs; no external image or shell-escape dependency is required. The source uses standard AMS, Latin Modern, microtype, TikZ and hyperref packages. It was built with pdfTeX 1.40.29 / TeX Live 2026. A disabled-shell-escape notice from epstopdf is harmless; no image conversion is requested.

## Verification programs and exact witnesses

Read anc/README.md. The standard replay comprises 23 programs: 18 floating-point diagnostics, two exact rational scalar checks and three outward-rounded interval certificate programs. Python 3.12 and the pinned libraries in anc/requirements.txt are required. The certificate proofs, including full operator error estimates, are in the paper. The ancillary programs are supporting reproducibility material, not a mathematical supplement.

The source archive is self-contained. No companion paper source, private repository, external mathematical appendix, or pre-generated bibliography is required. LICENSES.md states the manuscript and code licenses.
