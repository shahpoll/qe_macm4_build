# Troubleshooting FAQ

See also [`docs/Troubleshooting.md`](../Troubleshooting.md) for the full in-repo version.

## `make[1]: laxlib*.fh ... no input files`

- Symptom: preprocessing fails while building `LAXlib`.
- Fix: rerun `./configure` with `CPP="gcc -E"` (or export `CPP` globally) so that CLT 16 uses GNU cpp semantics.

## `wxml.f90: DP_XML has no implicit type`

- Symptom: CMake/OpenBLAS build fails compiling `upflib/wxml.f90`.
- Cause: QE 7.4.1 expects FoX sources when `__XML_STANDALONE` is off; cloning FoX or removing the `__XML_STANDALONE` define resolves it.

## `pp.x`/`pw.x` cannot find pseudopotentials

- Ensure pseudopotentials live in `./pp` (e.g., `Si.pbe-n-rrkjus_psl.1.0.0.UPF`).
- The provided inputs assume relative paths (`pseudo_dir='./pp'`).

## `projwfc.x` missing

- By default `make pw` does not build post-processing tools. Run `make -C PP/src projwfc.x dos.x bands.x` in the QE source tree.

## Where is Matplotlib?

- `pip --user` installs to `~/Library/Python/3.13/bin`. Add this directory to `PATH` if `python3 scripts/plot_si_*.py` cannot find `matplotlib`.

## GPU acceleration?

- Not on Apple Silicon. QE supports CUDA/ROCm only; Metal GPUs are not yet supported upstream.
