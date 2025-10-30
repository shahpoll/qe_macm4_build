# Troubleshooting (macOS 15.3, Apple Silicon)
- `make[1]: laxlib*.fh ... no input files`: rerun configure with `CPP="gcc -E"`. Ref: QE user guide, “Installation tricks and problems”.
- MPI quirks: prefer Homebrew OpenMPI bottles for Sequoia; if runtime issues persist, document and consider MPICH after removing OpenMPI.
- Linker proof: keep `otool -L` for `pw.x` to show OpenBLAS or Accelerate and MPI.
- If PP download is blocked, fetch manually from the PSLibrary element page and retry the Si run.
