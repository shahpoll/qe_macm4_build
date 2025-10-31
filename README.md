# QE on Apple Silicon (M4)

> Reproducible Quantum ESPRESSO 7.4.1 workflow on macOS 15.3 (Sequoia) using Homebrew GCC + OpenMPI + veclibfort (Accelerate).

![Si band structure](analysis/Si/plots/si_band_structure.png)

## Contents

1. [Quick start](#quick-start)
2. [Environment snapshots](#environment-snapshots)
3. [Silicon example](#silicon-example)
4. [Helper scripts](#helper-scripts)
5. [Docs & Wiki](#docs--wiki)

## Quick start

```sh
git clone https://github.com/shahpoll/qe_macm4_build.git
cd qe_macm4_build
brew install gcc open-mpi cmake veclibfort openblas wget python
./docs/AppleSilicon_QE_Guide.md   # follow step-by-step playbook
```

- **Accelerate toolchain** already configured under `artifacts/q-e-qe-7.4.1`.
- **SCF sanity check**: `mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.scf.in`.
- **Report bundle**: `qe_macm4_attempt_bundle_v3.zip` collects logs, docs, plots.

## Environment snapshots

- `logs/system.txt` — OS/kernel/date, CLT path, PATH content.
- `logs/versions.txt` — Homebrew packages, `gfortran`, `mpirun`, `cmake` versions.

## Silicon example

Artifacts reside under [`analysis/Si`](analysis/Si/README.md):

| File | Description |
|------|-------------|
| `inputs/Si.scf.in` | Baseline SCF input (Γ-point).
| `inputs/Si.nscf.in` | 12×12×12 NSCF grid for DOS/PDOS.
| `analysis/Si/logs/*.txt` | Full QE transcripts (SCF → NSCF → bands → DOS → projwfc).
| `analysis/Si/data/silicon.dos` | Total DOS sampled from -15 to 20 eV.
| `analysis/Si/data/si_band_summary.txt` | Fermi level and indirect/direct gaps.
| `analysis/Si/plots/si_total_dos.png` | DOS (Accelerate build).
| `analysis/Si/plots/si_pdos.png` | s/p-projected DOS per atom.

> ![Si total DOS](analysis/Si/plots/si_total_dos.png)

## Helper scripts

Located in [`scripts/`](scripts):

- `plot_si_bands.py` — render Γ–X–W–K–Γ–L band structure.
- `plot_si_dos.py` — total DOS visualization (`silicon.dos`).
- `plot_si_pdos.py` — overlay s/p PDOS from `projwfc.x` outputs.
- `analyze_si_bandgap.py` — compute indirect/direct gaps from band data.

All scripts use Matplotlib; install once via `python3 -m pip install --user matplotlib`.

## Docs & Wiki

- `docs/AppleSilicon_QE_Guide.md` — full setup guide (macOS Sequoia + CLT 16).
- `docs/Guide.md` — quick checklist.
- `docs/Si_Worklog.md` — chronological command log with macOS vs Ubuntu notes.
- `docs/Troubleshooting.md` — known issues (CPP override, FoX, Metal GPU status).

**Wiki suggestions:**

- Environment notes (OS/CLT quirks, veclibfort vs OpenBLAS).
- Troubleshooting FAQ mirrored from docs.
- Workflow recipes (Accelerate path today, OpenBLAS/CMake when fixed).
- Results gallery (embed band/DOS plots, band-gap summary).
