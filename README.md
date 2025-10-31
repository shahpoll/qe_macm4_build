# qe_macm4_build

Reproducible Quantum ESPRESSO 7.4.1 builds on an Apple Silicon Mac mini (M4). The repository captures:

- Toolchain captures (`logs/versions.txt`, `logs/system.txt`).
- Accelerate-based configure build logs and working binaries (`artifacts/q-e-qe-7.4.1`).
- Silicon validation workflow under `analysis/Si` with inputs, raw outputs, and plots.
- Guides:
  - `docs/AppleSilicon_QE_Guide.md` — end-to-end setup playbook.
  - `docs/Guide.md` — quick-start checklist.
  - `docs/Si_Worklog.md` — chronological notes with commands and caveats.
- Helper scripts live under `scripts/` (band plotting, DOS/PDOS visuals, band-gap summary).

Start with `docs/AppleSilicon_QE_Guide.md` to reproduce the exact environment. The Silicon example doubles as a regression test and a template for future case studies.
