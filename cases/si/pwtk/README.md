# Silicon (Si) via PWTK — QE 7.4.1 / macOS 15.3

This folder mirrors the manual QE workflow using PWTK scripts.

## Artifacts

- `scripts/si_workflow.pwtk` — orchestrates SCF → NSCF → BANDS → DOS → PDOS.
- `logs/pwtk_run.log` — stdout from the PWTK execution.
- QE inputs/outputs:
  - `pw.si.scf.in`, `pw.si.scf.out`, ...
  - `dos.si.in`, `dos.si.out`, etc.
- `data/` — consolidated results (`silicon.dos`, `silicon.bands.dat`, `silicon_pdos.*`, `si_band_summary.txt`).
- `plots/` — visuals regenerated from the data (band structure, total DOS, PDOS).

## Workflow reference

1. Ensure Tcl 8.6 and PWTK 3.2 are on the PATH:
   ```sh
   export PATH="/opt/homebrew/opt/tcl-tk@8/bin:$PATH"
   export PATH="/Users/pollob/qe_macm4_build/pwtk-3.2:$PATH"
   ```
2. Run the PWTK script:
   ```sh
   cd analysis_pwtk/Si
   pwtk scripts/si_workflow.pwtk | tee logs/pwtk_run.log
   ```
3. Generate plots and summaries:
   ```sh
   cd /Users/pollob/qe_macm4_build
   python3 scripts/plot_si_bands.py --pwtk
   python3 scripts/plot_si_dos.py   --pwtk
   python3 scripts/plot_si_pdos.py  --pwtk
   python3 scripts/analyze_si_bandgap.py --pwtk
   ```

Results can now be compared directly with the manual QE run under `analysis/Si/`.
