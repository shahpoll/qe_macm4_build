# Si Workflow via PWTK (attempt)

## 2025-11-01
- Downloaded `pwtk-3.2.tar.gz` and extracted under `/Users/pollob/qe_macm4_build/pwtk-3.2`.
- Created PWTK config `~/.pwtk/pwtk.tcl` pointing to QE bins (`artifacts/q-e-qe-7.4.1/bin`) and pseudopotentials (`pp/`).
- Installed Homebrew `tcl-tk` (9.0.2) to provide a modern Tcl interpreter.
- Downloaded `tcllib-1.21` for optional dependencies.
- Authored `analysis_pwtk/Si/scripts/si_workflow.pwtk` mirroring the SCF → NSCF → BANDS → DOS → PDOS pipeline.
- Repeatedly attempted to execute the script:
  ```sh
  cd analysis_pwtk/Si
  PATH="/Users/pollob/qe_macm4_build/pwtk-3.2:$PATH" \
    ../../pwtk-3.2/pwtk scripts/si_workflow.pwtk
  ```
- Execution aborts during PWTK initialization with missing Tcl commands (`print`, `varvalue`, `time2ms`, `getExecutable`, ...). These utilities are normally provided by bundled Tcl modules; the clean 3.2 tarball does not ship precompiled tcllib on macOS, and the interpreter exits before loading PWTK’s helper libraries.
- Given the chained failures, the PWTK-based workflow could not be executed in the current macOS environment without substantial patching of PWTK internals.

## Next steps
1. Obtain a prebuilt Tcl/Tk distribution that already bundles `tcllib`, or install the PWTK-recommended packages (`tcllib`, `tcl-tclreadline`, `tcl-thread`) from source and verify they register commands like `print`, `varvalue`, `time2ms`, and `getExecutable` before PWTK initialization.
2. Re-run `scripts/si_workflow.pwtk` once PWTK starts cleanly; capture output in `analysis_pwtk/Si/logs/pwtk_run.log`.
3. Mirror the outputs/plots into `analysis_pwtk/Si/data` and `analysis_pwtk/Si/plots` for comparison with the hand-driven QE workflow.

Current state: **blocked** on missing Tcl procedures during PWTK start-up.
