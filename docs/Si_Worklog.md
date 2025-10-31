# Silicon Worklog (QE 7.4.1 on macOS 15.3 / M4)

## Environment reminders
- Toolchain: Homebrew GCC 15.2.0, OpenMPI 5.0.8, veclibfort 0.4.3 (Accelerate), OpenBLAS 0.3.30.
- QE source: `artifacts/q-e-qe-7.4.1` (configure build).
- Pseudopotential: `pp/Si.pbe-n-rrkjus_psl.1.0.0.UPF` (PSLibrary).

## Command history
Timestamps correspond to the most recent run captured in the log files.

1. **Self-consistent field (SCF) baseline** — 30 Oct 2025 15:21 local time
   ```sh
   mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.scf.in \
     | tee logs/si_scf_attemptB2.txt
   ```
   - Output copied to `analysis/Si/logs/si_scf_attemptB2.txt`.
   - Produces `tmp/silicon.save`, energy convergence in five iterations, Fermi level 6.4576 eV.

2. **Band-path SCF (restarts from save)** — 31 Oct 2025 15:07 local time
   ```sh
   mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.bands.in \
     | tee analysis/Si/logs/si_bands_pw.txt
   ```
   - Reuses charge density from the prior SCF (`restart_mode='from_scratch'` with saved state).
   - Traverses Γ–X–W–K–Γ–L (6 segments, 40 points each).

3. **bands.x post-processing** — 31 Oct 2025 15:08 local time
   ```sh
   artifacts/q-e-qe-7.4.1/bin/bands.x -in inputs/Si.bands_post.in \
     | tee analysis/Si/logs/si_bands_post.txt
   ```
   - Writes eigenvalues to `silicon.bands.dat` and `silicon.bands.dat.gnu`.
   - Files archived in `analysis/Si/data/` for plotting.

4. **Band plot generation** — 31 Oct 2025 15:09 local time
   ```sh
   python3 -m pip install --user matplotlib  # one-time dependency
   python3 scripts/plot_si_bands.py
   ```
   - Generates `analysis/Si/plots/si_band_structure.png` with Γ/X/W/K/Γ/L ticks and Fermi level overlay.

5. **Dense k-point NSCF run (for DOS)** — 31 Oct 2025 19:05 local time
   ```sh
   mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.nscf.in \
     | tee analysis/Si/logs/si_nscf_pw.txt
   ```
   - Uses a 12×12×12 Monkhorst–Pack mesh, 16 bands, fixed occupations.
   - Updates `tmp/silicon.save` with eigenvalues suitable for post-processing.

6. **Total DOS via dos.x** — 31 Oct 2025 19:06 local time
   ```sh
   artifacts/q-e-qe-7.4.1/bin/dos.x -in inputs/Si.dos.in \
     | tee analysis/Si/logs/si_dos.txt
   ```
   - Produces `analysis/Si/data/silicon.dos` with energies, DOS, and integrated DOS.

7. **DOS plot generation** — 31 Oct 2025 19:07 local time
   ```sh
   python3 scripts/plot_si_dos.py
   ```
   - Creates `analysis/Si/plots/si_total_dos.png` with Fermi-level marker.

8. **Orbital projections (`projwfc.x`)** — 31 Oct 2025 19:32 local time
   ```sh
   artifacts/q-e-qe-7.4.1/bin/projwfc.x -in inputs/Si.projwfc.in \
     | tee analysis/Si/logs/si_projwfc.txt
   ```
   - Generates Lowdin charges and PDOS files (`analysis/Si/data/silicon.pdos_*`).

9. **PDOS plot generation** — 31 Oct 2025 19:33 local time
   ```sh
   python3 scripts/plot_si_pdos.py
   ```
   - Produces `analysis/Si/plots/si_pdos.png` (total DOS plus s/p components).

10. **Band-gap summary extraction** — 31 Oct 2025 19:34 local time
    ```sh
    python3 scripts/analyze_si_bandgap.py
    ```
    - Writes `analysis/Si/data/si_band_summary.txt` and prints indirect/direct gap metrics.

## Cautions and notes
- The configure build expects Apple’s `cpp` to handle `.F90`; when it mangles `laxlib_*.h`, rerun `./configure` with `CPP="gcc -E"` (already baked into Attempt B2).
- `bands.x` is not built by default in a minimal `make pw`; we invoked `make -C PP/src bands.x` to populate `bin/bands.x`.
- Matplotlib installs under `~/Library/Python/3.13/bin`; add this directory to `PATH` if you want `f2py`/`ttx` on the CLI.
- Floating-point underflow warnings are expected for ultrasoft PSPs in QE and can be ignored unless they grow during convergence.

## Platform perspective
- **macOS Sequoia vs older macOS:** The CLT 16 toolchain defaults to `/usr/bin/cpp`; add `CPP="gcc -E"` for configure builds. On Monterey/Ventura the flag is redundant but benign.
- **macOS vs Ubuntu:** Accelerate (via `veclibfort`) replaces OpenBLAS/MKL and requires no `LD_LIBRARY_PATH` tuning. Ubuntu users typically stick with OpenBLAS and GNU `cpp`.
- **Binary toolchain:** OpenMPI and dependencies live under `/opt/homebrew` with wrapper scripts; this differs from the `/usr/lib` layout on Linux distros.
- **GPU acceleration:** QE exposes CUDA/ROCm only. Apple Metal GPUs (M1–M4) remain unsupported, so all runs here are CPU-only.
