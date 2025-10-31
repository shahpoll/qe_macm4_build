# Workflow Recipes

This page mirrors the commands logged in `docs/Si_Worklog.md`, focusing on the silicon validation path.

## 1. Configure build (Accelerate + MPI)

```sh
cd artifacts/q-e-qe-7.4.1
./configure MPIF90=mpif90 CC=mpicc CPP="gcc -E" \
  BLAS_LIBS="-L$(brew --prefix veclibfort)/lib -lvecLibFort -framework Accelerate" \
  LAPACK_LIBS="-L$(brew --prefix veclibfort)/lib -lvecLibFort -framework Accelerate"
make -j pw
make -C PP/src bands.x dos.x projwfc.x
cd ../../
```

> `CPP="gcc -E"` is the Sequoia-specific fix. The final `make -C PP/src …` pulls in the post-processing binaries used later.

## 2. Silicon baseline (SCF)

```sh
mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.scf.in \
  | tee analysis/Si/logs/si_scf_attemptB2.txt
```

Outputs: `tmp/silicon.save`, plus the transcript in `analysis/Si/logs/si_scf_attemptB2.txt`.

## 3. Band structure path

```sh
mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.bands.in \
  | tee analysis/Si/logs/si_bands_pw.txt
artifacts/q-e-qe-7.4.1/bin/bands.x -in inputs/Si.bands_post.in \
  | tee analysis/Si/logs/si_bands_post.txt
python3 scripts/plot_si_bands.py
```

Artifacts: `analysis/Si/data/silicon.bands.dat`, `.gnu`, and the plot `analysis/Si/plots/si_band_structure.png`.

## 4. Dense DOS sampling

```sh
mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.nscf.in \
  | tee analysis/Si/logs/si_nscf_pw.txt
artifacts/q-e-qe-7.4.1/bin/dos.x -in inputs/Si.dos.in \
  | tee analysis/Si/logs/si_dos.txt
python3 scripts/plot_si_dos.py
```

Artifacts: `analysis/Si/data/silicon.dos` and `analysis/Si/plots/si_total_dos.png`.

## 5. Projected DOS (PDOS)

```sh
artifacts/q-e-qe-7.4.1/bin/projwfc.x -in inputs/Si.projwfc.in \
  | tee analysis/Si/logs/si_projwfc.txt
python3 scripts/plot_si_pdos.py
```

Artifacts: `analysis/Si/data/silicon.pdos_*` plus `analysis/Si/plots/si_pdos.png`.

## 6. Band-gap summary

```sh
python3 scripts/analyze_si_bandgap.py
```

Artifact: `analysis/Si/data/si_band_summary.txt`, reporting indirect/direct gaps using the band data and Fermi level from `silicon.dos`.

## Optional: CMake + OpenBLAS (experimental)

Uses the same source tree (`artifacts/q-e-qe-7.4.1`) with FoX cloned in `external/fox/`. Caveat: ensure `target_compile_definitions(qe_xml PRIVATE __XML_STANDALONE)` is removed if you rely on FoX.

```sh
cd artifacts/q-e-qe-7.4.1
OBLAS=$(brew --prefix openblas)
cmake -S . -B buildA \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DCMAKE_C_COMPILER=mpicc \
  -DQE_ENABLE_MPI=ON \
  -DQE_ENABLE_OPENMP=ON \
  -DBLAS_LIBRARIES="$OBLAS/lib/libopenblas.dylib" \
  -DLAPACK_LIBRARIES="$OBLAS/lib/libopenblas.dylib"
cmake --build buildA --target pw -j
```

> This configuration is left in the repo for parity with Linux builds but still needs the FoX toggle for wxml compilation.
