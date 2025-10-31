# Quantum ESPRESSO on Apple Silicon — Minimal Guide
1. Install Homebrew toolchain: `brew install gcc open-mpi cmake veclibfort openblas wget`.
2. Fetch QE 7.4.1 sources into `artifacts/q-e-qe-7.4.1` (e.g. `git clone --depth 1 --branch qe-7.4.1 https://gitlab.com/QEF/q-e.git q-e-qe-7.4.1`).
3. Configure with CMake + OpenBLAS:
   ```sh
   cd artifacts/q-e-qe-7.4.1
   OBLAS="$(brew --prefix openblas)"
   cmake -S . -B buildA \
     -DCMAKE_Fortran_COMPILER=mpif90 \
     -DCMAKE_C_COMPILER=mpicc \
     -DQE_ENABLE_MPI=ON \
     -DQE_ENABLE_OPENMP=ON \
     -DBLAS_LIBRARIES="$OBLAS/lib/libopenblas.dylib" \
     -DLAPACK_LIBRARIES="$OBLAS/lib/libopenblas.dylib"
   cmake --build buildA -j --target pw
   ```
4. Pull a UPF (example): `curl -L -o pp/Si.pbe-n-rrkjus_psl.1.0.0.UPF https://pseudopotentials.quantum-espresso.org/upf_files/Si.pbe-n-rrkjus_psl.1.0.0.UPF`.
5. Run a quick SCF test with the configure+Accelerate build: `mpirun -n 2 artifacts/q-e-qe-7.4.1/bin/pw.x -in inputs/Si.scf.in`.
6. Smoke-test the CMake binary (expect it to stop waiting for stdin if no `-in` flag is passed): `buildA/bin/pw.x -v`.
7. Optional fallback: `./configure MPIF90=mpif90 CC=mpicc BLAS_LIBS="-L$(brew --prefix veclibfort)/lib -lvecLibFort -framework Accelerate" LAPACK_LIBS="-L$(brew --prefix veclibfort)/lib -lvecLibFort -framework Accelerate"`; add `CPP="gcc -E"` if preprocessing fails on `laxlib_*.h`.

### Commands logged
```
cmake -S . -B buildA ...
cmake --build buildA -j --target pw
./configure MPIF90=mpif90 CC=mpicc ...
make -j pw
```
