# Thread 7 — Mac mini M4 build attempt (negative result if failed)

**BLAS/LAPACK chosen:** Accelerate+veclibfort (fallback Attempt A uses OpenBLAS)
**MPI:** Homebrew OpenMPI
**Rosetta used:** No

## Packages
cmake 4.1.2
gcc 15.2.0
open-mpi 5.0.8
openblas 0.3.30
veclibfort 0.4.3

## Attempt A — CMake + OpenBLAS
```
=== Attempt A: configure ===
-- The Fortran compiler identification is GNU 15.2.0
-- The C compiler identification is AppleClang 17.0.0.17000013
-- Checking whether Fortran compiler has -isysroot
-- Checking whether Fortran compiler has -isysroot - yes
-- Checking whether Fortran compiler supports OSX deployment target flag
-- Checking whether Fortran compiler supports OSX deployment target flag - yes
-- Detecting Fortran compiler ABI info
-- Detecting Fortran compiler ABI info - done
-- Check for working Fortran compiler: /opt/homebrew/bin/mpif90 - skipped
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /opt/homebrew/bin/mpicc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Setting build type to 'Release' as none was specified
-- Enable sanitizer QE_ENABLE_SANITIZER=none
-- C preprocessor used by qe_preprocess_source in qeHelpers.cmake: /usr/bin/cpp
-- Performing Test Fortran_ISYSTEM_SUPPORTED
-- Performing Test Fortran_ISYSTEM_SUPPORTED - Success
-- Found OpenMP_Fortran: -fopenmp (found version "4.5")
-- Found OpenMP: TRUE (found version "4.5") found components: Fortran
-- Found MPI_Fortran: /opt/homebrew/bin/mpif90 (found version "3.1")
-- Found MPI: TRUE (found version "3.1") found components: Fortran
-- Selected the Fortran 'mpi' module. QE_ENABLE_MPI_MODULE=ON
-- MPI settings used by CTest
     MPIEXEC_EXECUTABLE : /opt/homebrew/bin/mpiexec
     MPIEXEC_NUMPROC_FLAG : -n
     MPIEXEC_PREFLAGS : 
   Tests run as : /opt/homebrew/bin/mpiexec -n <NUM_PROCS>  <EXECUTABLE>
-- Found Git: /opt/homebrew/bin/git (found suitable version "2.50.1", minimum required is "2.13")
-- Source files are cloned from a git repository.
   sed supports -E
   Git branch: HEAD
   Git commit hash: 500de340b820e1cb8c05f2d8bb8fced102f377c1
-- Trying to find LAPACK from ARM Performance Library
-- Found BLAS: /opt/homebrew/opt/openblas/lib/libopenblas.dylib
-- Found LAPACK: /opt/homebrew/opt/openblas/lib/libopenblas.dylib
-- Looking for Fortran zhpev
-- Looking for Fortran zhpev - found
-- Installing Wannier90 via submodule
Submodule 'external/wannier90' (https://github.com/wannier-developers/wannier90.git) registered for path 'external/wannier90'
Cloning into '/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/external/wannier90'...
From https://github.com/wannier-developers/wannier90
 * branch            1d6b187374a2d50b509e5e79e2cab01a79ff7ce1 -> FETCH_HEAD
Submodule path 'external/wannier90': checked out '1d6b187374a2d50b509e5e79e2cab01a79ff7ce1'
-- Installing MBD via submodule
Submodule 'external/mbd' (https://github.com/libmbd/libmbd.git) registered for path 'external/mbd'
Cloning into '/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/external/mbd'...
Submodule path 'external/mbd': checked out '89a3cc199c0a200c9f0f688c3229ef6b9a8d63bd'
-- Found Git: /opt/homebrew/bin/git (found version "2.50.1")
-- Setting version tag to 0.13.0-2-g89a3cc1 from Git
-- Installing DeviceXlib via submodule
Submodule 'external/devxlib' (https://gitlab.com/max-centre/components/devicexlib.git) registered for path 'external/devxlib'
Cloning into '/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/external/devxlib'...
From https://gitlab.com/max-centre/components/devicexlib
 * branch            a6b89ef77b1ceda48e967921f1f5488d2df9226d -> FETCH_HEAD
Submodule path 'external/devxlib': checked out 'a6b89ef77b1ceda48e967921f1f5488d2df9226d'
-- Could NOT find VendorFFTW (missing: VendorFFTW_LIBRARIES VendorFFTW_INCLUDE_DIRS VendorFFTW_ID) 
-- Found PkgConfig: /opt/homebrew/bin/pkg-config (found version "2.5.1")
...
      |                                                                           1
Warning: More actual than formal arguments in procedure call at (1)
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/input.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/loc_scdm_k.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/electrons.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/punch.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/read_file_new.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/wfcinit.f90.o
[100%] Linking Fortran static library ../lib/libqe_pw.a
[100%] Built target qe_pw
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_ev_exe.dir/tools/ev.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_cell2ibrav_exe.dir/tools/cell2ibrav.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_scan_ibrav_exe.dir/tools/scan_ibrav.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_pwi2xsf_exe.dir/tools/pwi2xsf.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_kpoints_exe.dir/tools/kpoints.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_exe.dir/src/pwscf.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_ibrav2cell_exe.dir/tools/ibrav2cell.f90.o
[100%] Linking Fortran executable ../bin/pw.x
[100%] Linking Fortran executable ../bin/cell2ibrav.x
[100%] Linking Fortran executable ../bin/ibrav2cell.x
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_lax.a', '../lib/libqe_modules.a', '../lib/libqe_upflib.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
[100%] Linking Fortran executable ../bin/scan_ibrav.x
[100%] Linking Fortran executable ../bin/pwi2xsf.x
[100%] Linking Fortran executable ../bin/kpoints.x
[100%] Built target qe_pw_tools_cell2ibrav_exe
[100%] Built target qe_pw_tools_ibrav2cell_exe
[100%] Built target qe_pw_exe
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
[100%] Built target qe_pw_tools_pwi2xsf_exe
[100%] Built target qe_pw_tools_scan_ibrav_exe
[100%] Built target qe_pw_tools_kpoints_exe
[100%] Linking Fortran executable ../bin/ev.x
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
[100%] Built target qe_pw_tools_ev_exe
[100%] basic code for scf, structure optimization, MD
[100%] Built target pw
```

## Attempt B — configure + Accelerate
```
=== Attempt B: configure ===
-n directory LAXlib : ok
-n directory FFTXlib/src : ok
-n directory UtilXlib : ok
-n directory dft-d3 : ok
-n directory KS_Solvers/Davidson : ok
-n directory KS_Solvers/Davidson_RCI : ok
-n directory KS_Solvers/CG : ok
-n directory KS_Solvers/PPCG : ok
-n directory KS_Solvers/ParO : ok
-n directory KS_Solvers/DENSE : ok
-n directory KS_Solvers/RMM : ok
-n directory upflib : ok
-n directory XClib : ok
-n directory Modules : ok
-n directory LR_Modules : ok
-n directory PW/src : ok
-n directory CPV/src : ok
-n directory PW/tools : ok
-n directory PP/src : ok
-n directory PWCOND/src : ok
-n directory PHonon/Gamma : ok
-n directory PHonon/PH : ok
-n directory PHonon/FD : ok
-n directory HP/src : ok
-n directory atomic/src : ok
-n directory EPW/src : ok
-n directory EPW/ZG/src : ok
-n directory XSpectra/src : ok
-n directory NEB/src : ok
-n directory TDDFPT/src : ok
-n directory GWW/pw4gww : ok
-n directory GWW/gww : ok
-n directory GWW/head : ok
-n directory GWW/bse : ok
-n directory GWW/simple : ok
-n directory GWW/simple_bse : ok
-n directory GWW/simple_ip : ok
-n directory QEHeat/src : ok
-n directory KCW/src : ok
-n directory KCW/PP : ok

all dependencies updated successfully
checking build system type... aarch64-apple-darwin24.3.0
checking ARCH... mac686
checking setting AR... ... ar
checking setting ARFLAGS... ... ruv
checking for gfortran... gfortran
checking whether the Fortran compiler works... yes
checking for Fortran compiler default output file name... a.out
checking for suffix of executables... 
checking whether we are cross compiling... no
checking for suffix of object files... o
checking whether we are using the GNU Fortran compiler... yes
checking whether gfortran accepts -g... yes
checking for Fortran flag to compile .f90 files... none
checking for mpif90... mpif90
checking whether we are using the GNU Fortran compiler... yes
checking whether mpif90 accepts -g... yes
checking version of mpif90... gfortran 15.2
...
a - read_fhi.o
a - read_ncpp.o
a - read_ps.o
a - read_psml.o
a - read_upf_new.o
a - read_upf_v1.o
a - read_uspp.o
a - spinor.o
a - sph_ind.o
a - sph_bes.o
a - splinelib.o
a - simpsn.o
a - upf_auxtools.o
a - upf_const.o
a - upf_error.o
a - upf_invmat.o
a - upf_io.o
a - upf_ions.o
a - upf_kinds.o
a - upf_params.o
a - upf_parallel_include.o
a - upf_spinorb.o
a - upf_to_internal.o
a - upf_utils.o
a - uspp.o
a - uspp_param.o
a - write_upf_new.o
a - xmltools.o
a - ylmr2.o
a - ylmr2_gpu.o
a - dom.o
a - wxml.o
ranlib -c libupf.a    
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
gfortran -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include   -c mbd.F90
mpif90 -g -o virtual_v2.x virtual_v2.o atom.o atomic_number.o dylmr2.o gth.o paw_variables.o pseudo_types.o radial_grids.o read_cpmd.o read_fhi.o read_ncpp.o read_ps.o read_psml.o read_upf_new.o read_upf_v1.o read_uspp.o spinor.o sph_ind.o sph_bes.o splinelib.o simpsn.o upf_auxtools.o upf_const.o upf_error.o upf_invmat.o upf_io.o upf_ions.o upf_kinds.o upf_params.o upf_parallel_include.o upf_spinorb.o upf_to_internal.o upf_utils.o uspp.o uspp_param.o write_upf_new.o xmltools.o ylmr2.o ylmr2_gpu.o dom.o wxml.o -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
ar -r libmbd.a mbd.o mbd_constants.o mbd_coulomb.o mbd_damping.o mbd_defaults.o mbd_density.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_hamiltonian.o mbd_lapack.o mbd_linalg.o mbd_matrix.o mbd_methods.o mbd_rpa.o mbd_scs.o mbd_ts.o mbd_utils.o mbd_version.o mbd_vdw_param.o
ar: creating archive libmbd.a
```

## Errors
Attempt A (tail):
```
[ 91%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/mix_rho.f90.o
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/weights.f90.o
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/gcscf_input.f90.o
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/vhpsi_gpu.f90:159:22:

  153 |                       vns_d(1,1,na), ldimax, &
      |                      2
......
  159 |                       wfcU(1,offsetU(na)+1), 2*ldap, rtemp_d, ldimaxt, &
      |                      1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/vhpsi_gpu.f90:160:30:

  154 |                       proj_r(offsetU(na)+1,1), nwfcU, 0.0_dp, rtemp_d, ldimaxt )
      |                                                              2
......
  160 |                       1.0_dp, hpsi, 2*ldap )
      |                              1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/vhpsi.f90:382:28:

  379 |                             rvaux,ldimx, projauxr,ldimx, 0.0_dp, rtemp, ldimx)
      |                            2
......
  382 |                             wfcUaux, 2*np, rtemp, ldimx, &
      |                            1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/vhpsi.f90:383:36:

  379 |                             rvaux,ldimx, projauxr,ldimx, 0.0_dp, rtemp, ldimx)
      |                                                                 2
......
  383 |                             1.0_dp, hpsi, 2*ldap)
      |                                    1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/v_of_rho.f90.o
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/sum_band.f90.o
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/scissor.f90.o
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/force_hub.f90:2520:20:

 2413 |       CALL MYDGEMM( 'T','N',ldim, nbnd, 2*npw, 2.0_DP, dwfc, 2*npwx, spsi, &
      |                                                       2
......
 2520 |                     wfatdbeta, nwfcU, betapsi(1,nb_s), nh(nt), 1.0_dp, &
      |                    1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/force_hub.f90:2520:38:

 2413 |       CALL MYDGEMM( 'T','N',ldim, nbnd, 2*npw, 2.0_DP, dwfc, 2*npwx, spsi, &
      |                                                                     2
......
 2520 |                     wfatdbeta, nwfcU, betapsi(1,nb_s), nh(nt), 1.0_dp, &
      |                                      1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8)).
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/c_bands.f90.o
[ 94%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/run_driver.f90.o
[ 97%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/environ_pw_module.f90.o
[ 97%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/potinit.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/setlocal.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/forces.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/newd_gpu.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/newd.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/hinit0.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/addusdens.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/us_exx.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/print_clock_pw.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/s_psi.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/addusforce.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/s_1psi.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/addusstress.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/data_structure.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/hs_1psi.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/s_1psi_gpu.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/s_psi_acc.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/hs_1psi_gpu.f90.o
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/s_psi_acc.f90:296:55:

  281 |                                   qq_at(1,1,na), nhm, becpr(ofsbeta(na)+1,1),&
      |                                  2                     
......
  296 |           CALL MYDGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb, &
      |                                                       1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/s_psi_acc.f90:297:39:

  282 |                                   nkb, 0.0_dp, ps(ofsbeta(na)+1,1), nkb )
      |                                               2
......
  297 |                2 * lda, ps, nkb, 1.D0, spsi, 2*lda )
      |                                       1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/addusstress.f90:195:31:

  160 |                                qgm, 2*ngm_l, tbecsum, nij, 0.0_dp, aux2, 2*ngm_l )
      |                                             2
......
  195 |                                aux2, 2*ngm_l, 0.0_dp, fac, 3 )
      |                               1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/addusstress.f90:195:54:

  160 |                                qgm, 2*ngm_l, tbecsum, nij, 0.0_dp, aux2, 2*ngm_l )
      |                                                                   2
......
  195 |                                aux2, 2*ngm_l, 0.0_dp, fac, 3 )
      |                                                      1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (REAL(8)/COMPLEX(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/s_psi.f90:228:53:

  217 |                                   qq_at(1,1,na), nhm, becp%r(ofsbeta(na)+1,1),&
      |                                  2                   
......
  228 |           CALL DGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb, &
      |                                                     1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/s_psi.f90:229:39:

  218 |                                   nkb, 0.0_dp, ps(ofsbeta(na)+1,1), nkb )
      |                                               2
......
  229 |                2 * lda, ps, nkb, 1.D0, spsi, 2*lda )
      |                                       1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/move_ions.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/hinit1.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/init_run.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/exx.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/compute_becsum.f90.o
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/compute_becsum.f90:290:30:

  255 |                       1.0_dp, auxg1, nbnd_loc,    &
      |                              2
......
  290 |                       1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
      |                              1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/compute_becsum.f90:290:55:

  256 |                       auxg2, nbnd_loc, 0.0_dp, aux_gk, nhnt )
      |                      2                                 
......
  290 |                       1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
      |                                                       1
Warning: Type mismatch between actual argument at (1) and actual argument at (2) (COMPLEX(8)/REAL(8)).
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/pw_restart_new.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/clean_pw.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/h_psi.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/h_psi_gpu.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/loc_scdm.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/non_scf.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/run_pwscf.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/setup.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/summary.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/memory_report.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/pw2casino_write.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/stress.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/scale_h.f90.o
/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1/PW/src/stress.f90:108:75:

  108 |                          nspin, dfftp, g, alat, omega, sigmaxc, rho%kin_r )
      |                                                                           1
Warning: More actual than formal arguments in procedure call at (1)
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/input.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/loc_scdm_k.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/electrons.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/punch.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/read_file_new.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw.dir/src/wfcinit.f90.o
[100%] Linking Fortran static library ../lib/libqe_pw.a
[100%] Built target qe_pw
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_ev_exe.dir/tools/ev.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_cell2ibrav_exe.dir/tools/cell2ibrav.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_scan_ibrav_exe.dir/tools/scan_ibrav.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_pwi2xsf_exe.dir/tools/pwi2xsf.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_kpoints_exe.dir/tools/kpoints.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_exe.dir/src/pwscf.f90.o
[100%] Building Fortran object PW/CMakeFiles/qe_pw_tools_ibrav2cell_exe.dir/tools/ibrav2cell.f90.o
[100%] Linking Fortran executable ../bin/pw.x
[100%] Linking Fortran executable ../bin/cell2ibrav.x
[100%] Linking Fortran executable ../bin/ibrav2cell.x
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_lax.a', '../lib/libqe_modules.a', '../lib/libqe_upflib.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
[100%] Linking Fortran executable ../bin/scan_ibrav.x
[100%] Linking Fortran executable ../bin/pwi2xsf.x
[100%] Linking Fortran executable ../bin/kpoints.x
[100%] Built target qe_pw_tools_cell2ibrav_exe
[100%] Built target qe_pw_tools_ibrav2cell_exe
[100%] Built target qe_pw_exe
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
[100%] Built target qe_pw_tools_pwi2xsf_exe
[100%] Built target qe_pw_tools_scan_ibrav_exe
[100%] Built target qe_pw_tools_kpoints_exe
[100%] Linking Fortran executable ../bin/ev.x
ld: warning: ignoring duplicate libraries: '../lib/libqe_modules.a'
[100%] Built target qe_pw_tools_ev_exe
[100%] basic code for scf, structure optimization, MD
[100%] Built target pw
```
Attempt B (tail):
```
   INTEGER LAX_DESC_SIZE
   PARAMETER (LAX_DESC_SIZE=32)
   INTEGER LAX_DESC_IR
   PARAMETER (LAX_DESC_IR=1)
   INTEGER LAX_DESC_NR
   PARAMETER (LAX_DESC_NR=2)
   INTEGER LAX_DESC_IC
   PARAMETER (LAX_DESC_IC=3)
   INTEGER LAX_DESC_NC
   PARAMETER (LAX_DESC_NC=4)
   INTEGER LAX_DESC_NRCX
   PARAMETER (LAX_DESC_NRCX=5)
   INTEGER LAX_DESC_ACTIVE_NODE
   PARAMETER (LAX_DESC_ACTIVE_NODE=6)
   INTEGER LAX_DESC_N
   PARAMETER (LAX_DESC_N=7)
   INTEGER LAX_DESC_NX
   PARAMETER (LAX_DESC_NX=8)
   INTEGER LAX_DESC_NPR
   PARAMETER (LAX_DESC_NPR=9)
   INTEGER LAX_DESC_NPC
   PARAMETER (LAX_DESC_NPC=10)
   INTEGER LAX_DESC_MYR
   PARAMETER (LAX_DESC_MYR=11)
   INTEGER LAX_DESC_MYC
   PARAMETER (LAX_DESC_MYC=12)
   INTEGER LAX_DESC_COMM
   PARAMETER (LAX_DESC_COMM=13)
   INTEGER LAX_DESC_CNTX
   PARAMETER (LAX_DESC_CNTX=14)
   INTEGER LAX_DESC_MYPE
   PARAMETER (LAX_DESC_MYPE=15)
   INTEGER LAX_DESC_NRL
   PARAMETER (LAX_DESC_NRL=16)
   INTEGER LAX_DESC_NRLX
   PARAMETER (LAX_DESC_NRLX=17)
   !
   INTEGER LAX_STATUS_SIZE
   PARAMETER (LAX_STATUS_SIZE=32)
   INTEGER LAX_STATUS_NPROC        ! nproc_ortho
   PARAMETER (LAX_STATUS_NPROC=1)
   INTEGER LAX_STATUS_LEG          ! leg_ortho
   PARAMETER (LAX_STATUS_LEG=2)
   INTEGER LAX_STATUS_NP1          ! np_ortho(1)
   PARAMETER (LAX_STATUS_NP1=3)
   INTEGER LAX_STATUS_NP2          ! np_ortho(2)
   PARAMETER (LAX_STATUS_NP2=4)
   INTEGER LAX_STATUS_ME1          ! me_ortho(1)
   PARAMETER (LAX_STATUS_ME1=5)
   INTEGER LAX_STATUS_ME2          ! me_ortho(2)
   PARAMETER (LAX_STATUS_ME2=6)
   INTEGER LAX_STATUS_COMM         ! ortho_comm
   PARAMETER (LAX_STATUS_COMM=7)
   INTEGER LAX_STATUS_ROWCOMM      ! ortho_row_comm
   PARAMETER (LAX_STATUS_ROWCOMM=8)
   INTEGER LAX_STATUS_COLCOMM      ! ortho_col_comm
   PARAMETER (LAX_STATUS_COLCOMM=9)
   INTEGER LAX_STATUS_COMMID       ! ortho_comm_id
   PARAMETER (LAX_STATUS_COMMID=10)
   INTEGER LAX_STATUS_PARENTCOMM   ! ortho_parent_comm
   PARAMETER (LAX_STATUS_PARENTCOMM=11)
   INTEGER LAX_STATUS_ORTHOCNTX    ! ortho_cntx
   PARAMETER (LAX_STATUS_ORTHOCNTX=12)
   INTEGER LAX_STATUS_DISTDIAG     ! do_distr_diag_inside_bgrp
   PARAMETER (LAX_STATUS_DISTDIAG=13)

1 warning generated.
cc: error: no input files
cc: error: no input files
make[1]: *** [laxlib_hi.fh] Error 1
make[1]: *** Waiting for unfinished jobs....
make[1]: *** [laxlib_low.fh] Error 1
cc: error: no input files
cc: error: no input files
make[1]: *** [laxlib_kinds.fh] Error 1
make[1]: *** [laxlib.fh] Error 1
cc: error: no input files
make[1]: *** [laxlib_mid.fh] Error 1
cc: error: no input files
make[1]: *** [laxlib_param.fh] Error 1
make: *** [libla] Error 1
make: *** Waiting for unfinished jobs....
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_const.f90 -o upf_const.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c pseudo_types.f90 -o pseudo_types.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_invmat.f90 -o upf_invmat.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c paw_variables.f90 -o paw_variables.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c dylmr2.f90 -o dylmr2.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c splinelib.f90 -o splinelib.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c xmltools.f90 -o xmltools.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c spinor.f90 -o spinor.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c sph_ind.f90 -o sph_ind.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c simpsn.f90 -o simpsn.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c atomic_number.f90 -o atomic_number.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_error.f90 -o upf_error.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c radial_grids.f90 -o radial_grids.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c sph_bes.f90 -o sph_bes.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c ylmr2.f90 -o ylmr2.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c atom.f90 -o atom.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c uspp_param.f90 -o uspp_param.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c gth.f90 -o gth.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_cpmd.f90 -o read_cpmd.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_fhi.f90 -o read_fhi.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_ncpp.f90 -o read_ncpp.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_upf_v1.f90 -o read_upf_v1.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_uspp.f90 -o read_uspp.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_auxtools.f90 -o upf_auxtools.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_to_internal.f90 -o upf_to_internal.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c casino_pp.f90 -o casino_pp.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c atwfc_mod.f90 -o atwfc_mod.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c qrad_mod.f90 -o qrad_mod.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c rhoat_mod.f90 -o rhoat_mod.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c rhoc_mod.f90 -o rhoc_mod.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_spinorb.f90 -o upf_spinorb.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c init_us_0.f90 -o init_us_0.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c init_us_b0.f90 -o init_us_b0.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upf_ions.f90 -o upf_ions.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_psml.f90 -o read_psml.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_upf_new.f90 -o read_upf_new.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c write_upf_new.f90 -o write_upf_new.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c wxml.f90 -o wxml.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c uspp.f90 -o uspp.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c qvan2.f90 -o qvan2.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c dqvan2.f90 -o dqvan2.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c init_us_1.f90 -o init_us_1.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c beta_mod.f90 -o beta_mod.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c vloc_mod.f90 -o vloc_mod.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c gen_us_dj.f90 -o gen_us_dj.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c gen_us_dy.f90 -o gen_us_dy.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c init_us_2_acc.f90 -o init_us_2_acc.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c virtual_v2.f90 -o virtual_v2.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c upfconv.f90 -o upfconv.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c casino2upf.f90 -o casino2upf.o
gfortran -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include   -c mbd_methods.F90
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I../UtilXlib -I. -c read_ps.f90 -o read_ps.o
ar ruv libupf.a atwfc_mod.o beta_mod.o qrad_mod.o rhoat_mod.o rhoc_mod.o vloc_mod.o qvan2.o dqvan2.o gen_us_dj.o gen_us_dy.o init_us_0.o init_us_b0.o init_us_1.o init_us_2_acc.o atom.o atomic_number.o dylmr2.o gth.o paw_variables.o pseudo_types.o radial_grids.o read_cpmd.o read_fhi.o read_ncpp.o read_ps.o read_psml.o read_upf_new.o read_upf_v1.o read_uspp.o spinor.o sph_ind.o sph_bes.o splinelib.o simpsn.o upf_auxtools.o upf_const.o upf_error.o upf_invmat.o upf_io.o upf_ions.o upf_kinds.o upf_params.o upf_parallel_include.o upf_spinorb.o upf_to_internal.o upf_utils.o uspp.o uspp_param.o write_upf_new.o xmltools.o ylmr2.o ylmr2_gpu.o dom.o wxml.o       
mpif90 -g -o upfconv.x upfconv.o casino_pp.o atom.o atomic_number.o dylmr2.o gth.o paw_variables.o pseudo_types.o radial_grids.o read_cpmd.o read_fhi.o read_ncpp.o read_ps.o read_psml.o read_upf_new.o read_upf_v1.o read_uspp.o spinor.o sph_ind.o sph_bes.o splinelib.o simpsn.o upf_auxtools.o upf_const.o upf_error.o upf_invmat.o upf_io.o upf_ions.o upf_kinds.o upf_params.o upf_parallel_include.o upf_spinorb.o upf_to_internal.o upf_utils.o uspp.o uspp_param.o write_upf_new.o xmltools.o ylmr2.o ylmr2_gpu.o dom.o wxml.o -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate
mpif90 -g -o casino2upf.x casino2upf.o casino_pp.o atom.o atomic_number.o dylmr2.o gth.o paw_variables.o pseudo_types.o radial_grids.o read_cpmd.o read_fhi.o read_ncpp.o read_ps.o read_psml.o read_upf_new.o read_upf_v1.o read_uspp.o spinor.o sph_ind.o sph_bes.o splinelib.o simpsn.o upf_auxtools.o upf_const.o upf_error.o upf_invmat.o upf_io.o upf_ions.o upf_kinds.o upf_params.o upf_parallel_include.o upf_spinorb.o upf_to_internal.o upf_utils.o uspp.o uspp_param.o write_upf_new.o xmltools.o ylmr2.o ylmr2_gpu.o dom.o wxml.o -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate
ar: creating archive libupf.a
a - atwfc_mod.o
a - beta_mod.o
a - qrad_mod.o
a - rhoat_mod.o
a - rhoc_mod.o
a - vloc_mod.o
a - qvan2.o
a - dqvan2.o
a - gen_us_dj.o
a - gen_us_dy.o
a - init_us_0.o
a - init_us_b0.o
a - init_us_1.o
a - init_us_2_acc.o
a - atom.o
a - atomic_number.o
a - dylmr2.o
a - gth.o
a - paw_variables.o
a - pseudo_types.o
a - radial_grids.o
a - read_cpmd.o
a - read_fhi.o
a - read_ncpp.o
a - read_ps.o
a - read_psml.o
a - read_upf_new.o
a - read_upf_v1.o
a - read_uspp.o
a - spinor.o
a - sph_ind.o
a - sph_bes.o
a - splinelib.o
a - simpsn.o
a - upf_auxtools.o
a - upf_const.o
a - upf_error.o
a - upf_invmat.o
a - upf_io.o
a - upf_ions.o
a - upf_kinds.o
a - upf_params.o
a - upf_parallel_include.o
a - upf_spinorb.o
a - upf_to_internal.o
a - upf_utils.o
a - uspp.o
a - uspp_param.o
a - write_upf_new.o
a - xmltools.o
a - ylmr2.o
a - ylmr2_gpu.o
a - dom.o
a - wxml.o
ranlib -c libupf.a    
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
gfortran -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include   -c mbd.F90
mpif90 -g -o virtual_v2.x virtual_v2.o atom.o atomic_number.o dylmr2.o gth.o paw_variables.o pseudo_types.o radial_grids.o read_cpmd.o read_fhi.o read_ncpp.o read_ps.o read_psml.o read_upf_new.o read_upf_v1.o read_uspp.o spinor.o sph_ind.o sph_bes.o splinelib.o simpsn.o upf_auxtools.o upf_const.o upf_error.o upf_invmat.o upf_io.o upf_ions.o upf_kinds.o upf_params.o upf_parallel_include.o upf_spinorb.o upf_to_internal.o upf_utils.o uspp.o uspp_param.o write_upf_new.o xmltools.o ylmr2.o ylmr2_gpu.o dom.o wxml.o -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
ar -r libmbd.a mbd.o mbd_constants.o mbd_coulomb.o mbd_damping.o mbd_defaults.o mbd_density.o mbd_dipole.o mbd_formulas.o mbd_geom.o mbd_gradients.o mbd_hamiltonian.o mbd_lapack.o mbd_linalg.o mbd_matrix.o mbd_methods.o mbd_rpa.o mbd_scs.o mbd_ts.o mbd_utils.o mbd_version.o mbd_vdw_param.o
ar: creating archive libmbd.a
```

## Version captures
```
=== BREW VERSIONS ===
cmake 4.1.2
gcc 15.2.0
open-mpi 5.0.8
openblas 0.3.30
veclibfort 0.4.3
=== BREW CONFIG ===
HOMEBREW_VERSION: 4.6.18
ORIGIN: https://github.com/Homebrew/brew
HEAD: 00001bedf8eac442bdabd327148be0b97998b27b
Last commit: 6 days ago
Branch: stable
Core tap JSON: 25 Oct 17:10 UTC
Core cask tap JSON: 25 Oct 17:10 UTC
HOMEBREW_PREFIX: /opt/homebrew
HOMEBREW_CASK_OPTS: []
HOMEBREW_FORBID_PACKAGES_FROM_PATHS: set
HOMEBREW_MAKE_JOBS: 10
Homebrew Ruby: 3.4.5 => /opt/homebrew/Library/Homebrew/vendor/portable-ruby/3.4.5/bin/ruby
CPU: deca-core 64-bit arm_donan
Clang: 17.0.0 build 1700
Git: 2.50.1 => /opt/homebrew/bin/git
Curl: 8.7.1 => /usr/bin/curl
macOS: 15.3-arm64
CLT: 16.4.0.0.1.1747106510
Xcode: N/A
Rosetta 2: false
=== TOOL VERSIONS ===
GNU Fortran (Homebrew GCC 15.2.0) 15.2.0
Copyright (C) 2025 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

mpirun (Open MPI) 5.0.8

Report bugs to https://www.open-mpi.org/community/help/

cmake version 4.1.2

CMake suite maintained and supported by Kitware (kitware.com/cmake).
```

## Assessment
Build succeeded on Apple Silicon. Use this machine for light tests; keep heavy runs on the server.
## Attempt B2 (configure+Accelerate) result
ranlib -c libpw.a
if test -d PW ; then \
	( cd PW ; /Library/Developer/CommandLineTools/usr/bin/make TLDEPS= all || exit 1) ; fi
( cd src ; /Library/Developer/CommandLineTools/usr/bin/make all || exit 1 )
if test -n "" ; then \
	( cd ../.. ; /Library/Developer/CommandLineTools/usr/bin/make  || exit 1 ) ; fi
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../../dft-d3/ -c pwscf.f90 -o pwscf.o
mpif90 -g -o pw.x \
	   pwscf.o  libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
( cd ../../bin; ln -fs ../PW/src/pw.x . ; \
               ln -fs ../PW/src/pw.x dist.x ; ln -fs ../PW/src/pw.x manypw.x ; )
( cd tools ; /Library/Developer/CommandLineTools/usr/bin/make all || exit 1 )
if test -n "" ; then \
	( cd ../.. ; /Library/Developer/CommandLineTools/usr/bin/make  || exit 1 ) ; fi
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c ev.f90 -o ev.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c kpoints.f90 -o kpoints.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c pwi2xsf.f90 -o pwi2xsf.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c ibrav2cell.f90 -o ibrav2cell.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c cell2ibrav.f90 -o cell2ibrav.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c scan_ibrav.f90 -o scan_ibrav.o
mpif90 -O3 -g -fallow-argument-mismatch -cpp -D__FFTW -D__MPI -D__MPI_MODULE  -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -I. -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//include -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD -I/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//KS_Solvers  -I../src -I../../dft-d3/ -c rism1d.f90 -o rism1d.o
mpif90 -g -o rism1d.x \
		rism1d.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a   -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
mpif90 -g -o cell2ibrav.x \
		cell2ibrav.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
mpif90 -g -o ibrav2cell.x \
		ibrav2cell.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
( cd ../../bin ; ln -fs ../PW/tools/rism1d.x . )
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
( cd ../../bin ; ln -fs ../PW/tools/cell2ibrav.x . )
( cd ../../bin ; ln -fs ../PW/tools/ibrav2cell.x . )
mpif90 -g -o scan_ibrav.x \
		scan_ibrav.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
mpif90 -g -o pwi2xsf.x \
		pwi2xsf.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
mpif90 -g -o kpoints.x \
		kpoints.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
( cd ../../bin ; ln -fs ../PW/tools/scan_ibrav.x . )
( cd ../../bin ; ln -fs ../PW/tools/kpoints.x . )
( cd ../../bin ; ln -fs ../PW/tools/pwi2xsf.x . )
mpif90 -g -o ev.x \
		ev.o ../src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//Modules/libqemod.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//upflib/libupf.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//XClib/xc_lib.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//FFTXlib/src/libqefft.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//LAXlib/libqela.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//UtilXlib/libutil.a /Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//MBD/libmbd.a  -L/Users/pollob/qe_macm4_build/artifacts/q-e-qe-7.4.1//external/devxlib/src -ldevXlib  -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate   -L/opt/homebrew/opt/veclibfort/lib -lvecLibFort -framework Accelerate     
ld: warning: ignoring duplicate libraries: '-lvecLibFort'
( cd ../../bin ; ln -fs ../PW/tools/ev.x . )

## Linkage
### otool -L artifacts/q-e-qe-7.4.1/bin/pw.x
artifacts/q-e-qe-7.4.1/bin/pw.x:
	/opt/homebrew/opt/veclibfort/lib/libvecLibFort.dylib (compatibility version 0.0.0, current version 0.0.0)
	/opt/homebrew/opt/open-mpi/lib/libmpi_usempif08.40.dylib (compatibility version 81.0.0, current version 81.3.0)
	/opt/homebrew/opt/open-mpi/lib/libmpi_usempi_ignore_tkr.40.dylib (compatibility version 81.0.0, current version 81.1.0)
	/opt/homebrew/opt/open-mpi/lib/libmpi_mpifh.40.dylib (compatibility version 81.0.0, current version 81.1.0)
	/opt/homebrew/opt/open-mpi/lib/libmpi.40.dylib (compatibility version 81.0.0, current version 81.7.0)
	/opt/homebrew/opt/gcc/lib/gcc/current/libgfortran.5.dylib (compatibility version 6.0.0, current version 6.0.0)
	/System/Library/Frameworks/Accelerate.framework/Versions/A/Accelerate (compatibility version 1.0.0, current version 4.0.0)
	/opt/homebrew/opt/gcc/lib/gcc/current/libquadmath.0.dylib (compatibility version 1.0.0, current version 1.0.0)
	/usr/lib/libSystem.B.dylib (compatibility version 1.0.0, current version 1351.0.0)

## Si SCF quick run
PP missing: place Si.pbe-n-rrkjus_psl.1.0.0.UPF under ./pp and re-run. See docs/PP.md
