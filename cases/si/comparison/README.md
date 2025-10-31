# Silicon — Manual vs PWTK Comparison

Generated artefacts:

- `data/si_comparison.txt` — maximum absolute differences between manual and PWTK runs.
- `plots/si_manual_vs_pwtk.png` — overlay of band structure, total DOS, and PDOS.

Recreate with:
```sh
cd /Users/pollob/qe_macm4_build
python3 scripts/compare_si_runs.py
```

Notes:
- Band/DOS/PDOS curves are numerically identical up to floating-point noise (`<1e-10`).
- QE runtimes are not captured in current logs; add `/usr/bin/time` wrappers if you need benchmarking.
