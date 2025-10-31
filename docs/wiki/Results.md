# Results Gallery

Key artefacts from the silicon (Si) validation run.

## Band structure

![Si band structure](../Si/plots/si_band_structure.png)

- Path: Γ → X → W → K → Γ → L (40 points per segment).
- Fermi level (red dashed line via `plot_si_bands.py` script).

## Total DOS

![Si DOS](../Si/plots/si_total_dos.png)

- Derived from the 12×12×12 NSCF run + `dos.x`.
- Energy range: –15 to 20 eV with ΔE = 0.02 eV.

## Projected DOS (s/p channels)

![Si PDOS](../Si/plots/si_pdos.png)

- Aggregated s and p contributions for the two Si atoms.
- Generated with `projwfc.x` and plotted via `plot_si_pdos.py`.

## Band-gap summary (`si_band_summary.txt`)

```
# Silicon band summary (derived from silicon.bands.dat.gnu)
Fermi level (from silicon.dos)  :   6.2200 eV
Valence band maximum            :   6.2198 eV at k-index 0 (s = 0.0000)
Conduction band minimum         :   6.7897 eV at k-index 34 (s = 0.8500)
Indirect gap (Cmin - Vmax)      :   0.5699 eV
Direct gap (min over k)         :   2.5598 eV at k-index 0 (s = 0.0000)
```

These metrics are produced by `scripts/analyze_si_bandgap.py`.
