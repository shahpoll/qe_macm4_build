#!/usr/bin/env python3
"""Extract band gap information from silicon band structure output."""
from __future__ import annotations

import math
import pathlib
from typing import List, Tuple

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
BAND_FILE = ROOT / "analysis" / "Si" / "data" / "silicon.bands.dat.gnu"
DOS_FILE = ROOT / "analysis" / "Si" / "data" / "silicon.dos"
OUT_FILE = ROOT / "analysis" / "Si" / "data" / "si_band_summary.txt"


def read_fermi() -> float:
    with DOS_FILE.open("r", encoding="utf-8") as dos:
        line = dos.readline()
        if "EFermi" not in line:
            raise RuntimeError("Could not locate EFermi in silicon.dos header")
        return float(line.split("EFermi =")[1].split()[0])


def read_bands() -> Tuple[np.ndarray, np.ndarray]:
    bands: List[List[Tuple[float, float]]] = []
    current: List[Tuple[float, float]] = []
    with BAND_FILE.open("r", encoding="utf-8") as bf:
        for raw in bf:
            line = raw.strip()
            if not line:
                if current:
                    bands.append(current)
                    current = []
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            k = float(parts[0])
            e = float(parts[1])
            current.append((k, e))
    if current:
        bands.append(current)

    if not bands:
        raise RuntimeError("No band data parsed from silicon.bands.dat.gnu")

    nk = len(bands[0])
    for band in bands:
        if len(band) != nk:
            raise RuntimeError("Inconsistent k-point sampling across bands")
    k_path = np.array([pt[0] for pt in bands[0]])
    energies = np.array([[pt[1] for pt in band] for band in bands])
    return k_path, energies


def summarize() -> str:
    k_path, energies = read_bands()
    fermi = read_fermi()

    with np.errstate(invalid="ignore"):
        valence = np.where(energies <= fermi, energies, -math.inf)
        conduction = np.where(energies >= fermi, energies, math.inf)

    valence_max_by_k = valence.max(axis=0)
    conduction_min_by_k = conduction.min(axis=0)

    global_vmax = valence_max_by_k.max()
    global_cmin = conduction_min_by_k.min()

    indirect_gap = global_cmin - global_vmax
    k_vmax_idx = int(valence_max_by_k.argmax())
    k_cmin_idx = int(conduction_min_by_k.argmin())

    direct_deltas = conduction_min_by_k - valence_max_by_k
    direct_gap = direct_deltas.min()
    k_direct_idx = int(direct_deltas.argmin())

    lines = [
        "# Silicon band summary (derived from silicon.bands.dat.gnu)",
        f"Fermi level (from silicon.dos)  : {fermi:8.4f} eV",
        f"Valence band maximum            : {global_vmax:8.4f} eV at k-index {k_vmax_idx} (s = {k_path[k_vmax_idx]:.4f})",
        f"Conduction band minimum         : {global_cmin:8.4f} eV at k-index {k_cmin_idx} (s = {k_path[k_cmin_idx]:.4f})",
        f"Indirect gap (Cmin - Vmax)      : {indirect_gap:8.4f} eV",
        f"Direct gap (min over k)         : {direct_gap:8.4f} eV at k-index {k_direct_idx} (s = {k_path[k_direct_idx]:.4f})",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    summary = summarize()
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(summary, encoding="utf-8")
    print(summary)
if __name__ == "__main__":
    main()
