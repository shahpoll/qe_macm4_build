#!/usr/bin/env python3
"""
Plot projected DOS for Si (s and p channels per atom).

Inputs (from projwfc.x):
  analysis/Si/data/silicon.pdos_atm#1(Si)_wfc#1(s)
  analysis/Si/data/silicon.pdos_atm#1(Si)_wfc#2(p)
  analysis/Si/data/silicon.pdos_atm#2(Si)_wfc#1(s)
  analysis/Si/data/silicon.pdos_atm#2(Si)_wfc#2(p)
  analysis/Si/data/silicon.pdos_tot
Output:
  analysis/Si/plots/si_pdos.png
"""

from __future__ import annotations

import pathlib
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "analysis" / "Si" / "data"
OUT_PNG = ROOT / "analysis" / "Si" / "plots" / "si_pdos.png"


def read_pdos(path: pathlib.Path) -> Tuple[np.ndarray, np.ndarray]:
    energy = []
    density = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            energy.append(float(parts[0]))
            density.append(float(parts[1]))
    return np.array(energy), np.array(density)


def main() -> None:
    files = {
        "tot": DATA_DIR / "silicon.pdos_tot",
        "atm1_s": DATA_DIR / "silicon.pdos_atm#1(Si)_wfc#1(s)",
        "atm1_p": DATA_DIR / "silicon.pdos_atm#1(Si)_wfc#2(p)",
        "atm2_s": DATA_DIR / "silicon.pdos_atm#2(Si)_wfc#1(s)",
        "atm2_p": DATA_DIR / "silicon.pdos_atm#2(Si)_wfc#2(p)",
    }
    for name, path in files.items():
        if not path.exists():
            raise SystemExit(f"Missing PDOS file: {path} ({name})")

    e_tot, dos_tot = read_pdos(files["tot"])
    e1s, dos1s = read_pdos(files["atm1_s"])
    e1p, dos1p = read_pdos(files["atm1_p"])
    e2s, dos2s = read_pdos(files["atm2_s"])
    e2p, dos2p = read_pdos(files["atm2_p"])

    if not (np.allclose(e_tot, e1s) and np.allclose(e_tot, e1p)):
        raise SystemExit("Energy grids differ between PDOS files")

    s_total = dos1s + dos2s
    p_total = dos1p + dos2p

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.plot(e_tot, dos_tot, color="0.2", label="Total DOS")
    ax.plot(e_tot, s_total, color="tab:blue", label="Si s (total)")
    ax.plot(e_tot, p_total, color="tab:orange", label="Si p (total)")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.set_title("Silicon projected DOS — QE 7.4.1")
    ax.set_xlim(e_tot.min(), e_tot.max())
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    ax.grid(alpha=0.2, which="both", axis="both")
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG)
    plt.close(fig)


if __name__ == "__main__":
    main()
