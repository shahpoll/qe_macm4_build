#!/usr/bin/env python3
"""
Plot total DOS for silicon using Quantum ESPRESSO outputs.

Input: analysis/Si/data/silicon.dos
Output: analysis/Si/plots/si_total_dos.png
"""

from __future__ import annotations

import pathlib
import numpy as np
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOS_FILE = ROOT / "analysis" / "Si" / "data" / "silicon.dos"
OUT_PNG = ROOT / "analysis" / "Si" / "plots" / "si_total_dos.png"


def load_dos(path: pathlib.Path):
    energies = []
    dos = []
    fermi = None

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                if "EFermi" in line:
                    parts = line.split("EFermi =")
                    if len(parts) > 1:
                        fermi = float(parts[1].split()[0])
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            energies.append(float(cols[0]))
            dos.append(float(cols[1]))

    if not energies:
        raise SystemExit("No DOS data parsed.")
    return np.array(energies), np.array(dos), fermi


def main():
    energies, dos, fermi = load_dos(DOS_FILE)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.plot(energies, dos, color="tab:green")
    if fermi is not None:
        ax.axvline(fermi, color="tab:red", linestyle="--", linewidth=0.8, label="Fermi level")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.set_title("Silicon total DOS — QE 7.4.1, Apple Silicon")
    ax.set_xlim(energies[0], energies[-1])
    ax.set_ylim(bottom=0)
    if fermi is not None:
        ax.legend(frameon=False)
    ax.grid(alpha=0.2, which="both", axis="both")
    fig.tight_layout()
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG)
    plt.close(fig)


if __name__ == "__main__":
    main()
