#!/usr/bin/env python3
"""
Plot total DOS for silicon using Quantum ESPRESSO outputs.

Usage examples
--------------
Default (manual workflow)::

    python3 scripts/plot_si_dos.py

PWTK workflow::

    python3 scripts/plot_si_dos.py --pwtk

Custom base directory::

    python3 scripts/plot_si_dos.py --base /path/to/run
"""

from __future__ import annotations

import argparse
import pathlib
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_BASE = ROOT / "analysis" / "Si"
PT_BASE = ROOT / "analysis_pwtk" / "Si"


def compute_paths(base: pathlib.Path | None) -> Tuple[pathlib.Path, pathlib.Path]:
    base_dir = DEFAULT_BASE if base is None else base
    return base_dir / "data" / "silicon.dos", base_dir / "plots" / "si_total_dos.png"


def load_dos(path: pathlib.Path) -> Tuple[np.ndarray, np.ndarray, float | None]:
    energies: list[float] = []
    dos: list[float] = []
    fermi: float | None = None

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                if "EFermi" in line:
                    try:
                        fermi = float(line.split("EFermi =")[1].split()[0])
                    except (IndexError, ValueError):
                        pass
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            energies.append(float(cols[0]))
            dos.append(float(cols[1]))

    if not energies:
        raise SystemExit(f"No DOS data parsed from {path}")
    return np.array(energies), np.array(dos), fermi


def main(base: pathlib.Path | None = None) -> None:
    dos_file, out_png = compute_paths(base)
    energies, dos, fermi = load_dos(dos_file)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)
    ax.plot(energies, dos, color="tab:green")
    if fermi is not None:
        ax.axvline(fermi, color="tab:red", linestyle="--", linewidth=0.8, label="Fermi level")
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.set_title("Silicon total DOS — QE 7.4.1")
    ax.set_xlim(energies[0], energies[-1])
    ax.set_ylim(bottom=0)
    if fermi is not None:
        ax.legend(frameon=False)
    ax.grid(alpha=0.2, which="both", axis="both")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot silicon DOS")
    parser.add_argument("--base", type=pathlib.Path, help="Base directory (default analysis/Si)")
    parser.add_argument("--pwtk", action="store_true", help="Use analysis_pwtk/Si as base directory")
    args = parser.parse_args()

    base_dir = args.base
    if args.pwtk:
        base_dir = PT_BASE

    main(base_dir)
