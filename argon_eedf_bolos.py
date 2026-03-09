#!/usr/bin/env python3
"""Compute the Argon EEDF with BOLOS from an LXCat/BOLSIG+ cross-section file."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from bolos import grid, parser, solver


def parse_args() -> argparse.Namespace:
    cli = argparse.ArgumentParser(
        description=(
            "Load an Argon LXCat/BOLSIG+ cross-section file and "
            "compute the EEDF with BOLOS."
        )
    )
    cli.add_argument(
        "cross_sections",
        type=Path,
        help="LXCat/BOLSIG+ cross-section file",
    )
    cli.add_argument(
        "--species",
        default="Ar",
        help="Argon target name inside the cross-section file (default: Ar)",
    )
    cli.add_argument(
        "--en-td",
        type=float,
        default=120.0,
        help="Reduced electric field E/N in Td (default: 120)",
    )
    cli.add_argument(
        "--gas-temperature-k",
        type=float,
        default=300.0,
        help="Gas temperature in K (default: 300)",
    )
    cli.add_argument(
        "--cells",
        type=int,
        default=300,
        help="Number of energy-grid cells (default: 300)",
    )
    cli.add_argument(
        "--max-energy-ev",
        type=float,
        default=60.0,
        help="Maximum energy of the initial grid in eV (default: 60)",
    )
    cli.add_argument(
        "--initial-te-ev",
        type=float,
        default=2.0,
        help="Electron temperature of the initial Maxwellian guess in eV (default: 2.0)",
    )
    cli.add_argument(
        "--max-iterations",
        type=int,
        default=200,
        help="Maximum number of solver iterations (default: 200)",
    )
    cli.add_argument(
        "--rtol",
        type=float,
        default=1e-5,
        help="Relative tolerance for convergence (default: 1e-5)",
    )
    cli.add_argument(
        "--output",
        type=Path,
        default=Path("argon_eedf.csv"),
        help="Output CSV file (default: argon_eedf.csv)",
    )
    return cli.parse_args()


def build_solver(
    cross_section_path: Path,
    species: str,
    en_td: float,
    gas_temperature_k: float,
    cells: int,
    max_energy_ev: float,
) -> solver.BoltzmannSolver:
    energy_grid = grid.LinearGrid(0.0, max_energy_ev, cells)
    boltzmann = solver.BoltzmannSolver(energy_grid)

    with cross_section_path.open("r", encoding="utf-8") as handle:
        processes = parser.parse(handle)

    boltzmann.load_collisions(processes)

    if species not in boltzmann.target:
        available = ", ".join(sorted(boltzmann.target.keys()))
        raise ValueError(
            f"Target '{species}' was not found. Available targets: {available}"
        )

    boltzmann.target[species].density = 1.0
    boltzmann.kT = gas_temperature_k * solver.KB / solver.ELECTRONVOLT
    boltzmann.EN = en_td * solver.TOWNSEND
    boltzmann.init()
    return boltzmann


def solve_eedf(
    boltzmann: solver.BoltzmannSolver,
    initial_te_ev: float,
    max_iterations: int,
    rtol: float,
):
    # BOLOS needs an initial guess; the converged distribution is the Boltzmann EEDF.
    f0 = boltzmann.maxwell(initial_te_ev)
    coarse = boltzmann.converge(f0, maxn=max_iterations, rtol=rtol)

    mean_energy = boltzmann.mean_energy(coarse)
    refined_grid = grid.QuadraticGrid(0.0, max(15.0 * mean_energy, 5.0), boltzmann.n)
    previous_grid = boltzmann.grid

    boltzmann.grid = refined_grid
    boltzmann.init()

    interpolated = boltzmann.grid.interpolate(coarse, previous_grid)
    refined = boltzmann.converge(interpolated, maxn=max_iterations, rtol=rtol)
    return refined


def list_inelastic_rate_coefficients(boltzmann: solver.BoltzmannSolver, eedf) -> list[dict]:
    rates = []
    for target, process in boltzmann.iter_inelastic():
        rates.append(
            {
                "target": target,
                "process": str(process),
                "rate_coefficient_m3_s": boltzmann.rate(eedf, process),
            }
        )
    return rates


def compute_eedf(
    cross_section_path: Path,
    species: str = "Ar",
    en_td: float = 120.0,
    gas_temperature_k: float = 300.0,
    cells: int = 300,
    max_energy_ev: float = 60.0,
    initial_te_ev: float = 2.0,
    max_iterations: int = 200,
    rtol: float = 1e-5,
) -> dict:
    boltzmann = build_solver(
        cross_section_path=cross_section_path,
        species=species,
        en_td=en_td,
        gas_temperature_k=gas_temperature_k,
        cells=cells,
        max_energy_ev=max_energy_ev,
    )
    eedf = solve_eedf(
        boltzmann,
        initial_te_ev=initial_te_ev,
        max_iterations=max_iterations,
        rtol=rtol,
    )
    mobility_n_si = boltzmann.mobility(eedf)
    diffusion_n_si = boltzmann.diffusion(eedf)
    return {
        "solver": boltzmann,
        "energy_eV": boltzmann.cenergy,
        "eedf": eedf,
        "mean_energy_eV": boltzmann.mean_energy(eedf),
        "mobility_n_SI": mobility_n_si,
        "diffusion_n_SI": diffusion_n_si,
        "transport": {
            "mean_energy_eV": boltzmann.mean_energy(eedf),
            "mobility_n_SI": mobility_n_si,
            "diffusion_n_SI": diffusion_n_si,
        },
        "reaction_rates": list_inelastic_rate_coefficients(boltzmann, eedf),
        "species": species,
        "en_td": en_td,
        "gas_temperature_k": gas_temperature_k,
    }


def write_csv(output_path: Path, energies, eedf) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["energy_eV", "eedf"])
        writer.writerows(zip(energies, eedf))


def main() -> None:
    args = parse_args()

    result = compute_eedf(
        cross_section_path=args.cross_sections,
        species=args.species,
        en_td=args.en_td,
        gas_temperature_k=args.gas_temperature_k,
        cells=args.cells,
        max_energy_ev=args.max_energy_ev,
        initial_te_ev=args.initial_te_ev,
        max_iterations=args.max_iterations,
        rtol=args.rtol,
    )

    write_csv(args.output, result["energy_eV"], result["eedf"])

    print(f"species           : {args.species}")
    print(f"E/N [Td]          : {args.en_td}")
    print(f"gas temperature[K]: {args.gas_temperature_k}")
    print(f"mean energy [eV]  : {result['mean_energy_eV']:.6g}")
    print(f"mobility*n [SI]   : {result['mobility_n_SI']:.6g}")
    print(f"diffusion*n [SI]  : {result['diffusion_n_SI']:.6g}")
    print(f"saved CSV         : {args.output.resolve()}")


if __name__ == "__main__":
    main()
