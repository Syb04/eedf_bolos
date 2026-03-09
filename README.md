# Compute Argon EEDF with BOLOS

`argon_eedf_bolos.py` reads an Argon LXCat/BOLSIG+ cross-section file and solves for the EEDF with `BOLOS`.

`argon_eedf_argon_cross_section.ipynb` uses the local `argon cross section.txt` file, sweeps `E/N` from low to high values, computes the EEDF for each case, and visualizes the EEDF, transport coefficients, and reaction rates.

`argon_eedf_boltzmann_vs_maxwell.ipynb` compares the converged Boltzmann EEDF against a Maxwellian with the same mean electron energy and plots `f_Boltzmann / f_Maxwell`.

## Requirements

- `bolos`
- `ipywidgets`
- `matplotlib`
- `pandas`
- An Argon cross-section file in LXCat/BOLSIG+ format

## Script example

```powershell
python .\argon_eedf_bolos.py ".\argon cross section.txt" --species Ar --en-td 120 --output .\argon_eedf.csv
```

## Notebook

Open `argon_eedf_argon_cross_section.ipynb` in Jupyter or VS Code and run the cells from top to bottom.

## Output

- Console summary from the script:
  - mean electron energy
  - mobility
  - diffusion coefficient
- CSV columns:
  - `energy_eV`
  - `eedf`
- Notebook CSV exports:
  - transport summary vs `E/N`
  - EEDF vs energy for each `E/N`
  - reaction rate coefficients vs `E/N`
