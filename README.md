# Recombination study: figures, tables and scripts

This repository accompanies the manuscript on natural variation in wheat crossover
number in the Watkins landrace collection. It contains the figures, the tables and
the analysis scripts.

Raw genotype and phenotype data are not included; they are public and are identified
in the manuscript's data availability statement. The manuscript text is not included.

## Layout

| Directory | Content |
|---|---|
| `Scripts/` | R scripts for the analyses and figures, numbered in the order they are run |
| `Tables/manuscript/` | Tables 1-2 and S1-S21, numbered as in the paper |
| `Tables/analysis/` | Intermediate tables with the original working names (`Table22c_...`) |
| `Figures/main/` | Main figures, TIFF, 600 dpi |
| `Figures/supplementary/` | Supplementary figures, PNG, numbered as in the paper |
| `Figures/preview/` | Low-resolution PNG previews of the main figures |
| `Figures/legacy/` | Superseded and exploratory panels, not cited in the manuscript |

See `Figures/INDEX.md` and `Tables/INDEX.md` for the full file lists.

## Reproducing the analysis

Paths at the top of each script (`ROOT`) point to the working directory used in
the study. Change `ROOT` to the location of your local copy before running.
Scripts are numbered; those with a shared prefix (for example `250-`, `251-`,
`252-`) generate the main figures and expect the tables produced by the earlier
steps.

## Notes

- R 4.5.1 was used.
- `README_simulation.md` is the README of an earlier, unrelated simulation that
  was kept in this repository when it was reused for this study.
