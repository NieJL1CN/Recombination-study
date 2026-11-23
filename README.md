# Recombination study: figures, tables and scripts

This repository accompanies the manuscript on natural variation in wheat crossover
number in the Watkins landrace collection. It contains the figures, the derived
tables and the analysis scripts.

Raw genotype and phenotype data are not included here; they are public and are
identified in the manuscript's data availability statement. The manuscript text
is not included.

## Layout

| Directory | Content |
|---|---|
| `Scripts/` | R scripts for the analyses and figures (numbered in the order they are run) |
| `Tables/` | Derived tables referenced in the manuscript and its supplementary information (Table*.csv) |
| `Figures/` | Main figures (`main/`, TIFF) and supplementary figures (`supplementary/`, PNG) |

## Reproducing the analysis

Paths at the top of each script (`ROOT`) point to the working directory used in
the study. Change `ROOT` to the location of your local copy before running.
Scripts are numbered; those with a shared prefix (for example `250-`, `251-`,
`252-`) generate the main figures and expect the tables produced by the earlier
steps.

`Docs/SCRIPTS_TABLES_FIGURES.md` in the private working repository maps each
script to the tables and figures it produces; the naming convention is
`Scripts/<number>_<topic>.R` -> `Tables/Table<number>_<topic>.csv` ->
`Figures/Fig<n>_<panel>.png`.

## Notes

- R 4.5.1 was used.
- Figures are provided at publication resolution (TIFF, 600 dpi for the main
  figures; PNG for the supplementary figures).
- Two very large profile tables are included for completeness and are not needed
  to reproduce the figures.
