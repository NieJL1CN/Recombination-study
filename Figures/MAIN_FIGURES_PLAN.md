# Main figure plan (TAG submission)

TAG specifications: single column 89 mm, double column 183 mm, maximum height 247 mm; line art at 600-1200 dpi, combination 600 dpi, halftone 300 dpi; sans-serif fonts no smaller than 5 pt at print size (7 pt preferred); vector (PDF/EPS) preferred for line art.

## Audit of the current figure folder

- 148 image files.
- 6 exceed 183 mm in width; 1 is below 300 dpi.
- Most panels sit between 89 and 183 mm but are single panels, not composed figures.
- Six candidate panels were resized to the target width and saved in `Figures/main/`:
  - Fig1a_CMplot_density_AXIOM.1_183mm.png
  - Fig3h_LD_decay_89mm.png
  - Fig3i_CO_metaQTL_183mm.png
  - Fig4e_CO_QTL_trait_scan_89mm.png
  - FigS18_CO_landrace_vs_cultivar_89mm.png
  - FigS_fixation_ParW049_89mm.png

## Proposed five main figures

### Fig. 1. A dense consensus map and the crossover landscape (double column, 183 mm)
- (a) Linkage disequilibrium decay with the 15-Mb window indicated (source: Fig3h_LD_decay).
- (b) Marker density and collinearity of the consensus map (source: Fig1a_CMplot_density, Fig1d_collinearity).
- (c) Crossovers by chromosome region, R3 > R2b > R1 > R2a > C (new bar panel from Table18b).
- (d) Crossover number against chromosome length, r = 0.76 (source: Fig3e_CO_vs_length_chr).
- (e) Crossover number per RIL across the 18 families (source: Fig3b_CO_violin).

### Fig. 2. Nine crossover QTL (double column, 183 mm)
- (a) Genome-wide QTL map with the nine loci, LOD and direction (source: Fig3i_CO_metaQTL).
- (b) Effect size and R² for each QTL (new forest panel from Tables 1 and 6).
- (c) Allele means for the two homozygous classes at each peak (new panel from Table 6).

### Fig. 3. The QTL are independent of agronomic traits (single column, 89 mm)
- (a) LOD before and after adding heading date for the three testable families (Table 8).
- (b) P values for the six agronomic traits across the nine comparisons (Table 7).

### Fig. 4. Landrace and cultivar populations (double column, 183 mm)
- (a) Crossover number per RIL in 18 Watkins families and nine cultivar populations (source: FigS18).
- (b) Effect of the Paragon allele at 4A in nine independent Paragon populations, with pooled estimate (new forest panel from Table 16).
- (c) Repeatability: peaks at LOD >= 2 in shared families and their chromosome match to 2017 (Table 9).

### Fig. 5. What a recombination QTL is worth in breeding (double column, 183 mm)
- (a) Time to 50% TTdd against target-drag distance, with and without a QTL (Table 4).
- (b) Time to 50%, 99% and 100% under the three selection schemes (Table 13).
- (c) Probability of at least one TTdd in F2 against population size, for four recombination increases (Table 15).

## Status

- Panels available as separate images: Fig. 1a,b,d,e; Fig. 2a; Fig. 3 (needs composition); Fig. 4a; Fig. 5a,c.
- Panels to be drawn: Fig. 1c, Fig. 2b,c, Fig. 3a,b, Fig. 4b,c, Fig. 5b.
- Composition and final export at 183 or 89 mm, 600 dpi, fonts >= 7 pt, to be done with the academic-figure-skill scripts.
