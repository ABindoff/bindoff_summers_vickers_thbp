# Cognitive trajectories in the Tasmanian Healthy Brain Project: analysis code

Code for the analyses reported in the manuscript and its Supplementary Information.

## The data are synthetic

The study data cannot be shared. `data/synthetic/` holds a synthetic release made by
`thbp_synthesise_all_domains.R` (synthpop 1.9.3, CART, seed 2026). All four cognitive domains
were synthesised jointly from one wide file, conditional on baseline characteristics, and
records that replicated a unique real combination of group, sex, age and education were
removed. The release has 358 synthetic participants (the study has 416).

The synthetic data keep the broad structure of the study data, so the analyses run and give
results of a similar kind, but **they will not reproduce the published estimates**. The
point of this repository is that the analysis can be repeated and inspected, not that the
numbers can be matched.

## Files

| File | Contents |
|---|---|
| `analysis.R` | Main analysis: fits the four cognitive-trajectory models and prints the coefficients and evidence ratios reported in the manuscript, beside the published values |
| `smoothbp_supplementary_information2.Rmd` | Supplementary Information as submitted, reading the synthetic data |
| `data/synthetic/{epis,wm,ef,lang}.csv` | Episodic memory, working memory, executive function and language processing; one row per participant-wave |
| `thbp_synthesise_all_domains.R` | The script that made the synthetic data (needs the real data, so it will not run here) |

The synthetic files contain the composite domain scores, not the sub-tests; the composites
were built as described in the manuscript Methods.

## Running

Requires R (developed under 4.6.1) and `smoothbp` (≥ 0.2.8) with `posterior`. The
supplement also needs `tidyverse`, `gt`, `loo` and `brms`. Open the `.Rproj` file, or set
the working directory to this folder, then:

```r
source("analysis.R")
rmarkdown::render("smoothbp_supplementary_information2.Rmd", output_format = "html_document")
```

Each trajectory model runs 4 chains of 10,000 iterations, as in the manuscript, so allow
some time. `analysis.R` leaves the fitted models in `fits` for inspecting trace plots and
MCMC diagnostics (see the comments at the top of the script). The supplement's
default PDF output needs XeLaTeX and the Cambria font, hence the HTML output above.

## Running on the study data

Reviewers given the study data under agreement should place `episodic.rds`, `working.rds`,
`executive.rds` and `language.rds` in `data/`. When all four are present `analysis.R` loads
them unchanged in place of the synthetic data, and fits the models exactly as for the
manuscript (same data, same seeds). The script states which data it was run on.
These files are excluded by `.gitignore`.

## Differences from the submitted supplement

Only the data input has changed. Supplementary Tables 6, 8, 9 and 10 (transition sharpness,
transition-point prior bounds, transition point by group, and sub-test refits) come from
separate sensitivity runs that are not included, so those tables are omitted when rendered
here. The sub-test refits cannot be repeated on the synthetic data in any case. The prose in
the supplement is unchanged and describes the real-data results.
