# Robust Mahalanobis distances for multivariate functional data — code and results

Reproduces every figure in the paper and its supplement. The welding and
fertility examples each require one data file that is downloaded separately; see
[Data](#data).

## Two levels of reproduction

**Level 1 — all figures from the stored results (minutes).** Everything needed is
in `1_reproduce_figures/`: the scripts, and the result files they read
(`stored_results/`, ~130 MB). Level 1 runs in minutes on a laptop, 
no special hardware.

**Level 2 — the results from scratch (more than 10 CPU-years).**
`2_run_simulations/` holds the simulation scripts to reproduce the simulations. 
They are run setting by setting using the parameter block at the top of each file.
Computation time is several days on a large server.

Every simulation script writes one file per setting into
`2_run_simulations/raw_results/`, which ships empty: the tables built from those
files are already in `stored_results/`, and that is all the figure scripts read.
After re-running anything, rebuild them with `2_run_simulations/merge_results.R`.
This is the only link between the two levels, and the only script that writes to
`1_reproduce_figures/stored_results/`. All flattening and label preparation
happens there, so the figure scripts only load a table and plot it. For the
p = 3 and p = 50 settings it combines the four outlier types; for the
non-separable and heavy-tailed settings, each produced by a single run, it only
carries the file over. `build_shapley_subset.R` does the same for the Shapley
run, and `reduce_results.R` then applies the size reduction described below.

Each script opens its own cluster with `makeCluster(120)`, except the four
`simulation_main_p50_*_FIF.R` scripts, which use `makeCluster(70)`: MFIF holds its
forests in Python memory per worker, and at p = 50 one hundred workers exhausted
750 GB of RAM. Those scripts also release each forest once its scores are read
back (`gc()` and `gc_py$collect()` after every dictionary), so Python memory does
not accumulate across the design grid. **Set the worker count to match your
machine before re-running** — one line per script.

## Layout

```
functions/                     shared helper functions
LICENSE                        MIT for the code, CC BY 4.0 for the results
MFIF.py                        vendored MFIF implementation
MFIF_LICENSE.txt               its license, retained as required

1_reproduce_figures/           LEVEL 1 — minutes
  simulation/                  figure scripts + plots/
  fertilitydata/               real-data example + plots/ (data downloaded separately)
  nino_nina/                   real-data example + data/ + plots/
  welding/                     real-data example + plots/ (data downloaded separately)
  stored_results/              the result tables the figure scripts read

2_run_simulations/             LEVEL 2 — CPU-years
  simulation_main_*.R          the simulation scripts
  merge_results.R              rebuilds the merged tables from raw_results/
  build_shapley_subset.R       rebuilds the Shapley subset
  reduce_results.R             the size reduction described below
  raw_results/                 where the simulation scripts write; ships empty
```

## Which script makes which figure

All paths relative to `1_reproduce_figures/`.

| script | figures |
|---|---|
| `simulation/simulation_overview_plots.R` | outlier-type overview (simulates its own data) |
| `simulation/simulation_summary.R` | separable Gaussian settings, p = 3 and p = 50 |
| `simulation/simulation_summary_non_separable.R` | non-separable settings |
| `simulation/simulation_summary_t_distribution.R` | heavy-tailed settings |
| `simulation/shapley_simulation_summary.R` | Shapley localisation |
| `simulation/simulation_summary_non_separable_mmle_vs_sample.R` | separable vs sample covariance |
| `fertilitydata/fertilitydata.R` | fertility example |
| `nino_nina/nino_nina.R`, `nino_nina_map.R` | El Niño example |
| `welding/welding.R` | welding example |

## Requirements

Run `Rscript check_packages.R` first: it lists any missing packages and prints the
`install.packages()` call for them. It installs nothing.

R (tested with 4.4.1) and:

- `robustmatrix` (>= 0.1.5), which provides the MMCD estimators
- `covsep`, for the separability tests in the real-data examples
- `this.path`, used by every script to set its working directory to its own
  location, so they run unchanged from RStudio, `Rscript`, or `source()`
- `tidyverse`, `ggh4x`, `ggridges`, `ggnewscale`, `ggrepel`, `geomtextpath`,
  `gridExtra`, `latex2exp`, `lvplot`, `Polychrome`, `plot.matrix`, `tsutils`,
  `knitr`, `kableExtra` for the summaries and figures
- `fda`, `mrfDepth`, `fdaoutlier`, `pROC`, `robustbase` for the functional
  detectors compared against
- `abind`, `cellWise`, `MixMatrix`, `rSPDE`, `foreach`, `doParallel`, `doSNOW`
  to re-run the simulations
- `countrycode`, `zoo` for the real-data examples
- `sf`, `rnaturalearth` and `rnaturalearthdata` for the El Niño map
  (`nino_nina/nino_nina_map.R`). `rnaturalearthdata` is only a *suggested*
  dependency of `rnaturalearth`, so it is not pulled in automatically, but
  `ne_countries()` needs it. `rnaturalearth` also brings `terra`, and `rSPDE`
  brings `fmesher`; on Linux both need the system GDAL, GEOS and PROJ libraries.
- Python with `reticulate` for the MFIF comparisons; `MFIF.py` is vendored
  unmodified from <https://github.com/GuillaumeStaermanML/FIF> and is covered by
  its own license, reproduced in `MFIF_LICENSE.txt`

## Data

The simulation study is self-contained. Of the three real-data examples, one
carries its data, and two use data owned by others, which are downloaded from the
original source rather than redistributed here:

| example | data | included |
|---|---|---|
| El Niño | Niño 3.4 sea surface temperatures, NOAA ERSST v5 | yes, `nino_nina/data/` |
| Welding | dynamic resistance curves from resistance spot welding | no, `welding/data/README.md` |
| Fertility | age-specific fertility rates, Human Fertility Database | no, `fertilitydata/data/README.md` |

Each of the two `data/README.md` files names the single file the script needs,
where to download it, and the version the reported figures are based on. Both
examples need one file each and nothing else; their figures are included in the
respective `plots/` folders, and no other script depends on them.

The fertility rates are covered by a Human Fertility Database user agreement and
are downloaded with an individual account. The welding data are distributed by
the robust statistics group at KU Leuven together with the method we compare
against.

## Notes on the stored results

**Size reduction.** `reduce_results.R` shrinks the stored tables by dropping 
sixteen columns that no figure uses (intermediate covariance errors and the errors 
relative to the clean-data benchmark) and rounding the remaining numeric columns to 
6 significant digits. Design columns are untouched.
Run it after `merge_results.R`: a fresh simulation run writes the full column set,
so the merged tables inherit it. The shipped tables are already reduced, and the
script is idempotent, so running it again changes nothing.

**MFIF** yields an anomaly score but no cut-off, so it is evaluated by AUC alone.
To keep a common interface the simulation flags the `ceiling(n * eps)`
highest-scoring curves, which means MFIF is *given* the true number of outliers;
its `n_flagged`, `TP`, `FP`, `precision`, `recall` and `F.score` are therefore not
comparable with the other methods, so every figure reporting count-based scores
drops MFIF. 

For **p = 50** MFIF was run with 10 replications per configuration rather than
100, because of its runtime. It is excluded from the rank-based comparison, which
needs balanced replications, so it is kept in its own table
(`simulation_res_p50_FIF.RData`) instead of being merged into
`simulation_res_p50.RData`. At p = 3 it has the full 100 replications and is
merged into `simulation_res_p3.RData`.

**Raw-data methods** do not depend on the number of basis functions, so they are
computed once per configuration, in the iteration with the smallest basis count.
Their rows are present but empty in the other iterations; the summary scripts drop
them with `drop_na()`, which leaves a balanced design.

**Use of AI tools** AI coding assistants (Anthropic, Claude Opus 5) were used for
implementation support: restructuring the code, writing the scripts that assemble 
and reduce the stored results, and drafting this documentation. 
All ideas, designs and methodological choices are those of the
authors, who are responsible for the code and for the reported results.
