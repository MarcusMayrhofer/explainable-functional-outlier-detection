# Per-setting simulation output

This folder is where the scripts in `2_run_simulations/` write their results,
one file per setting, e.g.

```
simulation_res_p50_shift.RData
simulation_res_p50_shift_FIF.RData
```

It ships empty. The files are not redistributed because the merged tables built
from them are already included in `1_reproduce_figures/stored_results/`, which is
all the figure scripts read.

`merge_results.R` reads this folder and writes those merged tables, so it can
only be run after the corresponding settings have been re-computed. Keep this
folder in place: the simulation scripts write into it and `save()` fails if it
does not exist.
