# Resistance spot welding data

`welding.R` reads one file, which is not included here:

```
data/X_RSW.txt
```

## How to obtain it

The data are distributed by the robust statistics group at KU Leuven, on the
software page accompanying the paper that introduces the cellwise robust
covariance estimator we compare against:

<https://wis.kuleuven.be/statdatascience/robust/software>

Download `X_RSW.txt` from there and place it in this folder.

## What the file contains

A plain table of 115 rows, one per welded spot. The first column, `labels`, is
the reference outlier indicator; the remaining 750 columns are five dynamic
resistance curves of 150 time points each, concatenated end to end.
`welding.R` reshapes them into a $5 \times 150 \times 115$ array.

## Attribution

The dataset is used in

- Centofanti, F., Hubert, M. and Rousseeuw, P. J. (2025). Cellwise and casewise
  robust covariance in high dimensions. arXiv:2505.19925.
- Capezza, C., Centofanti, F., Lepore, A. and Palumbo, B. (2024). Robust
  multivariate functional control chart. *Technometrics*.

Please cite these when using the data.

Nothing else in the repository depends on this file: the figures of the welding
example are included in `../plots/`.
