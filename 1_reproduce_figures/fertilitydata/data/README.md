# Fertility data

`fertilitydata.R` reads one file, which is not included here:

```
data/asfrRR.txt
```

It comes from the Human Fertility Database and is redistributed under a user
agreement, so it has to be downloaded with your own account.

## How to obtain it

1. Register at <https://www.humanfertility.org> and log in.
2. Go to **Data** and download the zipped collection of all countries for
   *Age-specific fertility rates by calendar year and age (Lexis squares, ASFR,
   period)*.
3. Place `asfrRR.txt` in this folder.

The file the figures were produced from carries the header

```
Period fertility rates by calendar year and age (Lexis squares, age in completed years (ACY))
Last modified: 14/11/2023
```

and begins with two header lines followed by the columns `Code`, `Year`, `Age`,
`ASFR`. `fertilitydata.R` skips those two lines and keeps ages 15--45, so a later
release of the same file works as long as that layout is unchanged.

Nothing else in the repository depends on this file: the six figures of the
fertility example are included in `../plots/`, and every other figure reproduces
without it.
