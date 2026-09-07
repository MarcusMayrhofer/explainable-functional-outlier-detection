# El Niño–Southern Oscillation data

`data_ERSSTv5.txt` holds monthly sea surface temperatures and anomalies for the
four Niño regions, as published by the U.S. Climate Prediction Center (CPC) and
derived from the Extended Reconstructed Sea Surface Temperature, Version 5
(ERSST v5) dataset.

## Source

<https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii>

Retrieved 24 October 2025. The copy included here runs from January 1950 to
September 2025; the CPC updates the file monthly, so a later download extends it
and the number of complete June-to-May periods grows accordingly.

The region definitions used in `nino_nina_map.R` come from
<https://www.cpc.ncep.noaa.gov/data/indices/>.

## Format

Whitespace-separated, one header line, one row per month:

```
 YR   MON  NINO1+2  ANOM   NINO3    ANOM   NINO4    ANOM   NINO3.4  ANOM
1950   1   23.01   -1.55   23.56   -2.10   26.94   -1.38   24.55   -1.99
```

`YR` and `MON` are the calendar year and month; each region contributes an
absolute temperature in degrees Celsius and its anomaly against the 1991--2020
baseline. `nino_nina.R` reads the file with
`read.csv(..., sep = "", head = TRUE)`, which makes the header syntactic: the
four anomaly columns share the name `ANOM` and become `ANOM`, `ANOM.1`, `ANOM.2`,
`ANOM.3`, and `NINO1+2` becomes `NINO1.2`. The script relies on that ordering.

## Terms

The data are a product of the U.S. National Oceanic and Atmospheric
Administration and are in the public domain, and are therefore included here.
NOAA asks that the source be acknowledged; the paper cites ERSST v5 as
Huang et al. (2017).
