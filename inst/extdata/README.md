# Derived SIPP analysis matrix

`sipp_2014_panel.csv.gz` is the complete data matrix used by the SIPP example
in the package README and `inst/validation/longitudinal_clpm_example.R`. It is a
gzip-compressed CSV that can be read directly with `read.csv()`.

The matrix was derived from the U.S. Census Bureau's public-use 2014 Survey of
Income and Program Participation (SIPP) panel, waves 1--4. Those waves cover
the 2013--2016 reference years. From each public-use wave file, the derivation
used `SSUID`, `PNUM`, `MONTHCODE`, `TAGE`, `TPEARN`, and `TMWKHRS`; none of
those source columns or identifiers is included here. Survey weights are also
omitted because the current MCMSEM moment/SE pipeline is unweighted; this is a
methodological illustration, not a population-representative labor estimate.

The derivation:

1. selected December (`MONTHCODE == 12`) in each wave;
2. retained people aged 25--57 in wave 1, so the baseline cohort remained no
   older than 60 in wave 4;
3. treated non-positive earnings or usual weekly hours, and their paired value
   at that wave, as unavailable;
4. transformed earnings with the natural logarithm;
5. used the wave-1 complete cases to calculate a separate mean and standard
   deviation for log earnings and hours, then applied the same affine
   transformation `2 * (value - wave1_mean) / wave1_sd` at every wave; and
6. removed rows with no usable observation at any wave.

The reference values were 8.0239060027 and 0.9187296559 for the mean and
standard deviation of log earnings, and 40.7161889428 and 11.9133321986 for
hours. The resulting file has 24,505 rows and eight columns. Jointly observed
sample sizes are 22,049, 15,787, 12,190, and 10,446 in waves 1--4; 6,647 rows
are complete at all four waves. The compressed file's MD5 checksum is
`40890990d7bf747323627a66872e397e`.

Official sources:

- [SIPP 2014 Panel Data](https://www.census.gov/programs-surveys/sipp/data/datasets/2014-panel.html)
- [2014 Panel Wave 1 public-use data and documentation](https://www.census.gov/programs-surveys/sipp/data/datasets/2014-panel/wave-1.html)
- [2014 Panel Wave 4 public-use data and documentation](https://www.census.gov/programs-surveys/sipp/data/datasets/2014-panel/wave-4.html)
- [2014 Panel data dictionaries](https://www.census.gov/programs-surveys/sipp/tech-documentation/data-dictionaries/data-dictionaries-2014.html)

The original public-use wave files are intentionally not vendored in this
repository.
