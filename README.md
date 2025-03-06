# aquaGHG
## an R package for reproducible calculation of GHG fluxes from static or floating chamber measurements in aquatic ecosystems


## About the package

`aquaGHG` is an easy-to-use and generic R package with the objective to provide reproducible calculation of GHG fluxes from static or floating chamber measurements in aquatic ecosystems. 
It uses the `goFlux` R package (https://github.com/Qepanna/goFlux) to import raw measuements. Look up [this webpage](https://qepanna.quarto.pub/goflux/) for a
demonstration of the `goFlux` package usage.
Then, mainly provides two options:
- manually inspect the user's selection of incubation to select valid observations and proceed with subsequent flux calculation using both linear (LM) and non-linear (HM;[Hutchinson and Mosier, 1981](https://doi.org/10.2136/sssaj1981.03615995004500020017x)) flux
calculation methods.
- automatically processing of the data, from harmonized incubation timeseries to flux estimates using both LM and HM methods, without any manual step with several options regarding the automatic selection of valid data points.

## Import data using `goFlux`

## Manual inspection and selection of valid data points

## Automatic data selection

## Flag bubbling events when CH4 data is available

## Flux calculation

## Chosing betwen linear and non-linear methods

## Installation

To install a package from GitHub, one can use the package `devtools` (or
`remotes`) from the CRAN:

``` r
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
```

Then, install the `aquaGHG` package from GitHub. If it is not the first
time you install the package, it must first be detached before being
updated.

> The package is constantly being updated with new functions or
> de-bugging. To make sure that you are using the latest version,
> re-install the package every time you use it.

``` r
try(detach("package:aquaGHG", unload = TRUE), silent = TRUE)
devtools::install_github("camilleminaudo/aquaGHG")
```

**If prompted, it is recommended to update any pre-installed packages.**

## How to cite


