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
Refer to [this webpage](https://qepanna.quarto.pub/goflux/import.html) to learn how to import raw data using the goFlux package.

## Visualize incubation timeseries and save as pdf file
It can be covenient to first have look at the incubation timeseries. This can be easily achieved with the `plot.incubations()` function.
Plots can be saved with `gg_save_pdf()`.
``` r
p <- plot.incubations(mydata)
print(p)
# to save these plots in a dedicated path, do
gg_save_pdf(list = p, path = "mypath", filename = "myfilename.pdf")
```

## Selection of valid observations, trimming incubation timeseries
After import, the user can either define manually the start and end points of each measurement, or proceed with an automatised
selection of the observations with a set of different methods:

	### Manual inspection and selection of valid data points
	Just as in the `goFlux` R package, the **manual selection of the measurements** is based on `start.time`, provided separately 
	in an auxiliary file.
	The function `obs.win` splits the imported data into a list of data frame (divided by `UniqueID`) and creates an observation 
	window around the `start.time` to allow for a manual selection of the start and end points of each measurements, using the 
	function `click.peak2`.
	This manual selection is then used to trim the timeseries and calculate fluxes.
	
	### Automatic data selection
	The **automatic identification of the measurements** has two methods:
	- "trust.it.all": the user is willing to trust the incubation timeseries as much as possible, with the possibility to trim 
	the incubation timeseries by a defined time window (see `shoulder` parameter). This method allows for non-linear patterns in the timeseries.
	- "keep.it.linear": the method is first automatically identifying linear patterns after the beginning of the incubation, and discard the measurements
	when they diverge away from a linear pattern.

## Separate Diffusion from Ebullition in CH4 measurements
One typical pathway for CH4 in aquatic ecosystems is through ebullition, or bubbling. Here the function `flux.separator()` automatically
separates diffusive linear patterns from abrupt and highly non-linear ebullition patterns. This is done by first identifying the first
linear chunk in the incubation time series which is used to calculate diffusion. The total flux is calculated by calculating the total
concentration change inside the chamber over the incubation time, and ebullition is then calculated by difference:

$$\mathbf{Eqn~1}~~~~~~F_ebullition = F_total - F_diffusion$$

## Flux calculation
**Diffusive** fluxes are calculated with the `goFlux` functions. In brief, the function `goFlux` calculates fluxes from a variety of greenhouse
gases (CO<sub>2</sub>, CH<sub>4</sub>, N<sub>2</sub>O, NH<sub>3</sub>,
CO, and H<sub>2</sub>O) using both linear (LM) and non-linear (HM;
[Hutchinson and Mosier,
1981](https://doi.org/10.2136/sssaj1981.03615995004500020017x)) flux
calculation methods. 
The HM model for the chamber concentration $C_t$ at
time $t > 0$ after deployment is given by:

$$\mathbf{Eqn~2}~~~~~~C_t = \varphi~+~(C_0 - \varphi)e^{-~\kappa~t}$$

Where $\varphi$ is the assumed concentration of constant gas source
below the surface (also known as $C_i$), $C_0$ is the concentration in
the chamber at the moment of chamber closure and $\kappa$ (kappa)
determines the curvature of the model. A large kappa returns a strong
curvature.

$$\mathbf{Eqn~3}~~~~~~k.max = \frac{LM.flux}{MDF~\times~t}$$

Where $LM.flux$ and $MDF$ have the same units (nmol or
µmol·m<sup>-2</sup>·s<sup>-1</sup>) and $t$ is in seconds. Therefore,
the units of kappa-max is s<sup>-1</sup>. This limit of kappa-max is
included in the `goFlux` function, so that the non-linear flux estimate
cannot exceed this maximum curvature. See below for more details about
the minimal detectable flux (MDF).

All flux estimates, including the MDF, are multiplied by a $flux.term$
which is used to correct for water vapor inside the chamber, as well as
convert the units to obtain a term in nmol or
µmol·m<sup>-2</sup>·s<sup>-1</sup>:

$$\mathbf{Eqn~4}~~~~~~flux.term = \frac{(1 - H_2O)~V~P}{A~R~T}$$

Where $H_2O$ is the water vapor in mol·mol<sup>-1</sup>, $V$ is the
volume inside the chamber in liters, $P$ is the pressure in kPa, $A$ is
the surface area inside the chamber in m<sup>2</sup>, $R$ is the
universal gas constant in L·kPa·K<sup>-1</sup>·mol<sup>-1</sup>. Each
parameters are measured inside the chamber at $t = 0$.

More details and demonstrations about `goFlux` can be found [here](https://qepanna.quarto.pub/goflux/goFlux.html).

**Ebullition** is calculated by difference between the total estimated flux and the diffusion flux.



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


## Community Guidelines

Author: Camille Minaudo

The package is ready to use and fully functional, but errors may still
occur. To report any issues, suggest improvements or ask for new
features, [open an issue on
GitHub](https://github.com/camilleminaudo/aquaGHG/issues). Alternatively,
contact directly the maintainer, Camille Minaudo (<camille.minaudo@ub.edu>),
including script and raw data if necessary. Thank you for helping me in
the development of this tool! 🙏

## How to cite


