# radiation_code
Standalone [RegCM](https://github.com/ICTP/RegCM/) radiation code driven by
[ERA5](https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-complete?tab=overview)
atmospheric dataset on model level and the original reduced gaussian
[N640](https://confluence.ecmwf.int/display/EMOS/N640) grid.

# Description
The aim of this software is to provide a reference baseline for optimization.
The original code used in RegCM is coming from [NCAR CCSM4](https://www2.cesm.ucar.edu/models/ccsm4.0/cam/docs/description/cam4_desc.pdf) and can be found in
the [ccm-phys](https://github.com/graziano-giuliani/radiation_code/tree/main/ccm-phys)
directory for reference.

# RegCM modifications of CCM3 code

The RegCM is modified from the NCAR-CCSM by:

- Solar constant is function of time inline with [input4MIPS](https://aims2.llnl.gov/search/input4MIPs/)
by using the [SOLARIS-HEPPA](https://solarisheppa.geomar.de) dataset.
- GHG concentrations for CO2, NH4, N2O, CFC11, CFC11 are taken from the
input4MIPS CMIP SSP scenario data. In radiation\_code, the SSP370 scenario is
used.
- Solar zenith angle computation use code derived from [CESM](https://www.cesm.ucar.edu/)
[CLM](https://www.cesm.ucar.edu/models/clm) code.
- Aerosol optical properties are taken from [MACv2-SP](https://gmd.copernicus.org/articles/10/433/2017/)
Max Planck Institute Aerosol Climatology version 2 inline with CMIP6 input4MIPS.
- Saturation water vapor pressure and mixing ratio computation from pressure
and temperature use [Piotr J. Flatau, et al](https://journals.ametsoc.org/view/journals/apme/31/12/1520-0450_1992_031_1507_pftsvp_2_0_co_2.xml)
polynomial fit.
- Ice and cloud particles effective radii is derived respectively from
[Stengel](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2022GL102521)
and [Savijärvi](https://a.tellusjournals.se/articles/10.3402/tellusa.v50i1.14508)
- Effective cloud fraction for Longwave radiation is derived from
[Savijärvi](https://a.tellusjournals.se/articles/10.3402/tellusa.v50i1.14508)
- Cloud fraction profile at interfaces uses Maximum Random Overlap assumption.

# Data
Data to run the application needs to be downloaded from ECMWF. In the *data*
directory, the *mlevel* and *surface* directories must be populated. The user
is provided with two python scripts using the [CDS API](https://cds.climate.copernicus.eu/api-how-to)
that can be used with a *python3* to download the dataset. Static data populate
the other directories.

# Compilation
To compile the code, a Fortran 2003 compiler is required (tested with GNU
*gfortran9* on Ubuntu Linux), along with the [netCDF](https://www.unidata.ucar.edu/software/netcdf/)
Fortran library (tested with version 4.6.1) and the *make* program
(tested with GNU *make* 4.2.1). Type

`make`

and the program *radiation_code* is built and can be run. Only text output
is provided, any modification involves editing the Fortran soirce code and is
left to the user.

# Happy hacking.

Graziano.
