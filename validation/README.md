# Validation of OLJC and LOSH functions

As part of the function implementation process, each of the LJC and LOSH functions I am implementing throughout GSOC2020 will need to be validated on 1) a 'toy' dataset and 2) severla real-world datasets. These validations will confirm the behavior of the function and are necessary before making the functions widely available for public use.

## Validation plan

**LOSH**: 'toy' dataset will be the routinely used `neighborhoods.gpkg` and `listings.gpkg`. The 'real-world' datasets will compare the function's output from the `R` package `spdep::LOSH`.

**LJC**: 'toy' dataset will be the routinely used 4x4 spatial grid of `np.ones()` used in the validation of other `esda` functions. The 'real-world' datasets will compare the function's output from the GeoDa software.

## Folder organization 

There will be a separate notebook for each of the four statistics: LOSH, LJC univariate, LJC bivariate, and LJC multivariate.

- Validation_LOSH.ipynb
- Validation_LJC_univariate.ipynb
- Validation_LJC_bivariate.ipynb
- Validation_LJC_multivariate.ipynb