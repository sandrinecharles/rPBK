# pbtkDB

An R package for the inference of physiologically-based toxicokinetics model

# Development

## Package skeleton

```R
library("rstantools")
rstan_create_package(path = 'pbtkDB')
```

## Dev info

### Rstudio

Before 'Documentation' and 'Clear and Rebuild' and within 'Configure Build Tools...'
untick 'use devtools package function if available'.
Keep 'Generate Configurationwith Roxygen'

### little hack

- to load all internal function of a package during dev: `devtools::load_all()`

### A lighter package build

To make the package lighter, we have to remove the vignettes: see file `.Rbuildignore`

### Error to recompile during package dev

Sometimes, there is an Error to recompile during development after change of .stan files.
A solution is to remove the `rbioacc` folder in R repository of the win-library (see the path written in the error message).

An other solution is to build the package from the terminal using `R CMD -preclean INSTALL rbioacc` from parent directory of `rbioacc`.


## Dev tools

- S3 Object System: http://adv-r.had.co.nz/S3.html
- Google's R Style Guide: https://google.github.io/styleguide/Rguide.html
- testthat: https://github.com/r-lib/testthat
- covr: https://github.com/r-lib/covr
- to add a package `xxr`: `usethis::use_package("xxr")`
- to add data set `datar`: `usethis::use_data(datar)`



