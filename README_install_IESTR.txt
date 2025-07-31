Welcome in the IESTR package.

Before you start, make sure you installed clang on your computer so that Rcpp can work correctly.
Also, user need to set their working directory in the scripts.

You could install the package with the zip provided or install it from github with the commands:
library(devtools)
install_github("ArndCallebaut/IESTR")

Four scripts are providing material to do the figures of the article associated with the package, in "scripts_for_article" :
-building_example_functions.R contain function to create the theorical case study.
-plot_exemple_functions.R contain functions to plot informations and results of the package.
-SupplementaryMaterial.R contain additional plot, used for the figures of the paper (legends, ect)
-IESTR_main_script.R is the main script, and show how to use IESTR.

Note that the vignettes of the package use the same theorical study as a showcase of the package.

Several libraries are used for the example :
IESTR, dyplr, purrr, Matrix, grid, raster, lattice, ggplot2, viridis, hrbrthemes, methods.


