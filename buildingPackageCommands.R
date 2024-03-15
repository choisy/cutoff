sessionInfo()

detach("package:cutoff")
detach("package:cutoff2")
unloadNamespace("cutoff")
unloadNamespace("cutoff2")

setwd(here::here('cutoff2'))
roxygen2::roxygenise()

setwd(here::here())
devtools::build("cutoff2", vignettes=TRUE) # FALSE

devtools::check("cutoff2")

devtools::install("cutoff2", build_vignettes = TRUE) # FALSE

help.start()

setwd("cutoff2")
devtools::build_vignettes()
