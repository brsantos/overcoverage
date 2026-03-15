pkgname <- "overcoverage2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('overcoverage2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("model_BLB")
### * model_BLB

flush(stderr()); flush(stdout())

### Name: model_BLB
### Title: Run the BLB model
### Aliases: model_BLB

### ** Examples

# See tests/model_BLB_simulated.R for a simulated run



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
