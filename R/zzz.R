############ First.lib ###############

.First.lib <- function(lib, pkg){
   require(survival)
   library.dynam("frailtypack", pkg, lib)
}
############ End of .First.lib ###############

