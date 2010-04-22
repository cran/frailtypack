############ First.lib ###############

.onLoad <- function(lib, pkg){
   require(survival)
   library.dynam("frailtypack", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("frailtypack", libpath)

############ End of .First.lib ###############




