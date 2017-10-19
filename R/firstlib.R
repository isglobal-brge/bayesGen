
############ First.lib ###############

.onLoad <- function(lib, pkg){
   require(rjags)
   require(INLA)
   library.dynam("bayesGen", pkg, lib)
}

.onUnload <- function(libpath)
    library.dynam.unload("bayesGen", libpath)


############ End of .First.lib ###############


