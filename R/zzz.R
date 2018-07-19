# zzz.R - DESC
# /zzz.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

.onLoad <- function(libname, pkgname ) {
 if(.Platform$OS.type == "windows" & .Platform$r_arch != "x64")
 stop("This version of FLash on Windows can only run on 64-bit R.")
#	cat("The Creator Has a Mastertape\n")
}

