# zzz.R - DESC
# /zzz.R

# Copyright European Union, 2017
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

.onLoad <- function(libname, pkgname ) {
  if(.Platform$OS.type == "windows" & .Platform$r_arch != "i386")
	  stop("FLash on Windows can only run on 32-bit R. Please change and retry")
}

