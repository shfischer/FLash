# fwd.R
# FLash/R/fwd.R
# Copyright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, Cefas
# Last Change: 21 Jan 2010 09:51
# $Id: fwd.R 232 2009-04-30 15:44:58Z fscott $


setGeneric('validSRPar', function(object, ...)
		standardGeneric('validSRPar'))


.validSRPar<-function(object, sr, yrs=NULL, availability=NULL)
     {
     #### check that sr has dims in FLQuant
     if (!all(names(sr) %in% c("params","year","unit","season","area","iter")))
        stop("dims in sr not recognised")
        
     #### Check yrs
     if (!is.null(yrs)){
        yrs<-ac(yrs)
        if (!all(yrs %in% ac(dims(object)$minyear:dims(object)$maxyear)))
           stop("yrs exceed years in object")
     }else
        yrs<-ac(dims(object)$minyear:dims(object)$maxyear)

     if ("year" %in% names(sr)){
        if (!all(yrs[-length(yrs)] %in% dimnames(sr)$year)) ## cos of extra year needed in C++
           stop("yrs exceed years in sr")}
     else{
         dmns  <-dimnames(sr)

         params<-list(params=dmns$params,year=yrs)

         if ("iter"   %in% names(dmns)) params[["iter"]]  <-dmns$iter
         if ("unit"   %in% names(dmns)) params[["unit"]]  <-dmns$unit
         if ("area"   %in% names(dmns)) params[["area"]]  <-dmns$area
         if ("season" %in% names(dmns)) params[["season"]]<-dmns$season

         srTmp<-FLQuant(c(sr),        dimnames=dimnames(sr))
         res  <-FLQuant(as.numeric(0),dimnames=params)
         sr   <-sweep(res,(1:6)[c("params","year","unit","season","area","iter") %in% names(dimnames(sr))],srTmp,"+")
         }

     #### create FLQuant compatible FLPar
     sr <-FLPar(as.FLQuant(as.data.frame(sr)))

     #### Check availability
     #### Needed to distribute recruits if SR$area==1
     if (dims(object)$area>1){
        if (is.null(availability))
           stop("availability needs to be provided if multiple areas")

        if (!all(yrs %in% dimnames(availability)$year))
           stop("years in availability mismatch")
        }

     #### Check iters
     niter<-unique(c(dims(object)$iter,1))

     if (!max(dims(sr)$iter) %in% niter)
        stop("Iters in sr don't match those in object")

     if (!is.null(availability))
        if (!(dims(availability)$iter %in% niter))
           stop("Iters in availability don't match those in object")

     dmns<-list(params=dimnames(sr)$params,
                year  =yrs,
                unit  =dimnames(m(object))$unit,
                season=dimnames(m(object))$season,
                area  =dimnames(m(object))$area,
                iter  =dimnames(sr)$iter)

     #### check units, seasons and areas
     if (any(unlist(dims(        object)[c("season","area","unit")])>1) |
         any(unlist(dims(as.FLQuant(sr))[c("season","area","unit")])>1)){

        #### check units
        if (!all(dimnames(sr)$unit %in% dimnames(m(object))$unit))
           stop("unit in sr and object don't match")

        #### check season
        if (!all(dimnames(sr)$season %in% dimnames(m(object))$season))
           stop("season in sr and object don't match")

        #### check area
        if (!all(dimnames(sr)$area %in% dimnames(m(object))$area))
           stop("area in sr and object don't match")
        }

     res<-FLQuant(as.numeric(NA),dimnames=dmns)

     dm<-list(unit=dimnames(sr)$unit,season=dimnames(sr)$season,area=dimnames(sr)$area)
     if (dim(res)[3]==1 & dm$unit  =="unique") dm$unit  <-dimnames(res)$unit
     if (dim(res)[4]==1 & dm$season=="all")    dm$season<-dimnames(res)$season
     if (dim(res)[5]==1 & dm$area  =="unique") dm$area  <-dimnames(res)$area

##bug if sr doesn´t match control
     res[,dimnames(sr)$year,dm$unit,dm$season,dm$area,]<-as.FLQuant(sr)

     return(FLPar(res))
     }

setMethod('validSRPar', signature(object='FLStock'),
  function(object, sr, yrs=NULL, availability=NULL)
     {
     return(.validSRPar(object,sr,yrs,availability))
     })

setMethod('validSRPar', signature(object='FLBiol'),
  function(object, sr, yrs, availability)
     {
     return(.validSRPar(object,sr,yrs,availability))
     })

#setMethod('validSRPar', signature(object='FLBRP'),
#  function(object, sr, yrs)
#     {
#     return(.validSRPar(object,sr,yrs,availability(object)))
#     })

validSRRes<-function(object,sr,res){
   if (dims(res)$iter   != dims(object)$iter) stop("iters in residuals and object don't match")
   if (dims(res)$year   != dims(sr)$year)     stop("years in residuals and sr don't match")
   
   if (dims(res)$unit   != dims(sr)$units)    stop("unit in residuals and sr don't match")
   if (dims(res)$area   != dims(sr)$area)     stop("area in residuals and sr don't match")
   if (dims(res)$season != dims(sr)$season)   stop("season in residuals and sr don't match")
   }
    
