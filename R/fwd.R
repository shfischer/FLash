# fwd.R
# FLash/R/fwd.R
# Copyright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, Cefas
# Last Change: Thu Jul 04, 2013 at 11:21 AM +0200
# $Id: fwd.R 1797 2012-12-15 11:11:41Z lauriekell $


# fwd(FLStock, fwdControl, missing) {{{
# Dreadful hack to get around fwd with FLStock and fwdControl but neither argument is named
setMethod("fwd", signature(object="FLStock", fishery="fwdControl", control="missing"),
    function(object, fishery, sr=NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE, availability=NULL,maxF=2.0) {
        res <- fwd(object=object, control=fishery, sr = sr, sr.residuals= sr.residuals, sr.residuals.mult = sr.residuals.mult, availability=availability, maxF=maxF)
        return(res)
}) # }}}

# fwd(FLStock, missing, fwdControl) {{{
setMethod("fwd", signature(object="FLStock", fishery="missing", control="fwdControl"),
	function(object, control,
   sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE, availability=NULL,maxF=2.0)
    {

    if (is(sr,"FLBRP")) sr=list(params=params(sr),model=SRModelName(model(sr)))
    ## make sure slots have correct iters 
    if (is(sr,"FLSR")) nDim=dims(params(sr))$iter  else nDim=1
    if (!is.null(sr.residuals)) nDim=max(nDim, dims(sr.residuals)$iter, na.rm=TRUE)  
    if (nDim>1) m(object)=propagate(m(object),nDim)

    object<-CheckNor1(object)

    if (!(units(object@harvest)=="f"))
       stop("harvest slot has to have units of type 'f'")
    if (!validObject(control))
       stop("control not valid")

    yrs <- sort(as.numeric(as.character(unique(control@target[, "year"])))) 

    ## check years
    ## years in control have to be in order
    if (!all(yrs==sort(yrs)))
       stop("years in control not in order")
    ## no gaps in years
    if (length(min(yrs):max(yrs))!=length(unique(yrs)))
       stop("years in control not contiguous")
    ## years in control have to be in FLStock object
    if (!all(ac(yrs) %in% ac(dims(object)$minyear:(dims(object)$maxyear))))
       stop("years in control outside of those in stock object")
    ## Need year+1 in FLStock object
    if (max(yrs) == dims(object)$maxyear){
       endYr<-dims(object)$maxyear+1
       object <- FLCore::window(object, end=dims(object)$maxyear+1)
       #object =qapply(object, FLCore::window, end=dims(object)$maxyear+1)
       }
    else
       endYr<-NULL     
 
    if (is.null(availability)) availability<-sweep(stock.n(object),c(1:4,6),apply(stock.n(object),c(1:4,6), sum),"/")

    sr<-setSR(sr=sr, object=object, yrs=yrs, sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals.mult, availability=availability)
    
    ## check iters in control are '1 or n' and correct if necessary
    control@trgtArray <- chkTrgtArrayIters(object,control@trgtArray,sr)

    ## Season
    if (any(is.na(control@target$season)) & dims(object)$season==1)
       control@target$season<-1
    else if (any(is.na(control@target$season)) & dims(object)$season>1)
       stop("need to specific season in target")

    ## Unit
    if (any(is.na(control@target$unit)) & dims(object)$unit==1)
       control@target$unit<-1
    else if (any(is.na(control@target$unit)) & dims(object)$unit>1)
       stop("need to specific unit in target")

    ## Area
    if (any(is.na(control@target$area)) & dims(object)$area==1)
       control@target$area<-1
    else if (any(is.na(control@target$area)) & dims(object)$area>1)
       stop("need to specific area in target")

    control@target    <- chkTargetQuantity(control@target,object)

    stock.n(object)[1,ac(min(control@target[,"year"]))]<-NA

    ## Availability check
    if (dims(object)$area>1){
       if (is.null(availability)) stop("need to specify availability as areas>1")
       if (any(unlist(dims(availability))[c("age","min","max","unit","season","area","iter")]!=
               unlist(dims(m(object)))[   c("age","min","max","unit","season","area","iter")]))
          stop("dims mismatch in availability")
       if (!all(dimnames(availability)$year %in% (dimnames(m(object))$year)))
          stop("dim year mismatch in availability")
       }

   if (class(maxF)=="numeric" & length(maxF)==1) {
      dmns=dimnames(rec(object))
      dmns$age=1
      maxf=FLQuant(maxF, dimnames=dmns)}
   else maxf=maxF
    
	 x<-.Call("fwd_adolc_FLStock", object, matrixTarget(control@target),
  control@trgtArray, yrs, sr$model, sr$params, sr$residuals, sr$residuals.mult[[1]],
  availability,maxF=maxf, PACKAGE="FLash")

    if (is.numeric(x)) stop(x)

    units(x@harvest)<-"f"

    stock.n(x)[is.na(stock.n(x))]<-0.0
    catch(   x)<-computeCatch(   x)
    landings(x)<-computeLandings(x)
    discards(x)<-computeDiscards(x)
    stock(   x)<-computeStock(   x)

    #name@x<-name(object)
    #desc@x<-desc(object)
    if (!is.null(endYr)) x <- window(x, end=endYr-1)
      #object =qapply(x, FLCore::window, end=endYr-1)

    return(x)}) # }}}

# fwd(FLStock, missing, missing) {{{
setMethod("fwd", signature(object="FLStock", fishery="missing", control="missing"),
    function(object, sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE, availability=NULL,maxF=2.0,...)
    {

    # parse ... for ctrl
    args=list(...)

    # Hack for changing the argument name from 'ctrl' to 'control'
    if("ctrl" %in% names(args)) {
        res <- do.call("fwd", list(object=object, control=args[["ctrl"]], sr=sr, sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals.mult, availability=availability, maxF=maxF))
        return(res)
    }

    if (class(args[[1]])=="FLQuant"){
      control=args[[1]]
      quantity=names(args)[[1]]
    }
  
    control.=apply(control,1:5,mean,na.rm=TRUE)
   
    control.=cbind(quantity=quantity,as.data.frame(control.))

    names(control.)[seq(dim(control.)[2])[names(control.)=="data"]]="val"
    
    control.=fwdControl(control.)
    dmns=dimnames(control.@trgtArray)
    dmns$iter=dimnames(control)$iter
 
    control.@trgtArray=array(c(control),dim=unlist(lapply(dmns,length)),dimnames=dmns)
    control.@trgtArray[,c("min","max"),][]=NA
    
    nits=max(length(dmns$iter),  length(dimnames(sr.residuals)$iter))   
         
    if (nits>1 & dims(object)$iter==1){
       stock.n(object)=propagate(stock.n(object),nits)
       if (length(dimnames(sr.residuals)$iter)==1)
         sr.residuals=propagate(sr.residuals,nits)
       }
    res=fwd(object, control=control.,
            sr=sr,sr.residuals = sr.residuals, sr.residuals.mult = sr.residuals.mult, availability=availability,maxF=maxF)  
    
    return(res)}) # }}}

# fwd(FLStock, missing, FLQuants) {{{
setMethod("fwd", signature(object="FLStock", fishery="missing", control="FLQuants"),
    function(object, control, sr =NULL,
      sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
      availability=NULL,maxF=2.0){
    res=FLStocks(mlply(seq(length(control)),
          function(x,object,control,sr,sr.residuals,sr.residuals.mult,availability,maxF) {
            fwd(object,control=control[[x]],quantity=names(control)[x],
                sr=sr,sr.residuals=sr.residuals,sr.residuals.mult=sr.residuals.mult,
                availability=availability,maxF=maxF)},
                    object=object,control=control,
                    sr=sr,sr.residuals=sr.residuals,sr.residuals.mult=sr.residuals.mult,
                    availability=availability,maxF=maxF))                          
                                          
    return(res)}) # }}}

# fwd(FLStock, missing, FLQuant) {{{
setMethod("fwd", signature(object="FLStock", fishery="missing", control="FLQuant"),
    function(object, control,quantity,
               sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
               availability=NULL,maxF=2.0,...)
    {
    control.=apply(control,1:5,mean,na.rm=TRUE)
    control.=cbind(quantity=quantity,as.data.frame(control.,drop=T))
    names(control.)[seq(dim(control.)[2])[names(control.)=="data"]]="val"
    
    control.=fwdControl(control.)
    dmns=dimnames(control.@trgtArray)
    dmns$iter=dimnames(control)$iter
 
    control.@trgtArray=array(c(control.@trgtArray),dim=unlist(lapply(dmns,length)),dimnames=dmns)
    
    res=fwd(object, control=control.,
            sr=sr, sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals.mult, availability=availability, maxF=maxF)  
    
    return(res)}) # }}}

# fwd(FLStock, ANY, missing, ...) {{{

setMethod("fwd", signature(object="FLStock", fishery="ANY", control="missing"),
  function(object, fishery=missing, ..., sr=NULL,
    sr.residuals=FLQuant(1, dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE) {
    
    # PARSE ...
    args <- list(...)
    .qlevels <- quantityNms()
    # HACK: deal with f assigned to fishery, might fail
    if(!missing(fishery)) {

      if(!is(fishery, "FLQuant"))
        stop("targets can only be of class FLQuant if no fwdControl is provided")
      narg <- names(sys.calls()[[length(sys.calls())-1]])
      narg <- narg[!narg %in% c("", "object", "sr", "sr.residuals", "sr.residuals.mult",
        grep("^[f].*", .qlevels, value=TRUE, invert=TRUE))]
      args[[narg]] <- fishery
    }
    
    # Does ... exist?
    if(length(args) < 1)
      stop("No fwdControl provided and no FLQuant targets given, cannot do anything!")

    # NAMES in qlevels?
    if(any(!names(args) %in% .qlevels))
      stop(paste0("Names of input FLQuant(s) do not match current allowed targets: ",
            paste(.qlevels, collapse=", ")))

    args <- FLQuants(args)

    # COERCE to fwdControl
    control <- as(args, "fwdControl")
    
    return(fwd(object=object, control=control, sr=sr, sr.residuals=sr.residuals,
      sr.residuals.mult=sr.residuals.mult, availability=NULL,maxF=2.0))
  }
) # }}}
