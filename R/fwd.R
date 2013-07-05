# fwd.R
# FLash/R/fwd.R
# Copyright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Finlay Scott, Cefas
# Last Change: Thu Jul 04, 2013 at 11:21 AM +0200
# $Id: fwd.R 1797 2012-12-15 11:11:41Z lauriekell $

## fwd(FLStock)
if (!isGeneric("fwd"))
  setGeneric("fwd", function(object, ctrl, ...)
	  standardGeneric("fwd"))

setMethod("fwd", signature(object="FLStock",ctrl="fwdControl"),
	function(object, ctrl,
   sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
               availability=NULL,maxF=2.0)
    {      
    if (is(sr,"FLBRP")) sr=list(params=params(sr),model=SRModelName(model(sr)))
    ## make sure slots have correct iters 
    if (is(sr,"FLSR")) nDim=dims(params(sr))$iter  else nDim=1
    if (!is.null(sr.residuals)) nDim=max(nDim, dims(sr.residuals)$iter, na.rm=TRUE)  
    if (nDim>1) m(object)=propagate(m(object),nDim)

    object<-CheckNor1(object)

    if (!(units(object@harvest)=="f"))
       stop("harvest slot has to have units of type 'f'")
    if (!validObject(ctrl))
       stop("ctrl not valid")

    yrs <- sort(as.numeric(as.character(unique(ctrl@target[, "year"])))) 

    ## check years
    ## years in ctrl have to be in order
    if (!all(yrs==sort(yrs)))
       stop("years in ctrl not in order")
    ## no gaps in years
    if (length(min(yrs):max(yrs))!=length(unique(yrs)))
       stop("years in ctrl not contiguous")
    ## years in ctrl have to be in FLStock object
    if (!all(ac(yrs) %in% ac(dims(object)$minyear:(dims(object)$maxyear))))
       stop("years in ctrl outside of those in stock object")
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
    
    ## check iters in ctrl are '1 or n' and correct if necessary
    ctrl@trgtArray <- chkTrgtArrayIters(object,ctrl@trgtArray,sr)

    ## Season
    if (any(is.na(ctrl@target$season)) & dims(object)$season==1)
       ctrl@target$season<-1
    else if (any(is.na(ctrl@target$season)) & dims(object)$season>1)
       stop("need to specific season in target")

    ## Unit
    if (any(is.na(ctrl@target$unit)) & dims(object)$unit==1)
       ctrl@target$unit<-1
    else if (any(is.na(ctrl@target$unit)) & dims(object)$unit>1)
       stop("need to specific unit in target")

    ## Area
    if (any(is.na(ctrl@target$area)) & dims(object)$area==1)
       ctrl@target$area<-1
    else if (any(is.na(ctrl@target$area)) & dims(object)$area>1)
       stop("need to specific area in target")

    ctrl@target    <- chkTargetQuantity(ctrl@target,object)

    stock.n(object)[1,ac(min(ctrl@target[,"year"]))]<-NA

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
    
	 x<-.Call("fwd_adolc_FLStock", object, matrixTarget(ctrl@target), ctrl@trgtArray, yrs, sr$model, sr$params, sr$residuals, sr$residuals.mult[[1]], availability,maxF=maxf)

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

    return(x)}) 

# setMethod("fwd", signature(object="FLStock", ctrl="FLQuant"),
#     function(object, ctrl,
#                sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
#                availability=NULL,quantity="f",maxF=2.0)
#     {
#     ctrl.=apply(ctrl,1:5,mean,na.rm=TRUE)
#    
#     ctrl.=cbind(quantity=quantity,as.data.frame(ctrl.,drop=T))
#     names(ctrl.)[seq(dim(ctrl.)[2])[names(ctrl.)=="data"]]="val"
#     
#     ctrl.=fwdControl(ctrl.)
#     dmns=dimnames(ctrl.@trgtArray)
#     dm  =dim(ctrl.@trgtArray)
#     dmns$iter=dimnames(ctrl)$iter
#     dm[3]    =dim(ctrl)[3]
#  
#     ctrl.@trgtArray=array(c(ctrl.@trgtArray),dim=dm,dimnames=dmns)
# 
#     res=fwd(object,ctrl=ctrl.,
#             sr=sr,sr.residuals,sr.residuals.mult=sr.residuals.mult,
#                availability=availability,maxF=maxF)  
#     
#     return(res)})

setMethod("fwd", signature(object="FLStock", ctrl="missing"),
    function(object, ctrl,
               sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
               availability=NULL,maxF=2.0,...)
    {
    args=list(...)
    if (class(args[[1]])=="FLQuant"){
      ctrl=args[[1]]
      quantity=names(args)[[1]]}
  
    ctrl.=apply(ctrl,1:5,mean,na.rm=TRUE)
   
    ctrl.=cbind(quantity=quantity,as.data.frame(ctrl.,drop=T))
    names(ctrl.)[seq(dim(ctrl.)[2])[names(ctrl.)=="data"]]="val"
    
    ctrl.=fwdControl(ctrl.)
    dmns=dimnames(ctrl.@trgtArray)
    dmns$iter=dimnames(ctrl)$iter
 
    ctrl.@trgtArray=array(c(ctrl),dim=unlist(lapply(dmns,length)),dimnames=dmns)
    ctrl.@trgtArray[,c("min","max"),][]=NA
    
    nits=max(length(dmns$iter),  length(dimnames(sr.residuals)$iter))   
         
    if (nits>1 & dims(object)$iter==1){
       stock.n(object)=propagate(stock.n(object),nits)
       if (length(dimnames(sr.residuals)$iter)==1)
         sr.residuals=propagate(sr.residuals,nits)
       }
       
    res=fwd(object,ctrl=ctrl.,
            sr=sr,sr.residuals,sr.residuals.mult=sr.residuals.mult,
               availability=availability,maxF=maxF)  
    
    return(res)})

setMethod("fwd", signature(object="FLStock", ctrl="FLQuants"),
    function(object, ctrl,
               sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
               availability=NULL,maxF=2.0){    
    res=FLStocks(mlply(seq(length(ctrl)),
          function(x,object,ctrl,sr,sr.residuals,sr.residuals.mult,availability,maxF) {
            fwd(object,ctrl=ctrl[[x]],quantity=names(ctrl)[x],
                sr=sr,sr.residuals=sr.residuals,sr.residuals.mult=sr.residuals.mult,
                availability=availability,maxF=maxF)},
                    object=object,ctrl=ctrl,
                    sr=sr,sr.residuals=sr.residuals,sr.residuals.mult=sr.residuals.mult,
                    availability=availability,maxF=maxF))                          
                                          
    return(res)})

setMethod("fwd", signature(object="FLStock", ctrl="FLQuant"),
    function(object, ctrl,quantity,
               sr =NULL, sr.residuals=FLQuant(1,dimnames=dimnames(rec(object))), sr.residuals.mult=TRUE,
               availability=NULL,maxF=2.0,...)
    {    
    ctrl.=apply(ctrl,1:5,mean,na.rm=TRUE)
    ctrl.=cbind(quantity=quantity,as.data.frame(ctrl.,drop=T))
    names(ctrl.)[seq(dim(ctrl.)[2])[names(ctrl.)=="data"]]="val"
    
    ctrl.=fwdControl(ctrl.)
    dmns=dimnames(ctrl.@trgtArray)
    dmns$iter=dimnames(ctrl)$iter
 
    ctrl.@trgtArray=array(c(ctrl.@trgtArray),dim=unlist(lapply(dmns,length)),dimnames=dmns)
    
    res=fwd(object,ctrl=ctrl.,
            sr=sr,sr.residuals,sr.residuals.mult=sr.residuals.mult,
               availability=availability,maxF=maxF)  
    
    return(res)})


#************ FLBiol and fleet methods ***********************************************
# 
# setMethod("fwd", signature(object='FLBiol', fleets='FLCatch'),
#     function(object,fleets,ctrl,
#                 sr           =NULL,
#                 sr.residuals =NULL, sr.residuals.mult=TRUE){
#                 
#     object<-FLBiols(object)
#     if (length(object)>1) stop("Only implemented for 1 FLBiol")
#     fleets <- FLFleets(FLFleet(fleets))
#     fwd(object=object, fleets=fleets, ctrl=ctrl,sr=sr,sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals) 
#     })
# 
# setMethod("fwd", signature(object='FLBiol', fleets='FLFleet'),
#     function(object,fleets,ctrl,
#                 sr           =NULL,
#                 sr.residuals =NULL, sr.residuals.mult=TRUE)
# {                         
#     object<-FLBiols(object)
#     if (length(object)>1) stop("Only implemented for 1 FLBiol")
#     fwd(object = object, fleets = FLFleets(fleets), ctrl=ctrl, sr=sr,sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals) 
# }
# )
# 
# setMethod("fwd", signature(object='FLBiol', fleets='FLFleets'),
#     function(object,fleets,ctrl,
#                 sr           =NULL,
#                 sr.residuals =NULL, sr.residuals.mult=TRUE)
#    {                         
#     object<-FLBiols(object)
#     if (length(object)>1) stop("Only implemented for 1 FLBiol")
#     fwd(object = object, fleets = fleets(fleets), ctrl=ctrl, sr=sr,sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals) 
# 
# })
# 
# setMethod("fwd", signature(object='FLBiols', fleets='FLCatch'),
#     function(object,fleets,ctrl,
#                 sr           =NULL,
#                 sr.residuals =NULL, sr.residuals.mult=TRUE)
# {                         
#     if (length(object)>1) stop("Only implemented for 1 FLBiol")
#     fleets <- FLFleets(FLFleet(fleets))
#     fwd(object=object, fleets=fleets, ctrl=ctrl,sr=sr,sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals) 
# }
# )
# 
# setMethod("fwd", signature(object='FLBiols', fleets='FLFleet'),
#     function(object,fleets,ctrl,
#                 sr           =NULL,
#                 sr.residuals =NULL, sr.residuals.mult=TRUE)
# {                         
#     if (length(object)>1) stop("Only implemented for 1 FLBiol")
#     fwd(object = object, fleets = FLFleets(fleets), ctrl=ctrl, sr=sr,sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals) 
# }
# )
# 
# # The main one
# setMethod("fwd", signature(object='FLBiols', fleets='FLFleets'),
#     function(object,fleets,ctrl,
#                 sr           =NULL,
#                 sr.residuals =NULL, sr.residuals.mult=TRUE)
#    {                         
#     # Turn biol into FLBiols
# #    object<-FLBiols(object)
#     if (length(object)>1) stop("Only implemented for 1 FLBiol")
#     biol  <-CheckNor1(object)
# 
# #    ## Sort out fleets input into FLFleets 
# #    if (is(fleets,"FLCatch")){
# #        fleets<-FLFleets(FLFleet())
# #        fleets[[1]]@metiers[[1]]<-as(fleets,"FLMetier")}
# #   if (is(fleets,"FLFleet"))
# #      fleets<-FLFleets(fleets)
# #   #if (!is(fleets,"FLFleets"))
# #   #   stop("'fleets' not of type FLFleet(s)")     
#     if (length(fleets)>1) stop("Only implemented for 1 FLFleet")
#     fleets<-CheckNor1(fleets)
# 
#     yrs<-as.numeric(sort(unique(ctrl@target[,"year"])))
#     sr<-setSR(sr=sr, object=biol[[1]], yrs=yrs, sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals.mult)          
# 
#     if (!is(ctrl,"fwdControl")) stop("ctrl not a valid 'fwdControl' object")
# 
#     # check iters in residuals, biol and ctrl@trgtArray are 1 or n and correct trgtArray if necessary
# ctrl@trgtArray <- chkTrgtArrayIters(object,ctrl@trgtArray,sr)
# #    its<-(unique(c(length(dimnames(ctrl@trgtArray)$iters), dims(object[[1]])$iters, length(dimnames(sr.residuals)$iter))))
# #    if (length(its)>2 | (length(its)>1 & its[1]!=1)) stop("Iters not 1 or n")
# #    if (length(its)==2 & dimnames(ctrl@trgtArray)$iter == 1){
# #          dmns<-dimnames(ctrl@trgtArray)
# #          dmns$iters<-1:its[2]
# #          ctrl@trgtArray<-array(ctrl@trgtArray,dim=unlist(lapply(dmns,length)),dimnames=dmns)}
# 
# ctrl@target <- chkTargetQuantity(ctrl@target)
# #     if (!is(ctrl@target[,"quantity"],"factor"))
# #        ctrl@target[,"quantity"]<-factor(ctrl@target[,"quantity"],quantityNms())
# #     if (!all(as.character(ctrl@target[,"quantity"]) %in% quantityNms()))
# #         stop("invalid quantity in control target")
# 
#    if(!(length(slot(ctrl, 'effort'))>0) & length(fleets)==1){
#       ctrl@effArray  <-ctrl@trgtArray
#       ctrl@effArray[]<-NA
#       
#       ctrl@effort <-data.frame(year=ctrl@target[,"year"],min=NA,val=NA,max=NA,fleet=1,metier=1,
#                                rel.year=NA,rel.fleet=NA,rel.metier=NA,rel.bound=NA)
# 
#       }
#    if (!validObject(ctrl))
#        stop("ctrl not a valid 'fwdControl' object")
# 
# 
#    ##check dims of fleet and biol
#    dms<-chkDms(biol[[1]],fleets[[1]],yrs)
#    #if (class(dms)=="character") 
#    #   stop(dms)
#       
#    ## Check iters are 1 or n
#    iter<-unique(c(unlist(lapply(fleets, function(x) dims(x)$iters)),
#                   unlist(lapply(biol, function(x) dims(x)$iter))))
#    if ((length(iter) >2) | (length(iter)==2 & !(1 %in% iter)))
#       stop("Iter needs to be '1 or n'")
# 
#    if (!all(unlist(dims(biol[[1]])[  c("minyear","maxyear","season","area","iter")])==
#             unlist(dims(fleets[[1]])[c("minyear","maxyear","season","area","iter")])))
#       stop("year, season, area and iter have to match in FLBiol & FLFleet")
#    
#    ## needed cos m plays a big role in the C++ code
#    for (i in 1:length(biol))
#      if (dims(m(biol[[i]]))$iter==1 && dims(biol[[i]])$iter>dims(m(biol[[i]]))$iter)
#         m(biol[[i]])<-propagate(m(biol[[i]]),iter=dims(biol[[i]])$iter)
#    
#    ## FLSR    
# #   if (is(sr,"FLSR")){
# #      sr.model <-NULL
# #      sr.params<-NULL}
# #   else if (is(sr,"list") & all(c("model","params") %in% names(sr))){  
# #      sr.model <-sr$model
# #      sr.params<-sr$params
# #      sr       <-NULL}   
# #   else 
# #      stop("sr has to either be an FLSR object or a list with items 'model' & 'params'")
# #      
# #   if (!is.null(sr.residuals))
# #      if(dims(sr.residuals)$iter != dms["iter"])
# #         stop("iter in object and sr.residuals don't match")
#     #sr<-setSR(sr=sr, object = object, yrs=yrs, sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals.mult)          
# #browser()
# #    sr<-setSR(sr=sr, object=biol[[1]], yrs=yrs, sr.residuals=sr.residuals, sr.residuals.mult=sr.residuals.mult)          
#    #if (is.character(sr)) 
#    #   stop(sr)
# 
#    dms<-dms[c("minyear","maxyear","unit","season","area","iter")]
# 
#    for (i in 1:length(biol))
#       if ((yrs[length(yrs)] == dims(biol[[i]])$maxyear)+1  )
#          biol[[i]]<-window(biol[[i]],start=biol[[i]]@range["minyear"],end=yrs[length(yrs)]+1)
# 
#    res<-.Call("fwd_adolc_FLBiols", biol, fleets,
#                                    matrixTarget(ctrl@target), ctrl@trgtArray,
#                                    matrixEffort(ctrl@effort), ctrl@effArray,
#                                    yrs, dms, sr)
# 
#    names(res)<-c("landings.n","discards.n","effort","n","f","catch.n")
# 
#    return(res)
#    })
# 
# #***************************************************      
# 
# # This needs checking...
# #setMethod("fwd", signature(object="character"),
# #       function(object,fltNm,ctrl,
# #                sr           =NULL,
# #                sr.residuals =NULL, sr.residuals.mult=TRUE)
# #   {
# #
# #    biol<-get(object,envir = as.environment(parent.frame()))
# #    if (!is(biol,"FLBiol")) stop("named object needs to be of type FLBiol")
# #
# #    if (is(fltNm,"character"))
# #        fleet<-get(fltNm,envir = as.environment(parent.frame()))
# #    else stop("fltNm not character")
# #    if (!is(fleet,"FLFLeet")) stop("named object needs to be of type FLFleet")
# #   
# #   res<-fwd(biol,fleets,ctrl,sr,sr.residuals,sr.residuals.mult)
# #
# #   n(biol              )[,dimnames(res$n         )$year]<-res$n
# #   landings.n(fleet,1,1)[,dimnames(res$landings.n)$year]<-res$landings.n
# #   discards.n(fleet,1,1)[,dimnames(res$discards.n)$year]<-res$discards.n
# #   landings(  fleet,1,1)[,dimnames(res$landings.n)$year]<-apply(res$landings.n*landings.wt(fleet,1,1)[,dimnames(res$landings.n)$year],2:6,sum)
# #   discards(  fleet,1,1)[,dimnames(res$discards.n)$year]<-apply(res$discards.n*discards.wt(fleet,1,1)[,dimnames(res$discards.n)$year],2:6,sum)
# #   effort(fleet        )[,dimnames(res$effort    )$year]<-res$effort
# #
# #   name(    x)<-name(object)
# #   desc(    x)<-desc(object)
# #
# #   assign(object,biol, envir=as.environment(parent.frame()))
# #   assign(fltNm, fleet,envir=as.environment(parent.frame()))
# #
# #   invisible(res)
# #   })
# #
