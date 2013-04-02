# Rewrite fwdSetSRs to be methods rather than functions
# Still outputs the same thing (list of 4 elements, each of length number of SRs)

# Two dispatch methods.  One where SR is FLSR, the other SR is missing.
# If missing need to specify model and params.  Model can be text string, formula, or function (the sr function used by FLSR).
# Gives a total of four different methods of making the SR


# Will eventually be expanded so that the multiple SRs can be returned.  For the moment s a single SR is returned

if(!isGeneric('setSR'))
    setGeneric('setSR', function(sr, ...) standardGeneric('setSR'))

setMethod('setSR', signature(sr='FLSR'),
    function(sr,object,yrs,sr.residuals=NULL,sr.residuals.mult=TRUE,availability=NULL) {

    # Strip out sr model and params then call that method
    setSR(list(model = model(sr), params = params(sr)), object = object, yrs = yrs, sr.residuals = sr.residuals, sr.residuals.mult = sr.residuals.mult, availability=availability)
    }
)

setMethod('setSR', signature(sr='list'),
    function(sr,object,yrs,sr.residuals=NULL,sr.residuals.mult=TRUE,availability=NULL)
{
# ****** Check arguments are all present and of correct type *******
    if(!(is(sr,"list") & all(c("model","params") %in% names(sr))))
        stop("sr has to be a list with items 'model' & 'params'")
    model <- sr$model
    params <- sr$params

    if (!is(params,'FLPar') & !is(params,'FLQuant')) stop("params must be of type FLPar or FLQuant")
    # Sort out what type of object sr is; can be character string, formula or function
    if (!is(model,'formula') & !is(model,'character') & !is(model,'function'))  stop("model m ust be of type formula, character or function")
    if (!is(object,'FLStock') & !is(object,'FLBiol')) stop("object must be an FLStock or FLBiol")
    if (!is(yrs,'numeric') & !is(yrs,'character')) stop("yrs must be a numeric or character")
     

#****** Check yrs are in object year range *********
    if(is(yrs,'numeric')) yrs <- as(yrs,'character')
    if (!all(yrs  %in% dimnames(object@m)$year))
       stop("yrs range not found in object")

#****** Check and sort out residuals **********
    # Make an FLQuant of correct dimension with residuals
    if (!is(sr.residuals, "FLQuant") & !is.null(sr.residuals))
       stop("sr.residuals must be an FLQuant")
    dmns     <-dimnames(rec(object))
    # Add extra year - needed to smooth over bumps in fwd
    yrs   <-c(yrs,as.integer(yrs[length(yrs)])+1)
    dmns$year<-yrs
    # If no sr.residuals then thay are NULL and just a quant of 1s is passed
    residuals<-FLQuant(1,dimnames=dmns)
    resYrs<-yrs[yrs %in% dimnames(sr.residuals)$year]

    if(!is(sr.residuals,"NULL") & !any(dimnames(sr.residuals)$year %in% yrs[1:(length(yrs)-1)])) warning("Year range of residuals is not within yrs range")
    residuals[,resYrs]<-iter(sr.residuals[,resYrs], dimnames(residuals)$iter)

#****** Check and force iterations **********
    # Iters of stock numbers and residuals should be the same or one of them should be 1
    if (is(object,"FLStock") && (!all(dimnames(sr.residuals)$iter  %in% dimnames(object@stock.n)$iter) && !(dimnames(sr.residuals)$iter=="1" || dimnames(object@stock.n)$iter=="1")))
       stop("iters in sr.residuals do not those in object")
    if (is(object,"FLBiol") && (!all(dimnames(sr.residuals)$iter  %in% dimnames(object@n)$iter) && !(dimnames(sr.residuals)$iter=="1" || dimnames(object@n)$iter=="1")))
       stop("iters in sr.residuals do not those in object")

#****** Get model type ******************
    # make array of model type
    t.       <- vector(mode="numeric",length(yrs))
    names(t.)<-yrs
    if (is(model,"formula")) model<-SRModelName(model)
    t.[]     <-SRchar2code(model)
    model <-t.

#****** Check and force parameters for all years ********
    # Turn the FLPar or Quant into a Quant of right dimensions
    dmns <-list(params=SRParams(SRcode2char(ac(model[1]))),
                year  =yrs,
                unit  =dimnames(m(object))$unit,
                season=dimnames(m(object))$season,
                area  =dimnames(m(object))$area,
                iter  =dimnames(params)$iter)
    dmns.<-dmns
    dmns.[[2]]<-dmns[[2]][-length(dmns[[2]])]# dmns of original years

    # Trim off extra params if coming from FLSR (e.g. Ricker has extra param sigma)    
    params <- params[dmns$params]
    
    # Reorder if necessary
    if (length(length(dimnames(params)$params))>1)
      if (all(dimnames(sr)$params[1:2]==c("b","a")))
        params<-params[c("a","b")]

    if (!any(dimnames(params) %in% dmns.))
      stop("Dims for sr.params illegal")

    params<-validSRPar(object, sr=params, yrs=yrs, availability=availability)
    params<-FLQuant(params)

#****** Cobble together into output format **********
    # At the moment each element is only length 1.  Eventually, multiple SRs will be possible.
    res <- list(model=list(model),
        params=FLQuants(params),
        residuals=FLQuants(residuals),
        residuals.mult=list(sr.residuals.mult))

    return(res)})

SRchar2code<-function(strCode){
   res<-as.integer(switch(strCode, "mean"            = 1,
                                   "geomean"         = 1,
                                   "bevholt"         = 2,
                                   "ricker"          = 3,
                                   "segreg"          = 4,
                                   "shepherd"        = 5,
                                   "bevholt.d"       = 21,
                                   "bevholt.c.a"     = 22,
                                   "bevholt.c.b"     = 23,
                                   "bevholt.sv"      = 24,
                                   "bevholt.ndc"     = 25,
                                   "bevholt.ar1"     = 26,
                                   "ricker.d"        = 31,
                                   "ricker.c.a"      = 32,
                                   "ricker.c.b"      = 33,
                                   "ricker.sv"       = 34,
                                   "ricker.ndc"      = 35,
                                   "ricker.ar1"      = 36,
                                   default           = 0))

   return(res)
   }

SRcode2char<-function(strCode){
   res<-switch(strCode,  "1"    = "geomean",       
                         "2"    = "bevholt",       
                         "3"    = "ricker",        
                         "4"    = "segreg",        
                         "5"    = "shepherd",      
                         "21"   = "bevholt.d",     
                         "22"   = "bevholt.c.a",   
                         "23"   = "bevholt.c.b",   
                         "24"   = "bevholt.sv",    
                         "25"   = "bevholt.ndc",   
                         "26"   = "bevholt.ar1",   
                         "31"   = "ricker.d",      
                         "32"   = "ricker.c.a",    
                         "33"   = "ricker.c.b",    
                         "34"   = "ricker.sv",    
                         "35"   = "ricker.ndc",    
                         "36"   = "ricker.ar1",    
                         "0"    = default)         

   return(res)
   }

SRParams<-function(strCode){
   res<-as.integer(switch(strCode, "mean"            = 1,
                                   "geomean"         = 1,
                                   "bevholt"         = 2,
                                   "ricker"          = 2,
                                   "segreg"          = 2,
                                   "shepherd"        = 3,
                                   "bevholt.d"       = 3,
                                   "bevholt.c.a"     = 3,
                                   "bevholt.c.b"     = 3,
                                   "bevholt.sv"      = 2,
                                   "bevholt.ndc"     = 3,
                                   "bevholt.ar1"     = 2,
                                   "ricker.d"        = 3,
                                   "ricker.c.a"      = 3,
                                   "ricker.c.b"      = 3,
                                   "ricker.sv"       = 2,
                                   "ricker.ndc"      = 3,
                                   "ricker.ar1"      = 3,
                                   default           = 0))
                                   
                          
   return(c("a","b","c")[1:res])
   }   
