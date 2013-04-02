calcF<-function(m,catch,n)
   {
   its<-unique(c(dims(    m)$iter, 
                 dims(catch)$iter, 
                 dims(    n)$iter))

   if (dims(m)$iter==1) 
      m<-propagate(m,max(its))   
   if (dims(catch)$iter==1) 
      catch<-propagate(catch,max(its))
   if (dims(n)$iter==1) 
      n<-propagate(n,max(its))
   
   if (length(its)> 2) stop("iter mismatch")
   if (length(its)==2 & !any((1 %in% its))) stop("iters have to be 1 or n")
      
   res       <-.Call("CalcF",m,catch,n)
   units(res)<-"f"

   return(res)
   }

setGeneric('computeHarvest', function(object, ...)
		standardGeneric('computeHarvest'))

setMethod('computeHarvest', signature(object='FLStock'),
  function(object, catch)
     {
     if (names(dims(m(object)))[1]!="age") warning("quant dim not age, harvest only valid for age")
     
     res <-calcF(m(object),catch.n(object),stock.n(object))

     return(res)
     })     
     
setMethod('harvest', signature(object='FLBiol', catch='FLQuant'),
  function(object, catch)
     {
     if (names(dims(catch)[1])!="age") warning("quant dim not age, harvest only valid for age")

     res <-calcF(m(object),catch, n(object))

     return(res)
     })

setMethod('harvest', signature(object='FLBiol', catch='FLCatch'),
  function(object, catch)
     {
     if (names(dims(m(object)))[1]!="age") warning("quant dim not age, harvest only valid for age")

     res <-calcF(m(object),catch.n(catch),n(object))

     return(res)
     })

setMethod('harvest', signature(object='FLBiol', catch='FLMetier'),
  function(object, catch, spp)
     {
     if (names(dims(m(object)))[1]!="age") warning("quant dim not age, harvest only valid for age")

     res <-calcF(m(object),catch.n(catch)[[spp]],n(object))

     return(res)
     })

setMethod('harvest', signature(object='FLBiol', catch='FLMetiers'),
  function(object,catch,spp)
     {
     if (names(dims(m(object)))[1]!="age") warning("quant dim not age, harvest only valid for age")

     res <-calcF(m(object),catch.n(catch,spp),n(object))

     return(res)
     })

setMethod('harvest', signature(object='FLBiol', catch='FLFleet'),
  function(object, catch, spp, mtr)
     {
     if (names(dims(m(object)))[1]!="age") warning("quant dim not age, harvest only valid for age")

     res <-calcF(m(object),catch.n(catch, mtr, spp),n(object))

     return(res)
     })

