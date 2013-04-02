##Pseudocohort analysis
#Tim Earl, August 2010
#Adapted by José De Oliveira to work with FLR, and to account for dynamic plusgroup

#Reference:
# http://www.fao.org/docrep/w5449e/w5449e00.htm
# section 5.3

##Data

#ages  = 2:8
#years = 2000:2001
#catch = matrix(c(882,2720,1899,754,268,81,24,
#                  882,2720,1899,754,268,81,24
#                 ), ncol=length(ages), byrow=TRUE)
#dimnames(catch) = list(year=years,age=ages)
#m     = 0.2
#fages = 3:4
#iter = 20


#Catch should be a matrix with columns for each age class, and rows for each year.
#column (?and row?) names will be preserved by the function
#m is a value for natural mortality (or possibly a vector for each age)

setGeneric('pseudoCohort', function(object, ...)
  standardGeneric('pseudoCohort'))

setMethod('pseudoCohort', signature(object='FLStock'),
  function(object,nits=10,tol=10e-8,...){

    pGrp  =range(object)["plusgroup"]
    if (is.na(pGrp)) stop("Plus group not specified")

    object=setPlusGroup(object,pGrp)
    object=trim(object,age=range(object)["min"]    :range(object)["max"],
                       year=range(object)["minyear"]:range(object)["maxyear"])

    avcatch = apply(catch.n(object),c(1,6),mean,na.rm=TRUE)
    avm     = apply(m(object),c(1,6),mean)
    n       = avcatch*0
    f       = n
    fterm   = f[ac(pGrp)]
    fterm[] = 0.5

    h       = exp(avm/2)

    #Estimate numbers in oldest age group
    j<-0
    while (j<nits & fdif>tol){
      n[pGrp,]  =avcatch[pGrp,]*(fterm+avm[pGrp,])/(fterm*(1-exp(-fterm-avm[pGrp,])))
      n[pGrp-1,]=((n[pGrp,]-(n[pGrp,]/h[pGrp,]-avcatch[pGrp,])/h[pGrp,])*h[pGrp-1]+avcatch[pGrp-1,])*h[pGrp-1]

      for (i in rev(seq(along=iter(n,1)))[-c(1:2)])
         n[i] = (h[i]*n[i+1] + avcatch[i])*h[i]

      f[pGrp]            = fterm
      f[pGrp-1]          = log(n[pGrp-1])-log(n[pGrp]*(1-exp(-f[pGrp]-avm[pGrp])))-avm[pGrp-1] #applies if plusgroup
      f[-c(pGrp-1,pGrp)] = log(n[-c(pGrp-1,pGrp)])-log(n[-c(1,pGrp)])-avm[-c(pGrp-1,pGrp)]
      fdif               = sum(fterm-apply(f[pGrp-c(3:1),],6,mean))
      fterm              = apply(f[pGrp-c(3:1),],6,mean)

      j<-j+1}

    return(list(n=n,harvest=f,fdif=fdif,nIters=j))})

#pseuco(ple4[,1:10])
 

