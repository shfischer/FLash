##Truncated VPA analysis
#José De Oliveira


setGeneric('truncVPA', function(object, ...)
  standardGeneric('truncVPA'))

setMethod('truncVPA', signature(object='FLStock'),
  function(object,nits=10,...){
    pGrp  =range(object)["plusgroup"]
    if (is.na(pGrp)) stop("Plus group not specified")

    yrs  =ac(dims(object)$minyear:dims(object)$maxyear)
    ages =ac(dims(object)$min    :dims(object)$max)
    fages=range(object)["minfbar"]:range(object)["maxfbar"]

    object=setPlusGroup(object,pGrp)
    object=trim(object,age=range(object)["min"]    :range(object)["max"],
                       year=range(object)["minyear"]:range(object)["maxyear"])

    catch    = catch.n(object)
    m        = m(object)
    n        = catch*0
    f        = n
    A        = ages
    A1       = ages[length(ages)-1]
    aY       = ages[1:(length(ages)-2)]
    Y        = yrs[3]
    ftermA1  = f[A1]
    ftermY   = f[aY,Y]
    ftermA1[]= 0.5
    ftermY[] = 0.5
    h        = exp(m/2)

    j<-0
    while (j<nits){
      n[A1]  =catch[A1]*(ftermA1+m[A1])/(ftermA1*(1-exp(-ftermA1-m[A1])))
      n[aY,Y]=catch[aY,Y]*(ftermY+m[aY,Y])/(ftermY*(1-exp(-ftermY-m[aY,Y])))
      for (i in rev(1:(length(ages)-2)))
         n[i,yrs[-3]] = (h[i,yrs[-3]]*n[i+1,yrs[-1]] + catch[i,yrs[-3]])*h[i,yrs[-3]]

      f[A1]            = ftermA1
      f[aY,Y]          = ftermY
      f[aY,yrs[-3]]    =log(n[aY,yrs[-3]])-log(n[ages[2:(length(ages)-1)],yrs[-1]])-m[aY,yrs[-3]]
      fdif             = sum(ftermY-apply(f[aY,c(1,2)],c(1,6),mean)) + sum(ftermA1[,c(1,2)]-apply(f[rev(aY)[c(1,2)],c(1,2)],c(2,6),mean)) + sum(ftermA1[,3]-apply(apply(f[aY,c(1,2)],c(1,6),mean)[rev(aY)[c(1,2)]],6,mean))
      ftermY[]         = apply(f[aY,c(1,2)],c(1,6),mean)
      ftermA1[,c(1,2)] = apply(f[rev(aY)[c(1,2)],c(1,2)],c(2,6),mean)
      ftermA1[,3]      = apply(ftermY[rev(aY)[c(1,2)]],6,mean)}

    f[A]       =f[A1]
    n[A]       =catch[A]*(f[A]+m[A])/(f[A]*(1-exp(-f[A]-m[A])))
    N          =n[,Y]
    N[ages[-1]]=n[ages[-length(ages)],Y]*exp(-f[ages[-length(ages)],Y]-m[ages[-length(ages)],Y])
    N[A]       =N[A] + n[A,Y]*exp(-f[A,Y]-m[A,Y])
    N[1]       =exp(apply(log(n[1,]),6,sum)/3)
    F          =apply(f,c(1,6),mean)
    fbar       = apply(F[ac(fages)],6,mean)
    sel        = sweep(F,6,fbar,"/")

    return(list(n=n,harvest=f,nits=j))})
