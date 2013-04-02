# fwdStock<-function(obj,fbar,sr,sr.residuals=NULL,distribution=NULL){
#    if (is(sr,"FLBRP")) sr=list(params=params(sr),model=SRModelName(model(sr)))
#  
#    if (is.null(sr.residuals))
#       sr.residuals<-FLQuant(1,dimnames=dimnames(fbar))
#     
#    ## make sure slots have correct iters      
#    nDim=max(dims(sr)$iter, dims(sr.residuals)$iter, na.rm=TRUE)  
#    if (nDim>1) m(obj)=propagate(m(obj),nDim)
#    obj<-CheckNor1(obj)
# 
#  
#    #### check dims
#    if (!all(dimnames(fbar)$area   %in% dimnames(m(obj))$area))
#       stop("Areas in FBar and obj don't match")
#    if (!all(dimnames(fbar)$season %in% dimnames(m(obj))$season))
#       stop("Seasons in FBar and obj don't match")
#    if (!all(dimnames(fbar)$unit   %in% dimnames(m(obj))$unit))
#       stop("Units in FBar and obj don't match")
#    if (!all(dimnames(fbar)$year   %in% dimnames(m(obj))$year))
#       stop("Years in FBar and obj don't match")
#       
#    yrs <-dimnames(fbar)$year
# 
#    if (is.null(distribution)){
#       distribution<-sweep(stock.n(obj)[,yrs],c(1:4,6),apply(stock.n(obj)[,yrs],c(1:4,6),"sum"),"/")
#       if (dims(fbar)$maxyear==dims(obj)$maxyear){
#          distribution<-window(distribution,end=dims(obj)$maxyear+1)
#          distribution[,ac(dims(obj)$maxyear+1)]<-distribution[,ac(dims(obj)$maxyear)]
#          }}
#          
#    extendFlag<-FALSE
#    if (dims(fbar)$maxyear==dims(obj)$maxyear){
#       extendFlag<-TRUE
#       obj<-window(obj,end=dims(obj)$maxyear+1)
#       }
# 
#    ##Fishing Mortality
#    fbRng             <-ac(range(obj,"minfbar"):range(obj,"maxfbar"))
#    scaling           <-fbar/apply(harvest(obj)[fbRng,yrs],2:6,mean)
#    harvest(obj)[,yrs]<-sweep(harvest(obj)[,yrs]@.Data,2:6,scaling@.Data,"*")
# 
#    stock.n(obj)[1,yrs]<-0
#    for (iY in yrs){
#       #cat("Year\t",iY,"\n")
#       for (iU in dimnames(m(obj))$unit){
#          #cat("Unit\t",iU,"\n")
#          for (iS in dimnames(m(obj))$season){
#            #cat("Season\t",iS,"\n")
# 
#            ## Recruitment
#            if (iS %in% dimnames(params(sr))$season){
#              recY<-ac(as.integer(iY)-range(obj,"min"))
#              ssb.<-apply(ssb(obj)[,recY,iU,iS],c(3,5:6),sum)
# 
#              stock.n(obj)[1,iY,iU,iS]<-distribution[1,iY,iU,iS] #predict(sr,ssb=ssb.)*sr.residuals[,iY,iU,iS]
#              }
# 
#            z                         <-m(obj)[,iY,iU,iS]+harvest(obj)[,iY,iU,iS]
#            ## Catches
#            catch.n(   obj)[,iY,iU,iS]<-stock.n(obj)[,iY,iU,iS]*harvest(obj)[,iY,iU,iS]/(z)*(1-exp(-z))
#            landings.n(obj)[,iY,iU,iS]<-catch.n(obj)[,iY,iU,iS]*landings.n(obj)[,iY,iU,iS]/(landings.n(obj)[,iY,iU,iS]+discards.n(obj)[,iY,iU,iS])
#            discards.n(obj)[,iY,iU,iS]<-catch.n(obj)[,iY,iU,iS]-discards.n(obj)[,iY,iU,iS]
# 
#            ## next season?
#            if (iS!=dimnames(m(obj))$season[dims(obj)$season]){
#              stock.n(obj)[,iY,iU,ac(as.integer(iS)+1)]<-stock.n(obj)[,iY,iU,iS]*exp(-z)
#            ##next year
#            }else if (as.integer(iY)<dims(obj)$maxyear){
#              stock.n(obj)[-1,ac(as.integer(iY)+1),iU,1]<-stock.n(obj)[-dims(obj)$max,iY,iU,iS]*exp(-z[-dims(obj)$max])
# 
#            ## plusgroup
#            if (!is.na(range(obj,"plusgroup")))
#              stock.n(obj)[dims(obj)$max,ac(as.integer(iY)+1),iU,1]<-stock.n(obj)[dims(obj)$max,ac(as.integer(iY)+1),iU,1]+stock.n(obj)[dims(obj)$max,iY,iU,iS]*exp(-z[ dims(obj)$max])
#            }
#          }
#        }
#      }
# 
#    if (extendFlag)
#       obj<-window(obj,end=dims(obj)$maxyear-1)
# 
#    return(obj)
#    }

