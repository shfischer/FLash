#define ADOLC TAPELESS
#include <adouble.h>
#include "fwdFLStock.h"
#include "fwd.h"
#include "float.h"

//#include <iostream>
//using namespace std;

fwdStk::fwdStk(void)
   {
   indepLastYr=true;
   }      

fwdStk::~fwdStk(void)
   {
   ;
   }      

double fwdStk::getVal(FLRConst_Target quantity,  int iyr, int iunit, int iseason, int iarea, int iter)
   {
   double value=0.0;

   switch(quantity){
       case FLRConst_F:
         value = stk.Fbar(           iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_Z:
         value = stk.Zbar(           iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_FLandings:
         value = stk.FbarLandings(   iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_FDiscards:
         value = stk.FbarDiscards(   iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_SSB:
         value = SSB(                iyr,iunit, iseason, iarea, iter);
         break;
       case FLRConst_Biomass:
         value = stk.computeStock(   iyr+1,iunit, iseason, iarea, iter);
         break;
       case FLRConst_MnSz:
         value = stk.computeMnSz(    iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_Catch:
         value = stk.computeCatch(   iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_Landings:
         value = stk.computeLandings(iyr,  iunit, iseason, iarea, iter);
         break;
       case FLRConst_Discards:
         value = stk.computeDiscards(iyr,  iunit, iseason, iarea, iter);
         break;
       }

   return value;
   }

void fwdStk::project(adouble *x, adouble *func, double *Trgt, int iTrgt, int nrow, double *Ary, int iter)
   {
   int relYr, relSn, relArea, relUnt;

   // Time step
   iTrgt--;

   int iyr   = (int)(Trgt)[iTrgt],
       iSn   =  relSn  = __max((int)(Trgt)[iTrgt+fwdTargetPos_season*nrow],1),
       iunit = relUnt  = __max((int)(Trgt)[iTrgt+fwdTargetPos_unit*nrow],  1),
       iarea = relArea = __max((int)(Trgt)[iTrgt+fwdTargetPos_area*nrow],  1);
       relYr = stk.minyr-1;

   if (!R_IsNA((Trgt)[iTrgt+fwdTargetPos_relyear*nrow]))
       relYr = (int)(Trgt)[iTrgt+fwdTargetPos_relyear*nrow];
   if (!R_IsNA((Trgt)[iTrgt+fwdTargetPos_relyear*nrow]))
       relSn = (int)(Trgt)[iTrgt+fwdTargetPos_relseason*nrow];
   if (!R_IsNA((Trgt)[iTrgt+fwdTargetPos_relyear*nrow]))
       relArea = (int)(Trgt)[iTrgt+fwdTargetPos_relarea*nrow];
   if (!R_IsNA((Trgt)[iTrgt+fwdTargetPos_relyear*nrow]))
       relUnt = (int)(Trgt)[iTrgt+fwdTargetPos_relunit*nrow];
   
   //Population model
   FLQuant_adolc ad_f(stk.harvest, iyr,   iyr,   iter);
   FLQuant_adolc ad_n(stk.stock_n, iyr+1, iyr+1, iter);

   // F
   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      ad_f(iage,iyr,iunit,iSn,iarea,iter) = stk.harvest(iage, iyr, iunit, iSn, iarea, iter)*x[0];

   // recruits
   int SSB_yr = __min(__max(iyr-stk.minquant,stk.minyr),stk.maxyr);
   if      (iSn==1) stk.stock_n(stk.minquant,iyr,iunit,iSn,iarea,iter) = SR.recruits(1,stk.SSB(SSB_yr,iunit,    iSn,iarea,iter),iyr,iunit,iSn,      iarea,iter);
   else if (iSn >1) stk.stock_n(stk.minquant,iyr,iunit,iSn,iarea,iter) = stk.stock_n(stk.minquant,iyr,iunit,iSn-1,iarea,iter)*exp(-stk.harvest(stk.minquant,iyr,iunit,iSn-1,iarea,iter)-stk.m(stk.minquant,iyr,iunit,iSn-1,iarea,iter))+
	                                                                     SR.recruits(1,stk.SSB(SSB_yr,iunit,iSn,iarea,iter),iyr,iunit,iSn,  iarea,iter);

   //Project
   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      {
      //numbers-at-age next season/year
      if (iSn<stk.nseasons)
         ad_n(iage,iyr,iunit,iSn+1,iarea,iter) = stk.stock_n(iage,iyr,iunit,iSn,iarea,iter)*
                                                 exp(-ad_f(  iage,iyr,iunit,iSn,iarea,iter)
                                                     -stk.m( iage,iyr,iunit,iSn,iarea,iter));
      else {
         if (iage < stk.maxquant)
            ad_n(iage+1,iyr+1,iunit,1,iarea,iter) = stk.stock_n(iage,iyr,iunit,iSn,iarea,iter)*
                                                      exp(-ad_f(  iage,iyr,iunit,iSn,iarea,iter)
                                                          -stk.m( iage,iyr,iunit,iSn,iarea,iter));

         if (iage == stk.plusgrp)
            ad_n(iage,  iyr+1,iunit,1,iarea,iter) += stk.stock_n(stk.maxquant,iyr,iunit,iSn,iarea,iter)*
                                                       exp(-ad_f(  stk.maxquant,iyr,iunit,iSn,iarea,iter)
                                                           -stk.m( stk.maxquant,iyr,iunit,iSn,iarea,iter));
         }
      }

    // Redistribute between areas if more than 1 area
    if (stk.nareas>1){
       if (iSn<stk.nseasons){        //within year
          for (int iage=stk.minquant; iage<=stk.maxquant; iage++){
             adouble sum=0.0;
             for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                sum+=ad_n(iage,iyr,iunit,iSn+1,jarea,iter);
          
		     for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                ad_n(iage,iyr,iunit,iSn+1,jarea,iter)=sum*avail(iage,iyr,iunit,iSn+1,jarea,iter);
             }	      
	      }
	   else if (avail.maxyr()>iyr) { //next year only if avial provided
          for (int iage=stk.minquant; iage<=stk.maxquant; iage++){
             adouble sum=0.0;
             for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                sum+=ad_n(iage,iyr+1,iunit,1,jarea,iter);
          
		     for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                ad_n(iage,iyr+1,iunit,1,jarea,iter)=sum*avail(iage,iyr+1,iunit,1,jarea,iter);
             }	      	      
	      }
	   }

	   
   //-------------------- Target Stuff ----------------------//
   //min & max bounds should only occur if a target calculated in a previous step for that year
   // target value relative to reference year
   //Targets
   double min = (Ary)[(iTrgt+0*nrow+3*nrow*(iter-1))];
   double val = (Ary)[(iTrgt+1*nrow+3*nrow*(iter-1))];
   double max = (Ary)[(iTrgt+2*nrow+3*nrow*(iter-1))];

   FLRConst_Target quantity = (FLRConst_Target)(int)(Trgt)[iTrgt + fwdTargetPos_quantity*nrow];

   if (relYr>=stk.minyr)
      {
      double RelVal=getVal(quantity, relYr, relUnt=1, relSn=1, relArea=1, iter);

      min=min*RelVal;
      max=max*RelVal;
      val=val*RelVal;
      }

   if (R_IsNA(val)) // max and min bounds
      {
      val=getVal(quantity, iyr, iunit, iSn, iarea, iter);

      if (!R_IsNA(min) && val<min) val = min;
      if (!R_IsNA(max) && val>max) val = max;
      }

   //NLSE's
   switch (quantity)
      {
      case FLRConst_F:
         func[0] = Fbar(            ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_Z:
         func[0] = Zbar(            ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_FLandings:
         func[0] = FbarLandings(    ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_FDiscards:
         func[0] = FbarDiscards(    ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_SSB:
         func[0] = SSB(             ad_n,ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_Biomass:
         func[0] = computeStock(    ad_n,ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_Catch:
         func[0] = computeCatch(    ad_f,      iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_Landings:
         func[0] = computeLandings( ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_Discards:
         func[0] = computeDiscards( ad_f, iyr, iunit, iSn, iarea, iter) - val;
         break;
      case FLRConst_MnSz:
         func[0] = MnSz(            ad_n,      iyr, iunit, iSn, iarea, iter) - val;
         break;
      default:
         func[0] = 0.0;
         break;
      }
   
  double value = func[0].value();
  }

void fwdStk::project(double *x, int iyr, int iunit, int iseason, int iarea, int iter, bool OnlyReplaceNA, bool OnlyCalcN)
   {
   // set F
   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      stk.harvest(iage,iyr,iunit,iseason,iarea,iter) = __max(0.0,stk.harvest(iage, iyr, iunit, iseason, iarea, iter)*x[0]);

   double _fbar = stk.Fbar(iyr,iunit,iseason,iarea,iter);
   
   if (_fbar>MaxF(1,iyr,iunit,iseason,iarea,iter))
  	  for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
   	     stk.harvest(iage,iyr,iunit,iseason,iarea,iter)=stk.harvest(iage,iyr,iunit,iseason,iarea,iter)*MaxF(1,iyr,iunit,iseason,iarea,iter)/_fbar;
	  
   // recruits
   int SSB_yr = __min(__max(iyr-stk.minquant,stk.minyr),stk.maxyr);
   if (iseason==1) 
       {
       if ((OnlyReplaceNA && R_IsNA(stk.stock_n(stk.minquant,iyr,iunit,iseason,iarea,iter))) || !OnlyReplaceNA)
	      stk.stock_n(stk.minquant,iyr,iunit,iseason,iarea,iter) = SR.recruits(1,stk.SSB(SSB_yr,iunit,iseason,iarea,iter),iyr,iunit,iseason,iarea,iter);
       }
   else if (iseason >1) stk.stock_n(stk.minquant,iyr,iunit,iseason,iarea,iter) = stk.stock_n(stk.minquant,iyr,iunit,iseason-1,iarea,iter) * exp(stk.harvest(stk.minquant,iyr,iunit,iseason-1,iarea,iter)-stk.m(stk.minquant,iyr,iunit,iseason-1,iarea,iter))+
	                                                                             SR.recruits(1,stk.SSB(SSB_yr,iunit,iseason,iarea,iter),iyr,iunit,iseason,  iarea,iter);
   // Projection
   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      {       
      //numbers-at-age next season/year
      if (iseason<stk.nseasons)
         stk.stock_n(iage,iyr,iunit,iseason+1,iarea,iter) = stk.stock_n(      iage,iyr,iunit,iseason,iarea,iter)
		                                                    *exp(-stk.harvest(iage,iyr,iunit,iseason,iarea,iter)
															     -stk.m(      iage,iyr,iunit,iseason,iarea,iter));
      else
         {
         if (iage < stk.maxquant)
            stk.stock_n(iage+1,iyr+1,iunit,1,iarea,iter)  = stk.stock_n(       iage,iyr,iunit,iseason,iarea,iter)
			                                                 *exp(-stk.harvest(iage,iyr,iunit,iseason,iarea,iter)
																  -stk.m(      iage,iyr,iunit,iseason,iarea,iter));
         if (iage == stk.plusgrp)
            stk.stock_n(iage,  iyr+1,iunit,1,iarea,iter) += stk.stock_n(              stk.maxquant,iyr,iunit,iseason,iarea,iter)
			                                                        *exp(-stk.harvest(stk.maxquant,iyr,iunit,iseason,iarea,iter)
																	     -stk.m(      stk.maxquant,iyr,iunit,iseason,iarea,iter));
         }          
	  }

    // Redistribute between areas if more than 1 area
    if (stk.nareas>1){
       if (iseason<stk.nseasons){        //within year
          for (int iage=stk.minquant; iage<=stk.maxquant; iage++){
             double sum=0.0;
             for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                sum+=stk.stock_n(iage,iyr,iunit,iseason+1,jarea,iter);
          
		     for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                stk.stock_n(iage,iyr,iunit,iseason+1,jarea,iter)=sum*avail(iage,iyr,iunit,iseason+1,jarea,iter);
             }	      
	      }
	   else if (avail.maxyr()>iyr) { //next year only if avial provided
          for (int iage=stk.minquant; iage<=stk.maxquant; iage++){
             double sum=0.0;
             for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                sum+=stk.stock_n(iage,iyr+1,iunit,1,jarea,iter);
          
		     for (int jarea=stk.nareas; jarea<=stk.nareas; jarea++)
                stk.stock_n(iage,iyr+1,iunit,1,jarea,iter)=sum*avail(iage,iyr+1,iunit,1,jarea,iter);
             }	      	      
	      }
	   }

      if (!OnlyReplaceNA || (OnlyReplaceNA && R_IsNA(stk.stock_n(stk.minquant,iyr+1,iunit,iseason,iarea,iter))))    
         if (SR.recruits(1,stk.SSB(SSB_yr,iunit,iseason,iarea,iter),iyr+1,iunit,iseason,iarea,iter)>0) 
            stk.stock_n(stk.minquant,iyr+1,iunit,iseason,iarea,iter) = SR.recruits(1,stk.SSB(SSB_yr,iunit,iseason,iarea,iter),iyr+1,iunit,iseason,iarea,iter);

	  if (!OnlyCalcN)
         for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
            {
            double z    =  stk.m(         iage, iyr, iunit, iseason, iarea, iter) + stk.harvest(   iage, iyr, iunit, iseason, iarea, iter),
                   ctch =  stk.discards_n(iage, iyr, iunit, iseason, iarea, iter) + stk.landings_n(iage, iyr, iunit, iseason, iarea, iter);
                 
             stk.discards_n( iage, iyr, iunit, iseason, iarea, iter)=stk.discards_n( iage, iyr, iunit, iseason, iarea, iter)/ctch;
             stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)=stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)/ctch;
                  
             stk.catch_n( iage, iyr, iunit, iseason, iarea, iter) =
                stk.stock_n( iage, iyr, iunit, iseason, iarea, iter)*
                stk.harvest( iage, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z));

             stk.discards_n( iage, iyr, iunit, iseason, iarea, iter)=stk.discards_n( iage, iyr, iunit, iseason, iarea, iter)*stk.catch_n( iage, iyr, iunit, iseason, iarea, iter);
             stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)=stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)*stk.catch_n( iage, iyr, iunit, iseason, iarea, iter);
             }
   } 

adouble fwdStk::computeStock(FLQuant_adolc &n, FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= stk.minquant; iage<= stk.maxquant; iage++)
      {
      double mass;
      adouble deadeds;
     
         int _age = __min(iage+1,__max(stk.plusgrp,stk.maxquant));
         
         deadeds= exp(-stk.m(iage, iyr, iunit, iseason, iarea, iter)
                      -f(    iage, iyr, iunit, iseason, iarea, iter));
         
         mass   = stk.stock_wt(_age, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter);
      
      val += stk.stock_n(iage,iyr,iunit,iseason,iarea,iter)*deadeds*mass;
      }

   val +=            n(stk.minquant, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter)*
          stk.stock_wt(stk.minquant, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter);
      
   return val;
   }  

adouble fwdStk::SSB(FLQuant_adolc &n, FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter) 
   {
   adouble val = 0.0, ssb = 0.0;

   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      {
      adouble survivors;
        
      // if spawning before any fishing then project to end of year
      // but natural mortality might still occur
      if (stk.harvest_spwn(iage,iyr,iunit,iseason,iarea,iter) == 0)
         {
         double mass = stk.stock_wt(iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter)*
                       stk.mat(     iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter);
      
         survivors = exp(-stk.m(     iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter)*
                          stk.m_spwn(iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter));
         
         val = n(iage,iyr+1,iunit,iseason,iarea,iter)*mass*survivors;
         }
      else
         {
         double mass = stk.stock_wt(iage, iyr, iunit, iseason, iarea, iter)*
                       stk.mat(     iage, iyr, iunit, iseason, iarea, iter);
        
         survivors = exp(-stk.m(iage, iyr, iunit, iseason, iarea, iter)*stk.m_spwn(      iage, iyr, iunit, iseason, iarea, iter)
                         -f(    iage, iyr, iunit, iseason, iarea, iter)*stk.harvest_spwn(iage, iyr, iunit, iseason, iarea, iter));
           
         val = stk.stock_n(iage,iyr,iunit,iseason,iarea,iter)*mass*survivors;
         }

      if (val>0.0) ssb +=val;
      }

   return ssb;
   }  
                             
adouble fwdStk::computeCatch(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble RtnVal = 0.0, 
           val    = 0.0;

   double t1, t2;

   for (int iage= stk.minquant; iage<= stk.maxquant; iage++)
      {
      adouble z = stk.m(iage, iyr, iunit, iseason, iarea, iter) + f(iage, iyr, iunit, iseason, iarea, iter);

      val = stk.stock_n( iage, iyr, iunit, iseason, iarea, iter)*
            f(           iage, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z))*
            stk.catch_wt(iage, iyr, iunit, iseason, iarea, iter);

t1 =val.value();
t2 =RtnVal.value();

      RtnVal += val;
      }

   return RtnVal;
   } 

adouble fwdStk::computeDiscards(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= stk.minquant; iage<= stk.maxquant; iage++)
      {
      adouble z = stk.m(iage, iyr, iunit, iseason, iarea, iter) + f(iage, iyr, iunit, iseason, iarea, iter);

      val += stk.discards_n( iage, iyr, iunit, iseason, iarea, iter)/
           (stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)+
             stk.discards_n( iage, iyr, iunit, iseason, iarea, iter))* 
             stk.stock_n(    iage, iyr, iunit, iseason, iarea, iter)*
             f(              iage, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z))*
             stk.discards_wt(iage, iyr, iunit, iseason, iarea, iter);
 
     }

   return val;
   } 

adouble fwdStk::computeLandings(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= stk.minquant; iage<= stk.maxquant; iage++)
      {
      adouble z = stk.m(iage, iyr, iunit, iseason, iarea, iter) + f(iage, iyr, iunit, iseason, iarea, iter);

      val += stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)/
           (stk.landings_n( iage, iyr, iunit, iseason, iarea, iter)+
             stk.discards_n( iage, iyr, iunit, iseason, iarea, iter))* 
             stk.stock_n(    iage, iyr, iunit, iseason, iarea, iter)*
             f(              iage, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z))*
             stk.landings_wt(iage, iyr, iunit, iseason, iarea, iter);

double t=val.value(); 
     }

   return val;
   } 

adouble fwdStk::Fbar(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage=stk.minfbar; iage<=stk.maxfbar; iage++)
      val += f(iage, iyr, iunit, iseason, iarea, iter);

   return val/(stk.maxfbar-stk.minfbar+1);
   }                               

double fwdStk::Fbar(int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   double val = 0.0;

   for (int iage=stk.minfbar; iage<=stk.maxfbar; iage++)
      val += stk.harvest(iage, iyr, iunit, iseason, iarea, iter);

   return val/(stk.maxfbar-stk.minfbar+1);
   }                               

adouble fwdStk::Zbar(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= stk.minfbar; iage<= stk.maxfbar; iage++)
      val += f(    iage, iyr, iunit, iseason, iarea, iter)+
             stk.m(iage, iyr, iunit, iseason, iarea, iter);

   return val/(stk.maxfbar-stk.minfbar+1);
   }                               

adouble fwdStk::FbarLandings(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;
double t;

   for (int iage= stk.minfbar; iage<= stk.maxfbar; iage++){
      val += f(             iage, iyr, iunit, iseason, iarea, iter)*
             stk.landings_n(iage, iyr, iunit, iseason, iarea, iter)/
             stk.catch_n(   iage, iyr, iunit, iseason, iarea, iter);
      
      t=val.value();
      }

   return val/(stk.maxfbar-stk.minfbar+1);
   }                               

adouble fwdStk::FbarDiscards(FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= stk.minfbar; iage<= stk.maxfbar; iage++)
      val += f(             iage, iyr, iunit, iseason, iarea, iter)*
             stk.discards_n(iage, iyr, iunit, iseason, iarea, iter)/
             stk.catch_n(   iage, iyr, iunit, iseason, iarea, iter);

   return val/(stk.maxfbar-stk.minfbar+1);
   }   

double fwdStk::SSB(int iyr, int iunit, int iseason, int iarea, int iter) 
   {
   double val = 0.0, ssb = 0.0;

   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      {
      double survivors;
        
      // if spawning before any fishing then project to end of year
      // but natural mortality might still occur
      if (stk.harvest_spwn(iage,iyr,iunit,iseason,iarea,iter) == 0)
         {
         double mass = stk.stock_wt(iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter)*
                       stk.mat(     iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter);
      
         survivors = exp(-stk.m(     iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter)*
                          stk.m_spwn(iage, __min(iyr+1,stk.maxyr), iunit, iseason, iarea, iter));
         
         val = stk.stock_n(iage,iyr+1,iunit,iseason,iarea,iter)*mass*survivors;
         }
      else
         {
         double mass = stk.stock_wt(iage, iyr, iunit, iseason, iarea, iter)*
                       stk.mat(     iage, iyr, iunit, iseason, iarea, iter);
        
         survivors = exp(-stk.m(      iage, iyr, iunit, iseason, iarea, iter)*stk.m_spwn(      iage, iyr, iunit, iseason, iarea, iter)
                         -stk.harvest(iage, iyr, iunit, iseason, iarea, iter)*stk.harvest_spwn(iage, iyr, iunit, iseason, iarea, iter));
           
         val = stk.stock_n(iage,iyr,iunit,iseason,iarea,iter)*mass*survivors;
         }

      if (val>0.0) ssb +=val;
      }

   return ssb;
   }  

adouble fwdStk::MnSz(FLQuant_adolc &n, int iyr, int iunit, int iseason, int iarea, int iter) 
   {
   adouble val = 0.0, mnsz = 0.0, sumN=0.0;

   for (int iage=stk.minquant; iage<=stk.maxquant; iage++)
      {
      val = n(iage,iyr+1,iunit,iseason,iarea,iter)*stk.stock_wt(iage,iyr+1,iunit,iseason,iarea,iter);

      if (val>0.0) {
          sumN +=n(iage,iyr+1,iunit,iseason,iarea,iter);
          mnsz +=val;
          }
      }


   return mnsz/sumN;}  

void fwdStk::InitAvail(SEXP x) 
   {
   avail.Init(x);     
   
   // ensure relative abundance
   int iIter, iAge, iYr, iUnit, iSeason, iArea;  

   for (iAge=avail.minquant(); iAge<=avail.maxquant(); iAge++)
     for (iYr=avail.minyr(); iYr<=avail.maxyr(); iYr++)
       for (iUnit=1; iUnit<=avail.nunits(); iUnit++)
         for (iSeason=1; iSeason<=avail.nseasons(); iSeason++)
           for (iIter=1; iIter<=avail.niters(); iIter++)
             {
             double sum = 0.0;
             for (iArea=1; iArea<=avail.nareas(); iArea++)
                sum += avail(iAge,iYr,iUnit,iSeason,iArea,iIter);
    
             for (iArea=1; iArea<=avail.nareas(); iArea++)
                avail(iAge,iYr,iUnit,iSeason,iArea,iIter) /= sum;
             }   
   }

SEXP fwdStk::Init(SEXP xStk, SEXP xYrs, SEXP xSRModel,SEXP xSRParam,SEXP xSRResiduals,SEXP xMult,SEXP xAvail, SEXP xMaxF)    
    {
    SEXP Err = PROTECT(NEW_NUMERIC(1)); 

    //Set Stock 
    stk.Init(xStk); 

    //Set SRR
    //REAL(Err)[0]=1.0;
    if (!SR.Init(1, xYrs)) {
       UNPROTECT(1);
       return Err;}

    //REAL(Err)[0]=2.0;
    if (!SR.Init(1, xSRModel, xSRParam, xSRResiduals, xMult))  {
       UNPROTECT(1);
       return Err;}
    
    avail.Init(xAvail);
    MaxF.Init(xMaxF);

    UNPROTECT(1);

    return stk.Return();
    }    

SEXP fwdStk::run(SEXP xTrgt, SEXP xAry)    
    {
    SEXP Err = PROTECT(NEW_NUMERIC(1)); 

    // check target object
    REAL(Err)[0]=3.0;
    if (!isMatrix(xTrgt) || !isNumeric(xTrgt)) {
       ;
       //UNPROTECT(1);
       //return Err;
       }

    // check target min/max/value object
    REAL(Err)[0]=4.0;
    if (!isArray(xAry) || !isNumeric(xAry)) {
       ;
       //UNPROTECT(1);
       //return Err;
       }

    SEXP TrgtDims = GET_DIM(xTrgt);

    REAL(Err)[0]=5.0*0 +(double)(LENGTH(TrgtDims)*1.0) + (double)(INTEGER(TrgtDims)[1]*100.0);
    if (LENGTH(TrgtDims) != 2 || INTEGER(TrgtDims)[1] != 15)  {
       ;
       //UNPROTECT(1);
       //return Err;
       }
  
    //SEXP AryDims = GET_DIM(xAry);

    //REAL(Err)[0]=6.0*0.0 +  (double) stk.niters;  //INTEGER(AryDims)[0]- INTEGER(TrgtDims)[0] + INTEGER(AryDims)[1]-3 + INTEGER(AryDims)[2] - 
                            
    //if (LENGTH(AryDims) != 3 || INTEGER(AryDims)[0] != INTEGER(TrgtDims)[0] || 
    //                            INTEGER(AryDims)[1] != 3                    ||
    //                            INTEGER(AryDims)[2] != (int) stk.niters)  {
       //UNPROTECT(1);
       //return Err;
    //   }

    double *Trgt = NUMERIC_POINTER(xTrgt);
    double *Ary  = NUMERIC_POINTER(xAry);

    //ADol-C stuff
    int i, n=1;

    double  *depen,    *indep,   *r, **jac;
    adouble *depen_ad, *indep_ad;
    
    depen    = new   double[n];
    indep    = new   double[n];
    r        = new   double[n];
    jac      = new  double*[n];
    depen_ad = new  adouble[n];
    indep_ad = new  adouble[n];
    
	for (i=0; i<n; i++){
       jac[i]   = new double[n];
	   indep[i] = 1.0;
	   }

    int iTrgt = 0, 
        iter  = 0,
        nrow  = (int)INTEGER(TrgtDims)[0];


    //get N at start of year
    for (iter=1; iter<=stk.niters; iter++)
      for (int iunit=1; iunit<=stk.nunits; iunit++)
        for (int iarea=1; iarea<=stk.nareas; iarea++)
          for (int iseason=stk.nseasons; iseason<=stk.nseasons; iseason++)
              project(indep, SR.minyear()-1, iunit, iseason, iarea, iter, TRUE, TRUE);

    int iTape = 1, _Tape;
    for (iTrgt=1; iTrgt<=(int)(INTEGER(TrgtDims)[0]); iTrgt++)
       for (iter=1; iter<=stk.niters; iter++)
          {
          FLRConst_Target quantity = (FLRConst_Target)(int)(Trgt)[iTrgt-1 + fwdTargetPos_quantity*nrow];
          
          _Tape = 0; //iTrgt % ++iTape + 1;
      
          for (i=0; i<n; i++)
			  if (indepLastYr)
                indep[i]=0.01/Fbar((int)(Trgt)[iTrgt-1+fwdTargetPos_year  *nrow], 
                                   (int)(Trgt)[iTrgt-1+fwdTargetPos_unit  *nrow], 
                                   (int)(Trgt)[iTrgt-1+fwdTargetPos_season*nrow], 
                                   (int)(Trgt)[iTrgt-1+fwdTargetPos_area  *nrow],  iter);
			  else
                indep[i]=0.1;

           // Taping the computation of the jacobian 
           trace_on(_Tape);

           // marking independent variables 
           for (i=0; i<n; i++)
              indep_ad[i] <<= indep[i];

          project(indep_ad,depen_ad,Trgt,iTrgt,nrow,Ary,iter);

          // marking dependent variables 
          for (i=0; i<n; i++)
             depen_ad[i] >>= depen[i];

          trace_off(_Tape);

          //jacobian(tag,n,n,indep,jac);
          r[0]=1.0;
          function(_Tape,n,n,indep,r);
          int NIters=0;
          while (norm(r,n) > 1e-12 && norm(indep,n) < 100 && NIters++<50)
              {
              function(_Tape,n,n,indep,r);

              jac_solv(_Tape,n,indep,r,0,2);

              for (i=0; i<n; i++)
                  indep[i] -= r[i];	   
              }         

		  project(indep, (int)(Trgt)[iTrgt-1+fwdTargetPos_year  *nrow], 
                         (int)(Trgt)[iTrgt-1+fwdTargetPos_unit  *nrow], 
                         (int)(Trgt)[iTrgt-1+fwdTargetPos_season*nrow], 
                         (int)(Trgt)[iTrgt-1+fwdTargetPos_area  *nrow],  iter);
          }

    delete[] depen;
    delete[] indep;
    delete[] r;
    delete[] depen_ad;
    delete[] indep_ad;
    
    for (i=0; i<n; i++)
       delete[] jac[i];
    delete[] jac;

    UNPROTECT(1);

    return stk.Return();
    }    

