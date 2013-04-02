#include "fwd.h"
#include "Rmath.h"
#include <Rdefines.h>
#include <Rinternals.h>

#define const_nan 0.0

double norm(double *x, int n)
  {
  int i;
  
  double out=0;
  for (i=0; i<n; i++)
    out += x[i]*x[i];

  out = sqrt(out);
  
  return out;
  }

FLRConst_Target get_target(int i)
   {
   switch(i) {
      case 1:  return FLRConst_SSB;       
      case 2:  return FLRConst_Biomass;    
      case 3:  return FLRConst_Catch;      
      case 4:  return FLRConst_Landings;   
      case 5:  return FLRConst_Discards;   
      case 6:  return FLRConst_F;          
      case 7:  return FLRConst_Z;          
      case 8:  return FLRConst_FLandings; 
      case 9:  return FLRConst_FDiscards; 
      case 10: return FLRConst_Effort;     
      case 11: return FLRConst_Costs;      
      case 12: return FLRConst_Revenue;    
      case 13: return FLRConst_Profit;     
      case 14: return FLRConst_MnSz;     
      case 0:  return FLRConst_None;       
      }

   return FLRConst_None;
   }

int fwd::NBiol(void)
   {
   return flc::_nstock;
   }
      
fwd::fwd(SEXP xStk, SEXP xYrs, SEXP xSRModel, SEXP xSRParam, SEXP xSRResiduals, SEXP xMult)
  {
  Ctrl.AllocFlag     = FALSE;
  Trgt.AllocDataFlag = FALSE;
  Trgt.AllocAryFlag  = FALSE;

  Init(xStk);

  //Set SRR
  InitSR(1, xYrs);
    
  InitSR(1, xSRModel, xSRParam, xSRResiduals, xMult);
  }

fwd::fwd(SEXP xBiol, SEXP xFleet, SEXP xYrs, SEXP xDim, SEXP xSRR)
  {
  Ctrl.AllocFlag     = FALSE;
  Trgt.AllocDataFlag = FALSE;
  Trgt.AllocAryFlag  = FALSE;

  InitBiolFleet(xBiol,xFleet,xDim);
  
  InitSR(NBiol(), xYrs);

  //Set SRR
  //if (NBiol() != NElemList(xSRR))
  //    return;
 
  for (int i=0; i<NBiol(); i++)
     {
     //Set SRR
     InitSR(i+1, PROTECT(VECTOR_ELT(xSRR, 0)),
                 PROTECT(VECTOR_ELT(xSRR, 1)),
                 PROTECT(VECTOR_ELT(xSRR, 2)),
                 PROTECT(VECTOR_ELT(xSRR, 3)));
  
     UNPROTECT(4);
     }
  }

fwd::fwd(SEXP xBiol, SEXP xFleet, SEXP xYrs, SEXP xDim, SEXP xSRModel, SEXP xSRParam, SEXP xSRResiduals, SEXP xMult)
  {
  Ctrl.AllocFlag     = FALSE;
  Trgt.AllocDataFlag = FALSE;
  Trgt.AllocAryFlag  = FALSE;

  InitBiolFleet(xBiol,xFleet,xDim);

  if (!InitSR(1, xYrs)) 
     return;

  //Set SRR
  bool t =FALSE;
  
  t = InitSR(1, xSRModel,     
                xSRParam,     
                xSRResiduals, 
                xSRModel);
  }
  
bool fwd::run(SEXP xTrgt, SEXP xAryTrgt, SEXP xCtrl, SEXP xAryCtrl, SEXP xYrs)
   {
   if (!Trgt.Init(xTrgt,xAryTrgt,niters())) return false;
   if (!Ctrl.Init(xCtrl,xAryTrgt,niters())) return false;
   
   MinProjYear = (int)REAL(xYrs)[0]; 
   MaxProjYear = (int)REAL(xYrs)[0];

   int i; 
   for (i=1; i<LENGTH(xYrs); i++)
      {
      if (MinProjYear>REAL(xYrs)[i]) MinProjYear=(int)REAL(xYrs)[i];
      if (MaxProjYear<REAL(xYrs)[i]) MaxProjYear=(int)REAL(xYrs)[i];
      }

   //ADol-C stuff
   int n=1;

   double  *depen,    *indep,   *r, **jac;
   adouble *depen_ad, *indep_ad;
   
   int iter  = 0;

   //get N at start of year
   double x=1.0;
   for (iter=1; iter<=niters(); iter++)
      project(&x, MinProjYear-1 ,iter, TRUE, TRUE);

   int tag = 0;
   for (int iYr=MinProjYear; iYr<=MaxProjYear; iYr++)
      {
      int _tag = 0; //(tag % 6 + 1);

      n = (int)Trgt.n(iYr);

      //n of independent variables MUST = n of equations
      if (Ctrl.n(iYr) == n)
         {
         depen    = new   double[n];
         indep    = new   double[n];
         r        = new   double[n];
         jac      = new  double*[n];
         depen_ad = new  adouble[n];
         indep_ad = new  adouble[n];
   
         for (i=0; i<n; i++)
            jac[i] = new double[n];

         for (iter=1; iter<=niters(); iter++)
            {
            // set independent variables to estimaate
            int j=0;
            for (int iFleet=1; iFleet<=nfleet(); iFleet++)
               for (int iMetier=1; iMetier<=nmetier(); iMetier++)
                  if (Ctrl.fit(iYr, iFleet, iMetier))
                     indep[j++] = 1.0;
      
            // Taping the computation of the jacobian 
            trace_on(_tag);

            // marking independent variables 
            for (i=0; i<n; i++)
               indep_ad[i] <<= indep[i];
   
            project(indep_ad,depen_ad,iYr,iter);  

            // marking dependent variables 
            for (i=0; i<n; i++)
               depen_ad[i] >>= depen[i];

            trace_off(_tag);

            //jacobian(tag,m,n,indep,jac);
            r[0]=1.0;
            function(_tag,n,n,indep,r);
            int NIters=0;
	         while (norm(r,n) > 1e-10 && NIters++<50)
	            {
	            jac_solv(_tag,n,indep,r,0,2);

	            for (i=0; i<n; i++)
		             indep[i] -= r[i];	   

	            function(_tag,n,n,indep,r);
               }         
       
            project(indep, iYr, iter);
            }

         delete[] depen;
         delete[] indep;
         delete[] r;
         delete[] depen_ad;
         delete[] indep_ad;
   
         for (i=0; i<n; i++)
            delete[] jac[i];
         delete[] jac;
         }  
 
      tag++;
      }

   return true;
   }

adouble fwd::computeCatch(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   int iAge, iFleet=1, iMetier=1;
   for (iAge= minage(ispp); iAge<= maxage(ispp); iAge++)
      {
      adouble z = ad_f(ispp, iAge, iyr, iunit, iseason, iarea, iter)+
                     m(ispp, iAge, iyr, iunit, iseason, iarea, iter);
     
       val += n(                        ispp, iAge, iyr, iunit, iseason, iarea, iter)*
              ad_f(                     ispp, iAge, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z))*
              catch_wt(iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter);
       }  

   return val;
   } 

adouble fwd::computeStock(FLQuant2_adolc &ad_n, FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= minage(ispp); iage<= maxage(ispp); iage++)
      {
      double mass;
      adouble deadeds;
     
      int _age = __min(iage+1,__max(plusgrp(ispp),maxage(ispp)));
         
      deadeds= exp(-m(ispp, iage, iyr, iunit, iseason, iarea, iter)
                   -ad_f(    ispp, iage, iyr, iunit, iseason, iarea, iter));
         
      mass   = stock_wt(ispp, _age, __min(iyr+1,maxyr(ispp)), iunit, iseason, iarea, iter);
      
      val += n(ispp, iage,iyr,iunit,iseason,iarea,iter)*deadeds*mass;
      }

   val +=            n(ispp, minage(ispp), iyr+1,                        iunit, iseason, iarea, iter)*
          stock_wt(ispp, minage(ispp), __min(iyr+1,maxyr(ispp)), iunit, iseason, iarea, iter);

   return val;
   }  

adouble fwd::SSB(FLQuant2_adolc &ad_n, FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter) 
   {
   adouble val = 0.0;
  
   for (int iage=minage(ispp); iage<=maxage(ispp); iage++)
      {
      adouble survivors;
        
      // if spawning before any fishing then project to end of year
      // but natural mortality might still occur
      if (f_spwn(ispp,iage,iyr,iunit,iseason,iarea,iter) == 0)
         {
         double mass = stock_wt(   ispp, iage, __min(iyr+1,maxyr(ispp)), iunit, iseason, iarea, iter)*
                       fec(        ispp, iage, __min(iyr+1,maxyr(ispp)), iunit, iseason, iarea, iter);
      
         survivors   = exp(-m(     ispp, iage, __min(iyr+1,maxyr(ispp)), iunit, iseason, iarea, iter)*
                            m_spwn(ispp, iage, __min(iyr+1,maxyr(ispp)), iunit, iseason, iarea, iter));
         
         val += ad_n(ispp,iage,iyr+1,iunit,iseason,iarea,iter)*mass*survivors;
         }
      else
         {
         double mass = stock_wt(ispp, iage, iyr, iunit, iseason, iarea, iter)*
                       fec(     ispp, iage, iyr, iunit, iseason, iarea, iter);
        
         survivors = exp(-m(      ispp, iage, iyr, iunit, iseason, iarea, iter)*m_spwn(      ispp, iage, iyr, iunit, iseason, iarea, iter)
                         -ad_f(ispp, iage, iyr, iunit, iseason, iarea, iter)*f_spwn(ispp, iage, iyr, iunit, iseason, iarea, iter));
           
         val += ad_n(ispp, iage,iyr,iunit,iseason,iarea,iter)*mass*survivors;
         }
      }
   return val;
   }  
                              
adouble fwd::computeDiscards(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   int iAge, iFleet=1, iMetier=1;
   for (iAge= minage(ispp); iAge<= maxage(ispp); iAge++)
      {
      adouble z = ad_f(ispp, iAge, iyr, iunit, iseason, iarea, iter)+
                     m(ispp, iAge, iyr, iunit, iseason, iarea, iter);
     
      adouble t = n(                            ispp, iAge, iyr, iunit, iseason, iarea, iter)*
                  ad_f(                         ispp, iAge, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z))*
                  catch_wt(    iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter);

      double ratio = discards_sel(iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter)/
                     catch_sel(   iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter);
     
      if (R_IsNA(ratio) && ratio<0.0) ratio=0.0;
      
      val += t*ratio;   
      }  
  
   return val;
   } 

adouble fwd::computeLandings(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   int iAge, iFleet=1, iMetier=1;
   for (iAge= minage(ispp); iAge<= maxage(ispp); iAge++)
      {
      adouble z =    m(ispp, iAge, iyr, iunit, iseason, iarea, iter)+
                  ad_f(ispp, iAge, iyr, iunit, iseason, iarea, iter);
     
      adouble t =    n(ispp, iAge, iyr, iunit, iseason, iarea, iter)*
                  ad_f(ispp, iAge, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z));

      double ratio =  landings_sel(iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter)/
                     (landings_sel(iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter)+
                      discards_sel(iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter));
     
      if (R_IsNA(ratio) && ratio<0.0) ratio=0.0;
      
      val += t*ratio*landings_wt(iFleet, iMetier, ispp, iAge, iyr, iunit, iseason, iarea, iter);   
      }  
  
   return val;
   } 

adouble fwd::Fbar(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= minfbar(ispp); iage<= maxfbar(ispp); iage++)
      val += ad_f(ispp, iage, iyr, iunit, iseason, iarea, iter);

   return val/(maxfbar(ispp)-minfbar(ispp)+1);
   }                               

adouble fwd::Zbar(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= minfbar(ispp); iage<= maxfbar(ispp); iage++)
      val += ad_f(ispp, iage, iyr, iunit, iseason, iarea, iter)+
             m(   ispp, iage, iyr, iunit, iseason, iarea, iter);

   return val/(maxfbar(ispp)-minfbar(ispp)+1);
   }                               

adouble fwd::FbarLandings(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= minfbar(ispp); iage<= maxfbar(ispp); iage++)
      val += ad_f(      ispp, iage, iyr, iunit, iseason, iarea, iter)*
             landings_n(ispp, iage, iyr, iunit, iseason, iarea, iter)/
             catch_n(   ispp, iage, iyr, iunit, iseason, iarea, iter);

   return val/(maxfbar(ispp)-minfbar(ispp)+1);
   }                               

adouble fwd::FbarDiscards(FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter)      
   {
   adouble val = 0.0;

   for (int iage= minfbar(ispp); iage<= maxfbar(ispp); iage++)
      val += ad_f(      ispp, iage, iyr, iunit, iseason, iarea, iter)*
             discards_n(ispp, iage, iyr, iunit, iseason, iarea, iter)/
             catch_n(   ispp, iage, iyr, iunit, iseason, iarea, iter);

   return val/(maxfbar(ispp)-minfbar(ispp)+1);
   }   

double fwd::getVal(FLRConst_Target quantity, int ispp,  int iyr, int iunit, int iseason, int iarea, int iter)
   {   
   switch(quantity)
    {
    case FLRConst_F:
      return flc::Fbar(           ispp, iyr,  iunit, iseason, iarea, iter);
  	case FLRConst_FLandings:
      return flc::FbarLandings(   ispp, iyr,  iunit, iseason, iarea, iter);
   	case FLRConst_FDiscards:
      return flc::FbarDiscards(   ispp, iyr,  iunit, iseason, iarea, iter);
   	case FLRConst_SSB:
      return flc::SSB(            ispp, iyr-1,iunit, iseason, iarea, iter);
   	case FLRConst_Biomass:
      return flc::computeStock(   ispp, iyr,  iunit, iseason, iarea, iter);
    case FLRConst_Catch:
      return flc::computeCatch(   ispp, iyr,  iunit, iseason, iarea, iter);
    case FLRConst_Landings:
      return flc::computeLandings(ispp, iyr,  iunit, iseason, iarea, iter);
    case FLRConst_Discards:
      return flc::computeDiscards(ispp, iyr,  iunit, iseason, iarea, iter);
 	  default:
	    return 0.0;
	    break;
       }

   return 0.0;
   }

void fwd::project(double *x, int iyr, int iter, bool OnlyReplaceNA, bool OnlyCalcN)
   {
   int iunit  =1,
	   iseason=1,
       iarea  =1;

    for (int ispp=1; ispp<=NBiol(); ispp++)
       {
       for (int iage=minage(ispp); iage<=maxage(ispp); iage++)
          {
          int i = 0;
          f(ispp, iage,iyr,iunit,iseason,iarea,iter) = 0.0;
          for (int iFleet=1; iFleet<=nfleet(); iFleet++)
             for (int iMetier=1; iMetier<=nmetier(); iMetier++)
                if (!Ctrl.fit(iyr, iFleet, iMetier))
                   {
                   f(ispp, iage,iyr,iunit,iseason,iarea,iter) += __max(0.0,catch_sel(iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter)*                                                          
                                                                           catch_q(  iFleet, iMetier, ispp, 1,    iyr, iunit, iseason, iarea, iter)*
                                                                           effshare( iFleet, iMetier,       1,    iyr, iunit, iseason, iarea, iter)*
                                                                           effort(   iFleet,                1,    iyr, iunit, iseason, iarea, iter));                                                          

                   }
                else
                   {
                   if (!R_IsNA(Ctrl.min(iyr, iFleet, iMetier)))
                      x[i]=__max(x[i],Ctrl.min(iyr, iFleet, iMetier));
                   if (!R_IsNA(Ctrl.max(iyr, iFleet, iMetier)))
                      x[i]=__min(x[i],Ctrl.max(iyr, iFleet, iMetier));

                   f(ispp, iage,iyr,iunit,iseason,iarea,iter) += __max(0.0,catch_sel(iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter)*                                                          
                                                                           catch_q(  iFleet, iMetier, ispp, 1,    iyr, iunit, iseason, iarea, iter)*
                                                                           effshare( iFleet, iMetier,       1,    iyr, iunit, iseason, iarea, iter)*
                                                                           x[i]);

                   effort(iFleet,1,iyr,iunit,iseason,iarea,iter) = x[i];
                   i++;
                   }                                                         
            
          //numbers-at-age next year
          if (iage < maxage(ispp))
             n(ispp, iage+1,iyr+1,iunit,iseason,iarea,iter)  = n(ispp, iage,    iyr,iunit,iseason,iarea,iter)*exp(-f(ispp,     iage,iyr,iunit,iseason,iarea,iter)-m(ispp, iage,    iyr,iunit,iseason,iarea,iter));
          if (iage == plusgrp(ispp))
             n(ispp, iage,  iyr+1,iunit,iseason,iarea,iter) += n(ispp, maxage(ispp),iyr,iunit,iseason,iarea,iter)*exp(-f(ispp, maxage(ispp),iyr,iunit,iseason,iarea,iter)-m(ispp, maxage(ispp),iyr,iunit,iseason,iarea,iter));
          }

       int rec_yr = __min(__max(iyr-minage(ispp)+1,minyr(ispp)),maxyr(ispp));
    
       if (!OnlyReplaceNA || (OnlyReplaceNA && R_IsNA(n(ispp, minage(ispp),iyr+1,iunit,iseason,iarea,iter))))    
          n(ispp, minage(ispp),iyr+1,iunit,iseason,iarea,iter) = _sr.recruits(ispp,iyr,flc::SSB(ispp,rec_yr,iunit,iseason,iarea,iter),iter);

       if (!OnlyCalcN)
          for (int iage=minage(ispp); iage<=maxage(ispp); iage++)
             {
             double z = m(ispp, iage, iyr, iunit, iseason, iarea, iter) + f(ispp, iage, iyr, iunit, iseason, iarea, iter);

             for (int iFleet=1; iFleet<=nfleet(); iFleet++)
                for (int iMetier=1; iMetier<=nmetier(); iMetier++)
                   {
                   double c_, d_, l_;

                   c_ = 
                       n(                         ispp, iage, iyr, iunit, iseason, iarea, iter)*
                       catch_sel(iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter)*                                                          
                       catch_q(  iFleet, iMetier, ispp, 1,    iyr, iunit, iseason, iarea, iter)*
                       effshare( iFleet, iMetier,       1,    iyr, iunit, iseason, iarea, iter)*
                       effort(   iFleet,                1,    iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z));

                   d_=c_*discards_sel(iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter)/
                         catch_sel(   iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter);
           
                   l_=c_*landings_sel(iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter)/
                         catch_sel(   iFleet, iMetier, ispp, iage, iyr, iunit, iseason, iarea, iter);
             
                   catch_n(iFleet,iMetier,ispp,iage,iyr,iunit,iseason,iarea,iter)   =c_;
                   discards_n(iFleet,iMetier,ispp,iage,iyr,iunit,iseason,iarea,iter)=d_;
                   landings_n(iFleet,iMetier,ispp,iage,iyr,iunit,iseason,iarea,iter)=l_;
                   }
             }
       }
 } 

void fwd::project(adouble *x, adouble *func, int iYr, int iter)
   {
   int iunit  =1,
       iseason=1,
       iarea  =1;

   FLQuant2_adolc ad_f(f, NBiol(), iYr,   iYr,   iter);
   FLQuant2_adolc ad_n(n, NBiol(), iYr+1, iYr+1, iter);

   for (int ispp=1; ispp<=NBiol(); ispp++)
      {
      //-------------------- Control Stuff ----------------------//  
      for (int iage=minage(ispp); iage<=maxage(ispp); iage++)
         {
         ad_f(ispp,iage,iYr,iunit,iseason,iarea,iter) = 0.0;
      
         int i=0;         
         for (int iFleet=1; iFleet<=nfleet(); iFleet++)
            for (int iMetier=1; iMetier<=nmetier(); iMetier++)
               if (Ctrl.fix(iYr, iFleet, iMetier)) 
                  ad_f(ispp,iage,iYr,iunit,iseason,iarea,iter) = __max(0.0,catch_sel(iFleet, iMetier, ispp, iage, iYr, iunit, iseason, iarea, iter)*                                                          
                                                                           catch_q(  iFleet, iMetier, ispp, 1,    iYr, iunit, iseason, iarea, iter)*
                                                                           effshare( iFleet, iMetier,       1,    iYr, iunit, iseason, iarea, iter)*
                                                                           effort(   iFleet,                1,    iYr, iunit, iseason, iarea, iter));                                                          
               else
                  {
                  ad_f(ispp,iage,iYr,iunit,iseason,iarea,iter) = catch_sel(iFleet, iMetier, ispp, iage, iYr, iunit, iseason, iarea, iter)*
                                                                 catch_q(  iFleet, iMetier, ispp, 1,    iYr, iunit, iseason, iarea, iter)*
                                                                 effshare( iFleet, iMetier,       1,    iYr, iunit, iseason, iarea, iter)*
                                                                 x[i++];                                                          
                  
                  }

         //numbers-at-age next year
         if (iage < maxage(ispp)){
            ad_n(ispp,iage+1,iYr+1,iunit,iseason,iarea,iter)  = n(ispp,iage,iYr,iunit,iseason,iarea,iter)*exp(-ad_f(ispp,        iage,iYr,iunit,iseason,iarea,iter)-m(ispp,        iage,iYr,iunit,iseason,iarea,iter));
            }
         else if (iage == plusgrp(ispp))
            ad_n(ispp,iage,  iYr+1,iunit,iseason,iarea,iter) += n(ispp,iage,iYr,iunit,iseason,iarea,iter)*exp(-ad_f(ispp,maxage(ispp),iYr,iunit,iseason,iarea,iter)-m(ispp,maxage(ispp),iYr,iunit,iseason,iarea,iter));
         }
   
      int rec_yr = __min(__max(iYr-minage(ispp)+1,minyr(ispp)),maxyr(ispp));
    
      ad_n(ispp, minage(ispp),iYr+1,iunit,iseason,iarea,iter) = _sr.recruits(ispp,iYr,flc::SSB(ispp, rec_yr,iunit,iseason,iarea,iter),iter);
      }  
   
      //-------------------- Target Stuff ----------------------//  
      for (int iTrgt=1; iTrgt<=Trgt.n(iYr); iTrgt++)
         {
         FLRConst_Target quantity = Trgt.quantity(iYr,iTrgt);
         double rel               = Trgt.rel(     iYr,iTrgt);
         double min               = Trgt.min(     iYr,iTrgt,iter);
         double val               = Trgt.val(     iYr,iTrgt,iter);
         double max               = Trgt.max(     iYr,iTrgt,iter);
         double ispp              = Trgt.spp(     iYr,iTrgt);
 
         //min & max bounds should only occur if a target calculated in a previous step for that year 
         // target value relative to reference year
         if (!R_IsNA(rel) && rel<=iYr && rel>=minyr(ispp))
            { 
            double RelVal    =getVal(quantity, ispp, (int)rel, iunit, iseason, iarea, iter),
                   CurrentVal=getVal(quantity, ispp,      iYr, iunit, iseason, iarea, iter);

            if      (!R_IsNA(min) && CurrentVal<min*RelVal) 
               val = min*RelVal;
            else if (!R_IsNA(max) && CurrentVal>max*RelVal) 
               val = max*RelVal;
            else val = RelVal*val;
            }
         else if (R_IsNA(val)) // absolute target with max and min bounds 
            { 
            double CurrentVal=getVal(quantity, ispp, iYr, iunit, iseason, iarea, iter);

            if      (!R_IsNA(min)) val = (CurrentVal>min ? CurrentVal : min);
            else if (!R_IsNA(max)) val = (CurrentVal<max ? CurrentVal : max);
            else                   val =  CurrentVal; // this shouldn't happen   
            } 
      
         //NLSE's
         switch (quantity)
             {
   	       case FLRConst_F:
               func[iTrgt-1] = Fbar(                ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
   	       case FLRConst_Z:
	            func[iTrgt-1] = Zbar(                ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
   	       case FLRConst_FLandings:
	            func[iTrgt-1] = FbarLandings(        ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
   	       case FLRConst_FDiscards:
	            func[iTrgt-1] = FbarDiscards(        ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
    	       case FLRConst_SSB:
	            func[iTrgt-1] = SSB(            ad_n,ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
   	       case FLRConst_Biomass:
	            func[iTrgt-1] = computeStock(   ad_n,ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
	          case FLRConst_Catch:
	            func[iTrgt-1] = computeCatch(        ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
	          case FLRConst_Landings:
	            func[iTrgt-1] = computeLandings(     ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
	          case FLRConst_Discards:
	            func[iTrgt-1] = computeDiscards(     ad_f, ispp, iYr, iunit, iseason, iarea, iter) - val;
	          break;
 	          default:
	            func[iTrgt-1] = 1.0;
	          break;
             }
         }
   } 

bool control::Init(SEXP xCtrl, SEXP xAry, int niters)
    {
    SEXP dims;

    // Get target object
    if (!isMatrix(xCtrl) || !isNumeric(xCtrl))
      return false;

    dims = GET_DIM(xCtrl);

    if (LENGTH(dims) != 2 || INTEGER(dims)[1] != 10) 
       return false;

    double *Ctrl = NUMERIC_POINTER(xCtrl);
    
    nfleet =0;
    nmetier=0;

    minyear=(int)Ctrl[0+fwdControlPos_year*(int)INTEGER(dims)[0]];
    maxyear=(int)Ctrl[0+fwdControlPos_year*(int)INTEGER(dims)[0]];
    
   //year fleet metier  value rel.year rel.fleet rel.metier min max rel.bound
   int i = 0, j = 0, k =0;
   for (i=0; i<(int)INTEGER(dims)[0]; i++)
      {
      minyear=__min(minyear, (int)Ctrl[i+fwdControlPos_year*(int)INTEGER(dims)[0]]);
      maxyear=__max(maxyear, (int)Ctrl[i+fwdControlPos_year*(int)INTEGER(dims)[0]]);

      nmetier=__max(nmetier,(int)Ctrl[i+fwdControlPos_metier*(int)INTEGER(dims)[0]]); 
      nfleet =__max(nfleet, (int)Ctrl[i+fwdControlPos_fleet *(int)INTEGER(dims)[0]]); 
      }

   data = new double***[maxyear-minyear+1] - minyear;
   for (i = minyear; i<=maxyear; i++)
      {
      (data)[i] = new double**[nfleet] - 1;
      for (j=1; j<=nfleet; j++)
 	      {
         (data)[i][j] = new double*[nmetier] - 1; 
         for (k=1; k<=nmetier; k++)
// 	         (data)[i][j][k] = new double[10-3] - 3; 
 	         (data)[i][j][k] = new double[10] - 0; 
         }
      }
 
    int iyr,   
        ifleet, 
        imetier;

    for (i=0; i<(int)INTEGER(dims)[0]; i++)
       {
       iyr    = (int)Ctrl[i+fwdControlPos_year  *(int)INTEGER(dims)[0]]; 
       ifleet = (int)Ctrl[i+fwdControlPos_fleet *(int)INTEGER(dims)[0]]; 
       imetier= (int)Ctrl[i+fwdControlPos_metier*(int)INTEGER(dims)[0]]; 
  
       data[iyr][ifleet][imetier][fwdControlPos_val  ]    = Ctrl[i+fwdControlPos_val      *(int)INTEGER(dims)[0]];
       data[iyr][ifleet][imetier][fwdControlPos_relyear]  = Ctrl[i+fwdControlPos_relyear  *(int)INTEGER(dims)[0]];
       data[iyr][ifleet][imetier][fwdControlPos_relfleet] = Ctrl[i+fwdControlPos_relfleet *(int)INTEGER(dims)[0]];
       data[iyr][ifleet][imetier][fwdControlPos_relmetier]= Ctrl[i+fwdControlPos_relmetier*(int)INTEGER(dims)[0]];
       data[iyr][ifleet][imetier][fwdControlPos_min]      = Ctrl[i+fwdControlPos_min      *(int)INTEGER(dims)[0]];
       data[iyr][ifleet][imetier][fwdControlPos_max]      = Ctrl[i+fwdControlPos_max      *(int)INTEGER(dims)[0]];
       data[iyr][ifleet][imetier][fwdControlPos_relbound] = Ctrl[i+fwdControlPos_relbound *(int)INTEGER(dims)[0]];
       }

   AllocFlag = true;
   
   return(AllocFlag);
   }

control::control(void)
   {
   AllocFlag = false;
   
   minyear = 0;
   maxyear = 0;
   nfleet  = 0;
   nmetier = 0;
   }

bool control::unalloc(void)
   {
   if (!AllocFlag) return false;

   for (int i = minyear; i <= maxyear; i++) 
      {
      for (int j = 1; j <= nfleet; j++) 
   	   {
	      for (int k = 1; k <= nmetier; k++) 
   	       delete [] (data[i][j][k]+0);
   
   	   delete [] (data[i][j]+1);
   	   }
	   delete [] (data[i]+1);
      }               
   delete [] (data+minyear);
   
   AllocFlag = false;

   return true;
   }

control::~control(void)
   {
   unalloc();
   }

double control::value(int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return R_NaN;
   else
      return data[year][fleet][metier][fwdControlPos_val];
   }

double control::min(  int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return R_NaN;
   else
      return data[year][fleet][metier][fwdControlPos_min];
   }

double control::max(  int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return R_NaN;
   else
      return (data)[year][fleet][metier][fwdControlPos_max];
   }

int control::rel_year(  int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return (int)R_NaN;
   else
      return (int)data[year][fleet][metier][fwdControlPos_relyear];
   }

int control::rel_fleet( int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return (int)R_NaN;
   else
      return (int)data[year][fleet][metier][fwdControlPos_relfleet];
   }

int control::rel_metier(int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return (int)R_NaN;
   else
      return (int)data[year][fleet][metier][fwdControlPos_relmetier];
   }

bool control::rel_bound( int year, int fleet, int metier)
   {
   if (year>maxyear || year<minyear || fleet<1 || fleet>nfleet || metier<1 || metier>nmetier)
      return (bool)R_NaN;
   else
      return (bool)data[year][fleet][metier][fwdControlPos_relbound];
   }

bool target::Init(SEXP xTrgt, SEXP xAry, int niters)
    {
    SEXP TrgtDims;

    // Get Target object
    if (!isArray(xTrgt) || !isNumeric(xTrgt))
      return false;

    // Get Array object
    if (!isArray(xAry) || !isNumeric(xAry))
      return false;

    //Target data frame
    TrgtDims = GET_DIM(xTrgt);

    if (LENGTH(TrgtDims) != 2 || INTEGER(TrgtDims)[1] != 9) 
       return false;

    double *Trgt = NUMERIC_POINTER(xTrgt);
    
    _minyear=(int)Trgt[0+fwdTargetPos_year*(int)INTEGER(TrgtDims)[0]];
    _maxyear=(int)Trgt[0+fwdTargetPos_year*(int)INTEGER(TrgtDims)[0]];
    
    //year fleet metier  value rel.year rel.fleet rel.metier min max rel.bound
    int i = 0, iYr = 0;
    for (i=0; i<(int)INTEGER(TrgtDims)[0]; i++)
       {
       _minyear=__min(_minyear, (int)Trgt[i+fwdTargetPos_year*(int)INTEGER(TrgtDims)[0]]);
       _maxyear=__max(_maxyear, (int)Trgt[i+fwdTargetPos_year*(int)INTEGER(TrgtDims)[0]]);
       }
    
   _n  = new int[_maxyear-_minyear+1] - _minyear;
   for (iYr=_minyear; iYr<=_maxyear; iYr++)
      _n[iYr]=0;

   for (i=0; i<(int)INTEGER(TrgtDims)[0]; i++)
      {
      iYr=(int)Trgt[i+fwdTargetPos_year*(int)INTEGER(TrgtDims)[0]];
      _n[iYr]++;
      }
   
   data = new double**[_maxyear-_minyear+1] - _minyear;
   for (iYr = _minyear; iYr<=_maxyear; iYr++)
      {
      (data)[iYr] = new double*[_n[iYr]] - 1;
      for (i = 1; i<=_n[iYr]; i++)
         (data)[iYr][i] = new double[9-1] - 1;
      }

    int *in;
   
    in = new int[_maxyear-_minyear+1] - _minyear;
   
    for (iYr = _minyear; iYr<=_maxyear; iYr++)
       in[iYr]=0;

   _nrows = INTEGER(TrgtDims)[0];

    for (i=0; i<_nrows; i++)
       {
       iYr    = (int)Trgt[i+fwdTargetPos_year  *(int)INTEGER(TrgtDims)[0]]; 
      
       in[iYr]++;
  
       data[iYr][in[iYr]][fwdTargetPos_quantity] = Trgt[i+fwdTargetPos_quantity *(int)INTEGER(TrgtDims)[0]];
       data[iYr][in[iYr]][fwdTargetPos_min]      = Trgt[i+fwdTargetPos_min      *(int)INTEGER(TrgtDims)[0]];
       data[iYr][in[iYr]][fwdTargetPos_val]      = Trgt[i+fwdTargetPos_val      *(int)INTEGER(TrgtDims)[0]];
       data[iYr][in[iYr]][fwdTargetPos_max]      = Trgt[i+fwdTargetPos_max      *(int)INTEGER(TrgtDims)[0]];

       data[iYr][in[iYr]][fwdTargetPos_relyear]  = Trgt[i+fwdTargetPos_relyear  *(int)INTEGER(TrgtDims)[0]];
       data[iYr][in[iYr]][fwdTargetPos_spp]      = Trgt[i+fwdTargetPos_spp      *(int)INTEGER(TrgtDims)[0]];
       data[iYr][in[iYr]][fwdTargetPos_fleet]    = Trgt[i+fwdTargetPos_fleet    *(int)INTEGER(TrgtDims)[0]];
       data[iYr][in[iYr]][fwdTargetPos_metier]   = Trgt[i+fwdTargetPos_metier   *(int)INTEGER(TrgtDims)[0]];
       }
   
   delete [] (in + _minyear);

   AllocDataFlag = true;
   

   //Array
   SEXP AryDims = GET_DIM(xAry);

   if ((LENGTH( AryDims)    != 3)                    || 
       (INTEGER(AryDims)[0] != INTEGER(TrgtDims)[0]) || 
       (INTEGER(AryDims)[1] != 3)                    || 
       (INTEGER(AryDims)[2] != niters)) 
       return false;

   double *Ary = NUMERIC_POINTER(xAry);

   _min = new double**[_maxyear-_minyear+1] - _minyear;
   _val = new double**[_maxyear-_minyear+1] - _minyear;
   _max = new double**[_maxyear-_minyear+1] - _minyear;
   for (iYr = _minyear; iYr<=_maxyear; iYr++)
      {
      _min[iYr] = new double*[(int)n(iYr)] - 1;
      _val[iYr] = new double*[(int)n(iYr)] - 1;
      _max[iYr] = new double*[(int)n(iYr)] - 1;
      for (int iTrgt = 1; iTrgt<=n(iYr); iTrgt++)
         {
         _min[iYr][iTrgt] = new double[niters] - 1;
         _val[iYr][iTrgt] = new double[niters] - 1;
         _max[iYr][iTrgt] = new double[niters] - 1;
         }
      }

    i = 0;
    /*  
    int iter;
    for (iYr=_minyear; iYr<=_maxyear; iYr++)
       for (int iTrgt=1; iTrgt<=n(iYr); iTrgt++)
          for (iter=1; iter<=niters; iter++)
             _val[iYr][iTrgt][iter] = (Ary)[i++]; //(i-1+1*_nrows+3*_nrows*(iter-1))];

    for (iYr=_minyear; iYr<=_maxyear; iYr++)
       for (int iTrgt=1; iTrgt<=n(iYr); iTrgt++)
          for (iter=1; iter<=niters; iter++)
             _min[iYr][iTrgt][iter] = (Ary)[i++]; //(i-1+1*_nrows+3*_nrows*(iter-1))];

    for (iYr=_minyear; iYr<=_maxyear; iYr++)
       for (int iTrgt=1; iTrgt<=n(iYr); iTrgt++)
          for (iter=1; iter<=niters; iter++)
             _max[iYr][iTrgt][iter] = (Ary)[i++]; //(i-1+2*_nrows+3*_nrows*(iter-1))];
*/

    for (int iter=1; iter<=niters; iter++)
       {     
       for (iYr=_minyear; iYr<=_maxyear; iYr++)
          for (int iTrgt=1; iTrgt<=n(iYr); iTrgt++)
            _min[iYr][iTrgt][iter] = (Ary)[i++]; 
       for (iYr=_minyear; iYr<=_maxyear; iYr++)
          for (int iTrgt=1; iTrgt<=n(iYr); iTrgt++)
            _val[iYr][iTrgt][iter] = (Ary)[i++]; 
       for (iYr=_minyear; iYr<=_maxyear; iYr++)
          for (int iTrgt=1; iTrgt<=n(iYr); iTrgt++)
            _max[iYr][iTrgt][iter] = (Ary)[i++]; 
       }


   AllocAryFlag = true;

   return(AllocDataFlag && AllocAryFlag);
   }

target::target(void)
   {
   AllocDataFlag = 
   AllocAryFlag  = false;
   
   _minyear = 0;
   _maxyear = 0;
   }

void target::unalloc(void)
   {
   if (AllocDataFlag) 
      {
      for (int i = _minyear; i <= _maxyear; i++) 
         {
         for (int j = 1; j <= _n[i]; j++) 
            delete [] (data[i][j]+1);
   
         delete [] (data[i]+1);
         }               
      delete [] (data+_minyear);

      AllocDataFlag = false;
      }

   if (AllocAryFlag) 
      {
      for (int i = _minyear; i <= _maxyear; i++) 
         {
         for (int j = 1; j <= _n[i]; j++) 
            { 
            delete [] (_min[i][j]+1);
            delete [] (_max[i][j]+1);
            delete [] (_val[i][j]+1);
            }

         delete [] (_min[i]+1);
         delete [] (_max[i]+1);
         delete [] (_val[i]+1);
         }               
      delete [] (_min+_minyear);
      delete [] (_max+_minyear);
      delete [] (_val+_minyear);

      delete [] (  _n+_minyear);
   
	  AllocDataFlag = false;
      }

   }

target::~target(void)
   {
   unalloc();
   }

int target::spp(   int year, int i)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return (int)R_NaN;
   else
      return (int)data[year][i][fwdTargetPos_spp];
   }

int target::fleet( int year, int i)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return (int)R_NaN;
   else
      return (int)data[year][i][fwdTargetPos_fleet];
   }

int target::metier(int year, int i)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return (int)R_NaN;
   else
      return (int)data[year][i][fwdTargetPos_metier];
   }

double target::val(int year, int i, int iter)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return R_NaN;
   else if (iter<=0 || !AllocAryFlag)
      return data[year][i][fwdTargetPos_val];
   else
      return _val[year][i][iter];
   }

double target::min(int year, int i, int iter)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return R_NaN;
   else if (iter<=0 || !AllocAryFlag)
      return data[year][i][fwdTargetPos_min];
   else
      return _min[year][i][iter];
   }

double target::max(int year, int i, int iter)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return R_NaN;
   else if (iter<=0 || !AllocAryFlag)
      return data[year][i][fwdTargetPos_max];
   else
      return _max[year][i][iter];
   }
   
FLRConst_Target target::quantity(int year, int i)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return (FLRConst_Target)(int)R_NaN;
   else
      return (FLRConst_Target)(int)data[year][i][fwdTargetPos_quantity];
   }
   
double target::rel(int year, int i)
   {
   if (year>_maxyear || year<_minyear || i<1 || i>_n[year])
      return R_NaN;
   else
      return data[year][i][fwdTargetPos_relyear];
   }

bool control::fix(int year, int fleet, int metier)
   {
   if (!R_IsNA(value(year,fleet,metier)))
      return true;   
   else
      return false;
   }

bool control::fit(int year, int fleet, int metier)
   {
   if (R_IsNA(value(year,fleet,metier)))
      return true;   
   else
      return false;
   }

int control::n(int year)
   {
   int _n=0;
   for (int iFleet=1; iFleet<=nfleet; iFleet++)
      for (int iMetier=1; iMetier<=nmetier; iMetier++)
         if (fit(year, iFleet, iMetier))
             _n++;

   return(_n);
   }
