#include "FLash-VPA.h"

adouble func(adouble& f,adouble m, adouble c, adouble n)
    {
    //deriv((y ~ catch*(f+m)/(n*(1-exp(-f-m)))-f), c("f"), func = TRUE)
    adouble expr2 = c * (f + m),
            expr5 = exp(-f - m),
            expr7 = n * (1.0 - expr5),
            value = expr2/expr7 - f,
            grad  = c/expr7 - expr2 * (n * expr5)/pow(expr7,2.0) - 1.0;

    f-=value/grad;

    return(value);
    }
 
FLashVPA::FLashVPA(SEXP x)
   {
   InitFlag    = false;
   
   if (isFLStock(x) && !InitFlag)
      Init(x);
   }

bool FLashVPA::isFLashVPA(SEXP x)
   {
   const char *s = CHAR(STRING_ELT(GET_CLASS(x), 0));

   return strcmp(s, "FLashVPA")==0;
   }

void FLashVPA::Init(SEXP x)
   {

   SEXP      range = GET_SLOT(x, install("range"));
   SEXP RangeNames = GET_NAMES(range);
	  
   int n = length(RangeNames);

   minage   = maxage = 
   minyr    = maxyr  = 
   nunits   = 
   nseasons = 
   nareas   =
   niters   = 1;

   int i;

   for (i=0; i<n; i++)
      {
      const char *s = CHAR(VECTOR_ELT(RangeNames, i));

      if (      strcmp(s, "min")==0 || strcmp(s, "minage")==0 || strcmp(s, "minquant")==0)
         minage     = (short)(REAL(range)[i]);
      else  if (strcmp(s, "max")==0 || strcmp(s, "maxage")==0 || strcmp(s, "maxquant")==0)
         maxage     = (short)(REAL(range)[i]);
      else  if (strcmp(s, "plusgroup")==0)
         if (R_IsNA(REAL(range)[i]))
            FlagPlusGrp = false;
         else 
            {
            plusgrp     = (short)(REAL(range)[i]);
            FlagPlusGrp = true;
            }
        else  if (strcmp(s, "minyear")==0)
           minyr    = (short)(REAL(range)[i]);
        else  if (strcmp(s, "maxyear")==0)
           maxyr    = (short)(REAL(range)[i]);
      }

   Catch.Init(GET_SLOT(x, install("catch.n")));        
   M.Init(    GET_SLOT(x, install("m")));               
   F.Init(    GET_SLOT(x, install("harvest")));               
   N.Init(    GET_SLOT(x, install("stock.n")));               
   
   niters = M.niters();

   N_ad.Init(    N,     N.minyr(),     N.maxyr(),     N.niters());
   F_ad.Init(    F,     F.minyr(),     F.maxyr(),     F.niters());
   Catch_ad.Init(Catch, Catch.minyr(), Catch.maxyr(), Catch.niters());
   }

FLashVPA::~FLashVPA(void)      
   {
   ; //unalloc();
   }                               

SEXP FLashVPA::Return(void)
   {
   SEXP x, Range;

   PROTECT(x = NEW_OBJECT(MAKE_CLASS("FLAssess")));
   Range     = PROTECT(NEW_NUMERIC(5)); 
   
   REAL(Range)[0] = minage;
   REAL(Range)[1] = maxage;
   REAL(Range)[2] = plusgrp;
   REAL(Range)[3] = minyr;
   REAL(Range)[4] = maxyr;
       
   SET_SLOT(x, install("stock.n"), N_ad.Return());
   SET_SLOT(x, install("harvest"), F_ad.Return());
   SET_SLOT(x, install("catch.n"), Catch_ad.Return());
   
   UNPROTECT(2);

   return x;
   }

bool FLashVPA::SepVPA(SEXP xControl, SEXP xRefHarvest)
    {
    int iage, iyr, iyrcls, Iters =0, N_age, _maxage;
                                    
    adouble F1, Z, Z1, CumZ,
            Error,
            Error_reduction;
;   
    adouble *Selectivity, *OverallF, *CohortSize, **Residuals, **Observed, **H,
            *RsdlyrTotals, *RsdlageTotals;

    if (FlagPlusGrp) _maxage = maxage-1; else _maxage = maxage;

   Selectivity   = new adouble[_maxage -minage+1] - minage;
   RsdlageTotals = new adouble[_maxage -minage+1] - minage;
   OverallF      = new adouble[_maxage -minage+1] - minage;
   RsdlyrTotals  = new adouble[_maxage -minage+1] - minage;
   CohortSize    = new adouble[(maxyr-minage) - (minyr-_maxage) + 1]  - (minyr-_maxage);
    
   Residuals     =new adouble*[(_maxage -minage+1)] - minage;
   Observed      =new adouble*[(_maxage -minage+1)] - minage;
   H             =new adouble*[(_maxage -minage+1)] - minage;
    
   for (int iter=1; iter<=niters; iter++)
      {
      for (iage=minage; iage <= _maxage; iage++)
        {
        Residuals[iage] =new adouble[(maxyr -minyr +1)] - minyr;
        Observed[iage]  =new adouble[(maxyr -minyr +1)] - minyr;
        H[iage]         =new adouble[(maxyr -minyr +1)] - minyr;
        }

      int    refage = (short) INTEGER(GET_SLOT(xControl,install("sep.age")))[0];
      adouble refsep = _max(REAL(GET_SLOT(xControl,install("sep.sel")))[0],0.0);
      adouble refF   = _max(REAL(xRefHarvest)[0],0.0);

      Selectivity[refage] = refsep;
      Selectivity[_maxage] = 1.0;
      OverallF[maxyr]     = refF;

      //Initialise
      for (iage=minage; iage < _maxage; iage++)
          Selectivity[iage] = Selectivity[refage];
      for (iyr=minyr; iyr < maxyr; iyr++)
          OverallF[iyr] = OverallF[maxyr];

      //Calc Log Catch Ratios Residuals
      for (iage=minage, Error = 0.0f; iage < _maxage; iage++)
         for (iyr=minyr; iyr < maxyr; iyr++)
            {
            Observed[iage][iyr]   = log(Catch_ad(iage+1, iyr+1)/Catch_ad(iage,   iyr,1,1,1,iter));

            F_ad(iage,iyr,1,1,1,iter) = OverallF[iyr]*Selectivity[iage];
            F1          = OverallF[iyr+1]*Selectivity[iage+1];
            Z           = F_ad(iage,iyr,1,1,1,iter) + M(iage,   iyr,1,1,1,iter);
            Z1          = F1 + M(iage+1, iyr+1);

            Residuals[iage][iyr]  = Observed[iage][iyr]
                                      - log(F1*Z*(1-exp(-Z1))*exp(-Z)/(F_ad(iage,iyr,1,1,1,iter)*Z1*(1-exp(-Z))));

            Error += Residuals[iage][iyr]*Residuals[iage][iyr];
            }

      //Newton Rhapson to Solve for Overall F and Selectivity
      do  {
          //Sum Residuals by iyrs over iages
          for (iyr=minyr; iyr < maxyr; iyr++)
              for (iage=minage, RsdlyrTotals[iyr] = 0.0f; iage < _maxage; iage++)
                  RsdlyrTotals[iyr] += Residuals[iage][iyr];

          //Sum Residuals by iages over iyrs
          for (iage=minage; iage < _maxage; iage++)
              for (iyr=minyr, RsdlageTotals[iage] = 0.0f; iyr < maxyr; iyr++)
                  RsdlageTotals[iage] += Residuals[iage][iyr];

          //Calc New Overall F
          for (iyr=minyr; iyr < maxyr; iyr++)
              OverallF[iyr] /= exp(RsdlyrTotals[iyr] / (2.0f*(_maxage - minage + 1)));

          //Calc New Selectivity
          for (iage=minage; iage < _maxage; iage++)
              if (iage != refage)
                 Selectivity[iage] /= exp(RsdlageTotals[iage] / (2.0f*(maxyr - minyr + 1)));

          //Calc Log Catch Ratios Residuals
          for (iage=minage, Error_reduction = Error, Error = 0.0f; iage < _maxage; iage++)
              for (iyr=minyr; iyr < maxyr; iyr++)
                  {
                  F_ad(iage,iyr,1,1,1,iter) = OverallF[iyr]*Selectivity[iage];
                  F1          = OverallF[iyr+1]*Selectivity[iage+1];
                  Z           = F_ad(iage,iyr,1,1,1,iter) + M(iage,   iyr,1,1,1,iter);
                  Z1          = F1 + M(iage+1, iyr+1);

                  Residuals[iage][iyr]  = Observed[iage][iyr]
                                            - log(F1*Z*(1-exp(-Z1))*exp(-Z)/(F_ad(iage,iyr,1,1,1,iter)*Z1*(1-exp(-Z))));

                  Error += Residuals[iage][iyr]*Residuals[iage][iyr];
                  }

          Error_reduction -=Error;
          Iters++;
          }
      while (Error_reduction >= 1.0e-20f && Iters < VPA_ITS);

      //Calc Cohort Size
      for (iyrcls  = minyr - _maxage; iyrcls <= maxyr  - minage; iyrcls++)
          {
          CohortSize[iyrcls] = 0.0f;
          N_age              = 0;
          CumZ               = 0.0f;

          for (iage  = _max(minage, minyr - iyrcls); iage <= _min(_maxage, maxyr - iyrcls); iage++)
              {
              Z = OverallF[iyrcls + iage]*Selectivity[iage] + M(iage, iyrcls + iage);
              H[iage][iyrcls+iage] =  OverallF[iyrcls + iage]*Selectivity[iage]/Z
                                          *(1.0f-exp(-1*Z))*exp(-CumZ);


              //OverallF[iyrcls + iage] = 1.0;
              //Selectivity[iage]       = 1.0;

              //H[iage][iyrcls+iage] = 1.0;

              CohortSize[iyrcls] += log(Catch_ad(iage, iyrcls+iage)) - log(H[iage][iyrcls+iage]);
              CumZ += Z;
              N_age++;
              }
           CohortSize[iyrcls] = exp(CohortSize[iyrcls] / N_age);
           }

       //Calculate Fs and Ns
       for (iyrcls  = minyr - _maxage; iyrcls <= maxyr - minage;  iyrcls++)
          for (iage  = _max(minage, minyr - iyrcls); iage <= _min(_maxage, maxyr - iyrcls); iage++)
              {
              F_ad(iage,iyrcls+iage) = OverallF[iyrcls + iage]*Selectivity[iage];
              if (iage == _max(minage, minyr - iyrcls))
                 N_ad(iage,iyrcls+iage) = CohortSize[iyrcls];
              else
                 N_ad(iage,iyrcls+iage) = (N_ad(iage-1,iyrcls+iage-1)*exp(-F_ad(iage-1,iyrcls+iage-1)-M(iage-1, iyrcls+iage-1)));
              }

       if (FlagPlusGrp)
          for (iyr=minyr; iyr <= maxyr; iyr++)
              {
              F_ad(maxage, iyr,1,1,1,iter) = F_ad(_maxage, iyr,1,1,1,iter);
              N_ad(maxage, iyr,1,1,1,iter) = Catch_ad(maxage,iyr,1,1,1,iter)*F_ad(maxage,iyr,1,1,1,iter)/((F_ad(maxage,iyr,1,1,1,iter)+M(maxage,iyr,1,1,1,iter)*exp(1.0-exp(-F_ad(maxage,iyr,1,1,1,iter)-M(maxage, iyr,1,1,1,iter)))));
              }

        for (iyr  = minyr; iyr <= maxyr;  iyr++)
           for (iage  = minage; iage <= maxage; iage++)
              Catch_ad(iage,iyr,1,1,1,iter) = N_ad(iage,iyr,1,1,1,iter)*F_ad(iage,iyr,1,1,1,iter)/(F_ad(iage,iyr,1,1,1,iter)+M(iage, iyr,1,1,1,iter))*(1.0-exp(-F_ad(iage,iyr,1,1,1,iter)-M(iage, iyr,1,1,1,iter)));
       }

      delete [] (Selectivity   +minage);
      delete [] (OverallF      +minyr);
      delete [] (RsdlageTotals +minage);
      delete [] (RsdlyrTotals  +minyr);
      delete [] (CohortSize    +(minyr-_maxage));
                                    
      for (iage=minage; iage <= _maxage; iage++)
         {
         delete [] (Residuals[iage] + minyr);
         delete [] (Observed[ iage] + minyr);
         delete [] (H[        iage] + minyr);
         }
      delete [] (Residuals+minage);
      delete [] (Observed +minage);
      delete [] (H        +minage);
      
      return TRUE;
    }                         

void FLashVPA::FratioFunc(adouble FRatio, adouble F, double M, double M2, adouble C, adouble C2, adouble N, adouble *value, adouble *grad)
   {
   //deriv(~(F*(exp(F+M)-1.0)/(F+M))*(N-C2*(FRatio*F+M2)*(exp(-F*FRatio-M2))/(F*FRatio*(1-exp(-F*FRatio-M2)))) - C, "F")

    adouble expr1  = F + M;
    adouble expr2  = exp(expr1);
    adouble expr3  = expr2 - 1.0;
    adouble expr4  = F * expr3;
    adouble expr5  = expr4/expr1;
    adouble expr8  = C2 * (FRatio * F + M2);
    adouble expr12 = exp(-F * FRatio - M2);
    adouble expr13 = expr8 * expr12;
    adouble expr14 = F * FRatio;
    adouble expr15 = 1.0 - expr12;
    adouble expr16 = expr14 * expr15;
    adouble expr18 = N - expr13/expr16;
    adouble expr30 = expr12 * FRatio;
    *value  = expr5 * expr18 - C;
    *grad   = ((expr3 + F * expr2)/expr1 - expr4/(expr1*expr1)) * 
                expr18 - expr5 * ((C2 * FRatio * expr12 - expr8 * 
                expr30)/expr16 - expr13 * (FRatio * expr15 + expr14 * 
                expr30)/(expr16*expr16));
   }


bool FLashVPA::VPA(SEXP xFitPlusGroup, SEXP xFRatio, SEXP xFRatioFlag)
   {
   bool FlagFitPlusGrp;
  
   FLBool   FRatioFlag(xFRatioFlag);
   FLVector FRatio(    xFRatio);
   FlagFitPlusGrp =  (bool)LOGICAL(xFitPlusGroup)[0];

   short iter, iyr;
   
   adouble **FRatio_ad;

   FRatio_ad = new adouble*[niters] - 1;
    
   for (iter=1; iter<=niters; iter++)
      FRatio_ad[iter] =new adouble[(maxyr -minyr +1)] - minyr;
   
   for (iter=1; iter<=niters; iter++)
      for (iyr = minyr; iyr<=maxyr; iyr++)
         FRatio_ad[iter][iyr] = FRatio(iyr);

   for (iter=1; iter<=niters; iter++)
      RunVPA(iter,FRatio_ad[iter],FRatioFlag,FlagFitPlusGrp);
        
   for (iter=1; iter<=niters; iter++)
      delete [] (FRatio_ad[iter] + minyr);
   delete [] (FRatio_ad + 1);
    
   return TRUE;   
   }   

void FLashVPA::RunVPA(int iter, adouble *FRatio_ad, FLBool FRatioFlag, bool FlagFitPlusGrp)
   {   
   adouble f, dfdx, 
           Z, value, grad,         
           t, t1;

   int iage, iyr, iyrcls, Iters;

   //last year
   for (iage = maxage; iage >= minage; iage--)
     {
     Z  = F_ad(iage, maxyr,1,1,1,iter)+M(iage, maxyr,1,1,1,iter);
     t  = N_ad(iage, maxyr,1,1,1,iter);
     t1 = Catch_ad(iage, maxyr,1,1,1,iter);
     N_ad(iage, maxyr,1,1,1,iter) = Catch_ad(iage, maxyr,1,1,1,iter)*Z/(F_ad(iage, maxyr,1,1,1,iter)*(1.0-exp(-Z)));
     }

   for (iyr = maxyr; iyr >= minyr; iyr--)
     // Given terminal Fs, don't fit plus group
     if (!FRatioFlag(iyr) && !FlagFitPlusGrp)
        {
        Z                         = F_ad(maxage, iyr,1,1,1,iter)+M(maxage, iyr,1,1,1,iter);
        N_ad(maxage, iyr,1,1,1,iter) = Catch_ad(maxage, iyr,1,1,1,iter)*Z/(F_ad(maxage, iyr,1,1,1,iter)*(1.0-exp(-Z)));

        // Plus group F = last true age F or
        if (FlagPlusGrp)
           {
           F_ad(maxage-1, iyr,1,1,1,iter) = F_ad(maxage, iyr,1,1,1,iter);
           Z                           = F_ad(maxage-1, iyr,1,1,1,iter)+M(maxage-1, iyr,1,1,1,iter);
           N_ad(maxage-1, iyr,1,1,1,iter) = Catch_ad(maxage-1, iyr,1,1,1,iter)*Z/(F_ad(maxage-1, iyr,1,1,1,iter)*(1.0-exp(-Z)));
           }
        }
     // Terminal ages fitted using F ratio, with plusgroup
     else if (FRatioFlag(iyr) && FlagPlusGrp && iyr>minyr)
        {
        Z                             = F_ad(maxage, iyr-1,1,1,1,iter) + M(maxage, iyr-1,1,1,1,iter);
        N_ad(maxage,   iyr-1,1,1,1,iter) = Catch_ad(maxage,iyr-1,1,1,1,iter)*Z/(F_ad(maxage,iyr-1,1,1,1,iter)*(1.0-exp(-Z)));
        F_ad(maxage-1, iyr-1,1,1,1,iter) = F_ad(maxage-1, iyr,1,1,1,iter);

        Iters = 0;
        do
           {
           Iters++;

           FratioFunc(FRatio_ad[iyr], F_ad(maxage-1,iyr-1,1,1,1,iter), M(maxage-1,iyr-1,1,1,1,iter), M(maxage,iyr-1,1,1,1,iter), Catch_ad(maxage-1,iyr-1,1,1,1,iter), Catch_ad(maxage,iyr-1,1,1,1,iter), N_ad(maxage,iyr,1,1,1,iter), &value, &grad);

           //Newton Rhapson
           F_ad(maxage-1, iyr-1,1,1,1,iter) = F_ad(maxage-1, iyr-1,1,1,1,iter) - value/grad;
           Z                  = F_ad(maxage-1, iyr-1,1,1,1,iter) + M(maxage-1, iyr-1,1,1,1,iter);
           N_ad(maxage-1, iyr-1,1,1,1,iter) = Catch_ad(maxage-1, iyr-1,1,1,1,iter)*Z/(F_ad(maxage-1, iyr-1,1,1,1,iter)*(1.0-exp(-Z)));
           }
        while (fabs(value) >= VPA_TOL && Iters <= VPA_ITS);


        F_ad(maxage,   iyr-1,1,1,1,iter) = F_ad(maxage-1, iyr-1,1,1,1,iter)*FRatio_ad[iyr];
        Z                             = F_ad(maxage, iyr-1,1,1,1,iter) + M(maxage, iyr-1,1,1,1,iter);
        N_ad(maxage, iyr-1,1,1,1,iter)   = Catch_ad(maxage, iyr-1,1,1,1,iter)*Z/(F_ad(maxage, iyr-1,1,1,1,iter)*(1.0-exp(-Z)));
        }
     // Given terminal Fs, fitted plusgroup
     else if (!FRatioFlag(iyr) && FlagPlusGrp && iyr>minyr)
        {
        Iters = 0;

        adouble expr1, expr2, expr3, expr4, Z, PlusGrpN;

        Z                             = F_ad(maxage, iyr-1,1,1,1,iter) + M(maxage, iyr-1,1,1,1,iter);
        N_ad(maxage,   iyr-1,1,1,1,iter) = Catch_ad(maxage,iyr-1,1,1,1,iter)*Z/(F_ad(maxage,iyr-1,1,1,1,iter)*(1.0-exp(-Z)));
        F_ad(maxage-1, iyr-1,1,1,1,iter) = F_ad(maxage-1, iyr,1,1,1,iter);

        do
           {
           Iters++;

           PlusGrpN =N_ad(maxage, iyr,1,1,1,iter) - N_ad(maxage, iyr-1,1,1,1,iter)*exp(-F_ad(maxage,iyr-1,1,1,1,iter)-M(maxage, iyr-1,1,1,1,iter));

           //deriv(~ Fval*(exp(Fval+Mval)-1.0)/(Fval+Mval)*PlusGrpN-Cval, "Fval")
           expr1 = F_ad(maxage-1,iyr-1,1,1,1,iter) + M(maxage-1,iyr-1,1,1,1,iter);
           expr2 = exp(expr1);
           expr3 = expr2 - 1.0;
           expr4 = F_ad(maxage-1,iyr-1,1,1,1,iter) * expr3;
           value = expr4/expr1 * PlusGrpN - Catch_ad(maxage-1,iyr-1,1,1,1,iter);
           grad  = ((expr3 + F_ad(maxage-1,iyr-1,1,1,1,iter) * expr2)/expr1 - expr4/(expr1*expr1)) * PlusGrpN;

           //Newton Rhapson
           F_ad(maxage-1, iyr-1,1,1,1,iter) = (F_ad(maxage-1, iyr-1,1,1,1,iter) - value/grad);
           Z                    = F_ad(maxage-1, iyr-1,1,1,1,iter) + M(maxage-1, iyr-1,1,1,1,iter);
           N_ad(maxage-1, iyr-1,1,1,1,iter) = Catch_ad(maxage-1, iyr-1,1,1,1,iter)*Z/(F_ad(maxage-1, iyr-1,1,1,1,iter)*(1.0-exp(-Z)));
           }
        while (fabs(value) >= VPA_TOL && Iters <= VPA_ITS);
        }
     // Terminal ages fitted using F ratio, plusgroup not fitted
     else if (FRatioFlag(iyr) && FlagPlusGrp && !FlagFitPlusGrp && iyr>minyr)
        {
        ;
        }

   //Now go for the year classes
   for (iyrcls = maxyr - minage - 1; iyrcls >= minyr - maxage - 1; iyrcls--)
     for (iage = _min(maxage - (FlagPlusGrp /*&& !FRatioFlag(iyrcls+iage)*/ ? 2 : 1), maxyr - iyrcls - 1); iage >= _max(minage, minyr - iyrcls); iage--)
        {
        iyr = iyrcls+iage;

        Iters = 0;

        adouble t = N_ad(iage+1,iyr+1,1,1,1,iter)*exp(M(iage, iyr,1,1,1,iter));
        N_ad(iage,iyr,1,1,1,iter) = fmax(0.0, t + Catch_ad(iage, iyr,1,1,1,iter) * exp(0.5*M(iage, iyr,1,1,1,iter)));

        do
           {
           Iters++;
           //do Newton Raphson to estimate N
           f    = Calcf(M(iage, iyr,1,1,1,iter), Catch_ad(iage, iyr,1,1,1,iter), N_ad(iage,iyr,1,1,1,iter), N_ad(iage+1,iyr+1,1,1,1,iter));
           dfdx = Calcdfdx(M(iage, iyr,1,1,1,iter), N_ad(iage,iyr,1,1,1,iter), N_ad(iage+1,iyr+1,1,1,1,iter));

           //calc N
           N_ad(iage,iyr,1,1,1,iter) = NewtonRhapson(N_ad(iage,iyr,1,1,1,iter), f, dfdx);
           }
        while (fabs(f) >= VPA_TOL && Iters <= VPA_ITS);

        //calc F at iage
        F_ad(iage,iyr,1,1,1,iter) =  fmax(0.0,-log(N_ad(iage+1,iyr+1,1,1,1,iter)/N_ad(iage,iyr,1,1,1,iter)) - M(iage, iyr,1,1,1,iter));
        }
   }
