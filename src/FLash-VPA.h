#ifndef _INC_FLashVPA
#define _INC_FLashVPA

#include "flc_adolc.h"
      
#define VPA_TOL     1e-20
#define SEPVPA_TOL  1e-40
#define VPA_ITS     200     

#define _max(a,b) ((a)>(b)?(a):(b))
#define _min(a,b) ((a)<(b)?(a):(b))

class FLashVPA
{
public:        
   bool InitFlag;

   FLashVPA(SEXP);     
  ~FLashVPA(void);
                                                        
   void Init(SEXP);

   bool isFLashVPA(SEXP);

   void RunVPA(int iter, adouble *FRatio, FLBool FRatioFlag, bool FlagFitPlusGrp);
   bool VPA(SEXP,SEXP,SEXP);
   bool SepVPA(SEXP, SEXP);

   FLQuant_adolc N_ad, 
                 F_ad, 
                 Catch_ad;

   SEXP Return(void);

protected: 
    int minage, maxage,  plusgrp,  
        minyr,  maxyr,
        nunits, nseasons, nareas, niters;

   bool FlagPlusGrp;

   FLQuant  Catch,
            M,
            N,
            F;   

   inline adouble Calcf(double M, adouble Catch, adouble N, adouble N1)  {return (Catch - (1.0 - M/(log(N)-log(N1)))*(N-N1));}
   inline adouble Calcdfdx(double M, adouble N, adouble N1)              {return  (-1.0 -((log(N)-log(N1))*M -M*(N-N1)/N)/log(2*N/N1));}
   inline adouble NewtonRhapson(adouble x, adouble f, adouble dfdx)      {return (x - f/dfdx);}

   void FratioFunc(adouble FRatio, adouble F, double M, double M2, adouble C, adouble C2, adouble N, adouble *value, adouble *grad);

};


#endif /* _INC_FLashVPA */
