#include <math.h>
#include "fwd.h"
#include "fwdFLStock.h"
#include "FLash-VPA.h"

#ifdef WIN32
   #define SEXPDLLExport __declspec(dllexport) SEXP __cdecl    
#else
   #define SEXPDLLExport SEXP    
#endif

extern "C" SEXPDLLExport CalcF(SEXP xM, SEXP xCatch, SEXP xN)
   {
   FLQuant m(xM);
   FLQuant c(xCatch);
   FLQuant n(xN);
   FLQuant f(m.minquant(),m.maxquant(),m.minyr(),m.maxyr(),m.nunits(),m.nseasons(),m.nareas(),m.niters(),0.0);

   double _f = 0.1;

   int iAge, iYear, iUnit, iSeason, iArea, iIter;

   for (iIter = 1; iIter<=m.niters(); iIter++)
	   for (iArea = 1; iArea <= m.nareas(); iArea++)
		   for (iSeason = 1; iSeason <= m.nseasons(); iSeason++)
    		 for (iUnit = 1; iUnit <= m.nunits(); iUnit++)
    	 		 for (iYear = m.minyr(); iYear <= m.maxyr(); iYear++)
		   		   for (iAge = m.minquant(); iAge <= m.maxquant(); iAge++)
				    	 {
               _f = 0.1;
               int iter=0;
               while (fabs(F_func(&_f,m(iAge,iYear,iUnit,iSeason,iArea,iIter),c(iAge,iYear,iUnit,iSeason,iArea,iIter),n(iAge,iYear,iUnit,iSeason,iArea,iIter)))>1e-12 && iter++<50)
                 ;

					     f(iAge,iYear,iUnit,iSeason,iArea,iIter)=_f;
               }					
   return f.Return();
   }

extern "C" SEXPDLLExport VPA_ad(SEXP xStock, SEXP xFitPlusGroup, SEXP xFratio, SEXP xFratioFlag)
   {
   //Input VPA Control
   FLashVPA VPA(xStock);       
   
   //Run VPA
   VPA.VPA(xFitPlusGroup,xFratio,xFratioFlag);
                 
   //Display Ns and Fs 
   return VPA.Return();  
   }                 

extern "C" SEXPDLLExport SepVPA_ad(SEXP xStock, SEXP xControl, SEXP xRefHarvest)
   {
   //Input VPA Control
   FLashVPA VPA(xStock);       
   
   //Run Seperable VPA
   VPA.SepVPA(xControl, xRefHarvest);
                 
   //Display Ns and Fs 
   return VPA.Return();  
   }                 


// Uses fwdFLStock
extern "C" SEXPDLLExport fwd_adolc_FLStock(SEXP xStk,SEXP xTrgt,SEXP xAry,SEXP xYrs,SEXP xSRModel,SEXP xSRParam,SEXP xSRResiduals,SEXP xMult,SEXP xAvail,SEXP xMaxF) 
    {
	  fwdStk fwd;

    fwd.Init(xStk, xYrs, xSRModel, xSRParam, xSRResiduals, xMult, xAvail,xMaxF);    

    //MaxFBar = 2.0;

	  return fwd.run(xTrgt, xAry);}


extern "C" SEXPDLLExport fwd_adolc_FLBiol(SEXP xBiols, SEXP xFleets, SEXP xTrgt, SEXP xAryTrgt, SEXP xCtrl, SEXP xAryCtrl, SEXP xYrs, SEXP xDims, SEXP xSRModel,SEXP xSRParam,SEXP xSRResiduals,SEXP xMult,SEXP xAvail)    
    {
    SEXP ReturnVal = R_NilValue;
    SEXP Rerror = NEW_LOGICAL(0);
    
   int nparam = NElemList(xSRParam);
   if (nparam<1) 
      return Rerror;

   SEXP SRParam = PROTECT(VECTOR_ELT(xSRParam, 0));
   
   fwd fwd(xBiols, xFleets, xYrs, xDims, xSRModel, SRParam, xSRResiduals, xMult);

   fwd.run(xTrgt, xCtrl, xAryTrgt, xAryCtrl, xYrs);

   PROTECT(ReturnVal = allocVector(VECSXP,6));

    //SET_VECTOR_ELT(ReturnVal, 0, fwd.ReturnBiol());
    //SET_VECTOR_ELT(ReturnVal, 1, fwd.ReturnFleet());
    //SET_VECTOR_ELT(ReturnVal, 0, fwd.ReturnStock());

    SET_VECTOR_ELT(ReturnVal, 0, fwd.landings_n.Return(1,1,1));
    SET_VECTOR_ELT(ReturnVal, 1, fwd.discards_n.Return(1,1,1));
    SET_VECTOR_ELT(ReturnVal, 2, fwd.effort.Return(1));
 
    SET_VECTOR_ELT(ReturnVal, 3, fwd.n.Return(1));
    SET_VECTOR_ELT(ReturnVal, 4, fwd.f.Return(1));
    SET_VECTOR_ELT(ReturnVal, 5, fwd.catch_n.Return(1,1,1));
    
    UNPROTECT(2);

    return ReturnVal;
    }


// Uses fwd
extern "C" SEXPDLLExport _fwd_adolc_FLStock(SEXP xFLStock,SEXP xTrgt,SEXP xAry,SEXP xCtrl,SEXP xYrs,SEXP xSRModel,SEXP xSRParam,SEXP xSRResiduals,SEXP xMult,SEXP xAvail) 
    {
    SEXP ReturnVal = R_NilValue;
    
    fwd fwd(xFLStock, xYrs, xSRModel, xSRParam, xSRResiduals, xMult);

    //fwd.run(xTrgt, xCtrl);

    PROTECT(ReturnVal = allocVector(VECSXP,fwd.nstock()));

    for (int i=0; i<fwd.nstock(); i++)
       SET_VECTOR_ELT(ReturnVal, i, fwd.ReturnStock(i+1));
    
    return ReturnVal;
    }

extern "C" SEXPDLLExport fwd_adolc_FLBiols(SEXP xBiols, SEXP xFleets, SEXP xTrgt, SEXP xAryTrgt, SEXP xCtrl, SEXP xAryCtrl, SEXP xYrs, SEXP xDims, SEXP xSRR)    
   {
   SEXP ReturnVal = R_NilValue;
   SEXP Rerror = NEW_LOGICAL(0);

   if (NElemList(xBiols)<0)
      return Rerror;

 //  if (NElemList(xBiols)!=NElemList(xSRR))
 //     return false;

   if (NElemList(xFleets)<0)
      return Rerror;

   fwd fwd(xBiols, xFleets, xYrs, xDims, xSRR); //xSRModel, xSRParam, xSRResiduals, xMult);
	
   fwd.run(xTrgt, xAryTrgt, xCtrl, xAryCtrl, xYrs);

   PROTECT(ReturnVal = allocVector(VECSXP,6));

   //SET_VECTOR_ELT(ReturnVal, 0, fwd.ReturnBiol());
   //SET_VECTOR_ELT(ReturnVal, 1, fwd.ReturnFleet());
   //SET_VECTOR_ELT(ReturnVal, 0, fwd.ReturnStock());

   PROTECT(ReturnVal = allocVector(VECSXP,6));
  
   SEXP catch_n = R_NilValue,
        f       = R_NilValue,
        n       = R_NilValue;
  
   PROTECT(catch_n = allocVector(VECSXP,2));
   PROTECT(f       = allocVector(VECSXP,2));
   PROTECT(n       = allocVector(VECSXP,2));

   SET_VECTOR_ELT(ReturnVal, 0, fwd.landings_n.Return(1,1,1));
   SET_VECTOR_ELT(ReturnVal, 1, fwd.discards_n.Return(1,1,1));
   SET_VECTOR_ELT(ReturnVal, 2, fwd.effort.Return(1));
   SET_VECTOR_ELT(ReturnVal, 3, fwd.n.Return(1));
   SET_VECTOR_ELT(ReturnVal, 4, fwd.f.Return(1));
   SET_VECTOR_ELT(ReturnVal, 5, fwd.catch_n.Return(1,1,1));
    
   UNPROTECT(5);

   return ReturnVal;
   }

// Test functions for flc
extern "C" SEXPDLLExport flc_FLStock(SEXP xStock)
   {
   SEXP ReturnObject = R_NilValue;
  
   //FLStock_pointer stock(xStock);
   FLStock stock(xStock);

   for (int i=stock.m.minquant(); i<=stock.m.maxquant(); i++)
      stock.m(i,stock.m.minyr()) =2.0;

   return ReturnObject;
   //return stock.Return();
   }

extern "C" SEXPDLLExport TestFLBiolFLFleet(SEXP xBiol, SEXP xFleet, SEXP xDim)
   {
   SEXP ReturnVal = R_NilValue;
  
   flc bf;
   
   bf.InitBiolFleet(xBiol, xFleet, xDim);

   PROTECT(ReturnVal = allocVector(VECSXP,2));

   SET_VECTOR_ELT(ReturnVal, 0, bf.ReturnBiol());
   SET_VECTOR_ELT(ReturnVal, 1, bf.ReturnFleet());
   //SET_VECTOR_ELT(ReturnVal, 2, bf.ReturnStock());
    
   return ReturnVal;
   }

// Adapt stuff
extern "C" SEXPDLLExport AdaptFuncAD(SEXP xStock, SEXP xFitPlusGroup, SEXP xFratio, SEXP xFratioFlag, SEXP xQ, SEXP xIndex)
   {
   //Input VPA Control
   FLashVPA VPA(xStock);
       
   FLQuant  q(xQ); 
   FLQuant  index(xIndex);
 
   FLQuant_adolc q_ad;
   q_ad.Init(q,q.minyr(),q.maxyr(),q.niters());

   //Run VPA
   VPA.VPA(xFitPlusGroup,xFratio,xFratioFlag);

   //calculate objective function
   adouble ss_ad = 0.0;
   int iage, iyr;   
   for (iage=q.minquant(); iage<=q.maxquant(); iage++)
     for (iyr=index.minyr(); iyr<=index.maxyr(); iyr++)
       {
       adouble _N = VPA.N_ad(iage,iyr)*q_ad(iage,q.minyr(),1,1,1,1); 
       ss_ad += pow((index(iage,iyr)-_N)/_N,2.0);
       } 
   //end calculate objective function

   ss_ad /= (index.minquant()-index.minquant()+1)*(index.maxquant()-index.minquant()+1);
   ss_ad  = log(ss_ad);

   SEXP RtnVal     = PROTECT(NEW_NUMERIC(1)); 
   REAL(RtnVal)[0] = ss_ad.value();
 
   UNPROTECT(1);

   return RtnVal;
	}

extern "C" SEXPDLLExport AdaptGrad(SEXP xStock, SEXP xFitPlusGroup, SEXP xFratio, SEXP xFratioFlag, SEXP xQ, SEXP xIndex)
   {
   //Input VPA Control
   FLashVPA VPA(xStock);
       
   FLQuant  q(xQ); 
   FLQuant  index(xIndex);
 
   FLQuant_adolc q_ad;
   q_ad.Init(q,q.minyr(),q.maxyr(),index.niters());

   int n;
   double  *x, *g;

   n    = q.maxquant()-q.minquant()+1; 
   x    = new  double[n];
   g    = new  double[n];

   //initialise tape 
   int tag = 0, i, iage, iyr;   
   trace_on(tag); // tag = 1, keep = 0 by default

      //set independent variables
      for (iage=q.minquant(), i=0; iage<=q.maxquant(); iage++, i++)
         {
         x[i] = q(iage,q.minyr(),1,1,1,1); 
         q_ad(iage,q.minyr(),1,1,1,1) <<= x[i]; 
         }  
      //end set independent variables
 
      //Run VPA
      VPA.VPA(xFitPlusGroup,xFratio,xFratioFlag);
                 
      //calculate objective function
      double ss = 0.0;
      adouble ss_ad = 0.0;
      for (iage=q.minquant(); iage<=q.maxquant(); iage++)
           for (iage=q.minquant(), i=0; iage<=q.maxquant(); iage++, i++)
           {
           adouble _N = VPA.N_ad(iage,iyr)*q_ad(iage,q.minyr(),1,1,1,1); 
           ss_ad += pow((index(iage,iyr)-_N)/_N,2.0);
           } 

      ss_ad /= (index.minquant()-index.minquant()+1)*(index.maxquant()-index.minquant()+1);
      ss_ad  = log(ss_ad);

      ss_ad >>= ss;
     //end calculate objective function
    
   trace_off();
   //end initialise tape 
 
   //get gradient
   gradient(tag,n,x,g);
   FLVector Grad;
   Grad.Init(index.minquant(),index.maxquant(),0.0);
   for (iage=index.minquant(), i=0; iage<=index.maxquant(); iage++, i++)
      Grad(i) = g[i];
   //end get gradient
   
   //clean up
   delete[] x;    
   delete[] g;    
   
   return Grad.Return();
	
   SEXP RtnVal     = PROTECT(NEW_NUMERIC(1)); 
   REAL(RtnVal)[0] = ss_ad.value();
 
   UNPROTECT(1);

   return RtnVal;
   }

extern "C" SEXPDLLExport AdaptSetupTape(SEXP xStock, SEXP xFitPlusGroup, SEXP xFratio, SEXP xFratioFlag, SEXP xQ, SEXP xIndex)
   {
   SEXP RtnVal = R_NilValue;

   //Input VPA Control
   FLashVPA VPA(xStock);
       
   FLQuant  q(xQ); 
   FLQuant  index(xIndex);
 
   FLQuant_adolc q_ad;
   q_ad.Init(q,q.minyr(),q.maxyr(),index.niters());

   int n;
   double  *x;
   
   n    = q.maxquant()-q.minquant()+1; 
   x    = new  double[n];
   
   //initialise tape 
   int tag = 0, i, iage, iyr;   
   trace_on(tag); // tag = 1, keep = 0 by default

      //set independent variables
      for (iage=q.minquant(), i=0; iage<=q.maxquant(); iage++, i++)
         {
         x[i] = q(iage,q.minyr(),1,1,1,1); 
         q_ad(iage,q.minyr(),1,1,1,1) <<= x[i]; 
         }  
      //end set independent variables
 
      //Run VPA
      VPA.VPA(xFitPlusGroup,xFratio,xFratioFlag);
                 
      //calculate objective function
      double ss = 0.0;
      adouble ss_ad = 0.0;
      for (iyr=index.minyr(); iyr<=index.maxyr(); iyr++)
         for (iage=q.minquant(); iage<=q.maxquant(); iage++)
           {
           adouble _N = VPA.N_ad(iage,iyr)*q_ad(iage,q.minyr(),1,1,1,1); 
           ss_ad += pow((index(iage,iyr)-_N)/_N,2.0);
           } 

      ss_ad /= (index.minquant()-index.minquant()+1)*(index.maxquant()-index.minquant()+1);
      ss_ad  = log(ss_ad);

      ss_ad >>= ss;
     //end calculate objective function
    
   trace_off();
   //end initialise tape 
   
   //clean up
   delete[] x;    
   
   return RtnVal;
   }

extern "C" SEXPDLLExport AdaptTapeFunc(SEXP xQ)
   {
   SEXP RtnVal = R_NilValue;
   
   FLQuant  q(xQ); 
 
   int n, i, iage, tag=0;
   double *x, *y;

   n = q.maxquant()-q.minquant()+1; 
   x = new  double[n];
   y = new  double[1];
   
   for (iage=q.minquant(), i=0; iage<=q.maxquant(); iage++, i++)
      x[i] = q(iage,q.minyr(),1,1,1,1); 
 
   //get function
   function(tag,1,n,x,y);
   //end function
 
   RtnVal          = PROTECT(NEW_NUMERIC(1)); 
   REAL(RtnVal)[0] = y[0];
 
   //clean up
   delete[] x;    
   delete[] y;    
   	
   UNPROTECT(1);

   return RtnVal;
   }

extern "C" SEXPDLLExport AdaptTapeGrad(SEXP xQ)
   {
   SEXP RtnVal = R_NilValue;
   
   FLQuant  q(xQ); 

   int n, i, iage, tag=0;
   double *x, *g;
 
   n = q.maxquant()-q.minquant()+1; 
   x = new  double[n];
   g = new  double[n];
   
   for (iage=q.minquant(), i=0; iage<=q.maxquant(); iage++, i++)
      x[i] = q(iage,q.minyr(),1,1,1,1); 
 
   //get gradient
   gradient(tag,n,x,g);
   FLVector Grad;
   Grad.Init(q.minquant(),q.maxquant(),0.0);
   for (iage=q.minquant(), i=0; iage<=q.maxquant(); iage++, i++)
      Grad(i) = g[i];
   //end get gradient
   
   //clean up
   delete[] x;    
   delete[] g;    
   
   return Grad.Return();
	
   return RtnVal;}
