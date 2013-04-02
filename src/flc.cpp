/** $Id: flc.cpp 547 2010-03-13 12:26:11Z lauriekell $ **/
#include "flc.h"

double t1, t2, t3, t4, t5, t6;

void sv2ab(FLRConstSRR sr_model, double ***sr_params, int niters, int nunits)
  {
 	double steepness,vbiomass,spr0,a,b;

  if      (sr_model==FLRConst_BevHoltSV)
     sr_model=FLRConst_BevHolt;
  else if (sr_model==FLRConst_RickerSV)
     sr_model=FLRConst_Ricker;
  else
     return;  

	for (int i = 1; i <= niters; i++)
    for (int j = 2; j <= nunits; j++)
      {
      steepness=sr_params[i][j][1];
	    vbiomass =sr_params[i][j][2];
	    spr0     =sr_params[i][j][3];

      if (sr_model==FLRConst_BevHolt)
        {
        a = vbiomass*4.0*steepness/(spr0*(5.0*steepness-1.0));
	      b = a*spr0*(1.0/steepness - 1.0)/4.0;
	      }
      else if (sr_model==FLRConst_Ricker)
	      {
        b = log(5.0*steepness)/(vbiomass*0.8);
        a = exp(b*vbiomass)/spr0;
        }
             
     sr_params[i][j][1] = a;
     sr_params[i][j][2] = b;
	   }
   }

//deriv((y ~ catch*(f+m)/(n*(1-exp(-f-m)))-f), c("f"), func = TRUE)
double F_func(double* f,double m, double c, double n)
    {
    double expr2 = c * (*f + m);
    double expr5 = exp(-*f - m);
    double expr7 = n * (1.0 - expr5);
    double value = expr2/expr7 - *f;
    double grad  = c/expr7 - expr2 * (n * expr5)/pow(expr7,2.0) - 1.0;

    *f-=value/grad;

    return(value);
    }
 
FLRConstSRR get_sr(int i)
   {
   switch(i) {
      case 1:  return FLRConst_Mean;      
      case 2:  return FLRConst_BevHolt;    
      case 3:  return FLRConst_Ricker;    
      case 4:  return FLRConst_SegReg;   
      }

   return FLRConst_Mean;
   }

int get_sr(FLRConstSRR i)
   {
   switch(i) {
      case FLRConst_Mean:    return 1;
      case FLRConst_BevHolt: return 2;
      case FLRConst_Ricker:  return 3;
      case FLRConst_SegReg:  return 4;
      }

   return 1;
   }

sr::sr(void)
	{
	nstock = 0;
	}
	      
sr::sr(int _nstock)  
	{
	nstock = _nstock;
	}

bool sr::Init(int _nstock, SEXP xyrs)      
	{
   if (!isNumeric(xyrs)) 
      return false;
 
   int _minyr = (int)REAL(xyrs)[0], 
       _maxyr = (int)REAL(xyrs)[0];
 
   for (int i=1; i<LENGTH(xyrs); i++){
      if (_minyr>REAL(xyrs)[i]) _minyr=(int)REAL(xyrs)[i];
      if (_maxyr<REAL(xyrs)[i]) _maxyr=(int)REAL(xyrs)[i];}

   residuals.alloc_n7(_nstock);      
   param.alloc_n7(    _nstock);      

   return Init(_nstock, _minyr, _maxyr);
   }

bool sr::Init(int _nstock, int _minyr, int _maxyr)      
	{
	if (_nstock<1) 
		return false;

	if (nstock>0) 
		unalloc();

	nstock  = _nstock;
   _minyear = _minyr;
   _maxyear = _maxyr;
   
   residuals.alloc_n7(nstock);
   
   residuals_mult = new bool[nstock] - 1;
   model          = new FLRConstSRR*[nstock] - 1;
   
   for (int i=1; i<=nstock; i++)
	  {
     model[i]          = new FLRConstSRR[_maxyear-_minyear+1] - _minyear; 
     residuals_mult[i] = true;
     }

	return true;
   }
	      
bool sr::Init(int istock, SEXP xmodel, SEXP xparam, SEXP xresiduals, SEXP xmult)  
   {
   if (nstock<1 || istock>nstock || istock<0) 
     return false;

   //parameters
   param.Init(istock, PROTECT(VECTOR_ELT(xparam, istock-1)));
   
   // model
   SEXP vmodel, t;
   int  nmodel;
   
   t= PROTECT(VECTOR_ELT(xmodel, istock-1));

   if (!isNumeric(t)){
      UNPROTECT(2);
      return false;} 
   
   PROTECT(vmodel = AS_NUMERIC(t));
   nmodel         = LENGTH(vmodel);
   
   if (nmodel<(_maxyear-_minyear+1)){
      UNPROTECT(3);
      return false;} 
   
   double *dmodel = NUMERIC_POINTER(vmodel); 
   
   int i=0;
   for (int iyr=_minyear; iyr<=_maxyear; iyr++)
      model[istock][iyr]=(FLRConstSRR)((int)dmodel[i++]);
      

   // residuals
   residuals.Init(istock, PROTECT(VECTOR_ELT(xresiduals, istock-1)));

   if (residuals.minyr(istock) > _minyear || residuals.maxyr(istock) < _maxyear){
      UNPROTECT(1); 
      return false;}

   //multipicative residuals?
   SEXP vmult;   
   PROTECT(vmult = AS_LOGICAL(xmult));
   
   int *dmult = LOGICAL_POINTER(vmult);
   residuals_mult[istock] = (dmult[0]==1 ? true : false);
 	
   UNPROTECT(5); 
   
   return true;
   }
	      
sr::~sr(void)     
	{
	unalloc();
	}
	      	      	      	               
void sr::unalloc(void)
	{
	if (nstock<1) 
		return;
   
   delete [] (model         +1);
   delete [] (residuals_mult+1);
   }	      

double sr::recruits(int istock, double ssb, int iyr, int iunit, int iseason, int iarea, int iter)
   {
   double returnval=0.0,
          residual =0.0;
 
   int _yr = __max(__min(iyr, _maxyear),_minyear);

   //SSB as a function of SPR
   switch((int)model[istock][_yr]) 
      {
      case FLRConst_BevHolt: 
         returnval = param(istock,1,_yr,iunit,iseason,iarea,iter)*ssb/(ssb+param(istock,2,_yr,iunit,iseason,iarea,iter));
      break;
         
      case FLRConst_Ricker:
         returnval = param(istock,1,_yr,iunit,iseason,iarea,iter)*ssb*exp(-param(istock,2,_yr,iunit,iseason,iarea,iter)*ssb);
      break;
      
      case FLRConst_Shepherd:
         returnval = param(istock,1,_yr,iunit,iseason,iarea,iter) * ssb/pow(param(istock,2,_yr,iunit,iseason,iarea,iter) + ssb,param(istock,3,_yr,iunit,iseason,iarea,iter));
      break;
        
      case FLRConst_SegReg:
		  returnval = (ssb <= param(istock,2,_yr,iunit,iseason,iarea,iter) ?  ssb*param(istock,1,_yr,iunit,iseason,iarea,iter) : param(istock,1,_yr,iunit,iseason,iarea,iter)*param(istock,2,_yr,iunit,iseason,iarea,iter));
      break;
  
      case FLRConst_Mean: default:
         returnval = param(istock,1,_yr,iunit,iseason,iarea,iter);
      break;
      }

   if (residuals_mult[istock])
      {
      residual=(iyr<=_maxyear && iyr>=_minyear ? residuals(istock,residuals.minquant(istock),iyr,iunit,iseason,iarea,iter) : 1.0);
      returnval = returnval*residual;
      }
   else
      {
      residual=(iyr<=_maxyear && iyr>=_minyear ? residuals(istock,residuals.minquant(istock),iyr,iunit,iseason,iarea,iter) : 0.0);
      returnval = returnval+residual;
      }

   return(R_IsNA(returnval) ? 0 : returnval);
   } 

void flc::CalcF(void)
   {
   for (int i=1; i<=nstock(); i++)
     CalcF(i);
   }

void flc::CalcF(int istock)
   {
   if (istock<1 || istock>nstock())
      return;

   int iAge, iYear, iUnit, iSeason, iArea, iIter;

   for (iIter = 1; iIter<=niters(istock); iIter++)
	   for (iArea = 1; iArea <= nareas(istock); iArea++)
		  for (iSeason = 1; iSeason <= nseasons(istock); iSeason++)
    		 for (iUnit = 1; iUnit <= nunits(istock); iUnit++)
      		for (iYear = minyr(istock); iYear <= maxyr(istock); iYear++)
		   	  for (iAge = minage(istock); iAge <= maxage(istock); iAge++)
                 {
                 double _catch_n = 0.0;
                 double _f       = 0.1;
                 for (int i = 1; i <= nfleet(); i++)
                    for (int j = 1; j <= nmetier(); j++)
                       _catch_n += catch_n(istock, i,j,iAge,iYear,iUnit,iSeason,iArea,iIter); 
  
                 int j = 0;
				  	  while (fabs(F_func(&_f,m(istock,iAge,iYear,iUnit,iSeason,iArea,iIter), _catch_n, n(istock,iAge,iYear,iUnit,iSeason,iArea,iIter)))>1e-12 && j++<20)
							f(istock,iAge,iYear,iUnit,iSeason,iArea,iIter)=_f;

                 }
   }

flc::flc(void)
   {
   nfleet()  =
   nmetier() =
   nstock()    = 0;
   }      

flc::~flc(void)
  {
  unalloc_dims_range();
  }      

void flc::Init(SEXP x)
   {
   if (isFLStock( x)) InitStock(  x); else
   if (isFLStocks(x)) InitStocks( x); else
   if (isFLBiol(  x)) InitBiol(   x); else
   if (isFLBiols( x)) InitBiols(  x); else
   if (isFLFleet( x)) InitFleet(  x); else
   if (isFLFleets(x)) InitFleets( x);
   }

void flc::InitStock(SEXP x)
   {
   if (!isFLStock(x)) return;

   nfleet()  = 
   nmetier() = 
   nstock()    = 1;   
   
   //alloc_range();

   InitStock(1,1,1,x);

   minage(1)   = minage(  1,1,1);
   maxage(1)   = maxage(  1,1,1);
   plusgrp(1)  = plusgrp( 1,1,1);
   minfbar(1)  = minage(1);
   maxfbar(1)  = maxage(1);
   minyr(1)    = minyr(   1,1,1);
   maxyr(1)    = maxyr(   1,1,1);
   nunits(1)   = nunits(  1,1,1);
   nseasons(1) = nseasons(1,1,1);
   nareas(1)   = nareas(  1,1,1);
   niters(1)   = niters(  1,1,1);

   //setRange();

   int istock=1, ifleet=1, imetier=1;
   // estimate selection pattern and set effort to fbar and q to 1
   for (int iyr=minyr(istock); iyr<=maxyr(istock); iyr++)
      for (int iunit=1; iunit<=nunits(istock); iunit++)
         for (int iseason=1; iseason<=nseasons(istock); iseason++)
            for (int iarea=1; iarea<=nareas(istock); iarea++)
               for (int iter=1; iter<=niters(istock); iter++)
                  {
                  int iage;

                    effort(  ifleet,          1, iyr, iunit, iseason, iarea, iter) = 0.0;
                    effshare(ifleet, imetier, 1, iyr, iunit, iseason, iarea, iter) = 1.0;
                    for (iage = minage(istock); iage<=maxage(istock); iage++) 
                       catch_q( ifleet, imetier, istock, 1, iyr, iunit, iseason, iarea, iter) = 1.0;
                    

                  for (iage = minfbar(istock); iage<=maxfbar(istock); iage++)
                     effort(ifleet, imetier, iyr, iunit, iseason, iarea, iter) += f(istock, iage, iyr, iunit, iseason, iarea, iter);

                  effort(ifleet, imetier, iyr, iunit, iseason, iarea, iter) /= (maxfbar(istock)-minfbar(istock)+1);

                  for (iage = minage(istock); iage<=maxage(istock); iage++)
                       catch_sel(ifleet, imetier, istock, iage, iyr, iunit, iseason, iarea, iter) = f(istock, iage, iyr, iunit, iseason, iarea, iter)/effort(ifleet, 1, iyr, iunit, iseason, iarea, iter);
                  }
   }

void flc::InitStocks(SEXP x)
   {
   if (!isFLStocks(x)) return;

   nfleet()  = NElemList(x);
   nstock()    = nfleet();
   nmetier() = 1;

   int i,j;
/*   
   alloc_dims_fleet(nfleet());
   for (i=1; i<=nfleet(); i++)
      {
      alloc_dims_metiers(i);
      for (j=1; j<=nmetier(); j++)
         alloc_dims_stock(i, j, nstock());
       }
*/
   for (i=1; i<=nstock(); i++)
      InitStock(i, 1, i, VECTOR_ELT(x, i));
   
//   setRange();
   }
      
void flc::InitStock(int ifleet, int imetier, int istock, SEXP x)
   {
   if (!isFLStock( x))                 return;
   if (ifleet <1 || ifleet >nfleet() ) return;
   if (imetier<1 || imetier>nmetier()) return;
   if (istock <1 || istock >nstock() ) return;
   
   fl_minage[ ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[0];
   fl_maxage[ ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[1];
   fl_plusgrp[ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[2];
   fl_minyr[  ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[3];
   fl_maxyr[  ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[4];

   n.Init(       istock, GET_SLOT(x, install("stock.n"))); 
   stock_wt.Init(istock, GET_SLOT(x, install("stock.wt"))); 
   m.Init(       istock, GET_SLOT(x, install("m"))); 
   fec.Init(     istock, GET_SLOT(x, install("mat"))); 
   f.Init(       istock, GET_SLOT(x, install("harvest")));
   f_spwn.Init(  istock, GET_SLOT(x, install("harvest.spwn")));
   m_spwn.Init(  istock, GET_SLOT(x, install("m.spwn")));

   catch_.Init(      ifleet, imetier, istock, GET_SLOT(x, install("catch")));  
   catch_n.Init(     ifleet, imetier, istock, GET_SLOT(x, install("catch.n")));   
   catch_wt.Init(    ifleet, imetier, istock, GET_SLOT(x, install("catch.wt")));   
   catch_q.Init(     ifleet, imetier, istock, GET_SLOT(x, install("stock.n")));   
   landings.Init(    ifleet, imetier, istock, GET_SLOT(x, install("landings")));   
   landings_n.Init(  ifleet, imetier, istock, GET_SLOT(x, install("landings.n")));   
   landings_wt.Init( ifleet, imetier, istock, GET_SLOT(x, install("landings.wt")));   
   discards.Init(    ifleet, imetier, istock, GET_SLOT(x, install("discards")));   
   discards_n.Init(  ifleet, imetier, istock, GET_SLOT(x, install("discards.n")));   
   discards_wt.Init( ifleet, imetier, istock, GET_SLOT(x, install("discards.wt")));   

   catch_sel.Init(   ifleet, imetier, istock, GET_SLOT(x, install("harvest")));

   //check other dims and then set 'em   
   fl_nunits[  ifleet][imetier][istock] = 1;
   fl_nseasons[ifleet][imetier][istock] = 1;
   fl_nareas[  ifleet][imetier][istock] = 1;
   fl_niters[  ifleet][imetier][istock] = 1;

   effort.Init(   ifleet, 1, 1, minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);   
   fcost.Init(    ifleet, 1, 1, minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);   
   capacity.Init( ifleet, 1, 1, minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);  
   crewshare.Init(ifleet, 1, 1, minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);  

   effshare.Init(    ifleet, imetier, 1,            1,            minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);  
   vcost.Init(       ifleet, imetier, 1,            1,            minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);  
  
   landings_sel.Init(ifleet, imetier, istock, minage(ifleet,imetier,istock), maxage(ifleet,imetier,istock), minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);  
   discards_sel.Init(ifleet, imetier, istock, minage(ifleet,imetier,istock), maxage(ifleet,imetier,istock), minyr(ifleet,imetier,istock), maxyr(ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);  
   price.Init(       ifleet, imetier, istock, minage(ifleet,imetier,istock), maxage(ifleet,imetier,istock), minyr(ifleet,imetier,istock), maxyr( ifleet,imetier,istock), nunits(ifleet,imetier,istock), nareas(ifleet,imetier,istock), nseasons(ifleet,imetier,istock), niters(ifleet,imetier,istock), R_NaN);      
   }

bool flc::InitSR(int nstock, SEXP yrs) 
   {
   return _sr.Init(nstock, yrs);
   }

bool flc::InitSR(int istock, SEXP xmodel, SEXP  xparam, SEXP  xresiduals, SEXP xmult) 
   {
   if (istock<0 || istock>nstock())
      return false;

   return(_sr.Init(istock, xmodel, xparam, xresiduals, xmult));
   }

 SEXP flc::ReturnRange(int istock)
   {
   SEXP Range = R_NilValue;

   if (istock<1 || istock>nstock())
      return Range;

   Range          = PROTECT(NEW_NUMERIC(7)); 
   REAL(Range)[0] = minage(istock);
   REAL(Range)[1] = maxage(istock);
   REAL(Range)[2] = plusgrp(istock);
   REAL(Range)[3] = minyr(istock);
   REAL(Range)[4] = maxyr(istock);
   REAL(Range)[5] = minfbar(istock);
   REAL(Range)[6] = maxfbar(istock);

   SEXP names;
    
   PROTECT(names = allocVector(STRSXP, 7));
   SET_STRING_ELT(names, 0, mkChar("min"));
   SET_STRING_ELT(names, 1, mkChar("max"));
   SET_STRING_ELT(names, 2, mkChar("plusgroup"));
   SET_STRING_ELT(names, 3, mkChar("minyear"));
   SET_STRING_ELT(names, 4, mkChar("maxyear"));
   SET_STRING_ELT(names, 5, mkChar("minfbar"));
   SET_STRING_ELT(names, 6, mkChar("maxfbar"));
   setAttrib(Range, R_NamesSymbol, names);
   
   return Range;
   }

SEXP flc::ReturnRangeBiol(int istock)
   {
   SEXP Range = R_NilValue;

   if (istock<1 || istock>nstock())
      return Range;

   Range          = PROTECT(NEW_NUMERIC(5)); 
   REAL(Range)[0] = minage(istock);
   REAL(Range)[1] = maxage(istock);
   REAL(Range)[2] = plusgrp(istock);
   REAL(Range)[3] = minyr(istock);
   REAL(Range)[4] = maxyr(istock);

   SEXP names;
    
   PROTECT(names = allocVector(STRSXP, 5));
   SET_STRING_ELT(names, 0, mkChar("min"));
   SET_STRING_ELT(names, 1, mkChar("max"));
   SET_STRING_ELT(names, 2, mkChar("plusgroup"));
   SET_STRING_ELT(names, 3, mkChar("minyear"));
   SET_STRING_ELT(names, 4, mkChar("maxyear"));

   setAttrib(Range, R_NamesSymbol, names);
   
   return Range;
   }

SEXP flc::ReturnStock(int istock)
   {
   SEXP ReturnObj = R_NilValue;
   
   if (istock<1    || istock>nstock())
      return  ReturnObj;
   
   ReturnObj  = PROTECT(NEW_OBJECT(MAKE_CLASS("FLStock")));
   
   SET_SLOT(ReturnObj, install("range"),        ReturnRange(istock));

   SET_SLOT(ReturnObj, install("stock.n"),      n.Return(           istock));
   SET_SLOT(ReturnObj, install("stock.wt"),     stock_wt.Return(    istock));
   SET_SLOT(ReturnObj, install("m"),            m.Return(           istock));
   SET_SLOT(ReturnObj, install("mat"),          fec.Return(         istock));
   SET_SLOT(ReturnObj, install("harvest"),      f.Return(           istock));
   SET_SLOT(ReturnObj, install("harvest.spwn"), f_spwn.Return(      istock));
   SET_SLOT(ReturnObj, install("m.spwn"),       m_spwn.Return(      istock));
 
   FLQuant flq_stock(   1, 1, minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_catch(   1, 1, minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_landings(1, 1, minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_discards(1, 1, minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
  
   FLQuant flq_catch_wt(   minage(), maxage(), minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_landings_wt(minage(), maxage(), minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_discards_wt(minage(), maxage(), minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
  
   FLQuant flq_catch_n(   minage(), maxage(), minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_landings_n(minage(), maxage(), minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);
   FLQuant flq_discards_n(minage(), maxage(), minyr(), maxyr(), nunits(), nareas(), nseasons(), niters(), 0.0);

   for (int iyear=minyr(); iyear<=maxyr(); iyear++)
      for (int iunit=1; iunit<=nunits(); iunit++)
         for (int iseason=1; iseason<=nseasons(); iseason++)
            for (int iarea=1; iarea<=nareas(); iarea++)
               for (int iter=1; iter<=niters(); iter++)
                  {
                  flq_stock(   1, iyear, iunit, iseason, iarea, iter) = biomass(n,         stock_wt,    istock, iyear, iunit, iseason, iarea, iter);
                  flq_catch(   1, iyear, iunit, iseason, iarea, iter) = biomass(catch_n,   catch_wt,    istock, iyear, iunit, iseason, iarea, iter);
                  flq_landings(1, iyear, iunit, iseason, iarea, iter) = biomass(landings_n,landings_wt, istock, iyear, iunit, iseason, iarea, iter);
                  flq_discards(1, iyear, iunit, iseason, iarea, iter) = biomass(discards_n,discards_wt, istock, iyear, iunit, iseason, iarea, iter);
                  for (int iage=minage(); iage<=maxage(); iage++)
                     {
                     flq_catch_wt(   iage, iyear, iunit, iseason, iarea, iter) = wt(catch_wt,    catch_n,    istock, iage, iyear, iunit, iseason, iarea, iter);
                     flq_landings_wt(iage, iyear, iunit, iseason, iarea, iter) = wt(landings_wt, landings_n, istock, iage, iyear, iunit, iseason, iarea, iter);
                     flq_discards_wt(iage, iyear, iunit, iseason, iarea, iter) = wt(discards_wt, discards_n, istock, iage, iyear, iunit, iseason, iarea, iter);
 
                     flq_catch_n(   iage, iyear, iunit, iseason, iarea, iter) = sum(catch_n,    istock, iage, iyear, iunit, iseason, iarea, iter);
                     flq_landings_n(iage, iyear, iunit, iseason, iarea, iter) = sum(landings_n, istock, iage, iyear, iunit, iseason, iarea, iter);
                     flq_discards_n(iage, iyear, iunit, iseason, iarea, iter) = sum(discards_n, istock, iage, iyear, iunit, iseason, iarea, iter);
                     }
                  }

   SET_SLOT(ReturnObj, install("stock"),        flq_catch.Return());
   SET_SLOT(ReturnObj, install("catch"),        flq_catch.Return());
   SET_SLOT(ReturnObj, install("catch.n"),      flq_catch_n.Return());
   SET_SLOT(ReturnObj, install("catch.wt"),     flq_catch_wt.Return());
   SET_SLOT(ReturnObj, install("landings"),     flq_landings.Return());
   SET_SLOT(ReturnObj, install("landings.n"),   flq_landings_n.Return());
   SET_SLOT(ReturnObj, install("landings.wt"),  flq_landings_wt.Return());
   SET_SLOT(ReturnObj, install("discards"),     flq_discards.Return());
   SET_SLOT(ReturnObj, install("discards.n"),   flq_discards_n.Return());
   SET_SLOT(ReturnObj, install("discards.wt"),  flq_discards_wt.Return());
  
   UNPROTECT(2);

   return  ReturnObj;
   }

SEXP flc::ReturnBiol(int istock)
   {
   SEXP ReturnObj = R_NilValue;
   
   if (istock<1 || istock>nstock())
      return  ReturnObj;
   
   ReturnObj = PROTECT(NEW_OBJECT(MAKE_CLASS("FLBiol")));
   
   SET_SLOT(ReturnObj, install("range"),        ReturnRangeBiol(istock));

   SET_SLOT(ReturnObj, install("n"),    n.Return(           istock));
   SET_SLOT(ReturnObj, install("wt"),   stock_wt.Return(    istock));
   SET_SLOT(ReturnObj, install("m"),    m.Return(           istock));
   SET_SLOT(ReturnObj, install("fec"),  fec.Return(         istock));
   SET_SLOT(ReturnObj, install("spwn"), m_spwn.Return(      istock)); 
  
   UNPROTECT(2);

   return  ReturnObj;
   }

SEXP flc::ReturnFleet(int ifleet)
   {
   SEXP ReturnObj = R_NilValue;
   
   if (ifleet<1 || ifleet>nfleet()) return ReturnObj;

   SEXP fleet, Range;

   PROTECT(fleet  = NEW_OBJECT(MAKE_CLASS("FLFleet")));
   PROTECT(Range  = PROTECT(NEW_NUMERIC(5))); 
   
   REAL(Range)[0] = minage();
   REAL(Range)[1] = maxage();
   //REAL(Range)[2] = plusgrp();
   REAL(Range)[3] = minyr();
   REAL(Range)[4] = maxyr();
   
   SET_SLOT(fleet, install("range"),      Range);         
   SET_SLOT(fleet, install("effort"),     effort.Return(ifleet));         
   SET_SLOT(fleet, install("capacity"),   capacity.Return(ifleet));         
   SET_SLOT(fleet, install("crewshare"),  crewshare.Return(ifleet));         
   SET_SLOT(fleet, install("fcost"),      fcost.Return(ifleet));         
   
   SEXP metiers      = NEW_OBJECT(MAKE_CLASS("FLMetiers"));
   SEXP metiers_data = allocVector(VECSXP, nmetier());        
   
   for (int i=1; i<=nmetier(); i++)
      {
      SEXP metier = NEW_OBJECT(MAKE_CLASS("FLMetier"));

      SET_SLOT(metier, install("effshare"), effshare.Return(ifleet,i));         
      SET_SLOT(metier, install("vcost"),    vcost.Return(ifleet,i));         
 
      SEXP catches    = NEW_OBJECT(MAKE_CLASS("FLCatches"));
      SEXP catch_data = allocVector(VECSXP, nstock());        
   
      for (int j=1; j<=nstock(); j++)
         {
         SEXP _catch = NEW_OBJECT(MAKE_CLASS("FLCatch"));
      
         //SET_SLOT(_catch, install("catch"),        catch_.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("catch.n"),      catch_n.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("catch.wt"),     catch_wt.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("catch.sel"),    catch_sel.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("catch.q"),      catch_q.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("landings"),     landings.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("landings.n"),   landings_n.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("landings.wt"),  landings_wt.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("landings.sel"), landings_sel.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("discards"),     discards.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("discards.n"),   discards_n.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("discards.wt"),  discards_wt.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("discards.sel"), discards_sel.Return(ifleet,i,j));
         //SET_SLOT(_catch, install("price"),        price.Return(ifleet,i,j));

         SET_VECTOR_ELT(catch_data, j-1, _catch); 
         }
      catches = R_do_slot_assign(catches, install(".Data"), catch_data);
    
      SET_SLOT(metier, install("catches"), catches);
      
      SET_VECTOR_ELT(metiers_data, i-1, metier); 
      }
   metiers = R_do_slot_assign(metiers, install(".Data"), metiers_data);
      
   SET_SLOT(fleet, install("metiers"), metiers);
   
   return fleet;
   }

SEXP flc::ReturnFleets(void)
   {
   SEXP ReturnObject = PROTECT(NEW_OBJECT(MAKE_CLASS("FLFleets")));
   SEXP flfs         = allocVector(VECSXP,nfleet());
   
   int i, j, k;
   for (i=1; i<=nfleet(); i++)
      {
      SEXP flf = allocVector(VECSXP,nmetier());
      SEXP flms = PROTECT(NEW_OBJECT(MAKE_CLASS("FLMetiers")));
      
      for (j=1; j<=nmetier(); j++)
         {
         SEXP flm = allocVector(VECSXP,nstock());
         SEXP flcs   = PROTECT(NEW_OBJECT(MAKE_CLASS("FLCatches")));
            
         for (k=1; k<=nstock(); k++)
            {
            SEXP Range  = PROTECT(NEW_NUMERIC(5));   
            SEXP flc = PROTECT(NEW_OBJECT(MAKE_CLASS("FLCatch")));
            
            REAL(Range)[0] = minage(i,j,k);
            REAL(Range)[1] = maxage(i,j,k);
            REAL(Range)[2] = plusgrp( i,j,k);
            REAL(Range)[3] = minyr(   i,j,k);
            REAL(Range)[4] = maxyr(   i,j,k);
       
            SET_SLOT(flc, install("range"), Range);

            SET_SLOT(flc, install("catch"),       catch_.Return(      i,j,k));
            SET_SLOT(flc, install("catch.n"),     catch_n.Return(     i,j,k));
            SET_SLOT(flc, install("catch.wt"),    catch_wt.Return(    i,j,k));
            SET_SLOT(flc, install("catch.q"),     catch_q.Return(     i,j,k));
            SET_SLOT(flc, install("catch.sel"),   catch_sel.Return(   i,j,k));
            SET_SLOT(flc, install("landings"),    landings.Return(    i,j,k));
            SET_SLOT(flc, install("landings.n"),  landings_n.Return(  i,j,k));
            SET_SLOT(flc, install("landings.wt"), landings_wt.Return( i,j,k));
            SET_SLOT(flc, install("landings.sel"),landings_sel.Return(i,j,k)); 
            SET_SLOT(flc, install("discards"),    discards.Return(    i,j,k));
            SET_SLOT(flc, install("discards.n"),  discards_n.Return(  i,j,k));
            SET_SLOT(flc, install("discards.wt"), discards_wt.Return( i,j,k));
            SET_SLOT(flc, install("discards.sel"),discards_sel.Return(i,j,k)); 
            SET_SLOT(flc, install("price"),       price.Return(       i,j,k));

            SET_VECTOR_ELT(flcs, k-1, flc);
   
            UNPROTECT(2);
            }

         SET_SLOT(flm, install("effshare"), effshare.Return(i,j));
         SET_SLOT(flm, install("vcos"),     vcost.Return(    i,j));
         
         SET_VECTOR_ELT(flms, j-1, flm);
         }
      
      SET_SLOT(flf, install("fcost"),     fcost.Return(    i));
      SET_SLOT(flf, install("capacity"),  capacity.Return( i));
      SET_SLOT(flf, install("crewshare"), crewshare.Return(i));

      SET_VECTOR_ELT(flfs, i-1, flf);
      }

   //SET_SLOT(ReturnObject, install(".Data"), flfs);
   //ReturnObject = 
   R_do_slot_assign(ReturnObject, install(".Data"), flfs);
             
   UNPROTECT(1);

   return ReturnObject;
   }

double flc::wt(FLQuant4 &wt,FLQuant4 &n, int istock, int iage, int iyear, int iunit, int iseason, int iarea, int iter)
   {
   double _wt=0, _sum=0;

   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      for (int imetier=1; imetier<=nmetier(); imetier++)
         {
         _sum+=n( ifleet,imetier,istock,iage,iyear,iunit,iseason,iarea,iter);
         _wt +=wt(ifleet,imetier,istock,iage,iyear,iunit,iseason,iarea,iter)*n(ifleet,imetier,istock,iage,iyear,iunit,iseason,iarea,iter);
         }

   if (_sum>0)
      return (_wt/_sum); 
   else
      return (0.0); 
   }

double flc::biomass(FLQuant4 &wt,FLQuant4 &n, int istock, int iyear, int iunit, int iseason, int iarea, int iter)
   {
   double _sum=0;

   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      for (int imetier=1; imetier<=nmetier(); imetier++)
          for (int iage=minage(ifleet,imetier,istock); iage<=maxage(ifleet,imetier,istock); iage++)
            {
            double temp=wt(ifleet,imetier,istock,iage,iyear,iunit,iseason,iarea,iter)*n(ifleet,imetier,istock,iage,iyear,iunit,iseason,iarea,iter);
            if (temp>0.0 && !R_IsNA(temp)) _sum+=temp;
            }

   return (_sum); 
   }

double flc::biomass(FLQuant2 &wt,FLQuant2 &n, int istock, int iyear, int iunit, int iseason, int iarea, int iter)
   {
   double _sum=0;

   for (int iage=minage(istock); iage<=maxage(istock); iage++)
      {
      double temp=wt(istock,iage,iyear,iunit,iseason,iarea,iter)*n(istock,iage,iyear,iunit,iseason,iarea,iter);
      if (temp>0.0 && !R_IsNA(temp)) _sum+=temp;
      }

   return (_sum); 
   }

double flc::sum(FLQuant4 &n, int istock, int iage, int iyear, int iunit, int iseason, int iarea, int iter)
   {
   double _sum=0;

   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      for (int imetier=1; imetier<=nmetier(); imetier++)
         {
         double temp =n(ifleet,imetier,istock,iage,iyear,iunit,iseason,iarea,iter);
         if (temp>0.0 && !R_IsNA(temp)) _sum+=temp;
         }
   
   return (_sum); 
   }

bool flc::InitBiolFleet(SEXP xBiol, SEXP xFleet, SEXP xDim)
   {
   if (!isFLBiols(xBiol) && !isFLFleets(xFleet)) 
     return false;
   
   SEXP fleet   = PROTECT(VECTOR_ELT(xFleet, 0));
   SEXP metiers = PROTECT(GET_SLOT(fleet, install("metiers")));

   nstock() = NElemList(xBiol);
   nfleet() = NElemList(xFleet);
   nmetier()= NElemList(metiers);
   
   if (nfleet() < 1 || nmetier() < 1 || nstock() < 1)
      {
      UNPROTECT(2);

      return false;
      }

   SEXP dim;
   PROTECT(dim = AS_NUMERIC(xDim));
   if (LENGTH(dim)<6)
      {
      UNPROTECT(3);

      return false;
      }
        
   minyr()    = (int)REAL(dim)[0];
   maxyr()    = (int)REAL(dim)[1];
   nunits()   = (int)REAL(dim)[2];
   nseasons() = (int)REAL(dim)[3];
   nareas()   = (int)REAL(dim)[4];
   niters()   = (int)REAL(dim)[5];
   int flag   = (int)REAL(dim)[6];
  
   alloc_dims_range();

   InitBiols(xBiol); 
   InitFleets(xFleet);

   CalcF(1);

   UNPROTECT(3);

   return true;
   }

void flc::InitFleet(SEXP x)
   {
   if (!isFLFleet(x)) return;

   nfleet() = 1;
   effort.Init(   nfleet(), GET_SLOT(x, install("effort"))); 
   capacity.Init( nfleet(), GET_SLOT(x, install("capacity")));   
   crewshare.Init(nfleet(), GET_SLOT(x, install("crewshare")));
   fcost.Init(    nfleet(), GET_SLOT(x, install("fcost"))); 

   SEXP metiers = PROTECT(GET_SLOT(x, install("metiers")));
   int iProtect = 1;
   
   nmetier() = NElemList(metiers);

   if (nmetier() < 1)
      { 
      UNPROTECT(iProtect);
      return;
      }

   catch_.alloc(      nfleet(), nmetier(), nstock());
   catch_n.alloc(     nfleet(), nmetier(), nstock());
   catch_wt.alloc(    nfleet(), nmetier(), nstock());
   catch_sel.alloc(   nfleet(), nmetier(), nstock());
   catch_q.alloc(     nfleet(), nmetier(), nstock());
   landings.alloc(    nfleet(), nmetier(), nstock());
   landings_n.alloc(  nfleet(), nmetier(), nstock());
   landings_wt.alloc( nfleet(), nmetier(), nstock());
   landings_sel.alloc(nfleet(), nmetier(), nstock());
   discards.alloc(    nfleet(), nmetier(), nstock());
   discards_n.alloc(  nfleet(), nmetier(), nstock());
   discards_wt.alloc( nfleet(), nmetier(), nstock());
   discards_sel.alloc(nfleet(), nmetier(), nstock());
   price.alloc(       nfleet(), nmetier(), nstock());

   for (int i=1; i<=nmetier(); i++)
      {
      SEXP metier  = PROTECT(VECTOR_ELT(metiers, i-1));
      iProtect++;

      effshare.Init(1, i, GET_SLOT(metier, install("effshare"))); 
      vcost.Init(   1, i, GET_SLOT(metier, install("vcost"))); 
         
      SEXP catches = PROTECT(GET_SLOT(metier, install("catches")));
      iProtect++;
      if (NElemList(catches) != nstock()) 
         { 
         UNPROTECT(iProtect);
         return;
         }

      for (int j=1; j<=nstock(); j++)
         {
         SEXP _catch = PROTECT(VECTOR_ELT(catches, i-1));
         iProtect++;

         catch_.Init(      1, i, j, GET_SLOT(_catch, install("catch"))); 
         catch_n.Init(     1, i, j, GET_SLOT(_catch, install("catch.n"))); 
         catch_wt.Init(    1, i, j, GET_SLOT(_catch, install("catch.wt"))); 
         catch_sel.Init(   1, i, j, GET_SLOT(_catch, install("catch.sel"))); 
         catch_q.Init(     1, i, j, GET_SLOT(_catch, install("catch.q"))); 
         landings.Init(    1, i, j, GET_SLOT(_catch, install("landings"))); 
         landings_n.Init(  1, i, j, GET_SLOT(_catch, install("landings.n"))); 
         landings_wt.Init( 1, i, j, GET_SLOT(_catch, install("landings.wt"))); 
         landings_sel.Init(1, i, j, GET_SLOT(_catch, install("landings.sel"))); 
         discards.Init(    1, i, j, GET_SLOT(_catch, install("discards"))); 
         discards_n.Init(  1, i, j, GET_SLOT(_catch, install("discards.n"))); 
         discards_wt.Init( 1, i, j, GET_SLOT(_catch, install("discards.wt"))); 
         discards_sel.Init(1, i, j, GET_SLOT(_catch, install("discards.sel"))); 
         price.Init(       1, i, j, GET_SLOT(_catch, install("price"))); 
         }
      }

   UNPROTECT(iProtect);
   }

void flc::InitBiol(SEXP x)
   {
   if (!isFLBiol(x)) return;

   nfleet()  = 
   nmetier() = 
   nstock()    = 1;   
   
//   alloc_dims();

   InitBiol(1,1,1,x);

   minage(1)   = minage(  1,1,1);
   maxage(1)   = maxage(  1,1,1);
   plusgrp(1)  = plusgrp( 1,1,1);
   minfbar(1)  = minage(  1);
   maxfbar(1)  = maxage(  1);
   minyr(1)    = minyr(   1,1,1);
   maxyr(1)    = maxyr(   1,1,1);
   nunits(1)   = nunits(  1,1,1);
   nseasons(1) = nseasons(1,1,1);
   nareas(1)   = nareas(  1,1,1);
   niters(1)   = niters(  1,1,1);

//   setRange();
 
   return;
   }

void flc::InitBiol(int ifleet, int imetier, int istock, SEXP x)
   {
   if (!isFLBiol( x))                  return;
   if (ifleet <1 || ifleet >nfleet())  return;
   if (imetier<1 || imetier>nmetier()) return;
   if (istock   <1 || istock   >nstock())    return;
   
   fl_minage[ ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[0];
   fl_maxage[ ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[1];
   fl_plusgrp[ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[2];
   fl_minyr[  ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[3];
   fl_maxyr[  ifleet][imetier][istock] = (int)REAL(GET_SLOT(x, install("range")))[4];

   //f.Init(       istock, GET_SLOT(x, install("harvest")));
   n.Init(       istock, GET_SLOT(x, install("n"))); 
   stock_wt.Init(istock, GET_SLOT(x, install("wt"))); 
   m.Init(       istock, GET_SLOT(x, install("m"))); 
   fec.Init(     istock, GET_SLOT(x, install("fec"))); 
   f_spwn.Init(  istock, GET_SLOT(x, install("spwn")));
   m_spwn.Init(  istock, GET_SLOT(x, install("spwn")));
   f.Init(       istock, GET_SLOT(x, install("m"))); 

   //set other dims   
   fl_nunits[  ifleet][imetier][istock] = 1;
   fl_nseasons[ifleet][imetier][istock] = 1;
   fl_nareas[  ifleet][imetier][istock] = 1;
   fl_niters[  ifleet][imetier][istock] = n.niters(istock);
   }
           
/*
void flc::alloc_dims(void)
   {
   alloc_dims_nstock(nstock());
   alloc_dims_fleets();
   
   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      {
      alloc_dims_metiers(ifleet);
      
      for (int imetier=1; imetier<=nmetier(); imetier++)
         alloc_dims_stock(ifleet, imetier, nstock());
      }
   }

void flc::unalloc_dims(void)
   {
   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      {
      unalloc_dims_metiers(ifleet);
      
      for (int imetier=1; imetier<=nmetier(); imetier++)
         unalloc_dims_stock(ifleet, imetier);
      }

   unalloc_dims_fleets();  
   unalloc_dims_nstock();
   }

void flc::alloc_dims_nstock(int __nstock)
   {
   if (__nstock<1) 
      return;

   if (nstock()>1) 
      unalloc_dims_nstock();

   nstock() = __nstock;

   stock_minage   = new int[nstock()] - 1;
   stock_maxage   = new int[nstock()] - 1;
   stock_plusgrp  = new int[nstock()] - 1;
   stock_minfbar  = new int[nstock()] - 1;
   stock_maxfbar  = new int[nstock()] - 1;
   stock_minyr    = new int[nstock()] - 1;
   stock_maxyr    = new int[nstock()] - 1;
   stock_nunits   = new int[nstock()] - 1;
   stock_nseasons = new int[nstock()] - 1;
   stock_nareas   = new int[nstock()] - 1;
   stock_niters   = new int[nstock()] - 1;

   n.alloc_n7(nstock());
   m.alloc_n7(nstock());
   f.alloc_n7(nstock());
   stock_wt.alloc_n7(nstock());
   fec.alloc_n7(nstock());
   m_spwn.alloc_n7(nstock());
   f_spwn.alloc_n7(nstock());
   }

void flc::unalloc_dims_nstock(void)
   {
   delete[]  (stock_minage    + 1);
   delete[]  (stock_maxage    + 1);
   delete[]  (stock_plusgrp   + 1);
   delete[]  (stock_minfbar   + 1);
   delete[]  (stock_maxfbar   + 1);
   delete[]  (stock_minyr     + 1);
   delete[]  (stock_maxyr     + 1);
   delete[]  (stock_nunits    + 1);
   delete[]  (stock_nseasons  + 1);
   delete[]  (stock_nareas    + 1);
   delete[]  (stock_niters    + 1);

   n.unalloc_n7();
   m.unalloc_n7();
   f.unalloc_n7();
   stock_wt.unalloc_n7();
   fec.unalloc_n7();
   m_spwn.unalloc_n7();
   f_spwn.unalloc_n7();
   }

void flc::alloc_dims_fleet(int ifleet)
   {
   ;
   }

void flc::unalloc_dims_fleets(void)
   {
   ;
   }

void flc::alloc_dims_fleets(void)
   {
   fl_minage   = new int**[nfleet()] - 1;
   fl_maxage   = new int**[nfleet()] - 1;
   fl_plusgrp  = new int**[nfleet()] - 1;
   fl_minyr    = new int**[nfleet()] - 1;
   fl_maxyr    = new int**[nfleet()] - 1;
   fl_nunits   = new int**[nfleet()] - 1;
   fl_nseasons = new int**[nfleet()] - 1;
   fl_nareas   = new int**[nfleet()] - 1;
   fl_niters   = new int**[nfleet()] - 1;
   fl_InitFlag = new bool**[nfleet()] - 1;
   
   catch_.alloc_n9(      nfleet()); 
   catch_n.alloc_n9(     nfleet()); 
   catch_wt.alloc_n9(    nfleet()); 
   catch_q.alloc_n9(     nfleet()); 
   catch_sel.alloc_n9(   nfleet()); 
   landings.alloc_n9(    nfleet()); 
   landings_n.alloc_n9(  nfleet()); 
   landings_wt.alloc_n9( nfleet()); 
   landings_sel.alloc_n9(nfleet()); 
   discards.alloc_n9(    nfleet()); 
   discards_n.alloc_n9(  nfleet()); 
   discards_wt.alloc_n9( nfleet()); 
   discards_sel.alloc_n9(nfleet()); 
   price.alloc_n9(       nfleet());
   
   effshare.alloc_n8(nfleet());
   vcost.alloc_n8(nfleet());

   effort.alloc_n7(   nfleet());
   fcost.alloc_n7(    nfleet());
   capacity.alloc_n7( nfleet());
   crewshare.alloc_n7(nfleet());
   }

void flc::alloc_dims_metiers(int ifleet)
   {
   if (ifleet<1 || ifleet>nfleet())
      return;

   fl_minage[ifleet]   = new int*[nmetier()] - 1;
   fl_maxage[ifleet]   = new int*[nmetier()] - 1;
   fl_plusgrp[ifleet]  = new int*[nmetier()] - 1;
   fl_minyr[ifleet]    = new int*[nmetier()] - 1;
   fl_maxyr[ifleet]    = new int*[nmetier()] - 1;
   fl_nunits[ifleet]   = new int*[nmetier()] - 1;
   fl_nseasons[ifleet] = new int*[nmetier()] - 1;
   fl_nareas[ifleet]   = new int*[nmetier()] - 1;
   fl_niters[ifleet]   = new int*[nmetier()] - 1;
   fl_InitFlag[ifleet] = new bool*[nmetier()] - 1;
   
   effshare.alloc_n7(ifleet, nmetier());
   vcost.alloc_n7(   ifleet, nmetier());

   catch_.alloc_n8(      ifleet, nmetier()); 
   catch_n.alloc_n8(     ifleet, nmetier()); 
   catch_wt.alloc_n8(    ifleet, nmetier()); 
   catch_q.alloc_n8(     ifleet, nmetier()); 
   catch_sel.alloc_n8(   ifleet, nmetier()); 
   landings.alloc_n8(    ifleet, nmetier()); 
   landings_n.alloc_n8(  ifleet, nmetier()); 
   landings_wt.alloc_n8( ifleet, nmetier()); 
   landings_sel.alloc_n8(ifleet, nmetier()); 
   discards.alloc_n8(    ifleet, nmetier()); 
   discards_n.alloc_n8(  ifleet, nmetier()); 
   discards_wt.alloc_n8( ifleet, nmetier()); 
   discards_sel.alloc_n8(ifleet, nmetier()); 
   price.alloc_n8(       ifleet, nmetier());

   for (int imetier=1; imetier<=nmetier(); imetier++)
      {
      catch_.alloc_n7(      ifleet, imetier, nstock()); 
      catch_n.alloc_n7(     ifleet, imetier, nstock()); 
      catch_wt.alloc_n7(    ifleet, imetier, nstock()); 
      catch_q.alloc_n7(     ifleet, imetier, nstock()); 
      catch_sel.alloc_n7(   ifleet, imetier, nstock()); 
      landings.alloc_n7(    ifleet, imetier, nstock()); 
      landings_n.alloc_n7(  ifleet, imetier, nstock()); 
      landings_wt.alloc_n7( ifleet, imetier, nstock()); 
      landings_sel.alloc_n7(ifleet, imetier, nstock()); 
      discards.alloc_n7(    ifleet, imetier, nstock()); 
      discards_n.alloc_n7(  ifleet, imetier, nstock()); 
      discards_wt.alloc_n7( ifleet, imetier, nstock()); 
      discards_sel.alloc_n7(ifleet, imetier, nstock()); 
      price.alloc_n7(       ifleet, imetier, nstock());
      }
   }

void flc::unalloc_dims_metiers(int ifleet)
   {
   if (ifleet<1 || ifleet>nfleet())
      return;

   delete [] (fl_minage[ifleet]   + 1);
   delete [] (fl_maxage[ifleet]   + 1);
   delete [] (fl_plusgrp[ifleet]  + 1);
   delete [] (fl_minyr[ifleet]    + 1);
   delete [] (fl_maxyr[ifleet]    + 1);
   delete [] (fl_nunits[ifleet]   + 1);
   delete [] (fl_nseasons[ifleet] + 1);
   delete [] (fl_nareas[ifleet]   + 1);
   delete [] (fl_niters[ifleet]   + 1);
   delete [] (fl_InitFlag[ifleet] + 1);
   
   for (int imetier=1; imetier<=nmetier(); imetier++)
      {
      catch_.unalloc_n7(      ifleet, imetier); 
      catch_n.unalloc_n7(     ifleet, imetier); 
      catch_wt.unalloc_n7(    ifleet, imetier); 
      catch_q.unalloc_n7(     ifleet, imetier); 
      catch_sel.unalloc_n7(   ifleet, imetier); 
      landings.unalloc_n7(    ifleet, imetier); 
      landings_n.unalloc_n7(  ifleet, imetier); 
      landings_wt.unalloc_n7( ifleet, imetier); 
      landings_sel.unalloc_n7(ifleet, imetier); 
      discards.unalloc_n7(    ifleet, imetier); 
      discards_n.unalloc_n7(  ifleet, imetier); 
      discards_wt.unalloc_n7( ifleet, imetier); 
      discards_sel.unalloc_n7(ifleet, imetier); 
      price.unalloc_n7(       ifleet, imetier);
      }

   effshare.unalloc_n7(ifleet);
   vcost.unalloc_n7(   ifleet);

   catch_.unalloc_n8(      ifleet); 
   catch_n.unalloc_n8(     ifleet); 
   catch_wt.unalloc_n8(    ifleet); 
   catch_q.unalloc_n8(     ifleet); 
   catch_sel.unalloc_n8(   ifleet); 
   landings.unalloc_n8(    ifleet); 
   landings_n.unalloc_n8(  ifleet); 
   landings_wt.unalloc_n8( ifleet); 
   landings_sel.unalloc_n8(ifleet); 
   discards.unalloc_n8(    ifleet); 
   discards_n.unalloc_n8(  ifleet); 
   discards_wt.unalloc_n8( ifleet); 
   discards_sel.unalloc_n8(ifleet); 
   price.unalloc_n8(       ifleet);
   }

void flc::alloc_dims_stock(int ifleet, int imetier, int __nstock)
   {
   fl_minage[ifleet][imetier]   = new int[__nstock] - 1;
   fl_maxage[ifleet][imetier]   = new int[__nstock] - 1;
   fl_plusgrp[ifleet][imetier]  = new int[__nstock] - 1;
   fl_minyr[ifleet][imetier]    = new int[__nstock] - 1;
   fl_maxyr[ifleet][imetier]    = new int[__nstock] - 1;
   fl_nunits[ifleet][imetier]   = new int[__nstock] - 1;
   fl_nseasons[ifleet][imetier] = new int[__nstock] - 1;
   fl_nareas[ifleet][imetier]   = new int[__nstock] - 1;
   fl_niters[ifleet][imetier]   = new int[__nstock] - 1;
   fl_InitFlag[ifleet][imetier] = new bool[__nstock] - 1;

   catch_q.alloc_n7(     ifleet, imetier, __nstock);
   catch_.alloc_n7(      ifleet, imetier, __nstock); 
   catch_n.alloc_n7(     ifleet, imetier, __nstock); 
   catch_wt.alloc_n7(    ifleet, imetier, __nstock); 
   catch_q.alloc_n7(     ifleet, imetier, __nstock); 
   catch_sel.alloc_n7(   ifleet, imetier, __nstock); 
   landings.alloc_n7(    ifleet, imetier, __nstock); 
   landings_n.alloc_n7(  ifleet, imetier, __nstock); 
   landings_wt.alloc_n7( ifleet, imetier, __nstock); 
   landings_sel.alloc_n7(ifleet, imetier, __nstock); 
   discards.alloc_n7(    ifleet, imetier, __nstock); 
   discards_n.alloc_n7(  ifleet, imetier, __nstock); 
   discards_wt.alloc_n7( ifleet, imetier, __nstock); 
   discards_sel.alloc_n7(ifleet, imetier, __nstock); 
   price.alloc_n7(       ifleet, imetier, __nstock);
   }

void flc::unalloc_dims_fleet(void)
   {
   delete[]  (fl_minage   + 1);
   delete[]  (fl_maxage   + 1);
   delete[]  (fl_plusgrp  + 1);
   delete[]  (fl_minyr    + 1);
   delete[]  (fl_maxyr    + 1);
   delete[]  (fl_nunits   + 1);
   delete[]  (fl_nseasons + 1);
   delete[]  (fl_nareas   + 1);
   delete[]  (fl_niters   + 1);
   delete[]  (fl_InitFlag + 1);
   }

void flc::unalloc_dims_stock(int ifleet, int imetier)
   {
   delete[]  (fl_minage[ifleet][imetier]   + 1);
   delete[]  (fl_maxage[ifleet][imetier]   + 1);
   delete[]  (fl_plusgrp[ifleet][imetier]  + 1);
   delete[]  (fl_minyr[ifleet][imetier]    + 1);
   delete[]  (fl_maxyr[ifleet][imetier]    + 1);
   delete[]  (fl_nunits[ifleet][imetier]   + 1);
   delete[]  (fl_nseasons[ifleet][imetier] + 1);
   delete[]  (fl_nareas[ifleet][imetier]   + 1);
   delete[]  (fl_niters[ifleet][imetier]   + 1);
   delete[]  (fl_InitFlag[ifleet][imetier] + 1);
   }
*/

/*
void flc::setRange(void)
   {
   bool first=TRUE;
   
   for (int istock=1; istock<=nstock(); istock++)
      for (int ifleet=1; ifleet<=nfleet(); ifleet++)
         for (int imetier=1; imetier<=nmetier(); imetier++)
            if (InitFlag(ifleet,imetier,istock))
               if (first)
                  {
                  minage()   = minage(  ifleet,imetier,istock);
                  maxage()   = maxage(  ifleet,imetier,istock);
                  minyr()    = minyr(   ifleet,imetier,istock);
                  maxyr()    = maxyr(   ifleet,imetier,istock);
                  nunits()   = nunits(  ifleet,imetier,istock);
                  nseasons() = nseasons(ifleet,imetier,istock);
                  nareas()   = nareas(  ifleet,imetier,istock);
                  niters()   = niters(  ifleet,imetier,istock);
                  first = false;
                  }  
               else
                  {
                  minage()   = __max(minage(),   minage(  ifleet,imetier,istock));
                  maxage()   = __min(maxage(),   maxage(  ifleet,imetier,istock));
                  minyr()    = __max(minyr(),    minyr(   ifleet,imetier,istock));
                  maxyr()    = __min(maxyr(),    maxyr(   ifleet,imetier,istock));
                  nunits()   = __max(nunits(),   nunits(  ifleet,imetier,istock));
                  nseasons() = __max(nseasons(), nseasons(ifleet,imetier,istock));
                  nareas()   = __max(nareas(),   nareas(  ifleet,imetier,istock));
                  niters()   = __max(niters(),   niters(  ifleet,imetier,istock));
                  }  
      
   }
*/

int& flc::minage(void)
   {
   return fls_minage;
   }

int& flc::maxage(void)
   {
   return fls_maxage;
   }

int& flc::minyr(void)
   {
   return fls_minyr;
   }

int& flc::maxyr(void)
   {
   return fls_maxyr;
   }

int& flc::nunits(void)
   {
   return fls_nunits;
   }

int& flc::nseasons(void)
   {
   return fls_nseasons;
   }

int& flc::nareas(void)
   {
   return fls_nareas;
   }
 
int& flc::niters(void)
   {
   return fls_niters;
   }

int& flc::minage(int ifleet, int imetier, int istock)
   {
   return fl_minage[ifleet][imetier][istock];
   }

int& flc::maxage(int ifleet, int imetier, int istock)
   {
   return fl_maxage[ifleet][imetier][istock];
   }

int& flc::plusgrp(int ifleet, int imetier, int istock)
   {
   return fl_plusgrp[ifleet][imetier][istock];
   }

int& flc::minyr(int ifleet, int imetier, int istock)
   {
   return fl_minyr[ifleet][imetier][istock];
   }

int& flc::maxyr(int ifleet, int imetier, int istock)
   {
   return fl_maxyr[ifleet][imetier][istock];
   }

int& flc::nunits(int ifleet, int imetier, int istock)
   {
   return fl_nunits[ifleet][imetier][istock];
   }

int& flc::nseasons(int ifleet, int imetier, int istock)
   {
   return fl_nseasons[ifleet][imetier][istock];
   }

int& flc::nareas(int ifleet, int imetier, int istock)
   {
   return fl_nareas[ifleet][imetier][istock];
   }
 
int& flc::niters(int ifleet, int imetier, int istock)
   {
   return fl_niters[ifleet][imetier][istock];
   }

bool flc::InitFlag(int ifleet, int imetier, int istock)
   {
   return fl_InitFlag[ifleet][imetier][istock];
   }

int& flc::minage(int istock)
   {
   return stock_minage[istock];
   }

int& flc::maxage(int istock)
   {
   return stock_maxage[istock];
   }

int& flc::plusgrp(int istock)
   {
   return stock_plusgrp[istock];
   }

int& flc::minfbar(int istock)
   {
   return stock_minfbar[istock];
   }

int& flc::maxfbar(int istock)
   {
   return stock_maxfbar[istock];
   }

int& flc::minyr(int istock)
   {
   return stock_minyr[istock];
   }

int& flc::maxyr(int istock)
   {
   return stock_maxyr[istock];
   }

int& flc::nunits(int istock)
   {
   return stock_nunits[istock];
   }

int& flc::nseasons(int istock)
   {
   return stock_nseasons[istock];
   }

int& flc::nareas(int istock)
   {
   return stock_nareas[istock];
   }
 
int& flc::niters(int istock)
   {
   return stock_niters[istock];
   }

int& flc::nfleet(void)
   {
   return _nfleet;
   }

int& flc::nmetier(void)
   {
   return _nmetier;
   }

int& flc::nstock(void)
   {
   return _nstock;
   }

double flc::SSB(int istock, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
   double ssb=0.0, val;
   for (int iage=minage(istock); iage<=maxage(istock); iage++)
      {
      val =                n(istock, iage, iyr, iUnit, iSeason, iArea, iIter)*
                    stock_wt(istock, iage, iyr, iUnit, iSeason, iArea, iIter)*
                    fec(     istock, iage, iyr, iUnit, iSeason, iArea, iIter)*
              exp(-m(        istock, iage, iyr, iUnit, iSeason, iArea, iIter)*
                   m_spwn(   istock, iage, iyr, iUnit, iSeason, iArea, iIter)
                  -f(        istock, iage, iyr, iUnit, iSeason, iArea, iIter)*
                   f_spwn(   istock, iage, iyr, iUnit, iSeason, iArea, iIter));
     
      if (!R_IsNA(val) && val>0.0) ssb +=val;
      }

   return ssb;
   }  

double flc::computeStock(int iyr, int istock, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }

double flc::computeCatch(int istock, int iyr, int ifleet, int imetier, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }

double flc::partialF(int istock, int iFleet, int iMetier, int iAge, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
   if (istock   <1                           || istock    > nstock()                       ||
       iFleet <1                           || iFleet  > nfleet()                     ||
       iMetier<1                           || iMetier >nmetier()                     ||
       iAge   <minage(iFleet,iMetier,istock) || iAge    >maxage(  iFleet,iMetier,istock) ||
       iyr    <minyr( iFleet,iMetier,istock) || iyr     >maxyr(   iFleet,iMetier,istock) ||
       iUnit  <1                           || iUnit   >nunits(  iFleet,iMetier,istock) ||
       iSeason<1                           || iSeason >nseasons(iFleet,iMetier,istock) ||
       iArea  <1                           || iArea   >nareas(  iFleet,iMetier,istock) ||
       iIter  <1                           || iIter   >niters(  iFleet,iMetier,istock))
      return 0.0;
   else
      return __max(0.0,catch_sel(iFleet, iMetier, istock, iAge, iyr, iUnit, iSeason, iArea, iIter)*                                                          
                       catch_q(  iFleet, iMetier, istock, 1,    iyr, iUnit, iSeason, iArea, iIter)*
                       effshare( iFleet, iMetier,       1,    iyr, iUnit, iSeason, iArea, iIter)*
                       effort(   iFleet,                1,    iyr, iUnit, iSeason, iArea, iIter));                                                
   }

double flc::computeLandings(int iyr, int istock, int ifleet, int imetier, int iUnit, int iSeason, int iArea, int iIter)
{
return R_NaN;
}

double flc::computeDiscards(int iyr, int istock, int ifleet, int imetier, int iUnit, int iSeason, int iArea, int iIter)
{
return R_NaN;
}     

double flc::F(int istock, int iage, int iyr, int iUnit, int iSeason, int iArea, int iIter)
{
   double return_val=0.0;
   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      for (int imetier=1; imetier<=nmetier(); imetier++)
         {
         return_val +=  effort(   ifleet,                1,    iyr, iUnit, iSeason, iArea, iIter)
                       *effshare( ifleet, imetier,       1,    iyr, iUnit, iSeason, iArea, iIter)
                       *catch_q(  ifleet, imetier, istock, 1,    iyr, iUnit, iSeason, iArea, iIter) 
                       *catch_sel(ifleet, imetier, istock, iage, iyr, iUnit, iSeason, iArea, iIter);
         }
   return return_val;
}

double flc::Fbar(int istock, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
	double fbar=0.0;
	int iAge;
	for (iAge = minfbar(istock); iAge<= maxfbar(istock);iAge++)
		fbar += f(istock, iAge, iyr, iUnit,iSeason,iArea,iIter);
	return (fbar / (maxfbar(istock)-minfbar(istock)+1));
   }

double flc::FbarLandings(int istock, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
	double fbar=0.0;
	int iAge;
	for (iAge = minfbar(istock); iAge<= maxfbar(istock);iAge++)
		fbar += f(istock, iAge, iyr, iUnit,iSeason,iArea,iIter);
	return (fbar / (maxfbar(istock)-minfbar(istock)+1));
   }

double flc::FbarDiscards(int istock, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
	double fbar=0.0;
	int iAge;
	for (iAge = minfbar(istock); iAge<= maxfbar(istock);iAge++)
		fbar += FDiscards(istock, iAge, iyr, iUnit,iSeason,iArea,iIter);
	return (fbar / (maxfbar(istock)-minfbar(istock)+1));
   }

double flc::Zbar(int iyr, int istock, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }


double flc::FLandings(int istock, int iage, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
   double return_val=0.0;

   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      for (int imetier=1; imetier<=nmetier(); imetier++)
         return_val +=  effort(   ifleet,                1,    iyr, iUnit, iSeason, iArea, iIter)
                       *effshare( ifleet, imetier,       1,    iyr, iUnit, iSeason, iArea, iIter)
                       *catch_q(  ifleet, imetier, istock, 1,    iyr, iUnit, iSeason, iArea, iIter) 
                       *catch_sel(ifleet, imetier, istock, iage, iyr, iUnit, iSeason, iArea, iIter);

   return return_val;
   }

double flc::FDiscards(int istock, int iage, int iyr, int iUnit, int iSeason, int iArea, int iIter)
   {
   double return_val=0.0;

   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      for (int imetier=1; imetier<=nmetier(); imetier++)
         return_val +=  effort(   ifleet,                1,    iyr, iUnit, iSeason, iArea, iIter)
                       *effshare( ifleet, imetier,       1,    iyr, iUnit, iSeason, iArea, iIter)
                       *catch_q(  ifleet, imetier, istock, 1,    iyr, iUnit, iSeason, iArea, iIter) 
                       *catch_sel(ifleet, imetier, istock, iage, iyr, iUnit, iSeason, iArea, iIter);

   return return_val;
   }

double flc::Effort(int iyr, int ifleet, int imetier, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }

double flc::computeCosts(int iyr, int ifleet, int imetier, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }

double flc::computeRevenue(int iyr, int istock, int ifleet, int imetier, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }

double flc::computeProfit(int iyr, int ifleet, int iUnit, int iSeason, int iArea, int iIter)
   {
   return R_NaN;
   }

void flc::project(int istock, int iyr, int iter)
   {
   int iunit  =1,
       iseason=1,
       iarea  =1,
       iage;

   for (iage=minage(istock); iage<=maxage(istock); iage++)
      {
      f(istock, iage,iyr,iunit,iseason,iarea,iter) = 0.0;
      for (int ifleet=1; ifleet<=nfleet(); ifleet++)
         for (int imetier=1; imetier<=nmetier(); imetier++)
             f(istock, iage,iyr,iunit,iseason,iarea,iter) += __max(0.0,catch_sel(  ifleet, imetier, istock, iage, iyr, iunit, iseason, iarea, iter)*
                                                                        catch_q( ifleet, imetier, istock, 1,    iyr, iunit, iseason, iarea, iter)*
                                                                        effort(  ifleet,                1,    iyr, iunit, iseason, iarea, iter)*
                                                                        effshare(ifleet, imetier,       1,    iyr, iunit, iseason, iarea, iter));
      //numbers-at-age next year
      if (iage < maxage(istock))
         n(istock,iage+1,iyr+1,iunit,iseason,iarea,iter)  = n(istock, iage, iyr,iunit,iseason,iarea,iter)*exp(-f(istock,iage,iyr,iunit,iseason,iarea,iter)-m(istock,iage,iyr,iunit,iseason,iarea,iter));
      else if (iage == plusgrp(istock))
         n(istock,iage,  iyr+1,iunit,iseason,iarea,iter) += n(istock, iage, iyr,iunit,iseason,iarea,iter)*exp(-f(istock,maxage(istock),iyr,iunit,iseason,iarea,iter)-m(istock,maxage(istock),iyr,iunit,iseason,iarea,iter));
      }

   int rec_yr = __min(__max(iyr-minage(istock)+1,minyr(istock)),maxyr(istock));
      
   n(istock,minage(istock),iyr+1,iunit,iseason,iarea,iter) = _sr.recruits(istock,rec_yr,flc::SSB(istock,rec_yr,iunit,iseason,iarea,iter));

   for (iage=minage(istock); iage<=maxage(istock); iage++)
      {
      double z = m(istock, iage, iyr, iunit, iseason, iarea, iter) + f(istock, iage, iyr, iunit, iseason, iarea, iter);

      catch_n( istock, iage, iyr, iunit, iseason, iarea, iter) =
             n( istock, iage, iyr, iunit, iseason, iarea, iter)*
             f( istock, iage, iyr, iunit, iseason, iarea, iter)/z*(1-exp(-z));

      discards_n( istock, iage, iyr, iunit, iseason, iarea, iter)=catch_n(     istock, iage, iyr, iunit, iseason, iarea, iter)*
                                                                discards_sel(istock, iage, iyr, iunit, iseason, iarea, iter)/
                                                                catch_sel(   istock, iage, iyr, iunit, iseason, iarea, iter);
      
      landings_n( istock, iage, iyr, iunit, iseason, iarea, iter)=catch_n(     istock, iage, iyr, iunit, iseason, iarea, iter)*
                                                                landings_sel(istock, iage, iyr, iunit, iseason, iarea, iter)/
                                                                catch_sel(   istock, iage, iyr, iunit, iseason, iarea, iter);
      }
   } 

void flc::alloc_dims_range(void)
   {
   //stock
   stock_minage   = new int[nstock()] - 1;
   stock_maxage   = new int[nstock()] - 1;
   stock_plusgrp  = new int[nstock()] - 1;
   stock_minfbar  = new int[nstock()] - 1;
   stock_maxfbar  = new int[nstock()] - 1; 
   stock_minyr    = new int[nstock()] - 1; 
   stock_maxyr    = new int[nstock()] - 1;
   stock_nunits   = new int[nstock()] - 1;
   stock_nseasons = new int[nstock()] - 1;
   stock_nareas   = new int[nstock()] - 1;
   stock_niters   = new int[nstock()] - 1;

   //fleet, metier, stock
   fl_InitFlag = new bool**[nfleet()] - 1; 
   fl_minage   = new  int**[nfleet()] - 1; 
   fl_maxage   = new  int**[nfleet()] - 1; 
   fl_plusgrp  = new  int**[nfleet()] - 1; 
   fl_minyr    = new  int**[nfleet()] - 1; 
   fl_maxyr    = new  int**[nfleet()] - 1; 
   fl_nunits   = new  int**[nfleet()] - 1; 
   fl_nseasons = new  int**[nfleet()] - 1; 
   fl_nareas   = new  int**[nfleet()] - 1; 
   fl_niters   = new  int**[nfleet()] - 1; 
   for (int iFleet=1; iFleet<=nfleet(); iFleet++)
      {
      fl_InitFlag[iFleet] = new bool*[nmetier()] - 1; 
      fl_minage[iFleet]   = new  int*[nmetier()] - 1; 
      fl_maxage[iFleet]   = new  int*[nmetier()] - 1; 
      fl_plusgrp[iFleet]  = new  int*[nmetier()] - 1; 
      fl_minyr[iFleet]    = new  int*[nmetier()] - 1; 
      fl_maxyr[iFleet]    = new  int*[nmetier()] - 1; 
      fl_nunits[iFleet]   = new  int*[nmetier()] - 1; 
      fl_nseasons[iFleet] = new  int*[nmetier()] - 1; 
      fl_nareas[iFleet]   = new  int*[nmetier()] - 1; 
      fl_niters[iFleet]   = new  int*[nmetier()] - 1; 
      for (int iMetier=1; iMetier<=nmetier(); iMetier++)
         {
         fl_InitFlag[iFleet][iMetier] = new bool[nstock()] - 1; 
         fl_minage[iFleet][iMetier]   = new  int[nstock()] - 1; 
         fl_maxage[iFleet][iMetier]   = new  int[nstock()] - 1; 
         fl_plusgrp[iFleet][iMetier]  = new  int[nstock()] - 1; 
         fl_minyr[iFleet][iMetier]    = new  int[nstock()] - 1; 
         fl_maxyr[iFleet][iMetier]    = new  int[nstock()] - 1; 
         fl_nunits[iFleet][iMetier]   = new  int[nstock()] - 1; 
         fl_nseasons[iFleet][iMetier] = new  int[nstock()] - 1; 
         fl_nareas[iFleet][iMetier]   = new  int[nstock()] - 1; 
         fl_niters[iFleet][iMetier]   = new  int[nstock()] - 1; 
         for (int iStock=1; iStock<=nstock(); iStock++)
             fl_InitFlag[iFleet][iMetier][iStock] = FALSE;
         }
      }
   }

void flc::unalloc_dims_range(void)
   {
   if (nstock()<1 || nfleet()<1 || nmetier()<1)
      return;

   //fleet, metier, stock
   for (int iFleet=1; iFleet<=nfleet(); iFleet++)
      {
      for (int iMetier=1; iMetier<=nmetier(); iMetier++)
         {
         delete [] (fl_InitFlag[iFleet][iMetier] + 1); 
         delete [] (fl_minage[  iFleet][iMetier] + 1); 
         delete [] (fl_maxage[  iFleet][iMetier] + 1); 
         delete [] (fl_plusgrp[ iFleet][iMetier] + 1); 
         delete [] (fl_minyr[   iFleet][iMetier] + 1); 
         delete [] (fl_maxyr[   iFleet][iMetier] + 1); 
         delete [] (fl_nunits[  iFleet][iMetier] + 1); 
         delete [] (fl_nseasons[iFleet][iMetier] + 1); 
         delete [] (fl_nareas[  iFleet][iMetier] + 1); 
         delete [] (fl_niters[  iFleet][iMetier] + 1); 
         }

      delete [] (fl_InitFlag[iFleet] + 1); 
      delete [] (fl_minage[  iFleet] + 1);
      delete [] (fl_maxage[  iFleet] + 1);
      delete [] (fl_plusgrp[ iFleet] + 1);
      delete [] (fl_minyr[   iFleet] + 1);
      delete [] (fl_maxyr[   iFleet] + 1);
      delete [] (fl_nunits[  iFleet] + 1);
      delete [] (fl_nseasons[iFleet] + 1);
      delete [] (fl_nareas[  iFleet] + 1);
      delete [] (fl_niters[  iFleet] + 1);
      }

   delete [] (fl_InitFlag + 1); 
   delete [] (fl_minage   + 1);
   delete [] (fl_maxage   + 1);
   delete [] (fl_plusgrp  + 1);
   delete [] (fl_minyr    + 1);
   delete [] (fl_maxyr    + 1);
   delete [] (fl_nunits   + 1);
   delete [] (fl_nseasons + 1);
   delete [] (fl_nareas   + 1);
   delete [] (fl_niters   + 1);
   }

void flc::InitBiols(SEXP xBiol)
   {
   n.alloc_n7(       nstock());
   m.alloc_n7(       nstock());
   f.alloc_n7(       nstock());
   stock_wt.alloc_n7(nstock());
   fec.alloc_n7(     nstock());
   m_spwn.alloc_n7(  nstock());
   f_spwn.alloc_n7(  nstock());

   for (int iStock=1; iStock<=nstock(); iStock++)
      {
      SEXP x      = PROTECT(VECTOR_ELT(xBiol, iStock-1)),
           xRange = PROTECT(duplicate(GET_SLOT(x,  install("range"))));
        
      stock_minage[ iStock] = (int)REAL(xRange)[0];
      stock_maxage[ iStock] = (int)REAL(xRange)[1];
      stock_plusgrp[iStock] = (int)REAL(xRange)[2];
      stock_minyr[  iStock] = (int)REAL(xRange)[3];
      stock_maxyr[  iStock] = (int)REAL(xRange)[4];
   
      if (LENGTH(xRange) >= 6) 
         stock_minfbar[iStock] = (int)REAL(xRange)[5];
      if (LENGTH(xRange) >= 7) 
         stock_maxfbar[iStock] = (int)REAL(xRange)[6];
 
      SEXP xM    = PROTECT(duplicate(GET_SLOT(x,  install("m")))),
           Quant = PROTECT(duplicate(GET_SLOT(xM, install(".Data")))),
           dims  = GET_DIM(Quant);

      stock_nunits[  iStock] = INTEGER(dims)[2];
      stock_nseasons[iStock] = INTEGER(dims)[3];
      stock_nareas[  iStock] = INTEGER(dims)[4];
      stock_niters[  iStock] = INTEGER(dims)[5]; 
      UNPROTECT(4);
      
      n.Init(       iStock, GET_SLOT(x, install("n")));
      m.Init(       iStock, GET_SLOT(x, install("m")));
      f.Init(       iStock, GET_SLOT(x, install("m")));
      stock_wt.Init(iStock, GET_SLOT(x, install("wt")));
      fec.Init(     iStock, GET_SLOT(x, install("fec")));
      m_spwn.Init(  iStock, GET_SLOT(x, install("spwn")));
      f_spwn.Init(  iStock, GET_SLOT(x, install("spwn")));
      }
   }

void flc::InitFleets(SEXP xFleet)
   {
   effort.alloc_n7(   nfleet());
   fcost.alloc_n7(    nfleet());
   capacity.alloc_n7( nfleet());
   crewshare.alloc_n7(nfleet());

   effshare.Init(nfleet(), nmetier());
   vcost.Init(   nfleet(), nmetier());

   catch_.alloc_n9(      nfleet());
   catch_n.alloc_n9(     nfleet());
   catch_wt.alloc_n9(    nfleet());
   catch_q.alloc_n9(     nfleet());
   catch_sel.alloc_n9(   nfleet());

   landings.alloc_n9(    nfleet());
   landings_n.alloc_n9(  nfleet());
   landings_wt.alloc_n9( nfleet()); 
   landings_sel.alloc_n9(nfleet());  
   discards.alloc_n9(    nfleet());
   discards_n.alloc_n9(  nfleet());
   discards_wt.alloc_n9( nfleet()); 
   discards_sel.alloc_n9(nfleet());  
   price.alloc_n9(       nfleet());

   for (int ifleet=1; ifleet<=nfleet(); ifleet++)
      {
      catch_.alloc_n8(      ifleet, nmetier());
      catch_n.alloc_n8(     ifleet, nmetier());
      catch_wt.alloc_n8(    ifleet, nmetier());
      catch_q.alloc_n8(     ifleet, nmetier());
      catch_sel.alloc_n8(   ifleet, nmetier());
      landings.alloc_n8(    ifleet, nmetier());
      landings_n.alloc_n8(  ifleet, nmetier());
      landings_wt.alloc_n8( ifleet, nmetier()); 
      landings_sel.alloc_n8(ifleet, nmetier());  
      discards.alloc_n8(    ifleet, nmetier());
      discards_n.alloc_n8(  ifleet, nmetier());
      discards_wt.alloc_n8( ifleet, nmetier()); 
      discards_sel.alloc_n8(ifleet, nmetier());  
      price.alloc_n8(       ifleet, nmetier());

      for (int imetier=1; imetier<=nmetier(); imetier++)
         {
         catch_.alloc_n7(      ifleet, imetier, nstock());
         catch_n.alloc_n7(     ifleet, imetier, nstock());
         catch_wt.alloc_n7(    ifleet, imetier, nstock());
         catch_q.alloc_n7(     ifleet, imetier, nstock());
         catch_sel.alloc_n7(   ifleet, imetier, nstock());
         landings.alloc_n7(    ifleet, imetier, nstock());
         landings_n.alloc_n7(  ifleet, imetier, nstock());
         landings_wt.alloc_n7( ifleet, imetier, nstock()); 
         landings_sel.alloc_n7(ifleet, imetier, nstock());  
         discards.alloc_n7(    ifleet, imetier, nstock());
         discards_n.alloc_n7(  ifleet, imetier, nstock());
         discards_wt.alloc_n7( ifleet, imetier, nstock()); 
         discards_sel.alloc_n7(ifleet, imetier, nstock());  
         price.alloc_n7(       ifleet, imetier, nstock());
         }
      }

   for (int iFleet=1; iFleet<=nfleet(); iFleet++)
      {
      SEXP fleet   = PROTECT(VECTOR_ELT(xFleet, iFleet-1)),
           metiers = PROTECT(duplicate(GET_SLOT(fleet, install("metiers"))));
         
      effort.Init(   iFleet, GET_SLOT(fleet, install("effort")));
      fcost.Init(    iFleet, GET_SLOT(fleet, install("fcost")));
      capacity.Init( iFleet, GET_SLOT(fleet, install("capacity")));
      crewshare.Init(iFleet, GET_SLOT(fleet, install("crewshare")));
     
      for (int iMetier=1; iMetier<=nmetier(); iMetier++)
         {
         SEXP metier  = PROTECT(VECTOR_ELT(metiers, iMetier-1)),
              catches = PROTECT(duplicate(GET_SLOT(metier, install("catches"))));
      
         effshare.Init(iFleet, iMetier, GET_SLOT(metier, install("effshare")));
         vcost.Init(   iFleet, iMetier, GET_SLOT(metier, install("vcost")));

         for (int iStock=1; iStock<=nstock(); iStock++)
            {	
            SEXP _catch = PROTECT(VECTOR_ELT(catches, iStock-1));

            catch_q.Init(     iFleet, iMetier, iStock, GET_SLOT(_catch, install("catch.q")));  
            catch_.Init(      iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings"))); 
            catch_n.Init(     iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings.n")));  
            catch_wt.Init(    iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings.wt")));  
            catch_sel.Init(   iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings.sel")));  
            landings.Init(    iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings")));  
            landings_n.Init(  iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings.n")));  
            landings_wt.Init( iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings.wt")));  
            landings_sel.Init(iFleet, iMetier, iStock, GET_SLOT(_catch, install("landings.sel")));  
            discards.Init(    iFleet, iMetier, iStock, GET_SLOT(_catch, install("discards")));  
            discards_n.Init(  iFleet, iMetier, iStock, GET_SLOT(_catch, install("discards.n")));  
            discards_wt.Init( iFleet, iMetier, iStock, GET_SLOT(_catch, install("discards.wt")));  
            discards_sel.Init(iFleet, iMetier, iStock, GET_SLOT(_catch, install("discards.sel")));  
            price.Init(       iFleet, iMetier, iStock, GET_SLOT(_catch, install("price"))); 

            for(int iStock=1; iStock<=nstock(); iStock++)
              for (int iIter = 1; iIter<=niters(iStock); iIter++)
	             for (int iArea = 1; iArea <= nareas(iStock); iArea++)
		            for (int iSeason = 1; iSeason <= nseasons(iStock); iSeason++)
    		           for (int iUnit = 1; iUnit <= nunits(iStock); iUnit++)
      		          for (int iYr = minyr(iStock); iYr <= maxyr(iStock); iYr++)
		   	            for(int iAge=minage(iStock);iAge<=maxage(iStock);iAge++)
                           {
                           catch_(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter) +=
                                 discards(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter);

                           catch_n(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter) +=
                                 discards_n(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter);

                           catch_wt(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter)=
                                    (landings_sel(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter)*
                                     landings_wt( iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter)+
                                     discards_sel(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter)*
                                     discards_wt( iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter))/
                                    (landings_sel(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter)+
                                     discards_sel(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter));                           
                           catch_sel(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter) +=
                                 discards_sel(iFleet,iMetier,iStock,iAge,iYr,iUnit,iSeason,iArea,iIter);                           
                           }

            UNPROTECT(1);
            }
         UNPROTECT(2);
         }      
	  UNPROTECT(2);
      }
	}

bool flc::InitBiolsFleets(SEXP xBiols, SEXP xFleets, SEXP xDim)
   {
   if (!isFLBiols(xBiols) && !isFLFleets(xFleets)) 
     return false;
   
   nstock() = NElemList(xBiols);
   nfleet() = NElemList(xFleets);
//   nmetier()= NElemList(metiers);
 
   SEXP fleet   = PROTECT(VECTOR_ELT(xFleets, 0));
   SEXP metiers = PROTECT(GET_SLOT(fleet, install("metiers")));

   if (nfleet() < 1 || nmetier() < 1 || nstock() < 1)
      {
      UNPROTECT(2);

      return false;
      }

   SEXP dim;
   PROTECT(dim = AS_NUMERIC(xDim));
   if (LENGTH(dim)<6)
      {
      UNPROTECT(3);

      return false;
      }
        
   minyr()    = (int)REAL(dim)[0];
   maxyr()    = (int)REAL(dim)[1];
   nunits()   = (int)REAL(dim)[2];
   nseasons() = (int)REAL(dim)[3];
   nareas()   = (int)REAL(dim)[4];
   niters()   = (int)REAL(dim)[5];
   int flag   = (int)REAL(dim)[6];
  
   alloc_dims_range();

   InitBiols(xBiols); 
   InitFleets(xFleets);

   CalcF(1);

   UNPROTECT(3);

   return true;
   }
