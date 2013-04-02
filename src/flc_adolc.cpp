#include "flc_adolc.h"

FLQuant_adolc::FLQuant_adolc(void)
   {
   ; 
   }

FLQuant_adolc::FLQuant_adolc(FLQuant& _flq, int arg_minyr, int arg_maxyr, int _iter)
   {
   Init(_flq, arg_minyr, arg_maxyr, _iter); 
   }

void FLQuant_adolc::Init(FLQuant& _flq, int arg_minyr, int arg_maxyr, int _iter)
   {
   // set pointer to address of initialising quant
   flq = &_flq;
   _minyr  = arg_minyr;
   _maxyr  = __min(flq->maxyr(),arg_maxyr);
   _miniter= __max(1,_iter);
   _maxiter= __min(flq->niters(),_iter);

   int nvar = (_maxiter-_miniter+1)*flq->nareas()*flq->nseasons()*flq->nunits()*(_maxyr-_minyr+1)*(flq->maxquant()-flq->minquant()+1);
   data  = new adouble[nvar];
   
   // fill up array with values from initialising quant
   for(int iIter = _miniter; iIter <= _maxiter; iIter++)
	   for (int iArea = 1; iArea <= nareas(); iArea++)
	      for (int iSeason = 1; iSeason <= nseasons(); iSeason++)
		      for (int iUnit = 1; iUnit <= nunits(); iUnit++)
		         for (int iYear = _minyr; iYear <= _maxyr; iYear++)
			         for (int iAge = minquant(); iAge <= maxquant(); iAge++)
			            if (iYear>=flq->minyr() && iYear<=flq->maxyr()) 
							data[i(iAge,iYear,iUnit,iSeason,iArea,iIter)] = (*flq)(iAge,iYear,iUnit,iSeason,iArea,iIter);
                        else
						    data[i(iAge,iYear,iUnit,iSeason,iArea,iIter)] = 0.0;
   }

     
int FLQuant_adolc::i(int _age, int _yr, int _unit, int _season, int _area, int _iter)
   {
   int index = (((((_iter-_miniter)*nareas()+_area-1)*nseasons()+_season-1)*nunits()+_unit-1)*(maxyr()-minyr()+1)+_yr-minyr())*(maxquant()-minquant()+1)+_age-minquant();
   
   return index;
   }

adouble& FLQuant_adolc::operator () (int _age, int _yr, int _unit, int _season, int _area, int _iter)
   {
   return data[i(_age,_yr,_unit,_season,_area,_iter)];
   }

FLQuant_adolc::~FLQuant_adolc(void)
   {
   delete[] data;
   }

//Only return 5 dim array 
SEXP FLQuant_adolc::Return(void)      
    {
    SEXP Quant, v, 
         d1, d2, d3, d4, d5, d6, 
         dim, dimnames, names;    

    int i, iAge, iYear, iUnit, iArea, iSeason, iIter;

    //Create new S4 object    

    PROTECT(Quant = NEW_OBJECT(MAKE_CLASS("FLQuant")));

    //Create array for slot    
    //Set dimensions of array
    PROTECT(dim     = allocVector(INTSXP, 6));       
    INTEGER(dim)[0] = maxquant()-minquant() +1;
    INTEGER(dim)[1] = maxyr()   -minyr()    +1;
    INTEGER(dim)[2] = nunits(); 
    INTEGER(dim)[3] = nseasons(); 
    INTEGER(dim)[4] = nareas();
    INTEGER(dim)[5] = niters();

   //allocate memory
    PROTECT(v = Rf_allocArray(REALSXP, dim)); 
    
    //Create dimension names
    PROTECT(dimnames = allocVector(VECSXP, 6));
    
    PROTECT(d1 = allocVector(INTSXP, maxquant()-minquant() +1));
    for (iAge=minquant(),i=0; iAge<=maxquant(); iAge++, i++)
        INTEGER(d1)[i] = iAge; 
    SET_VECTOR_ELT(dimnames, 0, d1);
    
    PROTECT(d2 = allocVector(INTSXP, maxyr()-minyr()+1));
    for (iYear=minyr(), i=0; iYear<=maxyr(); iYear++, i++)
        INTEGER(d2)[i] = iYear; 
    SET_VECTOR_ELT(dimnames, 1, d2);
     
    if (nunits()==1)
       {
       PROTECT(d3 = allocVector(STRSXP, nunits()));
       SET_STRING_ELT(d3, 0, mkChar("unique"));
       }
    else
       {
       PROTECT(d3 = allocVector(INTSXP, nunits()));
       for (iUnit=1, i=0; iUnit<=nunits(); iUnit++, i++)
          INTEGER(d3)[i] = iUnit; 
       }
    SET_VECTOR_ELT(dimnames, 2, d3);
       
    if (nseasons()==1)
       {
       PROTECT(d4 = allocVector(STRSXP, nseasons()));
       SET_STRING_ELT(d4, 0, mkChar("all"));
       }
    else
       {
       PROTECT(d4 = allocVector(INTSXP, nseasons()));
       for (iSeason=1, i=0; iSeason<=nseasons(); iSeason++, i++)
          INTEGER(d4)[i] = iSeason; 
       }
    SET_VECTOR_ELT(dimnames, 3, d4);
    

    if (nareas()==1)
       {
       PROTECT(d5 = allocVector(STRSXP, nareas()));
       SET_STRING_ELT(d5, 0, mkChar("unique"));
       }
    else
       {
       PROTECT(d5 = allocVector(INTSXP, nareas()));
       for (iArea=1, i=0; iArea<=nareas(); iArea++, i++)
          INTEGER(d5)[i] = iArea; 
       }
    SET_VECTOR_ELT(dimnames, 4, d5);

    PROTECT(d6 = allocVector(INTSXP, niters()));
    for (iIter=1, i=0; iIter<=niters(); iIter++, i++)
        INTEGER(d6)[i] = iIter; 
    SET_VECTOR_ELT(dimnames, 5, d6);
    
    //Create names for dimensions
    PROTECT(names = allocVector(STRSXP, 6));
    SET_STRING_ELT(names, 0, mkChar("age"));
    SET_STRING_ELT(names, 1, mkChar("year"));
    SET_STRING_ELT(names, 2, mkChar("unit"));
    SET_STRING_ELT(names, 3, mkChar("season"));
    SET_STRING_ELT(names, 4, mkChar("area"));
    SET_STRING_ELT(names, 5, mkChar("iter")); 

    setAttrib(dimnames, R_NamesSymbol, names);
    setAttrib(v, R_DimNamesSymbol, dimnames);
   
    //Set data
    i=0;
    for(iIter = 1; iIter <= niters(); iIter++)
	    for (iArea = 1; iArea <= nareas(); iArea++)
	  	    for (iSeason = 1; iSeason <= nseasons(); iSeason++)
     		    for (iUnit = 1; iUnit <= nunits(); iUnit++)
	    		    for (iYear = minyr(); iYear <= maxyr(); iYear++)
			 		    for (iAge = minquant(); iAge <= maxquant(); iAge++)
			      			    REAL(v)[i++] = (double) data[FLQuant_adolc::i(iAge,iYear,iUnit,iSeason,iArea,iIter)].value(); 
                   
    //Set slot
    Quant = R_do_slot_assign(Quant, install(".Data"), v);

    UNPROTECT(11);
    
    return Quant;
    }

FLQuant2_adolc::FLQuant2_adolc(FLQuant2& _flq2, int _i7, int arg_minyr, int arg_maxyr, int _iter)
   {
   // set pointer to address of initialising Quant2
   flq2 = &_flq2;
   _minyr  = __max(flq2->minyr(_i7),arg_minyr);
   _maxyr  = __min(flq2->maxyr(_i7),arg_maxyr);
   _miniter= __max(1,_iter);
   _maxiter= __min(flq2->niters(_i7),_iter);
   _min7   = _i7;
   _max7   = _i7;

   int nvar = (_max7-_min7+1)*(_maxiter-_miniter+1)*flq2->nareas(_i7)*flq2->nseasons(_i7)*flq2->nunits(_i7)*(_maxyr-_minyr+1)*(flq2->maxquant(_i7)-flq2->minquant(_i7)+1);
   data  = new adouble[nvar];
   
   // fill up array with values from initialising Quant2
   for(int iIter = _miniter; iIter <= _maxiter; iIter++)
	   for (int iArea = 1; iArea <= nareas(_i7); iArea++)
	      for (int iSeason = 1; iSeason <= nseasons(_i7); iSeason++)
		      for (int iUnit = 1; iUnit <= nunits(_i7); iUnit++)
		         for (int iYear = _minyr; iYear <= _maxyr; iYear++)
			         for (int iAge = minquant(_i7); iAge <= maxquant(_i7); iAge++)
			            data[i(_i7,iAge,iYear,iUnit,iSeason,iArea,iIter)] = (*flq2)(_i7,iAge,iYear,iUnit,iSeason,iArea,iIter);
   }

     
int FLQuant2_adolc::i(int _i7, int _age, int _yr, int _unit, int _season, int _area, int _iter)
   {
   int index = ((((((_i7-_min7)*(_iter-_miniter)*nareas(_i7)+_area-1)*nseasons(_i7)+_season-1)*nunits(_i7)+_unit-1)*(maxyr()-minyr()+1)+_yr-minyr())*(maxquant(_i7)-minquant(_i7)+1)+_age-minquant(_i7))+_i7-_min7;
   
   return index;
   }

adouble& FLQuant2_adolc::operator () (int _i7, int _age, int _yr, int _unit, int _season, int _area, int _iter)
   {
   return data[i(_i7, _age,_yr,_unit,_season,_area,_iter)];
   }

FLQuant2_adolc::~FLQuant2_adolc(void)
   {
   delete[] data;
   }

