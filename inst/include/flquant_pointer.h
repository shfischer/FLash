#ifndef _INC_flquant_pointer
#define _INC_flquant_pointer

#include "FLCoreClasses.h"

extern int             outofbounds_int;
extern bool            outofbounds_bool;
extern double          outofbounds_double;

class FLQuant_pointer
{
public:        
   FLQuant_pointer(void);      
   FLQuant_pointer(SEXP);
  ~FLQuant_pointer(void);      

   void Init(SEXP);      
   SEXP Return(void);      

   inline int&  minquant()  { return flq_minquant; }
   inline int&  maxquant()  { return flq_maxquant; }
   inline int&  plusgrp()   { return flq_plusgrp;  }
   inline int&  minyr()     { return flq_minyr;    }
   inline int&  maxyr()     { return flq_maxyr;    }
   inline int&  nunits()    { return flq_nunits;   }
   inline int&  nseasons()  { return flq_nseasons; }
   inline int&  nareas()    { return flq_nareas;   }
   inline int&  niters()    { return flq_niters;   }
   inline bool& InitFlag()  { return flq_InitFlag; }

   int                i(int _age, int _yr, int _unit=1, int _season=1, int _area=1, int _iter=1);
   double& operator () (int _age, int _yr, int _unit=1, int _season=1, int _area=1, int _iter=1);
protected: 
   bool flq_InitFlag;

   int flq_minquant, flq_maxquant, flq_plusgrp,
       flq_minyr,    flq_maxyr,
       flq_nunits,
       flq_nseasons,
       flq_nareas,
       flq_niters;

   double *data;   
   };                  

#endif /* _INC_flquant_pointer */
