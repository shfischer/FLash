#ifndef _INC_FLCore_adolc
#define _INC_FLCore_adolc

#include <adolc.h>      // use of ALL ADOL-C interfaces
#include "FLCoreClasses.h"
 
class FLQuant_adolc
{
public:
   FLQuant_adolc(void);
   FLQuant_adolc(FLQuant&, int arg_minyr, int arg_maxyr, int _iter);
  ~FLQuant_adolc(void);

   void Init(FLQuant&, int arg_minyr, int arg_maxyr, int _iter);
  
   adouble& operator () (int _age, int _yr, int _unit=1, int _season=1, int _area=1, int _iter=1);

   inline int minquant() {return flq->minquant();}
   inline int maxquant() {return flq->maxquant();}
   inline int minyr()    {return _minyr;}
   inline int maxyr()    {return _maxyr;}
   inline int nunits()   {return flq->nunits();}
   inline int nareas()   {return flq->nareas();}
   inline int nseasons() {return flq->nseasons();}
   inline int niters()   {return __max(__min(flq->niters(),_maxiter),_miniter);}

   SEXP Return(void);      
protected:
   int i(int _age, int _yr, int _unit=1, int _season=1, int _area=1, int _iter=1);

   int _minyr,   _maxyr,
       _miniter, _maxiter;

   FLQuant *flq;
 
   adouble* data;
   };                  

class FLQuant2_adolc
{
public:
  ~FLQuant2_adolc(void);
   FLQuant2_adolc(FLQuant2 &, int _n7, int arg_minyr, int arg_maxyr, int _iter);
 
   adouble& operator () (int i7, int _age, int _yr, int _unit=1, int _season=1, int _area=1, int _iter=1);

   inline int minquant(int i7=1) {return flq2->minquant(i7);}
   inline int maxquant(int i7=1) {return flq2->maxquant(i7);}
   inline int minyr(void)        {return _minyr;}
   inline int maxyr(void)        {return _maxyr;}
   inline int nunits(  int i7=1) {return flq2->nunits(i7);}
   inline int nareas(  int i7=1) {return flq2->nareas(i7);}
   inline int nseasons(int i7=1) {return flq2->nseasons(i7);}
   inline int niters(  int i7=1) {return __max(__min(flq2->niters(i7),_maxiter),_miniter);}
   inline int min7(void)         {return _min7;}
   inline int max7(void)         {return _max7;}
 
protected:
   int i(int _q, int _age, int _yr, int _unit=1, int _season=1, int _area=1, int _iter=1);
   int _minyr,   _maxyr,
       _min7,    _max7,
       _miniter, _maxiter;

   FLQuant2 *flq2;
 
   adouble* data;
   };                  

#endif /* _INC_FLCore_adolc */

