#ifndef _INC_flc
#define _INC_flc

#include "FLCoreClasses.h"
#include <Rmath.h>

FLRConstSRR get_sr(int i);

int get_sr(FLRConstSRR i);

void sv2ab(FLRConstSRR model, double ***, int, int);

double F_func(double* f,double m, double c, double n);

class sr
{
public:        
   sr(void);      
   sr(int);
  ~sr(void);      

   inline int minyear() {return _minyear;} 
   inline int maxyear() {return _maxyear;} 
   
   bool   Init(int _nstock, SEXP xyrs);
   bool   Init(int _nstock, int _minyr, int _maxyr);     
   bool   Init(int istock, SEXP xmodel, SEXP xparam, SEXP xresiduals, SEXP xmult);
   bool   InitBiolsFleets(SEXP, SEXP, SEXP);
   
   double recruits(int istock, double ssb, int iyr, int iunit=1, int iseason=1, int iarea=1, int iter=1);

protected: 
   FLRConstSRR       **model;
   FLQuant2            residuals;
   FLQuant2            param;
   bool               *residuals_mult;

   int nstock,
	  _minyear, _maxyear;
  
   void unalloc(void);
   };                  

class flc
{
public:        
   flc(void);      
  ~flc(void);      

   void Init(      SEXP);      
   void InitStock( SEXP);      
   void InitStocks(SEXP);      
   void InitBiol(  SEXP);      
   void InitBiols( SEXP);      
   void InitFleet( SEXP);      
   void InitFleets(SEXP);    
  
   bool InitBiolFleet(  SEXP, SEXP, SEXP);
   bool InitBiolsFleets(SEXP, SEXP, SEXP);

   void InitStock(int, int, int, SEXP);      
   void InitBiol( int, int, int, SEXP);      

   bool InitSR(int, SEXP); 
   bool InitSR(int, SEXP, SEXP, SEXP, SEXP); 

   bool InitSR(SEXP); 
   bool InitSR(SEXP, SEXP, SEXP, SEXP); 

   void project(int,int,int);

   SEXP ReturnRange(int);
   SEXP ReturnRangeBiol(int);
   SEXP ReturnStock(int istock=1);      
   SEXP ReturnStocks(void);      
   SEXP ReturnBiol(  int istock=1);      
   SEXP ReturnBiols( void);      
   SEXP ReturnFleet(int ifleet=1);      
   SEXP ReturnFleets(void);      

  int  &nfleet(  void),
       &nmetier( void),
       &nstock(  void),
       
       //range
       &minage(  void), &maxage(void), &plusgrp(void),
       &minyr(   void), &maxyr( void),
       &nunits(  void),
       &nseasons(void),
       &nareas(  void),
       &niters(  void),

       //stock
       &minage(  int), &maxage(int), &plusgrp(int), &minfbar(int), &maxfbar(int), 
       &minyr(   int), &maxyr( int),
       &nunits(  int),
       &nseasons(int),
       &nareas(  int),
       &niters(  int),

       //fleet, metier, stock 
       &minage(  int,int,int), &maxage(int,int,int), &plusgrp(int,int,int),
       &minyr(   int,int,int), &maxyr( int,int,int),
       &nunits(  int,int,int),
       &nseasons(int,int,int),
       &nareas(  int,int,int),
       &niters(  int,int,int);

   bool InitFlag(int,int,int);

   sr _sr;

   FLQuant2 n,
            m,
            f,
            stock_wt,
            fec,
            m_spwn,
            f_spwn;   

   FLQuant2 effort,
            fcost,
            capacity,
            crewshare;

   FLQuant3 effshare,
            vcost;

   FLQuant4 catch_, 
            catch_n, 
            catch_wt, 
            catch_q, 
            catch_sel, 
            landings, 
            landings_n, 
            landings_wt, 
            landings_sel, 
            discards, 
            discards_n, 
            discards_wt, 
            discards_sel, 
            price;

   double SSB(            int istock, int iyr,                              int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double computeStock(   int istock, int iyr,                              int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double Zbar(           int istock, int iyr,                              int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double Fbar(           int istock, int iyr,                              int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double FbarLandings(   int istock, int iyr,                              int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double FbarDiscards(   int istock, int iyr,                              int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   
   double computeCatch(   int istock, int iyr, int ifleet=0, int iMetier=0, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double computeLandings(int istock, int iyr, int iFleet=0, int iMetier=0, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double computeDiscards(int istock, int iyr, int iFleet=0, int iMetier=0, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double computeRevenue( int istock, int iyr, int iFleet=0, int iMetier=0, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   
   void   CalcF(void);
   void   CalcF(int istock);
   double F(        int istock, int iage, int iyr, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double FLandings(int istock, int iage, int iyr, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double FDiscards(int istock, int iage, int iyr, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   
   double Effort(       int iyr, int iFleet=0, int iMetier=0, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double computeCosts( int iyr, int iFleet=0, int iMetier=0, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);
   double computeProfit(int iyr, int iFleet=0,                int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);

   double partialF(int istock, int ifleet, int imetier, int iage, int iyr, int iUnit=1, int iSeason=1, int iArea=1, int iIter=1);

protected: 
   bool ***fl_InitFlag;

   int _nfleet,
       _nmetier,
       _nstock,
       
       //range
       fls_minage, fls_maxage,
       fls_minyr,  fls_maxyr,
       fls_nunits,
       fls_nseasons,
       fls_nareas,
       fls_niters,

       //stock
       *stock_minage, *stock_maxage, *stock_plusgrp, *stock_minfbar, *stock_maxfbar, 
       *stock_minyr,  *stock_maxyr,
       *stock_nunits,
       *stock_nseasons,
       *stock_nareas,
       *stock_niters,

       //fleet, metier, stock
       ***fl_minage, ***fl_maxage, ***fl_plusgrp,
       ***fl_minyr,  ***fl_maxyr,
       ***fl_nunits,
       ***fl_nseasons,
       ***fl_nareas,
       ***fl_niters;
         
  double biomass(FLQuant2 &,FLQuant2 &, int,      int, int iunit=1, int iseason=1, int iarea=1, int iter=1);

  double sum(    FLQuant4 &,            int, int, int, int iunit=1, int iseason=1, int iarea=1, int iter=1);
  double wt(     FLQuant4 &,FLQuant4 &, int, int, int, int iunit=1, int iseason=1, int iarea=1, int iter=1);
  double biomass(FLQuant4 &,FLQuant4 &, int,      int, int iunit=1, int iseason=1, int iarea=1, int iter=1);

  void   alloc_dims_range(void);
  void unalloc_dims_range(void);

  void   alloc_dims_biol(int);
  void unalloc_dims_biol(int);

//  void alloc_dims_nfleets(void);
//  void alloc_dims_nmetiers(int);

//  void unalloc_dims(void);
//  void unalloc_dims_nstock(void);
//  void unalloc_dims_nfleets(void);
//  void unalloc_dims_nmetiers(int);

//  void alloc_dims_fleet(   int ifleet=1);
//  void alloc_dims_nmetiers(int ifleet=1);
//  void alloc_dims_stock(     int ifleet=1,  int iMetier=1, int nstock=1);

//  void unalloc_dims_fleet(void);
//  void unalloc_dims_metiers(int ifleet=1);
//  void unalloc_dims_stock(   int ifleet=1, int iMetier=1);
  };        

#endif /* _INC_flc */



