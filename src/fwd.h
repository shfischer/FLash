#ifndef _INC_fwd
#define _INC_fwd

#include "flc_adolc.h"
#include "flc.h"

//#define MaxFBar 2.0

double norm(double *, int);

//trgtNms    <-function() return(c("year","min","val","max","quantity","season","area","unit","spp","fleet","metier","rel.year","rel.season","rel.area","rel.unit"))
//effNms     <-function() return(c("year","min","val","max","fleet","metier","rel.year","rel.fleet","rel.metier","rel.bound"))
//quantityNms<-function() return(c("ssb","biomass","catch","landings","discards","f","z","f.landings","f.discards","effort","costs","revenue","profit","mnsz"))

typedef enum tag_fwdControlPos
	{
    fwdControlPos_year      = 0, 
    fwdControlPos_min       = 1,
    fwdControlPos_val       = 2,
    fwdControlPos_max       = 3, 
    fwdControlPos_season    = 4,
    fwdControlPos_area      = 5,
    fwdControlPos_fleet     = 6,
    fwdControlPos_metier    = 7,
    fwdControlPos_relyear   = 8,
    fwdControlPos_relseason = 9,
    fwdControlPos_relarea   = 10,
    fwdControlPos_relfleet  = 11, 
    fwdControlPos_relmetier = 12,
    fwdControlPos_relbound  = 13
	} fwdControlPos;

typedef enum tag_fwdTargetPos
	{
    fwdTargetPos_year     = 0,
    fwdTargetPos_min      = 1,
    fwdTargetPos_val      = 2,
    fwdTargetPos_max      = 3,
    fwdTargetPos_quantity = 4,
    fwdTargetPos_season   = 5,
    fwdTargetPos_area     = 6,
    fwdTargetPos_unit     = 7,
    fwdTargetPos_spp      = 8,
    fwdTargetPos_fleet    = 9,
    fwdTargetPos_metier   = 10,
    fwdTargetPos_relyear  = 11,
    fwdTargetPos_relseason= 12,
    fwdTargetPos_relarea  = 13,
    fwdTargetPos_relunit  = 14
	} fwdTargetPos;

class control
{
public:        
   control(void);
   ~control(void);

   bool AllocFlag;

   bool Init(SEXP xCtrl, SEXP xAry, int niters);

   double value(int year, int fleet, int metier);
   double min(  int year, int fleet, int metier);
   double max(  int year, int fleet, int metier);
   
   int    n(int year);

   int    rel_year(  int year, int fleet, int metier);
   int    rel_fleet( int year, int fleet, int metier);
   int    rel_metier(int year, int fleet, int metier);
   bool   rel_bound( int year, int fleet, int metier);

   bool   fix(int year, int fleet, int metier);
   bool   fit(int year, int fleet, int metier);
private:
   double ****data;

   int minyear, maxyear,
       nfleet, nmetier;
   
   bool unalloc(void);
};

//year value min max quantity rel spp fleet metier
class target
{
public:        
   target(void);
   ~target(void);

   bool AllocDataFlag, AllocAryFlag;

   bool Init(SEXP xTrgt, SEXP xAry, int niters);
   
   int spp(   int year, int i);
   int fleet( int year, int i);
   int metier(int year, int i);

   double val(int year, int i, int iter=0);
   double min(int year, int i, int iter=0);
   double max(int year, int i, int iter=0);
   
   inline double n(int year) {return _n[year];};
   inline double niters(void) {return _niters;};
   
   FLRConst_Target quantity(int year, int i);
   
   double rel(int year, int i);
private:
   double ***data,
          ***_min,
          ***_max,
          ***_val;

   int *_n,
        _minyear, _maxyear,
        _niters,
        _nrows;

   void unalloc(void);
};

class fwd : public flc
{
public:        
      fwd(SEXP xStk,                 SEXP xYrs,             SEXP xSRModel, SEXP xSRParam, SEXP xSRResiduals, SEXP xMult);
      fwd(SEXP xBiols, SEXP xFleets, SEXP xYrs, SEXP xDims, SEXP xSRR);
      fwd(SEXP xBiols, SEXP xFleets, SEXP xYrs, SEXP xDims, SEXP xSRModel, SEXP xSRParam, SEXP xSRResiduals, SEXP xMult);
      
      bool run(SEXP xTrgt);
      bool run(SEXP xTrgt, SEXP xCtrl, SEXP xAryTrgt, SEXP xAryCtrl, SEXP xYrs);

      adouble computeStock(   FLQuant2_adolc &ad_n, FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble SSB(            FLQuant2_adolc &ad_n, FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter); 
      adouble computeCatch(                         FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble computeDiscards(                      FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble computeLandings(                      FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble Fbar(                                 FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble Zbar(                                 FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble FbarLandings(                         FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      adouble FbarDiscards(                         FLQuant2_adolc &ad_f, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);      
      
      control Ctrl;
      target  Trgt;

      double getVal(FLRConst_Target quantity, int ispp, int iyr, int iunit, int iseason, int iarea, int iter);

      void    project(adouble *x, adouble *func, int iYr, int iter);
      void    project(double *x, int iYr, int iter, bool OnlyReplaceNA=FALSE, bool OnlyCalcN=FALSE);

      int NBiol(void);
      int MinProjYear, MaxProjYear;
      };

#endif /* _INC_fwd */



