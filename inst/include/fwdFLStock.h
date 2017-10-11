#ifndef _INC_fwdFLStock
#define _INC_fwdFLStock

#include "flc_adolc.h"
#include "flc.h"

double norm(double *, int);

class fwdStk 
{
public:     
	fwdStk(void);
	~fwdStk(void);

  SEXP Init(SEXP xStk, SEXP xYrs, SEXP xSRModel,SEXP xSRParam,SEXP xSRResiduals,SEXP xMult,SEXP xAvail,SEXP xMaxFBar);    

	void InitAvail(  SEXP x); 

	double getVal(FLRConst_Target quantity,  int iyr, int iunit, int iseason, int iarea, int iter);

	void project(double  *x,                 int iyr, int iunit, int iseason, int iarea, int iter, bool OnlyReplaceNA=FALSE, bool OnlyCalcN=FALSE);
    void project(adouble *x, adouble *func, double *Trgt, int iTrgt, int nrow, double *Ary, int iter);

	adouble computeStock(   FLQuant_adolc &n, FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter);      
	adouble SSB(            FLQuant_adolc &n, FLQuant_adolc &f, int iyr, int iunit, int iseason, int iarea, int iter);
	adouble computeCatch(   FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble computeDiscards(FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble computeLandings(FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble Zbar(           FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble Fbar(           FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble FbarLandings(   FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble FbarDiscards(   FLQuant_adolc &f,                   int iyr, int iunit, int iseason, int iarea, int iter);
	adouble MnSz(           FLQuant_adolc &n,                   int iyr, int iunit, int iseason, int iarea, int iter);

	double SSB(                                                 int iyr, int iunit, int iseason, int iarea, int iter);
	double Fbar(                                                int iyr, int iunit, int iseason, int iarea, int iter);

  SEXP run(SEXP xTrgt, SEXP xAry) ;   
protected:        
	FLStock stk;
	FLQuant avail;
  FLQuant MaxF;

	bool indepLastYr;

	sr SR;
    };

#endif /* _INC_fwdFLStock */



