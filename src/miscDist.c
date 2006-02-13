/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2002	S. M. Iacus
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Exports
 *	addDist(...)
 *  setDist(...)
 *
 * to be called as  .Call(.)  in ../R/miscDist.R
 */



#include <R.h>
#include <Rmath.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Complex.h>


/* IFS estimators */      
SEXP addDist(SEXP d, SEXP x, SEXP k, SEXP n);
SEXP setDist(SEXP d, SEXP x, SEXP k, SEXP n);


/* addDist:
   ------
   
   parameters:
   -----------
   
   d   : the pointer to the vector containing elements of the dist object
   x   : the vector of indexes to change in the original matrix
   k   : the constant to add vector-wise
   n   : the dimension of the dist object
*/



SEXP addDist(SEXP d, SEXP x, SEXP k, SEXP n)
{
	int i, j, pos;
	double *dObj, *kappa;
	int dSize, xSize;
	int *nSize, *xIdx;
	
	if(!isNumeric(d)) error("`d' must be numeric");
	if(!isInteger(x)) error("`x' must be an integer");
	if(!isNumeric(k)) error("`k' must be numeric");
	if(!isInteger(n)) error("`n' must be an integer");
	
	PROTECT(d = AS_NUMERIC(d));
	PROTECT(x = AS_INTEGER(x));
	PROTECT(k = AS_NUMERIC(k));
	PROTECT(n = AS_INTEGER(n));
	
	dSize = LENGTH(d); 
	xSize = LENGTH(x);
	
	xIdx = INTEGER_POINTER(x);
	dObj = NUMERIC_POINTER(d);
	nSize = INTEGER_POINTER(n);
	kappa = NUMERIC_POINTER(k);
	
	for(i=0; i< xSize; i++)
		if(xIdx[i] > *nSize)
			error("subscript out of bounds in `x'");
	
	for(i=0; i< xSize-1; i++)
		for(j=i+1; j<xSize; j++){
			pos = (xIdx[i]-1)*(2*(*nSize)-xIdx[i])/2 + xIdx[j]-xIdx[i];
			REAL(d)[pos-1] = dObj[pos-1] + *kappa;
		}
			UNPROTECT(4);
	return(d);
}

SEXP setDist(SEXP d, SEXP x, SEXP k, SEXP n)
{
	int i, j, pos;
	double *dObj, *kappa;
	int dSize, xSize;
	int *nSize, *xIdx;
	
	if(!isNumeric(d)) error("`d' must be numeric");
	if(!isInteger(x)) error("`x' must be an integer");
	if(!isNumeric(k)) error("`k' must be numeric");
	if(!isInteger(n)) error("`n' must be an integer");
	
	PROTECT(d = AS_NUMERIC(d));
	PROTECT(x = AS_INTEGER(x));
	PROTECT(k = AS_NUMERIC(k));
	PROTECT(n = AS_INTEGER(n));
	
	dSize = LENGTH(d); 
	xSize = LENGTH(x);

	xIdx = INTEGER_POINTER(x);
	nSize = INTEGER_POINTER(n);
	kappa = NUMERIC_POINTER(k);
	
	for(i=0; i< xSize; i++){
		if(xIdx[i] > *nSize)
			error("subscript out of bounds in `x'");
	}	

	for(i=0; i< xSize-1; i++)
		for(j=i+1; j<xSize; j++){
			pos = (xIdx[i]-1)*(2*(*nSize)-xIdx[i])/2 + xIdx[j]-xIdx[i];
			REAL(d)[pos-1] = *kappa;
		  }
	UNPROTECT(4);
	return(d);
}


static R_CMethodDef R_CDef[] = {
   {"addDist", (DL_FUNC)&addDist, 4},
   {"setDist", (DL_FUNC)&setDist, 4},
   {NULL, NULL, 0},
};

void
R_init_rrp(DllInfo *info)
{
    R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
}


