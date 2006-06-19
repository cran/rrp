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


/* dist object manipulation */      
SEXP addDist(SEXP d, SEXP x, SEXP k, SEXP n, SEXP rho);
SEXP setDist(SEXP d, SEXP x, SEXP k, SEXP n, SEXP rho);


/* addDist/setDist:
   ------
   
   parameters:
   -----------
   
   d   : the pointer to the vector containing elements of the dist object
   x   : a list of indexes to change in the original matrix
   k   : a vector of constants to add vector-wise
   n   : the dimension of the dist object
*/



SEXP addDist(SEXP d, SEXP x, SEXP k, SEXP n, SEXP rho)
{
	int i, j, h, pos;
	double *dObj, *kappa;
	int dSize, xSize;
	int *nSize, *xIdx;
	int n2,n1;
	
	if(!isNumeric(d)) error("`d' must be numeric");
	if(!isNumeric(k)) error("`k' must be numeric");
	if(!isInteger(n)) error("`n' must be an integer");
	
	PROTECT(d = AS_NUMERIC(d));
	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	PROTECT(n = AS_INTEGER(n));
	
	dSize = LENGTH(d); 
    n2 = length(eval(x, rho));
    if (n2 == NA_INTEGER)
     error("error");
	
	dObj = NUMERIC_POINTER(d);
	nSize = INTEGER_POINTER(n);
	kappa = NUMERIC_POINTER(k);

	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   pos = (xIdx[i]-1)*(2*(*nSize)-xIdx[i])/2 + xIdx[j]-xIdx[i];
	   REAL(d)[pos-1] = dObj[pos-1] + kappa[h];
	  }
	 }
	}
	UNPROTECT(4);
    return(d);
}

SEXP setDist(SEXP d, SEXP x, SEXP k, SEXP n, SEXP rho)
{
	int h, i, j, pos;
	double  *kappa;
	int dSize, xSize;
	int *nSize, *xIdx, n2, n1;
	
	if(!isNumeric(d)) error("`d' must be numeric");
	if(!isNumeric(k)) error("`k' must be numeric");
	if(!isInteger(n)) error("`n' must be an integer");
	
	PROTECT(d = AS_NUMERIC(d));
	PROTECT(x = AS_LIST(x));
	PROTECT(k = AS_NUMERIC(k));
	PROTECT(n = AS_INTEGER(n));
	
	dSize = LENGTH(d); 
    n2 = length(eval(x, rho));
    if (n2 == NA_INTEGER)

	xIdx = INTEGER_POINTER(x);
	nSize = INTEGER_POINTER(n);
	kappa = NUMERIC_POINTER(k);
			  
	for(h=0; h<n2; h++){
     n1 = LENGTH(VECTOR_ELT(x,h));
	 xIdx = INTEGER_POINTER(AS_INTEGER(VECTOR_ELT(x,h)));
	 for(i=0; i<n1-1; i++){
	  for(j=i+1; j<n1; j++){
	   pos = (xIdx[i]-1)*(2*(*nSize)-xIdx[i])/2 + xIdx[j]-xIdx[i];
	   REAL(d)[pos-1] = kappa[h];
	  }
	 }
	}

	UNPROTECT(4);
	return(d);
}


static R_CMethodDef R_CDef[] = {
   {"addDist", (DL_FUNC)&addDist, 5},
   {"setDist", (DL_FUNC)&setDist, 5},
   {NULL, NULL, 0},
};

void
R_init_rrp(DllInfo *info)
{
    R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
}


