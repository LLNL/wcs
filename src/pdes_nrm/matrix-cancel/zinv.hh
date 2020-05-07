/******************************************************************************
 *                                                                            *
 *    Copyright 2020   Lawrence Livermore National Security, LLC and other    *
 *    Whole Cell Simulator Project Developers. See the top-level COPYRIGHT    *
 *    file for details.                                                       *
 *                                                                            *
 *    SPDX-License-Identifier: MIT                                            *
 *                                                                            *
 ******************************************************************************/

#ifndef ZINV__
#define ZINV__

#include <cassert>

template<typename myuint>
myuint intinv(myuint b) {
  myuint t0 = 0,t = 1,q,r;
  myuint a = 0;

  if(b <= 1) return b;

  q = (~a) / b;
  if(b*q+b == 0) return 0;
  r = a - q*b;
  while(r > 0) {
    const myuint temp = t0 - q*t;
    t0 = t;
    t = temp;
    a = b;
    b = r;
    q = a/b;
    r = a - q*b;
  }
  if(b == 1) return t;
  else return 0;
}

template<typename myuint>
int matinv(int n,myuint X[],myuint Y[]) {
  int p[n],q[n];

  assert(((myuint) -1) > 0);

  for(int i = 0; i<n; i++)
    p[i] = q[i] = i;
    
  for(int i = 0; i<n; i++) {
    for(int j = 0; j<n; j++)
      Y[i*n+j] = 0;
    Y[i*n+i] = 1;
  }

  for(int i = 0; i<n; i++) {
    myuint xi = 0;
    /* Find pivot element */ {
      int j,k;
      for(j = i; j<n; j++)
	for(k = i; k<n; k++) {
	  xi = intinv(X[p[j]*n+q[k]]);
	  if(xi != 0) {
	    int t = p[i];
	    p[i] = p[j];
	    p[j] = t;
	    t = q[i];
	    q[i] = q[k];
	    q[k] = t;
	    goto found;
	  }
	}
    found:
      if(xi == 0) return 0;
    }
    
    /* Make diagonal element 1 */ {
      for(int j = i; j<n; j++)
	X[p[i]*n+q[j]] *= xi;
      for(int j = 0; j<n; j++)
	Y[p[i]*n+q[j]] *= xi;
    }

    /* Elimination */ {
      for(int j = i+1; j<n; j++) {
	myuint xji = X[p[j]*n+q[i]];
	for(int k = i; k<n; k++)
	  X[p[j]*n+q[k]] -= xji*X[p[i]*n+q[k]];
	for(int k = 0; k<n; k++)
	  Y[p[j]*n+k] -= xji*Y[p[i]*n+k];
      }
    }
  }

  /* Back substitution */ {
    for(int i = n-1; i>=0; i--) {
      for(int j = 0; j<i; j++) {
	myuint xji = X[p[j]*n+q[i]];
	for(int k = 0; k<n; k++)
	  Y[p[j]*n+k] -= xji*Y[p[i]*n+k];
	X[p[j]*n+q[i]] -= xji*X[p[i]*n+q[i]];
      }
    }
  }

  /* Undo permutation, and copy result to Y */ {
    for(int i = 0; i<n; i++)
      for(int j = 0; j<n; j++)
	X[q[i]*n+j] = Y[p[i]*n+j];
    for(int i = 0; i<n*n; i++)
      Y[i] = X[i];
  }
  return 1;
}

template<typename myuint>
void matmul(int n,myuint A[],myuint B[],myuint AB[]) {
  for(int i = 0; i<n; i++)
    for(int j = 0; j<n; j++) {
      myuint s = 0;
      for(int k = 0; k<n; k++)
	s = s + A[i*n+k]*B[k*n+j];
      AB[i*n+j] = s;
    }
}

#endif
