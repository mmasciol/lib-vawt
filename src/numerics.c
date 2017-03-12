/****************************************************************
 *   Copyright (C) 2016 mdm                                     *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#include <stdlib.h>
#include "numerics.h"
#include "sys.h"
#include "math.h"
#include "ac.h"


static inline int 
ij(const int row, const int col, const int rows)
{
    return row*rows + col;
}


ERROR_CODE 
powell_initialize(struct PowellCoeff* pc, char* msg, ERROR_CODE* ierr)
{
    pc->j = -9999;
    pc->n = -9999;
    pc->maxfev = -9999;
    pc->ml = -9999;
    pc->mu = -9999;
    pc->mode = -9999;
    pc->nprint = -9999;
    pc->info = -9999;
    pc->nfev = -9999;
    pc->ldfjac = -9999;
    pc->lr = -9999;
    pc->xtol = -999.9;
    pc->epsfcn = -999.9;
    pc->factor = -999.9;
    pc->diag = NULL;
    pc->r = NULL;
    pc->qtf = NULL;
    pc->fjac = NULL;
    pc->wa1 = NULL;
    pc->wa2 = NULL;
    pc->wa3 = NULL;
    pc->wa4 = NULL;
    return SAFE;
}

ERROR_CODE
powell_allocate(struct ActuatorCylinder* ac, const int nt, char* msg, ERROR_CODE* ierr)
{
    int j = 0;
    struct PowellCoeff* pc = NULL;
    ERROR_CODE success = SAFE;
    
    ac->pc = malloc(sizeof(struct PowellCoeff)); CHECK_MALLOC(ac->pc, FATAL_8, "<PowellCoeff>");            
    pc = ac->pc;
    success = powell_initialize(pc, msg, ierr); CHECKERRK(FATAL_5, "Failed to initialize Powell Hybrid parameters");

    pc->n = 2*nt;
    pc->ldfjac = pc->n;
    pc->lr = ceil((pc->n)*((pc->n)+1)/2);
    pc->xtol = 1E-3;//sqrt(__cminpack_func__(dpmpar)(1));            
    pc->maxfev = 200;
    pc->ml = (pc->n)-1;
    pc->mu = (pc->n)-1;
    pc->epsfcn = 0.;
    pc->mode = 1;
    pc->factor = 1;
    pc->nprint = 0;
    
    pc->diag = malloc(pc->n*sizeof(double*)); CHECK_MALLOC(pc->diag, FATAL_8, "diag");            
    pc->r = malloc(pc->lr*sizeof(double*)); CHECK_MALLOC(pc->r, FATAL_8, "r");
    pc->qtf = malloc(pc->n*sizeof(double*)); CHECK_MALLOC(pc->qtf, FATAL_8, "qtf");
    pc->fjac = malloc((pc->n)*(pc->n)*sizeof(double*)); CHECK_MALLOC(pc->fjac, FATAL_8, "fjac");
    pc->wa1 = malloc(pc->n*sizeof(double*)); CHECK_MALLOC(pc->wa1, FATAL_8, "wa1");
    pc->wa2 = malloc(pc->n*sizeof(double*)); CHECK_MALLOC(pc->wa2, FATAL_8, "wa2");
    pc->wa3 = malloc(pc->n*sizeof(double*)); CHECK_MALLOC(pc->wa3, FATAL_8, "wa3");
    pc->wa4 = malloc(pc->n*sizeof(double*)); CHECK_MALLOC(pc->wa4, FATAL_8, "wa4");
    for (j=1; j<=pc->n; j++) {
        pc->diag[j-1] = 1.;
    }    
    return SAFE;
CLEAN_UP:
    return FATAL;
}

ERROR_CODE
powell_free(struct PowellCoeff* pc, char* msg, ERROR_CODE* ierr)
{
    FREE_OBJ(pc->r);
    FREE_OBJ(pc->diag);
    FREE_OBJ(pc->qtf);
    FREE_OBJ(pc->fjac);
    FREE_OBJ(pc->wa1);
    FREE_OBJ(pc->wa2);
    FREE_OBJ(pc->wa3);
    FREE_OBJ(pc->wa4);
    return SAFE;
}


double 
min(double* a, const int n) 
{
    int i = 0;
    double min = a[0];

    for (i=0 ; i<n ; i++) {
        if (a[i]<min) {
            min = a[i];
        }
    }
    return min;
}


double 
max(double* a, const int n) 
{
    int i = 0;
    double max = a[0];

    for (i=0 ; i<n ; i++) {
        if (a[i]>max) {
            max = a[i];
        }
    }
    return max;
}



double*
linear_interpolation(double* alpha, const int n, double* in_alpha, double* X, const int m, char* msg, ERROR_CODE* ierr)
{
    double* arr = NULL;
    double x0 = 0.0;
    double x1 = 0.0;
    double y0 = 0.0;
    double y1 = 0.0;
    int i = 0;
    int j = 0;
    arr = malloc(n*sizeof(double));
    if (!arr) {
        CHECKERRQ(FATAL_8, "Linear Interpolation");
        return NULL;
    }

    for (i=0 ; i<n ; i++) { /* interates on alpha */
        for (j=0 ; j<m-1 ; j++) { /* interates on in_alpha/X */            
            if (in_alpha[j]<=alpha[i] && alpha[i]<=in_alpha[j+1]) {
                x0 = in_alpha[j];
                x1 = in_alpha[j+1];
                y0 = X[j];
                y1 = X[j+1];
                arr[i] = y1 + (alpha[i]-x1)*((y0-y1)/(x0-x1));
            }            
        }
    }
    return arr;
} 



ERROR_CODE
linear_bivariate_spline(double* fx, const double x, const double y, 
                        const double* alpha, const double* A, 
                        const int m, const int n, int* index_x, int* index_y, char* msg, ERROR_CODE* ierr)
{
    /* should be m x n matrices. m=len(alpha), n=len(Re tables). 
     * x ranges from 1 ... m and is alpha
     * y ranges from 1 ... n and is Re
     *
     */
    ERROR_CODE success = SAFE;
    int i = *index_x;
    int its = 0;
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    /* find x in alpha and start at index xi */
    /* @todo: this function does no bound checking. Throw and error if outside the array bounds */
    do {
        if (alpha[i-1]<=x && x<=alpha[i]) {
            *index_x = i;
            break;
        } 
        if (alpha[i]<=x) {
            i++;
        } else {
            i--;
        }
        its++;
        if (1000<its) {            
            CALL_CHECKERRK(FATAL_25, "Maximum number of iteration exceeded in linear bivariate spline interpolation.");
        }
    } while(1);
    x1 = alpha[i-1];
    x2 = alpha[i];
    y1 = A[i-1];
    y2 = A[i];
    *fx = y2 - (x2-x)*((y2-y1)/(x2-x1));   
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
linear_2D_spline(double* fx,
                 double* A,
                 int* xstart,  
                 const double x, 
                 const double* x_arr,
                 const int xlen,
                 int* ystart, 
                 const double y,                              
                 const double* y_arr,
                 const int ylen, 
                 char* msg, 
                 ERROR_CODE* ierr)
{
    /*
      fx = f(x,y)
      dimension(A) = m, n
      -- index_x and index_y are starting indices. We pass them back to the caller to give a good restart position
      -- 
    */
    ERROR_CODE success = SAFE;
    int its = 0;
    int i = *xstart;
    int j = *ystart;
    int i1 = 0;
    int i2 = 0;
    int j1 = 0;
    int j2 = 0;
    double x1 = 0.0;
    double x2 = 0.0;
    double y1 = 0.0;
    double y2 = 0.0;
    double Q11 = 0.0;
    double Q12 = 0.0;
    double Q21 = 0.0;
    double Q22 = 0.0;

    if (xlen<=i) {
        i = xlen-1;
    }

    if (ylen<=j) {
        j = ylen-1;
    }

    /* find x in x_arr and start at index istart */
    do {
        if (i==0) {
            i+=1;
            *xstart = i;
            i1 = i-1;
            i2 = i;
            break;
        } else if (x_arr[i-1]<=x && x<=x_arr[i]) {
            *xstart = i; 
            i1 = i-1;
            i2 = i;
            break;
        } 
        if (x_arr[i]<=x) {
            i++;
            its++;
        } else {
            i--;
            its++;
        }
        if (500<its) {
            CALL_CHECKERRK(FATAL_37,"Exceeded iteration count on x array <500>.");
        } else if (xlen<=i) {
            CALL_CHECKERRK(FATAL_37,"Exceeded array length on x array.");
        }
    } while(1);

    
    /* find y in y_arr and start at index jstart */
    its = 0;
    if (y_arr[ylen-1]<y) { /* this is here to complere a unit circle in case y_arr[len] < 2*PI */
        j1 = 0;
        j1 = ylen-1;
    } else {
        do {
            if (j==0) {
                j+=1;
                *ystart = j;
                j1 = j-1;
                j2 = j;
                break;
            } else if (y_arr[j-1]<=y && y<=y_arr[j]) {
                *ystart = j;
                j1 = j-1;
                j2 = j;
                break;
            } 
            if (y_arr[j]<=y) {
                j++;
                its++;
            } else {
                j--;
                its++;
            }
            if (500<its) {
                CALL_CHECKERRK(FATAL_37,"Exceeded iteration count on y array <500>.");
            } else if (ylen<=j) {
                CALL_CHECKERRK(FATAL_37,"Exceeded array length on y array.");
            }
        } while(1);
    }
    x1 = x_arr[i1];
    x2 = x_arr[i2];
    y1 = y_arr[j1];
    y2 = y_arr[j2];
    Q11 = A[ij(i1,j1,ylen)]; // (x1, y1);
    Q12 = A[ij(i1,j2,ylen)];   // (x1, y2);
    Q21 = A[ij(i2,j1,ylen)];   // (x2, y1);
    Q22 = A[ij(i2,j2,ylen)];     // (x2, y2);
    *fx = ((y2-y)/(y2-y1))*(((x2-x)/(x2-x1))*Q11 + ((x-x1)/(x2-x1))*Q21) 
        + ((y-y1)/(y2-y1))*(((x2-x)/(x2-x1))*Q12 + ((x-x1)/(x2-x1))*Q22);
    return SAFE;
CLEAN_UP:
    return FATAL;
}


int
get_index(const double val, double* arr, const int n)
{
    int i = 0;
    for (i=0 ; i<n ; i++) {
        if (arr[i]<=val && val<=arr[i+1]) {
            return i;
        }
    }
    return -1; /* error, check if i is negative upon return */
}


ERROR_CODE
integrate(double* sum, double* y, double* x, const int n, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    double add = 0.0;
    *sum = 0.0;

    for (i=0 ; i<n-1 ; i++) {
        add += ((x[i+1]-x[i])*(y[i+1]+y[i]));
    }
    *sum = 0.5*(add);
    
    return SAFE;
}
