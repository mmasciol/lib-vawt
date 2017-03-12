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


#ifndef _NUMERICS_H
#define _NUMERICS_H


#include "error.h"


struct ActuatorCylinder;


struct PowellCoeff {
    int j;
    int n;
    int maxfev;
    int ml;
    int mu;
    int mode;
    int nprint;
    int info;
    int nfev;
    int ldfjac;
    int lr;
    double xtol;
    double epsfcn;
    double factor;
    double* diag;
    double* r;
    double* qtf;
    double* fjac;
    double* wa1;
    double* wa2;
    double* wa3;
    double* wa4;
};

double min(double* a, const int n);
double max(double* a, const int n);
double* linear_interpolation(double* alpha, const int n, double* in_alpha, double* X, const int m, char* msg, ERROR_CODE* ierr);
ERROR_CODE linear_bivariate_spline(double* fx, const double x, const double y, const double* alpha, const double* A, const int m, const int n, int* index_x, int* index_y, char* msg, ERROR_CODE* ierr);
ERROR_CODE linear_2D_spline(double* fx,
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
                            ERROR_CODE* ierr);

ERROR_CODE integrate(double* sum, double* y, double* x, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE powell_allocate(struct ActuatorCylinder* ac, const int nt, char* msg, ERROR_CODE* ierr);
ERROR_CODE powell_free(struct PowellCoeff* pc, char* msg, ERROR_CODE* ierr);
int get_index(const double val, double* arr, const int n);
#endif /* _NUMERICS_H */
