/****************************************************************
 * Original work by A. Ning (c) 2016                            *
 * Modified by mdm 2016                                         *
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


#ifndef _AF_H
#define _AF_H

#include <stdio.h>
#include "error.h"
#include "./simclist/simclist.h"


struct bstrList;
struct ActuatorCylinder;
struct Domain;


// struct Polar {
//     double* alpha;
//     double* cl;
//     double* cd;
//     double* cm;
//     int len;
//     char type[32];
// };

struct Airfoil {
    double Re;   /* Reynold's number, this is the set from the aerodynamic input file */
    double aol;  /*< AOL - Zero Cn angle of attack (deg) */
    double cna;  /*< CNA - Cn slope for zero lift (dimensionless) */
    double cns;  /*< CNS - Cn extrapolated to value at positive stall angle of attack */
    double cnsl; /*< CNSL - Cn at stall value for negative angle of attack */
    double aod;  /*< AOD - Angle of attack for minimum CD (deg) */
    double cd0;  /*< CDO - Minimum CD value */
    double* alpha;
    double* cl;
    double* cd;
    double* cm;
    int na; /* @todo: alpha, cl, cd, and cm should be m x n matrices. m=len(alpha), n=len(Re tables). 
              *        See: createDataGrid in AirfoilPreppy.py.
              *        Note: a holder function af_set_grid is defined in af.c; to be populated later
              */
    char type[32];
    void* vtk_point;
};


struct ExtrapolatedFoil {
    int N;
    double* alpha;
    double* cl;
    double* cd;    
    double* alpha1;
    double* alpha2;
    double* alpha3;
    double* alpha4;
    double* alpha5;
    double* alpha6;
    double* alpha7;
    double* cd1;
    double* cd2;
    double* cd3;
    double* cd4;
    double* cd5;
    double* cd6;
    double* cd7;
    double* cl1;
    double* cl2;
    double* cl3;
    double* cl4;
    double* cl5;
    double* cl6;
    double* cl7;
    double cl_adj;
    double alpha_high;
    double alpha5_max;
    double cl_high;
    double cd_high;
    double cm_high;
    double alpha_low;
    double cl_low;
    double cd_low;
    double cdmax;
    double A;
    double B;
};

ERROR_CODE af_read_file(list_t* restrict af, char file_name[255], char type_name[32], char* msg, ERROR_CODE* ierr);
ERROR_CODE af_free(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_blend(list_t* restrict af, struct Airfoil* af1, struct Airfoil* af2, const double weight, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_initialize(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate(struct Airfoil* af, const double cd_max, const double aspect_ratio, const double cd_min, const int n_alpha, char* msg, ERROR_CODE* ierr);

void af_increase_array(double** data, const int count);
double* af_merge_aoa(double* arr1, const int n1, double* arr2, const int n2, int* len, char* msg, ERROR_CODE* ierr);

// ERROR_CODE af_get_coefficients(struct Airfoil* af, struct ActuatorCylinder* ac, const int n, const double* alpha, const double* Re, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_cl(double* value, struct Airfoil* af, const double alpha, const double Re, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_cd(double* value, struct Airfoil* af, const double alpha, const double Re, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_new_foil(list_t* restrict af, char* msg, ERROR_CODE* ierr);

ERROR_CODE af_initialize_extrapolated_foil(struct ExtrapolatedFoil* nf, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_free_extrapolated_foil(struct ExtrapolatedFoil* nf, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_viterna_extrapolation(double* alpha, const int n, const double cdmax, double* cd, double* cl, const double A, const double B, const double cl_adj, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_1(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_2(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_3(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_4(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_5(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_6(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_foil_7(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_extrapolate_pack(struct Airfoil* af, struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_Re(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_aol(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cna(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cns(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cnsl(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_aod(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cd0(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_alpha(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cl(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cd(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_set_cm(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr);


ERROR_CODE af_initialize_list(struct Domain* domain, char* msg, ERROR_CODE* ierr);
size_t af_size_meter(const void *el);


ERROR_CODE af_set_beddoes_coefficients(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_cna(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_cns(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_cnsl(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_cd0(struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE af_get_aol(struct Airfoil* af, char* msg, ERROR_CODE* ierr);



#endif /* _AF_H */
