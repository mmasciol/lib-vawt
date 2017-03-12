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


#ifndef _AC_H
#define _AC_H


#include <stdio.h>
#include <stdbool.h>
#include "error.h"
#include "./simclist/simclist.h"


struct Domain;
struct Airfoil;
struct BladeElement;


struct ACWork {
    double* w;         /* what we are solving for, w = [wx, wy]^T */
    double* wx;        /* of length netheta */
    double* wy;        /* of length netheta */
    double* Vt;        /* of length netheta */
    double* Vn;        /* of length netheta */
    double* cn;        /* of length netheta */
    double* ct;        /* of length netheta */
    double* Pn;        /* of length netheta */
    double* integrand; /* of length netheta */
    double* cl;        /* these are used to temporarily store interpolated cl and cd coefficients */
    double* cd;
    double* residual;
    double* W;     /* of length netheta */
    double* phi;   /* of length netheta */
    double* alpha; /* of length netheta */
    double* Re;    /* of length netheta */    
};


struct ACPerformance {
    double* fx;    /* of length netheta */
    double* fy;    /* of length netheta */
    double* fz;    /* of length netheta */
    double* Q;     /* torque contribution for all blades */
    double Cp;
    double Ct;
    double Cx;
    double Cy;
    double Cz;
    int n;
};


struct ActuatorCylinder {    
    struct PowellCoeff* pc;
    struct ACWork* work;
    struct ACPerformance* perf;
    struct Blade* blade;
    struct Environment* env;
    double* A;
    double* theta;    
    double* z;
    double* fr; /* nxm matrix, used by the solver to store -cn*q*chord */
    double* ft; /* nxm matrix, used by the solver to store ct*q*chord/cos(delta) */
    double* fz; /* nxm matrix, used by the solver to store -cn*q*(chord)*tan(delta) */
    double omega;
    double the0;
    double swept_area;
    int nt;     /* number of theta discretization */
    int nz;     /* length of vertcal discretization z */
    int nb;     /* number of blades */
    int index;  /* current index being solved */
    bool dstall;
    char* passed_msg;
    ERROR_CODE* passed_ierr;
};


size_t ac_size_meter(const void *el);
double ac_periodic_integrand(const int n, const double* fx, const double* theta, char* msg, ERROR_CODE* ierr);
int ac_func(void* ac_ptr, int n, const double* w, double* fvec, int iflag);
ERROR_CODE ac_velocity(struct ActuatorCylinder* ac, struct BladeElement* elem, const int index, const double Uinf, const double* w, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_error(struct ActuatorCylinder* ac, struct BladeElement* elem, const double* w, double* fvec, const int index, const double Uinf, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_initialize(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_perf_initialize(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_free(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_work_free(struct ACWork* work, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_perf_free(struct ACPerformance* perf, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_generate_matrix(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_generate_memory(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_allocate_perf(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_force(struct ActuatorCylinder* ac, /*const double u, const double omega,*/ char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_inflow_wind_velocity(struct ActuatorCylinder* ac, const double u, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_angular_velocity(struct ActuatorCylinder* ac, const double omega, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_z(struct ActuatorCylinder* ac, double* z, const int nz, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_airfoil_chord(struct ActuatorCylinder* ac, const double chord, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_n_theta(struct ActuatorCylinder* ac, const int nt, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_airfoil_radius(struct ActuatorCylinder* ac, const double r, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_airfoil_twist(struct ActuatorCylinder* ac, const double twist, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_airfoil_swept_offset(struct ActuatorCylinder* ac, const double twist, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
// ERROR_CODE ac_set_solidity_ratio(struct ActuatorCylinder* ac, const double solidity_ratio, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_airfoil_geometry(struct ActuatorCylinder* ac, struct Airfoil* af, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_power(struct ActuatorCylinder* ac, struct Blade* blade, /*const double theta, */ char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_initial_theta_offset(struct ActuatorCylinder* ac, const double the0, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_set_swept_area(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_allocate_A_matrix(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_allocate_theta(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_allocate_prerequisites_for_A_matrix(const int h, const int n, double** value, double*** Rx, double*** Ry, double*** Wx, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_delete_prerequisites_for_A_matrix(const int n, double** value, double*** Rx, double*** Ry, double*** Wx, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_wind_shear(struct ActuatorCylinder* ac, const double v_ref, const double ref_height, const double shear, double** V, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_evaluate(struct ActuatorCylinder* ac, const int index, const double V,  double* w, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_initialize_list(struct Domain* domain, char* msg, ERROR_CODE* ierr);
ERROR_CODE ac_new_turbine(list_t* restrict ac, char* msg, ERROR_CODE* ierr);


#endif /* _AC_H */
