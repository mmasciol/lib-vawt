/****************************************************************
 * Original work by NREL (c) 2013                               *
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


#ifndef _DSTALL_H
#define _DSTALL_H


#include <stdio.h>
#include <stdbool.h>
#include "error.h"


struct ActuatorCylinder;
struct Blade;
struct Airfoil;
struct BladeElement;


struct DStall {
    bool first_entry;
    bool flow_shift;
    bool vortex;
    bool separation;
    bool super_sonic;
    int na;      /* length of cl/cd array in airfoil structure */
    double CLA;
    double CDA;
    double CMA;

    double a1;
    double a2;
    double b1;
    double b2;

    double alpha_e;
    double dt;
    double Cn_cp;    
    double alpha_n;
    double ds;
    double fk;
    double Cni_q;    
    double fn_pp; 
    double ft_pp;
    double tv;   /* vortex shedding time constant */
    double tf;   /* separation time constant */
    double tv_l; /* Non-dimensional vortex transit time along airfoil chord */
    double c;    /* Speed of sound */
    double tp;   /* pressure lag time constant */
    double Cmi;
    double Cmq;
    double ct;    
    double cn;    
    double* fn;   /* length of af->len */
    double* ft;   /* length of af->len */

    /* previously arrays */
    double alpha_f[2]; 
    double alpha[2];
    double Cn_p[2]; 
    double Cn_pot[2]; 
    double fn_p[2]; 
    double ft_p[2]; 
    double Cv[2];
    double adot[2];
    double dn[2]; 
    double dq[2]; 
    double qx[2]; 
    double Xn[2]; 
    double Yn[2]; 
    double Dp[2]; 
    double dCn_p[2];
    double tau[2]; 
    double Df_n[2]; 
    double Cn_v[2]; 
    double Df_t[2]; 
    double dqp[2]; 
    double Df_alpha_f_e[2];
};


ERROR_CODE dstall_initialize(struct Blade* blade, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_free(struct DStall* d);
ERROR_CODE dstall_initialize_fn(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_initialize_ft(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_get_coefficients(struct BladeElement* elem, struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_beddoes(struct BladeElement* elem, const double ALPHA, const double U, double* cla, double* cda, char* msg, ERROR_CODE* ierr);

ERROR_CODE dstall_beddoes_init(struct BladeElement* elem, const double dt, double ALPHA, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_attach(struct DStall* d, struct Airfoil* af, const double chord, const double U, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_separate(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_reset(struct DStall* d);

double dstall_saturation(const double X, const double VAL, const double SLOPE );
int dstall_sign(const double x);
ERROR_CODE dstall_vortex(struct DStall* d, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_set_Cv(struct DStall* d, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_set_fn_p(struct DStall* d, const int idx, const double aol, const double coeff, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_set_ft_p(struct DStall* d, const int idx, const double aol, const double coeff, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_check_leading_edge_separation(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_normal_coefficient(struct DStall* d, const double tf_e, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_tangent_coefficient(struct DStall* d, const double tf_e, char* msg, ERROR_CODE* ierr);
ERROR_CODE dstall_pitch_moment(struct DStall* d, struct Airfoil* af, const double tf_e, char* msg, ERROR_CODE* ierr);

#endif /* _DSTALL_H */
