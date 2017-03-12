/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                    *
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


#ifndef _API_H
#define _API_H


#include "error.h"


#if defined(_MSC_VER)
#    define EXTERNCALL __declspec( dllexport )
#else
#    if defined(_MINGW)
#        define EXTERNCALL __declspec( dllexport )
#    else
#        define EXTERNCALL 
#    endif
#endif


struct Airfoil;
struct ActuatorCylinder;
struct Domain;
struct Environment;


EXTERNCALL struct Airfoil* C_new_airfoil(struct Domain* domain, char file_name[255], char type_name[32], char* msg, ERROR_CODE* ierr) ;
EXTERNCALL struct Domain* C_new_domain(char* msg, ERROR_CODE* ierr);
EXTERNCALL struct Airfoil* C_set_blend_airfoils(struct Domain* domain, struct Airfoil* af1, struct Airfoil* af2, const double weight, char type_name[32], char* msg, ERROR_CODE* ierr);
EXTERNCALL struct Environment* C_new_environment(struct Domain* domain, char* msg, ERROR_CODE* ierr);
EXTERNCALL struct ActuatorCylinder* C_new_actuator_cylinder(struct Domain* domain, const int nb, const int nt, double* z, const int nz, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_extrapolate_airfoil(struct Airfoil* af, const double cd_max, const double aspect_ratio, const double cd_min, const int n_alpha, char* msg, ERROR_CODE* ierr);
// EXTERNCALL void C_set_solidity_ratio(struct ActuatorCylinder* ac, const double solidity_ratio, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_assign_airfoil(struct ActuatorCylinder* ac, struct Airfoil* af, const int i_blade, const int index, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_r(struct ActuatorCylinder* ac, double* r, const int len_r, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_chord(struct ActuatorCylinder* ac, double* chord, const int len_chord, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_twist(struct ActuatorCylinder* ac, double* twist, const int len_twist, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_initial_theta(struct ActuatorCylinder* ac, const double the0, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_load_airfoil_file(struct Airfoil* af, char file_name[255], char type_name[32], char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_end_program(struct Domain* domain, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_shear(struct Environment* env, const double shear, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_rho(struct Environment* env, const double rho, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_mu(struct Environment* env, const double mu, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_ac_instantaneous_force(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_ac_performance(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_ac_solve(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_ac_commit(struct ActuatorCylinder* ac, struct Environment* env, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_wind_velocity(struct Environment* env, const double u, char* msg, ERROR_CODE* ierr);    
EXTERNCALL void C_set_angular_velocity(struct ActuatorCylinder* ac, const double omega, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_dynamic_stall(struct ActuatorCylinder* ac, bool flag, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_airfoil_chord(struct ActuatorCylinder* ac, const double chord, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_airfoil_twist(struct ActuatorCylinder* ac, const double twist, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_airfoil_radius(struct ActuatorCylinder* ac, const double r, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_airfoil_geometry(struct ActuatorCylinder* ac, struct Airfoil* af, const int ib, const int ie, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_ref_height(struct Environment* env, const double ref_h, char* msg, ERROR_CODE* ierr);
EXTERNCALL void C_set_airfoil_swept_offset(struct ActuatorCylinder* ac, const double swept, const int ib, const int ie, char* msg, ERROR_CODE* ierr);


#endif /* _API_H */
