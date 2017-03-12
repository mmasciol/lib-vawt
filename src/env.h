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


#ifndef _ENV_H
#define _ENV_H

#include <stdio.h>
#include "error.h"
#include "./simclist/simclist.h"


struct Domain;


struct Environment {
    double* V; /* length of nz */
    double ref_h; 
    double u_inf;
    double rho;
    double mu;
    double shear;
};


size_t env_size_meter(const void *el);
ERROR_CODE env_free(struct Environment* env, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_new_environment(list_t* restrict env, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_initialize_list(struct Domain* domain, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_initialize(struct Environment* env, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_set_rho(struct Environment* env, const double rho, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_set_mu(struct Environment* env, const double mu, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_set_shear(struct Environment* env, const double shear, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_set_ref_height(struct Environment* env, const double ref_h, char* msg, ERROR_CODE* ierr);
double env_get_rho(struct Environment* env, char* msg, ERROR_CODE* ierr);
double env_get_mu(struct Environment* env, char* msg, ERROR_CODE* ierr);
double env_get_shear(struct Environment* env, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_set_wind_velocity(struct Environment* env, const double u, char* msg, ERROR_CODE* ierr);
ERROR_CODE env_wind_shear(struct Environment* env, const double* z, const int nz, char* msg, ERROR_CODE* ierr);


#endif /* _ENV_H */
