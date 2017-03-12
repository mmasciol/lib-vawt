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
#include <math.h>
#include "sys.h"
#include "domain.h"
#include "env.h"


size_t
env_size_meter(const void *el) 
{
    return sizeof(struct Environment);
}


ERROR_CODE 
env_initialize_list(struct Domain* domain, char* msg, ERROR_CODE* ierr)
{
    list_init(&domain->env);  
    list_attributes_copy(&domain->env, env_size_meter, 1); 
    return SAFE;
}


ERROR_CODE
env_new_environment(list_t* restrict env, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    struct Environment new_env;
    struct Environment* env_ = NULL;

    success = env_initialize(&new_env, msg, ierr); CHECKERRK(FATAL_5,"Failed to initialize <Environment>");
    list_append(env, &new_env);   
    env_ = (struct Environment*)list_get_at(env, list_size(env)-1); CHECK_MALLOC(env_, FATAL_8, "<Environment>");
    return SAFE;
CLEAN_UP:
    return FATAL;

}

ERROR_CODE
env_initialize(struct Environment* env, char* msg, ERROR_CODE* ierr)
{
    env->V = NULL;
    env->ref_h = -999.9; 
    env->u_inf = -999.9;
    env->rho = -999.9; 
    env->mu = -999.9;
    env->shear = -999.9;
    return SAFE;
}


ERROR_CODE
env_free(struct Environment* env, char* msg, ERROR_CODE* ierr)
{
    FREE_OBJ(env->V);
    return SAFE;
}


ERROR_CODE
env_set_rho(struct Environment* env, const double rho, char* msg, ERROR_CODE* ierr)
{
    printf("%f",env->rho);
    printf("    >%f  \n",rho);
    env->rho = rho;
    return SAFE;
}


ERROR_CODE
env_set_ref_height(struct Environment* env, const double ref_h, char* msg, ERROR_CODE* ierr)
{
    env->ref_h = ref_h;
    return SAFE;
}


ERROR_CODE
env_set_inflow_velocity(struct Environment* env, const double u, char* msg, ERROR_CODE* ierr)
{
    env->u_inf = u;
    return SAFE;
}


ERROR_CODE
env_set_mu(struct Environment* env, const double mu, char* msg, ERROR_CODE* ierr)
{
    env->mu = mu;
    return SAFE;
}


ERROR_CODE
env_set_shear(struct Environment* env, const double shear, char* msg, ERROR_CODE* ierr)
{
    env->shear = shear;
    return SAFE;
}


double
env_get_rho(struct Environment* env, char* msg, ERROR_CODE* ierr)
{
    return env->rho;
}


double
env_get_mu(struct Environment* env, char* msg, ERROR_CODE* ierr)
{
    return env->mu;
}


double
env_get_shear(struct Environment* env, char* msg, ERROR_CODE* ierr)
{
    return env->shear;
}


ERROR_CODE
env_set_wind_velocity(struct Environment* env, const double u, char* msg, ERROR_CODE* ierr)
{
    env->u_inf = u;
    return SAFE;
}


ERROR_CODE 
env_wind_shear(struct Environment* env, const double* z, const int nz, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    const double ref_h = env->ref_h;
    const double shear = env->shear;
    if (!(env->V)) {
        env->V = malloc((nz)*sizeof(double)); CHECK_MALLOC((env->V), FATAL_8, "V");
    }
    
    for(i=0 ; i<nz ; i++) {
        env->V[i] = (env->u_inf)*pow((z[i]/ref_h), shear);
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}
