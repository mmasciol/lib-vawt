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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sys.h"
#include "blade.h"
#include "ac.h"
#include "numerics.h"
#include "dstall.h"


ERROR_CODE
blade_set_number(struct ActuatorCylinder* ac, const int nb, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    ERROR_CODE success = SAFE;
    double azimuth = -999.9;

    // if (ac->z==NULL) {
    //     printf("<%p> %d\n",(void*)ac->z, ac->n_theta);
    //     CALL_CHECKERRK(FATAL_5, "AC vertical discretization must be set before an airfoil can be assigned");
    // }
    ac->blade = malloc(nb*sizeof(struct Blade)); CHECK_MALLOC(ac->blade, FATAL_8, "<Blade>");
    for(i=0 ; i<nb ; i++) {
        success = blade_initialize(&ac->blade[i], msg, ierr); CHECKERRK(FATAL_5, "");
        azimuth = DEG2RAD*(360.0/nb)*i;
        success = blade_set_azimuth(&ac->blade[i], azimuth, msg, ierr); CHECKERRK(FATAL_32, "");
        ac->nb += 1;
    }
    if (ac->nb!=nb) {
        CALL_CHECKERRK(FATAL_9, "Number of blades allocated not consistent with the number requested");
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
blade_set_delta(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    // ERROR_CODE success = SAFE;
    int i = 0;
    int j = 0;
    int m = 0;
    const int nb = ac->nb;
    double slope = 0.0;
    double dz = 0.0;
    double dr = 0.0;
    double* z = ac->z;
    for(i=0 ; i<nb ; i++) {
        m = ac->blade[i].m;
        for(j=0 ; j<m ; j++) {
            if (j==0) {
                dz = z[j]-z[j+1];
                dr = (ac->blade[i].elem[j].r) - (ac->blade[i].elem[j+1].r);
                slope = dr/dz;
            } else if (j==m-1) {
                dz = z[j-1]-z[j];
                dr = (ac->blade[i].elem[j-1].r) - (ac->blade[i].elem[j].r);
                slope = dr/dz;
            } else {
                dz = z[j]-z[j+1];
                dr = (1/(2.0*dz))*(ac->blade[i].elem[j+1].r - ac->blade[i].elem[j-1].r);
                slope = dr/dz;
            }
            ac->blade[i].elem[j].delta = atan(-slope);
        }
    }
    return SAFE;
}


ERROR_CODE 
blade_set_azimuth(struct Blade* blade, const double azimuth, char* msg, ERROR_CODE* ierr)
{
    blade->azimuth = azimuth;
    return SAFE;
}


ERROR_CODE 
blade_set_rmax(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;
    int nb = ac->nb;
    int ne = 0;
    struct Blade* blade = NULL;
    for (i=0 ; i<nb ; i++) {
        blade = ac->blade;
        ne = blade->m;
        for (j=0 ; j<ne ; j++) {
            if (blade->rmax<=blade->elem[j].r) {
                blade->rmax = blade->elem[j].r;
            }
        }
    }
    return SAFE;
}


ERROR_CODE
blade_increment_azimuth(struct ActuatorCylinder* ac, const double dazimuth, char* msg, ERROR_CODE* ierr)
{
    const int nb = ac->nb;
    int i = 0;
    for(i=0 ; i<nb ; i++) {
        ac->blade[i].azimuth += dazimuth;
        if (ac->blade[i].azimuth<0.0) {
            ac->blade[i].azimuth += 2*PI;
        } else if (2*PI<ac->blade[i].azimuth) {
            ac->blade[i].azimuth -= 2*PI;
        }
    }
    return SAFE;
}

ERROR_CODE
blade_increment_element_azimuth(struct ActuatorCylinder* ac, const double theta, char* msg, ERROR_CODE* ierr)
{
    const int nb = ac->nb;
    int nz = 0.0;
    int i = 0;
    int j = 0;
    double azimuth = 0.0;
    for(i=0 ; i<nb ; i++) {
        nz = ac->blade[i].m;
        azimuth = ac->blade[i].azimuth;
        for(j=0 ; j<nz ; j++) {
            ac->blade[i].elem[j].the = theta + azimuth + ac->blade[i].elem[j].swept_offset;
            if (ac->blade[i].elem[j].the<0.0) {
                ac->blade[i].elem[j].the += 2*PI;
            } else if (2*PI<ac->blade[i].elem[j].the) {
                ac->blade[i].elem[j].the -= 2*PI;
            }                
        }
    }
    return SAFE;    
}

ERROR_CODE 
blade_initialize(struct Blade* blade, char* msg, ERROR_CODE* ierr)
{
    blade->elem = NULL;
    blade->azimuth = -999.9;
    blade->Rb = NULL; 
    blade->Tb = NULL; 
    blade->Zb = NULL; 
    blade->Ri = NULL; 
    blade->Ti = NULL;
    blade->Q = NULL;
    blade->Qavg = -999.9;
    blade->rmax = -999.9;
    blade->n = -9999;
    blade->m = -9999;
    blade->vtk = NULL;
    return SAFE;
}


ERROR_CODE 
blade_free(struct Blade* blade, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    const int nz = blade->m;
    for (i=0 ; i<nz ; i++) {
        if (blade->elem[i].d) {
            dstall_free(blade->elem[i].d);
            FREE_OBJ(blade->elem[i].d);
        }
        blade->elem[i].af = NULL;        
    }
    FREE_OBJ(blade->elem);    
    FREE_OBJ(blade->Rb);
    FREE_OBJ(blade->Tb);
    FREE_OBJ(blade->Zb);
    FREE_OBJ(blade->Ri);
    FREE_OBJ(blade->Ti);
    FREE_OBJ(blade->Q);
    return SAFE;
}


ERROR_CODE
blade_allocate_memory(struct Blade* blade, char* msg, ERROR_CODE* ierr)
{
    const int m = blade->m; /* len(z) */
    int i = 0;
    // printf("%d\n\n\n",m);
    blade->elem = malloc(m*sizeof(struct BladeElement)); CHECK_MALLOC(blade->elem, FATAL_8, "Blade <BladeElement>");    
    blade->Rb = malloc(m*sizeof(double)); CHECK_MALLOC(blade->Rb, FATAL_8, "Blade <Rb>");
    blade->Tb = malloc(m*sizeof(double)); CHECK_MALLOC(blade->Tb, FATAL_8, "Blade <Tb>");
    blade->Zb = malloc(m*sizeof(double)); CHECK_MALLOC(blade->Zb, FATAL_8, "Blade <Zb>");
    blade->Ri = malloc(m*sizeof(double)); CHECK_MALLOC(blade->Ri, FATAL_8, "Blade <Ri>");
    blade->Ti = malloc(m*sizeof(double)); CHECK_MALLOC(blade->Ti, FATAL_8, "Blade <Ti>");
    blade->Q = malloc(m*sizeof(double)); CHECK_MALLOC(blade->Q, FATAL_8, "Blade <Q>");
    for (i=0 ; i<m ; i++) {
        blade->elem[i].af = NULL;
        blade->elem[i].d = NULL;
        blade->elem[i].chord = 0.0;
        blade->elem[i].r = 0.0;
        blade->elem[i].W = 0.0;
        blade->elem[i].phi = 0.0;
        blade->elem[i].alpha = 0.0;
        blade->elem[i].Re = 0.0;
        blade->elem[i].twist = 0.0;
        blade->elem[i].swept_offset = 0.0;
        blade->elem[i].the = 0.0;
        blade->elem[i].delta = 0.0;
        blade->Rb[i] = 0.0; 
        blade->Tb[i] = 0.0; 
        blade->Zb[i] = 0.0; 
        blade->Ri[i] = 0.0; 
        blade->Ti[i] = 0.0; 
        blade->Q[i] = 0.0; 
    }
    // blade->x_theta = malloc(n*sizeof(double)); CHECK_MALLOC(ac->coeff->x_theta, FATAL_8, "x_theta");       
    // blade->y_theta = malloc(n*sizeof(double)); CHECK_MALLOC(ac->coeff->y_theta, FATAL_8, "y_theta");       
    // blade->z_theta = malloc(n*sizeof(double)); CHECK_MALLOC(ac->coeff->z_theta, FATAL_8, "z_theta");       
    return SAFE;
CLEAN_UP:
    return FATAL;
}

