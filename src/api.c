/****************************************************************
 *   Copyright (C) 2015 mdm                                     *
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
#include <stdbool.h>
#include <math.h>
#include "api.h"
#include "sys.h"
#include "error.h"
#include "blade.h"
#include "af.h"
#include "ac.h"
#include "domain.h"
#include "blade.h"
#include "numerics.h"
#include "dstall.h"
#include "env.h"

#include "./simclist/simclist.h"
#include "./bstring/bstrlib.h"

#if defined(WITH_VTK)
#    include "./vtk/render_api.h"
#endif

extern const char ERROR_STRING[][1024];


EXTERNCALL struct Airfoil*
C_get_airfoil(struct Domain* domain, const int index, char* msg, ERROR_CODE* ierr) 
{
    bstring b = NULL;
    ERROR_CODE success = SAFE;
    struct Airfoil* af = (struct Airfoil*)list_get_at(&domain->af, index);      
    if (!af) {
        b = bformat("Airfoil index <%d>", index);
        CALL_CHECKERRK(FATAL_28,(const char*)b->data);
        
    }
    return af;
CLEAN_UP:
    bdestroy(b);
    return NULL;
}


EXTERNCALL void
C_set_dynamic_stall(struct ActuatorCylinder* ac, bool flag, char* msg, ERROR_CODE* ierr)
{
    ac->dstall = flag;
    return;
}


EXTERNCALL struct Airfoil* 
C_new_airfoil(struct Domain* domain, char file_name[255], char type_name[32], char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    struct Airfoil* af_ = NULL;
    const int n = list_size(&domain->af);

    ierr_reset(msg, ierr);
    success = af_new_foil(&domain->af, msg, ierr); CHECKERRK(FATAL_5,"");
    success = af_read_file(&domain->af, file_name, type_name, msg, ierr); CHECKERRK(FATAL_5,"");
    af_ = (struct Airfoil*)list_get_at(&domain->af, n);
    success = af_set_beddoes_coefficients(af_,msg, ierr); CHECKERRK(FATAL_11,"Error setting dynamic stall coefficients");    
    if (!af_) {
        CALL_CHECKERRK(FATAL_28,"Could not allocate a new <Airfoil>");                    
    }
    return af_;
CLEAN_UP:
    return NULL;
}


EXTERNCALL struct ActuatorCylinder*
C_new_actuator_cylinder(struct Domain* domain, const int nb, const int nt, double* z, const int nz, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    struct ActuatorCylinder* ac_ = NULL;
    const int n = list_size(&domain->ac);

    ierr_reset(msg, ierr);
    success = ac_new_turbine(&domain->ac, msg, ierr); CHECKERRK(FATAL_5,"");
    ac_ = (struct ActuatorCylinder*)list_get_at(&domain->ac, n);
    if (!ac_) {
        CALL_CHECKERRK(FATAL_28,"Could not allocate a new <ActuatorCylinder>");                    
    }
    success = blade_set_number(ac_, nb, msg, ierr); CHECKERRK(FATAL_30, "");
    success = ac_set_z(ac_, z, nz, msg, ierr); CHECKERRK(FATAL_14, "Blade cross-sectional vertical height");
    success = ac_set_n_theta(ac_, nt, msg, ierr); CHECKERRK(FATAL_14, "Rotational discretization");
    return ac_;
CLEAN_UP:
    return NULL;
}


EXTERNCALL struct Environment*
C_new_environment(struct Domain* domain, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    struct Environment* env_ = NULL;
    const int n = list_size(&domain->env);

    ierr_reset(msg, ierr);
    success = env_new_environment(&domain->env, msg, ierr); CHECKERRK(FATAL_5,"");
    env_ = (struct Environment*)list_get_at(&domain->env, n);
    if (!env_) {
        CALL_CHECKERRK(FATAL_28,"Could not allocate a new <Environment>");                    
    }
    return env_;
CLEAN_UP:
    return NULL;
}


EXTERNCALL struct Airfoil*
C_set_blend_airfoils(struct Domain* domain, struct Airfoil* af1, struct Airfoil* af2, const double weight, char type_name[32], char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    struct Airfoil* af_ = NULL;
    const int n = list_size(&domain->af);
    
    ierr_reset(msg, ierr);
    success = af_new_foil(&domain->af, msg, ierr); CHECKERRK(FATAL_5,"");
    success = af_blend(&domain->af, af1, af2, weight, msg, ierr); CHECKERRK(FATAL_11,"");
    af_ = (struct Airfoil*)list_get_at(&domain->af, n);
    copy_target_string(af_->type, (unsigned char*)type_name);
    success = af_set_beddoes_coefficients(af_,msg, ierr); CHECKERRK(FATAL_11,"Error setting dynamic stall coefficients");    
    if (!af_) {
        CALL_CHECKERRK(FATAL_28,"Could not blend a new airfoil");                    
    }
    return af_;
CLEAN_UP:
    return NULL;
}


EXTERNCALL struct Domain*
C_new_domain(char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    struct Domain* domain = malloc(sizeof(struct Domain));   

    ierr_reset(msg, ierr); 
    if (!domain) {
        success = FATAL;
        CHECKERRK(FATAL_8,"<Domain>");
    }
    success = domain_initialize(domain, msg, ierr); 
    return domain;
CLEAN_UP:
    return NULL;
}


#if defined(WITH_VTK)
EXTERNCALL void*
C_set_renderer_engine(char* msg, ERROR_CODE* ierr)
{
    ierr_reset(msg, ierr);    
    return new_vtk();
}
#else
EXTERNCALL void*
C_set_renderer_engine(char* msg, ERROR_CODE* ierr)
{
    ierr_reset(msg, ierr);
    CHECKERRQ(WARNING_5,"");
    return NULL;
}
#endif


#if defined(WITH_VTK)
EXTERNCALL void 
C_vawt_render(void* vtk, struct Domain* domain, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;
    int nz = 0;
    int nb = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    struct ActuatorCylinder* ac_ = NULL;
    
    ierr_reset(msg, ierr);

    list_iterator_start(&domain->ac); 
    while (list_iterator_hasnext(&domain->ac)) { /* tell whether more values available */ 
        ac_ = (struct ActuatorCylinder*)list_iterator_next(&domain->ac);
        nb = ac_->nb;
        for (i=0 ; i<nb ; i++) {
            nz = ac_->nz;
            for (j=0 ; j<nz ; j++) {
                if (ac_->blade[i].vtk==NULL) {
                    ac_->blade[i].vtk = malloc(nz*sizeof(void*));
                }                
                x = 0.0;
                y = ac_->blade[i].elem[j].r;
                z = ac_->z[j];
                ac_->blade[i].vtk[j] = new_point();
                C_vtk_airfoil_hook(ac_->blade[i].vtk[j], vtk, ac_->blade[i].elem[j].chord);
                C_vtk_airfoil_rotate(ac_->blade[i].vtk[j], 0.0, 0.0, RAD2DEG*(ac_->blade[i].azimuth+ac_->blade[i].elem[j].swept_offset));
                C_vtk_airfoil_displace(ac_->blade[i].vtk[j], x, y, z);                
                C_vtk_airfoil_transform(ac_->blade[i].vtk[j]);
            }
        }
    }
    list_iterator_stop(&domain->ac); /* ending the iteration "session" */        

    vtk_start(vtk);    

    list_iterator_start(&domain->ac); 
    while (list_iterator_hasnext(&domain->ac)) { /* tell whether more values available */ 
        ac_ = (struct ActuatorCylinder*)list_iterator_next(&domain->ac);
        nb = ac_->nb;
        for (i=0 ; i<nb ; i++) {
            nz = ac_->nz;
            for (j=0 ; j<nz ; j++) {
                delete_point(ac_->blade[i].vtk[j]);
            }
            FREE_OBJ(ac_->blade[i].vtk);
        }
    }
    list_iterator_stop(&domain->ac); /* ending the iteration "session" */        

    delete_vtk(vtk);
    return;
}
#else
EXTERNCALL void 
C_vawt_render(void* vtk, struct Domain* domain, char* msg, ERROR_CODE* ierr)
{
    ierr_reset(msg, ierr);
    CHECKERRQ(WARNING_5,"");
    return;
}
#endif


EXTERNCALL void
C_set_initial_theta(struct ActuatorCylinder* ac, const double the0, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_initial_theta_offset(ac, the0, msg, ierr); CHECKERRK(FATAL_29,"");
    success = blade_increment_azimuth(ac, ac->the0, msg, ierr); CHECKERRK(FATAL_29,"");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void 
C_end_program(struct Domain* domain, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = domain_free(domain, msg, ierr); CHECKERRK(FATAL_23, "Domain is not allocated; could not deallocate/free.");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_ac_commit(struct ActuatorCylinder* ac, struct Environment* env, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    // int i = 0;
    ierr_reset(msg, ierr);
    success = blade_set_delta(ac, msg, ierr); CHECKERRK(FATAL_8, "");
    success = blade_set_rmax(ac, msg, ierr); CHECKERRK(FATAL_8, "");
    success = ac_set_swept_area(ac, msg, ierr); CHECKERRK(FATAL_5,"Swept area could not be calculated");
    success = ac_generate_matrix(ac, msg, ierr); CHECKERRK(FATAL_12,"");
    success = ac_generate_memory(ac, msg, ierr); CHECKERRK(FATAL_20,"");
    // if (ac->dstall) {
    //     for (i=0 ; i<ac->nb ; i++) {        
    //         success = dstall_initialize(&ac->blade[i], msg, ierr); CHECKERRK(FATAL_5,"Could not initialize dynamic stall parameters");
    //         // success = dstall_allocate(&ac->blade[i], msg, ierr); CHECKERRK(FATAL_8,"Could not allocate dynamic stall parameters");
    //     }
    // }
    success = ac_allocate_perf(ac, msg, ierr); CHECKERRK(FATAL_8, "Could not allocate Actuator Cylinder <ACPerformance>");
    ac->env = env;
    /* @todo: things to check for:
     *        1) number of airfoils assigned equals the numebr of blades
     *        2) z (and other arrays) must be equally spaced
     */
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_airfoil_chord(struct ActuatorCylinder* ac, const double chord, const int ib, const int ie, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_airfoil_chord(ac, chord, ib, ie, msg, ierr); CHECKERRK(FATAL_14, "Blade chord");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_airfoil_twist(struct ActuatorCylinder* ac, const double twist, const int ib, const int ie, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_airfoil_twist(ac, twist, ib, ie, msg, ierr); CHECKERRK(FATAL_14, "Blade twist");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_airfoil_swept_offset(struct ActuatorCylinder* ac, const double swept, const int ib, const int ie, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_airfoil_swept_offset(ac, swept, ib, ie, msg, ierr); CHECKERRK(FATAL_14, "Blade twist");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_airfoil_radius(struct ActuatorCylinder* ac, const double r, const int ib, const int ie, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_airfoil_radius(ac, r, ib, ie, msg, ierr); CHECKERRK(FATAL_14, "Blade cross-sectional radius");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_airfoil_geometry(struct ActuatorCylinder* ac, struct Airfoil* af, const int ib, const int ie, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    int i = 0;
    ierr_reset(msg, ierr);
    for (i=0 ; i<ac->nb ; i++) {
        if (!ac->blade[i].elem) {
            ac->blade[i].m = ac->nz;
            success = blade_allocate_memory(&ac->blade[i], msg, ierr); CHECKERRK(FATAL_5,"Could not allocate blade");            
        }
    }
    success = ac_set_airfoil_geometry(ac, af, ib, ie, msg, ierr); CHECKERRK(FATAL_14, "Aifoil could not be assigned to actuator cylinder");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_ac_solve(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    int i = 0;
    
    ierr_reset(msg, ierr);
    if (ac->dstall) {
        for (i=0 ; i<ac->nb ; i++) {        
            success = dstall_initialize(&ac->blade[i], msg, ierr); CHECKERRK(FATAL_5,"Could not initialize dynamic stall parameters");
        }
    }
    success = env_wind_shear(ac->env, ac->z, ac->nz, msg, ierr); CHECKERRK(FATAL_22,"");
    success = ac_force(ac, msg, ierr); CHECKERRK(FATAL_13, "");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_ac_performance(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;
    int k = 0;
    const int nz = ac->nz;
    const int nt = ac->nt;
    const int nb = ac->nb;
    ERROR_CODE success = SAFE;
    double* fr = NULL;
    double* ft = NULL;
    double* fz = NULL;
    double* Q = NULL;
    double P = -999.9;
    struct Environment* env = ac->env;    
    const double rho = env->rho;
    const double q = 0.5*rho*pow(env->u_inf,2);
    double theta = -999.9;
    double Q_bar = 0.0;
    double r = 0.0;
  
    fr = malloc(nt*sizeof(double)); CHECK_MALLOC(fr, FATAL_8, "AC performance calculation <fr>");
    ft = malloc(nt*sizeof(double)); CHECK_MALLOC(ft, FATAL_8, "AC performance calculation <ft>");
    fz = malloc(nt*sizeof(double)); CHECK_MALLOC(fz, FATAL_8, "AC performance calculation <fz>");
    Q = malloc(nt*sizeof(double)); CHECK_MALLOC(Q, FATAL_8, "AC performance calculation <Q>");

    ac->perf->Ct = 0.0;
    for (i=0 ; i<ac->perf->n ; i++) { 
        ac->perf->fx[i] = 0.0;
        ac->perf->fy[i] = 0.0;
        ac->perf->fz[i] = 0.0;
        ac->perf->Q[i] = 0.0;
    }        
 
    for (i=0 ; i<nb ; i++) {
        for (j=0 ; j<nt ; j++) { /* Rp has dimensions m=n_z, n=n_theta */
            success = blade_increment_element_azimuth(ac, ac->theta[j], msg, ierr);
            success = ac_power(ac, &ac->blade[i], /*theta,*/ msg, ierr); CHECKERRK(FATAL_21, "");
             
            /* integrate along the length of the blade */
            success = integrate(&fr[j], ac->blade[i].Ri, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<fr>");
            success = integrate(&ft[j], ac->blade[i].Ti, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<ft>");
            success = integrate(&fz[j], ac->blade[i].Zb, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<fz>");
 
            for (k=0 ; k<nz ; k++) {
                theta = ac->blade[i].elem[k].the;
                ac->perf->fx[j] += -fr[j]*sin(theta) - ft[j]*cos(theta); 
                ac->perf->fy[j] += fr[j]*cos(theta) - ft[j]*sin(theta);  
                ac->perf->fz[j] += fz[j];
            }

            for (k=0 ; k<nz ; k++) {
                r = ac->blade[i].elem[k].r;
                ac->blade[i].Q[k] = r*(ac->blade[i].Tb[k]);
                // sumx += ac->z[k]*(-(ac->blade[i].Ri[k])*(sin(ac->blade[i].elem[k].the)) - (ac->blade[i].Ti[k])*(cos(ac->blade[i].elem[k].the)));
                // sumy += ac->z[k]*((ac->blade[i].Ri[k])*(cos(ac->blade[i].elem[k].the)) - (ac->blade[i].Ti[k])*(sin(ac->blade[i].elem[k].the)));                
            }
            success = integrate(&Q[j], ac->blade[i].Q, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<Q>");
            ac->perf->Q[j] += Q[j];
        }
        success = integrate(&ac->blade[i].Qavg, Q, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Qbar>");
    }
    
    /* average Q */
    for (i=0 ; i<nb ; i++) {
        Q_bar += ac->blade[i].Qavg;
    }
     
    Q_bar /= nb;
    Q_bar *= nb/(2*PI);
    P = ac->omega*Q_bar;
    ac->perf->Cp = P/(q*(env->u_inf)*ac->swept_area); 

    for(i=0 ; i<nt ; i++) {
        ac->perf->fx[i] = ac->perf->fx[i]/(q*ac->swept_area*nz);
        ac->perf->fy[i] = ac->perf->fy[i]/(q*ac->swept_area*nz);
    }

    success = integrate(&ac->perf->Cx, ac->perf->fx, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Cx>");
    success = integrate(&ac->perf->Cy, ac->perf->fy, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Cy>");
    success = integrate(&ac->perf->Cz, ac->perf->fz, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Cz>");
    ac->perf->Ct = ac->perf->Cx/(2*PI);
    FREE_OBJ(fr);
    FREE_OBJ(ft);
    FREE_OBJ(fz);    
    FREE_OBJ(Q);
    return;
CLEAN_UP:
    FREE_OBJ(fr);
    FREE_OBJ(ft);
    FREE_OBJ(fz);
    FREE_OBJ(Q);
    return;
}

EXTERNCALL double
C_ac_get_stationary_performance(struct ActuatorCylinder* ac, const double the0, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;    
    const int nz = ac->nz;
    const int nb = ac->nb;
    ERROR_CODE success = SAFE;
    // double* Q = NULL;
    // double P = -999.9;
    // struct Environment* env = ac->env;    
    // const double rho = env->rho;
    // const double q = 0.5*rho*pow(env->u_inf,2);
    // double* x = NULL;
    // double* y = NULL;
    // double* z = NULL;
    // double theta = -999.9;
    // double Q_bar = 0.0;
    // double sumx = 0.0;
    // double sumy = 0.0;
    // double sumfx = 0.0;
    // double sumfy = 0.0;
    // double r = 0.0;

    

    double fr = 0.0;// malloc(nt*sizeof(double)); CHECK_MALLOC(fr, FATAL_8, "AC performance calculation <fr>");
    double ft = 0.0;// malloc(nt*sizeof(double)); CHECK_MALLOC(ft, FATAL_8, "AC performance calculation <ft>");
    double fz = 0.0;// malloc(nt*sizeof(double)); CHECK_MALLOC(fz, FATAL_8, "AC performance calculation <fz>");
    double Q = 0.0;// malloc(nt*sizeof(double)); CHECK_MALLOC(fz, FATAL_8, "AC performance calculation <fz>");
//     Q = malloc(nt*sizeof(double)); CHECK_MALLOC(Q, FATAL_8, "AC performance calculation <Q>");
//     x = malloc(nt*sizeof(double)); CHECK_MALLOC(x, FATAL_8, "AC performance calculation <x>");
//     y = malloc(nt*sizeof(double)); CHECK_MALLOC(y, FATAL_8, "AC performance calculation <y>");
//     z = malloc(nt*sizeof(double)); CHECK_MALLOC(z, FATAL_8, "AC performance calculation <z>");
//     // ac->perf->Ct = 0.0;
//     for (i=0 ; i<ac->perf->n ; i++) { 
//         ac->perf->fx[i] = 0.0;
//         ac->perf->fy[i] = 0.0;
//         ac->perf->fz[i] = 0.0;
//         ac->perf->Q[i] = 0.0;
//     }        
//
    ierr_reset(msg, ierr);
    success = env_wind_shear(ac->env, ac->z, ac->nz, msg, ierr); CHECKERRK(FATAL_22,"");
    for (i=0 ; i<nb ; i++) {
        fr = 0.0;
        ft = 0.0;
        fz = 0.0;

        for (j=0 ; j<nz ; j++) {
            success = blade_increment_element_azimuth(ac, the0, msg, ierr);
            /* @change */
            // double theta = the0 + (ac->blade[i].azimuth) + (ac->blade[i].elem[j].swept_offset);            
            // if (theta<0.0) {
            //     theta+=2*PI;
            // } else if (2*PI<theta) {
            //     theta-=2*PI;
            // }
            if (ac->blade[i].elem[j].the<0.0) {
                ac->blade[i].elem[j].the += 2*PI;
            } else if (2*PI<ac->blade[i].elem[j].the) {
                ac->blade[i].elem[j].the -= 2*PI;
            }
            double Uinf = ac->env->V[j];
            /* @change */
            // double Vt = Uinf*cos(theta);
            // double Vn = Uinf*sin(theta);
            double Vt = Uinf*cos(ac->blade[i].elem[j].the);
            double Vn = Uinf*sin(ac->blade[i].elem[j].the);
            double W = sqrt(pow(Vt,2) + pow(Vn,2));
            double q = 0.5*(ac->env->rho)*pow(W,2);
            double phi = atan2(Vn, Vt);
            double alpha = phi - ac->blade[i].elem[j].twist;
            // printf("%f  %f\n",Uinf, theta);
            if (alpha<-PI) {
                alpha+=PI;
            } else if (PI<alpha) {
                alpha-=PI;
            }

            double value = 0.0;
            /* i is supposed to be in nt */
            success = af_get_cl(&value, ac->blade[i].elem[j].af, alpha, ac->work->Re[i], msg, ierr); CHECKERRK(FATAL_26, "<af_get_cl>");
            double cl = value;        
            success = af_get_cd(&value, ac->blade[i].elem[j].af, alpha, ac->work->Re[i], msg, ierr); CHECKERRK(FATAL_26, "<af_get_cd>");
            double cd = value;
            double cn = cl*cos(phi) + cd*sin(phi);
            double ct = cl*sin(phi) - cd*cos(phi);
        
            const double chord = ac->blade[i].elem[j].chord; /* @todo: this is solved for a single blade, needs to change */
            const double delta = ac->blade[i].elem[j].delta;
            double r = ac->blade[i].elem[j].r;
            fr += -cn*q*(chord);
            ft += ct*q*(chord)/cos(delta);
            fz += -cn*q*(chord)*tan(delta);
            Q += ft*r;

//             success = ac_power_2(ac, &ac->blade[i], theta, msg, ierr); CHECKERRK(FATAL_21, "");
//             
//             /* integrate along the length of the blade */
//             success = integrate(&fr[j], ac->blade[i].Ri, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<fr>");
//             success = integrate(&ft[j], ac->blade[i].Ti, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<ft>");
//             success = integrate(&fz[j], ac->blade[i].Zb, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<fz>");
//             
//             ac->perf->fx[j] += -fr[j]*sin(theta) - ft[j]*cos(theta); /* @todo: assumption here is the blade is vertical and straight (not swept) */
//             ac->perf->fy[j] += fr[j]*cos(theta) - ft[j]*sin(theta);  /* @todo: assumption here is the blade is vertical and straight (not swept) */
//             ac->perf->fz[j] += fz[j];
//             for (k=0 ; k<nz ; k++) {
//                 r = ac->blade[i].elem[k].r;
//                 ac->blade[i].Q[k] = (r)*(ac->blade[i].Tb[k]);
//                 sumx += ac->z[k]*(-(ac->blade[i].Ri[k])*(sin(theta)) - (ac->blade[i].Ti[k])*(cos(theta)));
//                 sumy += ac->z[k]*((ac->blade[i].Ri[k])*(cos(theta)) - (ac->blade[i].Ti[k])*(sin(theta)));
//             }
//             success = integrate(&Q[j], ac->blade[i].Q, ac->z, nz, msg, ierr); CHECKERRK(FATAL_27, "<Q>");
//             ac->perf->Q[j] += Q[j];
        }
//         success = integrate(&ac->blade[i].Qavg, Q, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Qbar>");
    }
    return Q;
// 
//     /* average Q */
//     for (i=0 ; i<nb ; i++) {
//         Q_bar += ac->blade[i].Qavg;
//     }
//     
//     Q_bar /= nb;
//     Q_bar *= nb/(2*PI);
//     P = ac->omega*Q_bar;
//     ac->perf->Cp = P/(q*(env->u_inf)*ac->swept_area); 
// 
//     for (i=0 ; i<nt ; i++) {
//         x[i] = ac->perf->fx[i]/(q*ac->swept_area);
//         y[i] = ac->perf->fy[i]/(q*ac->swept_area);
//         z[i] = ac->perf->fz[i]/(q*ac->swept_area);
//         sumfx += ac->perf->fx[i];
//         sumfy += ac->perf->fy[i];
//     }
//     //printf("%f  %f\n",sumx/sumfx, sumy/sumfy);
//     //checkpoint();
//     success = integrate(&ac->perf->Cx, x, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Cx>");
//     success = integrate(&ac->perf->Cy, y, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Cy>");
//     success = integrate(&ac->perf->Cz, z, ac->theta, nt, msg, ierr); CHECKERRK(FATAL_27, "<Cz>");
//     //printf("%f  %f  %f\n",ac->perf->Cx, ac->perf->Cy, ac->perf->Cz);
//     ac->perf->Ct = sqrt(pow(ac->perf->Cx,2) + pow(ac->perf->Cy,2))/(2*PI);
// 
//     FREE_OBJ(fr);
//     FREE_OBJ(ft);
//     FREE_OBJ(fz);    
//     FREE_OBJ(Q);
//     FREE_OBJ(x);
//     FREE_OBJ(y);
//     FREE_OBJ(z);
CLEAN_UP:
//     FREE_OBJ(fr);
//     FREE_OBJ(ft);
//     FREE_OBJ(fz);
//     FREE_OBJ(Q);
//     FREE_OBJ(x);
//     FREE_OBJ(y);
//     FREE_OBJ(z);
    return -999.9;
}


EXTERNCALL void
C_ac_instantaneous_force(struct ActuatorCylinder* ac, /* theta should be an argument */ char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    int i = 0;
    int j = 0;
    // int k = 0;
    const int nt = ac->nt;
    const int nb = ac->nb;
    // const int m = ac->n_z;
    // double Qbar = 0.0;
    // double P = 0.0;
    // double q = 0.0;
    ierr_reset(msg, ierr);

    for (i=0 ; i<nb ; i++) {
        for (j=0 ; j<nt ; j++) { /* Rp has dimensions m=n_z, n=n_theta */           
            success = ac_power(ac, &ac->blade[i], /*ac->theta[j], */msg, ierr); CHECKERRK(FATAL_21, "");
            // ac->coeff->R_theta_1_blade[j] = 0.0;
            // ac->coeff->T_theta_1_blade[j] = 0.0;
            // ac->coeff->Z_theta_1_blade[j] = 0.0;
            // success = integrate(&ac->coeff->R_theta_1_blade[j], ac->blade[i].Ri, ac->z, m, msg, ierr); CHECKERRK(FATAL_27, "R_Theta");
            // success = integrate(&ac->coeff->T_theta_1_blade[j], ac->blade[i].Ti, ac->z, m, msg, ierr); CHECKERRK(FATAL_27, "T_Theta");
            // success = integrate(&ac->coeff->Z_theta_1_blade[j], ac->blade[i].Zb, ac->z, m, msg, ierr); CHECKERRK(FATAL_27, "Z_Theta");
            // for (k=0 ; k<m ; k++) {
            //     ac->T[k] = (ac->r[k])*ac->blade[i].Tb[k];
            // }
            // success = integrate(&ac->coeff->Q_theta_1_blade[j], ac->T, ac->z, m, msg, ierr); CHECKERRK(FATAL_27, "Q_Theta");
        }
    }
    
    // success = integrate(&Qbar, ac->coeff->Q_theta_1_blade, ac->theta, n, msg, ierr); CHECKERRK(FATAL_27, "Qbar");
    // Qbar = ac->nb/(2*PI)*Qbar;
    // P = ac->omega*Qbar;
    // q = 0.5*(ac->rho)*pow(ac->Uinf,2);
    // ac->CP = P/(q*(ac->Uinf)*30.0); // 30 should be the equvalent of self.S    
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_z(struct ActuatorCylinder* ac, double* z, const int len_z, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_z(ac, z, len_z, msg, ierr); CHECKERRK(FATAL_14, "Blade cross-sectional vertical height");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_rho(struct Environment* env, const double rho, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = env_set_rho(env, rho, msg, ierr); CHECKERRK(FATAL_14, "Atmospheric density (rho)");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_mu(struct Environment* env, const double mu, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = env_set_mu(env, mu, msg, ierr); CHECKERRK(FATAL_14, "Atmospheric dynamic viscosity (mu)");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_shear(struct Environment* env, const double shear, char* msg, ERROR_CODE* ierr) 
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = env_set_shear(env, shear, msg, ierr); CHECKERRK(FATAL_14, "Atmospheric shear ratio");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_extrapolate_airfoil(struct Airfoil* af, const double cd_max, const double aspect_ratio, const double cd_min, const int n_alpha, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = af_extrapolate(af, cd_max, aspect_ratio, cd_min, n_alpha, msg, ierr); CHECKERRK(FATAL_15,"");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_wind_velocity(struct Environment* env, const double u, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = env_set_wind_velocity(env, u, msg, ierr); CHECKERRK(FATAL_22,"");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_ref_height(struct Environment* env, const double ref_h, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = env_set_ref_height(env, ref_h, msg, ierr); CHECKERRK(FATAL_22,"");
    return;
CLEAN_UP:
    return;
}


EXTERNCALL void
C_set_angular_velocity(struct ActuatorCylinder* ac, const double omega, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    ierr_reset(msg, ierr);
    success = ac_set_angular_velocity(ac, omega, msg, ierr); CHECKERRK(FATAL_24,"");
    return;
CLEAN_UP:
    return;
}
