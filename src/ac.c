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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cminpack.h"
#include "sys.h"
#include "ac.h"
#include "af.h"
#include "numerics.h"
#include "domain.h"
#include "blade.h"
#include "dstall.h"
#include "env.h"
#include "./bstring/bstrlib.h"



static inline int 
ij(const int i, const int j, const int m)
{
    return i*m + j;
}


ERROR_CODE 
ac_initialize(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ac->pc = NULL;
    ac->blade = NULL;
    ac->env = NULL;
    ac->A = NULL;
    ac->theta = NULL;   
    ac->z = NULL;
    ac->fr = NULL;
    ac->ft = NULL;
    ac->fz = NULL;
    ac->omega = -999.9;
    ac->the0 = -999.9;
    ac->swept_area = -999.9;
    ac->nt = -9999;
    ac->nz = -9999;
    ac->nb = 0;
    ac->index = -9999;
    ac->dstall = false;
    ac->work = NULL;
    ac->perf = NULL;
    ac->passed_msg = NULL;
    ac->passed_ierr = NULL;
    return SAFE;
}


ERROR_CODE
ac_perf_initialize(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ac->perf->fx = NULL;
    ac->perf->fy = NULL;
    ac->perf->fz = NULL;
    ac->perf->Q = NULL;
    ac->perf->Cp = -9999;
    ac->perf->Ct = -9999;
    ac->perf->Cx = -9999;
    ac->perf->Cy = -9999;
    ac->perf->Cz = -9999;
    ac->perf->n = ac->nt;
    return SAFE;
}


size_t
ac_size_meter(const void *el) 
{
    return sizeof(struct ActuatorCylinder);
}


ERROR_CODE 
ac_initialize_list(struct Domain* domain, char* msg, ERROR_CODE* ierr)
{
    list_init(&domain->ac);  
    list_attributes_copy(&domain->ac, ac_size_meter, 1); 
    return SAFE;
}


ERROR_CODE
ac_new_turbine(list_t* restrict ac, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    struct ActuatorCylinder new_ac;
    struct ActuatorCylinder* ac_ = NULL;

    success = ac_initialize(&new_ac, msg, ierr); CHECKERRK(FATAL_5,"Failed to initialize <ActuatorCylinder>");
    list_append(ac, &new_ac);   
    ac_ = (struct ActuatorCylinder*)list_get_at(ac, list_size(ac)-1); CHECK_MALLOC(ac_, FATAL_8, "<ActuatorCylinder>");
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
ac_free(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    int i = 0;

    ac->passed_msg = NULL;
    ac->passed_ierr = NULL;
    ac->env = NULL;
    FREE_OBJ(ac->A);
    FREE_OBJ(ac->theta);   
    FREE_OBJ(ac->z);
    FREE_OBJ(ac->fr);
    FREE_OBJ(ac->ft);
    FREE_OBJ(ac->fz);

    if (ac->work) {
        success = ac_work_free(ac->work, msg, ierr); CHECKERRK(WARNING_2, "<ACWork>");                
    }
    FREE_OBJ(ac->work);

    for (i=0 ; i<ac->nb ; i++) {        
        success = blade_free(&ac->blade[i], msg, ierr); CHECKERRK(WARNING_2, "<Blade>");
    }
    FREE_OBJ(ac->blade);    

    if (ac->perf) {
        success = ac_perf_free(ac->perf, msg, ierr); CHECKERRK(WARNING_2, "<ACPerformance>");
    }
    FREE_OBJ(ac->perf);

    if (ac->pc) {
        success = powell_free(ac->pc, msg, ierr); CHECKERRK(WARNING_2, "<PowellCoeff>");
    }
    FREE_OBJ(ac->pc);

    return SAFE;
CLEAN_UP:
    return WARNING; 
}


ERROR_CODE
ac_work_free(struct ACWork* work, char* msg, ERROR_CODE* ierr)
{
    FREE_OBJ(work->w);
    FREE_OBJ(work->wx);
    FREE_OBJ(work->wy);
    FREE_OBJ(work->Vt);
    FREE_OBJ(work->Vn);
    FREE_OBJ(work->cn);
    FREE_OBJ(work->ct);
    FREE_OBJ(work->Pn);
    FREE_OBJ(work->integrand);
    FREE_OBJ(work->cl);
    FREE_OBJ(work->cd);
    FREE_OBJ(work->residual);
    FREE_OBJ(work->W);
    FREE_OBJ(work->phi);
    FREE_OBJ(work->alpha);
    FREE_OBJ(work->Re);
    return SAFE;
}


ERROR_CODE
ac_perf_free(struct ACPerformance* perf, char* msg, ERROR_CODE* ierr)
{
    FREE_OBJ(perf->fx);
    FREE_OBJ(perf->fy);
    FREE_OBJ(perf->fz);
    FREE_OBJ(perf->Q);    
    return SAFE;
}


ERROR_CODE 
ac_allocate_A_matrix(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;
    const int N = ac->nt;
    const int M = 2*(ac->nt);
    ac->A = malloc(M*N*sizeof(double)); CHECK_MALLOC(ac->A, FATAL_8, "A");
    for(i=0 ; i<M ; i++) {
        for(j=0 ; j<N ; j++) {            
            ac->A[ij(i,j,N)] = 0.0;
        }
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
ac_allocate_theta(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    const int nt = ac->nt;
    const double dtheta = 2.0*PI/nt; 

    ac->theta = malloc(nt*sizeof(double*)); CHECK_MALLOC(ac->theta, FATAL_8, "the");
    for(i=0 ; i<nt ; i++) {        
        ac->theta[i] = dtheta/2.0 + i*dtheta;        
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
ac_set_initial_theta_offset(struct ActuatorCylinder* ac, const double the0, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;    
    if (!ac) {
        CALL_CHECKERRK(FATAL_23,"Attempting to set intial theta offset. Actuator Cylinder must be created first");
    }
    if (the0<-PI || PI<the0) {
        CALL_CHECKERRK(FATAL_23,"Initial azimuth offset must be between -PI < azimuth < PI radians");
    }
    ac->the0 = the0;
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
ac_set_angular_velocity(struct ActuatorCylinder* ac, const double omega, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    if (!ac) {
        CALL_CHECKERRK(FATAL_23,"Attempting to set VAWT angular velocity. Actuator Cylinder must be created first");
    }
    ac->omega = omega;
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
ac_allocate_prerequisites_for_A_matrix(const int h, const int n, double** value, double*** Rx, double*** Ry, double*** Wx, char* msg, ERROR_CODE* ierr) 
{
    int i = 0;
    int j = 0;
    const double dtheta = 2.0*PI/n; 
    *Rx = malloc(n*sizeof(double*)); CHECK_MALLOC(*Rx, FATAL_8, "Rx");
    *Ry = malloc(n*sizeof(double*)); CHECK_MALLOC(*Ry, FATAL_8, "Ry");
    *Wx = malloc(n*sizeof(double*)); CHECK_MALLOC(*Wx, FATAL_8, "Wx");

    *value = malloc(h*sizeof(double)); CHECK_MALLOC(*value, FATAL_8, "value");
    for(i=0 ; i<h ; i++) {
        (*value)[i] = 0.0;
    }

    for(i=0 ; i<n ; i++) {
        (*Rx)[i] = NULL;
        (*Ry)[i] = NULL;
        (*Rx)[i] = malloc(n*sizeof(double)); CHECK_MALLOC((*Rx)[i], FATAL_8, "Rx[i]");
        (*Ry)[i] = malloc(n*sizeof(double)); CHECK_MALLOC((*Ry)[i], FATAL_8, "Ry[i]");
        (*Wx)[i] = malloc(n*sizeof(double)); CHECK_MALLOC((*Wx)[i], FATAL_8, "Wx[i]");
        for(j=0 ; j<n ; j++) {
            (*Rx)[i][j] = dtheta/2.0;
            (*Ry)[i][j] = 0.0;
            (*Wx)[i][j] = 0.0;
        }
    }

    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
ac_set_swept_area(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr){
    ERROR_CODE success = SAFE;
    double sum = 0.0;
    double* r = NULL;
    int i = 0;
    const int nz = ac->nz;

    r = malloc(nz*sizeof(double)); CHECK_MALLOC(r, FATAL_8, "<r>");
    for (i=0 ; i<nz ; i++) {
        r[i] = ac->blade[0].elem[i].r;
    }
    /* @todo: should check the swept area for all blade for consistency */    
    integrate(&sum, r, ac->z, ac->nz, msg, ierr); CHECKERRK(FATAL_27, "Could not calculate swept area.");
    ac->swept_area = 2*sum;
    FREE_OBJ(r);
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
ac_delete_prerequisites_for_A_matrix(const int n, double** value, double*** Rx, double*** Ry, double*** Wx, char* msg, ERROR_CODE* ierr) 
{
    int i = 0;
    for(i=0 ; i<n ; i++) {
        FREE_OBJ((*Rx)[i]);
        FREE_OBJ((*Wx)[i]);
        FREE_OBJ((*Ry)[i]);
    }
    FREE_OBJ(*Rx);
    FREE_OBJ(*Wx);
    FREE_OBJ(*Ry);    
    FREE_OBJ(*value);
    return SAFE;
}


ERROR_CODE 
ac_generate_memory(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    const double nt = ac->nt; 
    const int nz = ac->nz;
    int i = 0;
    int j = 0;    
    ac->work = malloc(sizeof(struct ACWork)); CHECK_MALLOC(ac->work, FATAL_8, "Actuatory Cylinder <ACWork>");
    ac->work->Vt = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->Vt, FATAL_8, "Actuatory Cylinder <Vt>");
    ac->work->Vn = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->Vn, FATAL_8, "Actuatory Cylinder <Vn>");
    ac->work->w = malloc(2*nt*sizeof(double)); CHECK_MALLOC(ac->work->w, FATAL_8, "Actuatory Cylinder <wx>");

    for(j=0 ; j<2*nt ; j++) {                
        ac->work->w[j] = 0.0;
    }

    ac->work->wx = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->wx, FATAL_8, "Actuatory Cylinder <wx>");
    ac->work->wy = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->wy, FATAL_8, "Actuatory Cylinder <wy>");
    ac->work->ct = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->ct, FATAL_8, "Actuatory Cylinder <ct>");
    ac->work->cn = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->cn, FATAL_8, "Actuatory Cylinder <cn>");
    ac->work->Pn = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->Pn, FATAL_8, "Actuatory Cylinder <Pn>");
    ac->work->integrand = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->integrand, FATAL_8, "Actuatory Cylinder <integrand>");
    ac->work->cl = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->cl, FATAL_8, "Actuatory Cylinder <cl>");
    ac->work->cd = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->cd, FATAL_8, "Actuatory Cylinder <cd>");
    ac->work->residual = malloc(2*nt*sizeof(double)); CHECK_MALLOC(ac->work->residual, FATAL_8, "Actuatory Cylinder <residual>");
    ac->work->W = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->W, FATAL_8, "Actuatory Cylinder <W>");
    ac->work->phi = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->phi, FATAL_8, "Actuatory Cylinder <phi>");
    ac->work->alpha = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->alpha, FATAL_8, "Actuatory Cylinder <phi>");
    ac->work->Re = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->work->Re, FATAL_8, "Actuatory Cylinder <Re>");
    ac->fr = malloc(nt*nz*sizeof(double)); CHECK_MALLOC(ac->fr, FATAL_8, "Actuatory Cylinder <fr>"); 
    ac->ft = malloc(nt*nz*sizeof(double)); CHECK_MALLOC(ac->ft, FATAL_8, "Actuatory Cylinder <ft>"); 
    ac->fz = malloc(nt*nz*sizeof(double)); CHECK_MALLOC(ac->fz, FATAL_8, "Actuatory Cylinder <fz>"); 
    // ac->Q = malloc(nt*sizeof(double)); CHECK_MALLOC(ac->Q, FATAL_8, "Actuatory Cylinder <Q>");
    for(i=0 ; i<nz ; i++) {
        for(j=0 ; j<nt ; j++) {
            ac->fr[ij(i,j,nt)] = 0.0;
            ac->ft[ij(i,j,nt)] = 0.0;
            ac->fz[ij(i,j,nt)] = 0.0;            
        }
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
ac_allocate_perf(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    const int n = ac->nt;
    ac->perf = malloc(sizeof(struct ACPerformance)); CHECK_MALLOC(ac->perf, FATAL_8, "Actuator Cylinder <ACPerformance>");
    success = ac_perf_initialize(ac, msg, ierr); CHECKERRK(FATAL_5,"Performance coefficients");
    ac->perf->fx = malloc(n*sizeof(double)); CHECK_MALLOC(ac->perf->fx, FATAL_8, "ACPerformance <fx>");       
    ac->perf->fy = malloc(n*sizeof(double)); CHECK_MALLOC(ac->perf->fy, FATAL_8, "ACPerformance <fy>");       
    ac->perf->fz = malloc(n*sizeof(double)); CHECK_MALLOC(ac->perf->fz, FATAL_8, "ACPerformance <fz>");       
    ac->perf->Q = malloc(n*sizeof(double)); CHECK_MALLOC(ac->perf->Q, FATAL_8, "ACPerformance <Q>");           
    return SAFE;
CLEAN_UP:
    return FATAL;
}

ERROR_CODE 
ac_generate_matrix(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    const int n = ac->nt;
    const double dtheta = 2.0*PI/n; 
    double** Rx = NULL;
    double** Ry = NULL;
    double** Wx = NULL;
    double* value = NULL;
    double den = 0.0;
    double v1 = 0.0;
    double v2 = 0.0;
    double phi = 0.0;
    double x = 0.0;
    double y = 0.0;
    double dx = 0.0;
    const int h = 100; /* number of grid spacing for integration */
    int i = 0;
    int j = 0;
    int k = 0;
    const int rows = (ac->nt);
    
    success = ac_allocate_A_matrix(ac, msg, ierr);
    success = ac_allocate_theta(ac, msg, ierr);
    success = ac_allocate_prerequisites_for_A_matrix(h, n, &value, &Rx, &Ry, &Wx, msg, ierr); CHECKERRK(FATAL_8,"Prerequisites");

    for(i=0 ; i<n ; i++) {
        phi = 0.0;
        x = -sin(ac->theta[i]);
        y = cos(ac->theta[i]);
        
        for(j=0 ; j<n ; j++) {
            phi = ac->theta[j] - dtheta/2.0;
            dx = ((ac->theta[j]+dtheta/2.0) - (ac->theta[j]-dtheta/2.0))/(double)(h-1);
            for(k=0 ; k<h ; k++) {
                v1 = x + sin(phi);                
                v2 = y - cos(phi);
                den = (pow(v1,2.0) + pow(v2,2.0)); 
                if (den<1E-6) {
                    value[k] = 0.0;
                } else {
                    value[k] = (v1*cos(phi) + v2*sin(phi))/den;
                }
                phi+=dx;
            }
            for(k=0 ; k<h-1 ; k++) {
                Ry[i][j] += (dx/2.0)*(value[k] + value[k+1]);
            }
        }
        
        if(i<n/2) {
            Rx[i][i] = PI*(1.0/(double)n - 1.0);
        } else {
            Rx[i][i] = PI*(1.0/(double)n + 1.0);
        }        
    }
    
    
    for(i=n/2 ; i<n ; i++) {
        Wx[i][n-(i+1)] = -1;
    }
    
    for(i=0 ; i<n ; i++) {
        for(j=0 ; j<n ; j++) {
            ac->A[ij(i,j,rows)] = Rx[i][j]/(2.0*PI) + Wx[i][j];
            ac->A[ij(i+n,j,rows)] = Ry[i][j]/(2.0*PI);            
        }
    }
    success = ac_delete_prerequisites_for_A_matrix(n, &value, &Rx, &Ry, &Wx, msg, ierr);
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
ac_set_z(struct ActuatorCylinder* ac, double* z, const int nz, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    ac->nz = nz;
    ac->z = malloc(nz*sizeof(double*)); CHECK_MALLOC(ac->z, FATAL_8, "z");
    for(i=0 ; i<nz ; i++) {
        ac->z[i] = z[i];
    }
    return SAFE;
CLEAN_UP:
    return FATAL;    
}


ERROR_CODE
ac_set_n_theta(struct ActuatorCylinder* ac, const int nt, char* msg, ERROR_CODE* ierr)
{
    ac->nt = nt;
    return SAFE;
}


ERROR_CODE
ac_set_airfoil_twist(struct ActuatorCylinder* ac, const double twist, const int ib, const int ie, char* msg, ERROR_CODE* ierr)
{
    const int n = ac->nz;
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    if (ib<0 || (ac->nb-1)<ib) {
        b = bformat("<%d> is not a valid blade number", ib);     
        CALL_CHECKERRK(FATAL_9, (const char*)b->data);        
    }
    if (n<=ie) {
        CALL_CHECKERRK(FATAL_14,"Number of airfoils must assigned to ActuatorCylinder must be the same length as spatial discretization.");
    }
    ac->blade[ib].elem[ie].twist = twist;
    return SAFE;
CLEAN_UP:
    bdestroy(b);
    return FATAL;        
}


ERROR_CODE
ac_set_airfoil_swept_offset(struct ActuatorCylinder* ac, const double swept, const int ib, const int ie, char* msg, ERROR_CODE* ierr)
{
    const int n = ac->nz;
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    if (ib<0 || (ac->nb-1)<ib) {
        b = bformat("<%d> is not a valid blade number", ib);     
        CALL_CHECKERRK(FATAL_9, (const char*)b->data);        
    }
    if (n<=ie) {
        CALL_CHECKERRK(FATAL_14,"Number of airfoils must assigned to ActuatorCylinder must be the same length as spatial discretization.");
    }
    ac->blade[ib].elem[ie].swept_offset = swept;
    return SAFE;
CLEAN_UP:
    bdestroy(b);
    return FATAL;        
}


ERROR_CODE
ac_set_airfoil_radius(struct ActuatorCylinder* ac, const double r, const int ib, const int ie, char* msg, ERROR_CODE* ierr)
{
    const int n = ac->nz;
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    if (ib<0 || (ac->nb-1)<ib) {
        b = bformat("<%d> is not a valid blade number", ib);     
        CALL_CHECKERRK(FATAL_9, (const char*)b->data);        
    }
    if (n<=ie) {
        CALL_CHECKERRK(FATAL_14,"Number of airfoils must assigned to ActuatorCylinder must be the same length as spatial discretization.");
    }
    ac->blade[ib].elem[ie].r = r;
    return SAFE;
CLEAN_UP:
    bdestroy(b);
    return FATAL;        
}


ERROR_CODE
ac_set_airfoil_geometry(struct ActuatorCylinder* ac, struct Airfoil* af, const int ib, const int ie, char* msg, ERROR_CODE* ierr)
{
    const int n = ac->nz;
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    if (ib<0 || (ac->nb-1)<ib) {
        b = bformat("<%d> is not a valid blade number", ib);     
        CALL_CHECKERRK(FATAL_9, (const char*)b->data);        
    }
    if (n<=ie) {
        CALL_CHECKERRK(FATAL_14,"Number of airfoils must assigned to ActuatorCylinder must be the same length as spatial discretization.");
    }
    ac->blade[ib].elem[ie].af = af;
    return SAFE;
CLEAN_UP:
    bdestroy(b);
    return FATAL;        
}


ERROR_CODE
ac_set_airfoil_chord(struct ActuatorCylinder* ac, const double chord, const int ib, const int ie, char* msg, ERROR_CODE* ierr)
{
    const int n = ac->nz;
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    if (ib<0 || (ac->nb-1)<ib) {
        b = bformat("<%d> is not a valid blade number", ib);     
        CALL_CHECKERRK(FATAL_9, (const char*)b->data);        
    }
    if (n<=ie) {
        CALL_CHECKERRK(FATAL_14,"Number of airfoils must assigned to ActuatorCylinder must be the same length as spatial discretization.");
    }
    ac->blade[ib].elem[ie].chord = chord;
    return SAFE;
CLEAN_UP:
    bdestroy(b);
    return FATAL;        
}


ERROR_CODE
ac_power(struct ActuatorCylinder* ac, struct Blade* blade, /*const double theta, */char* msg, ERROR_CODE* ierr)
/* @todo: theta is not a used input argument */
{
    ERROR_CODE success = SAFE;
    const int nt = ac->nt;
    const int nz = ac->nz;
    int j = 0;
    int istart = 5;
    int jstart = 5;
    const double dtheta = (2*PI/nt)/2.0;

    for (j=0 ; j<nz ; j++) { 
        blade->Rb[j] = 0.0;
        blade->Tb[j] = 0.0;
        blade->Zb[j] = 0.0;
        // /* @change */
        // the = theta + blade->elem[j].swept_offset;
        // if (the<0.0) {
        //     the+=2*PI;
        // } else if (2*PI<the) {
        //     the-=2*PI;
        // }                
        /* @todo: this 2D spline can produce inconsistent results. polynomial interpolation would be better. */
        // success = linear_2D_spline(&blade->Rb[j], ac->fr, &istart, ac->z[j], ac->z, ac->nz, &jstart, the, ac->theta, ac->nt, msg, ierr); CHECKERRK(FATAL_27, "2D spline iterating on <Rb>");
        // success = linear_2D_spline(&blade->Tb[j], ac->ft, &istart, ac->z[j], ac->z, ac->nz, &jstart, the, ac->theta, ac->nt, msg, ierr); CHECKERRK(FATAL_27, "2D spline iterating on <Tb>");
        // success = linear_2D_spline(&blade->Zb[j], ac->fz, &istart, ac->z[j], ac->z, ac->nz, &jstart, the, ac->theta, ac->nt, msg, ierr); CHECKERRK(FATAL_27, "2D spline iterating on <Zb>");
        success = linear_2D_spline(&blade->Rb[j], ac->fr, &istart, ac->z[j], ac->z, ac->nz, &jstart, blade->elem[j].the, ac->theta, ac->nt, msg, ierr); CHECKERRK(FATAL_27, "2D spline iterating on <Rb>");
        success = linear_2D_spline(&blade->Tb[j], ac->ft, &istart, ac->z[j], ac->z, ac->nz, &jstart, blade->elem[j].the, ac->theta, ac->nt, msg, ierr); CHECKERRK(FATAL_27, "2D spline iterating on <Tb>");
        success = linear_2D_spline(&blade->Zb[j], ac->fz, &istart, ac->z[j], ac->z, ac->nz, &jstart, blade->elem[j].the, ac->theta, ac->nt, msg, ierr); CHECKERRK(FATAL_27, "2D spline iterating on <Zb>");
        blade->Ri[j] = blade->Rb[j]*cos(dtheta) - blade->Tb[j]*sin(dtheta); /* @todo: check */
        blade->Ti[j] = blade->Rb[j]*sin(dtheta) + blade->Tb[j]*cos(dtheta); /* @todo: check */
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
ac_force(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    const int nz = ac->nz;
    const int nt = ac->nt;    
    int i = 0;
    int j = 0;
    double cn = 0.0;
    double ct = 0.0;
    double q = 0.0;
    double r = 0.0;
    double rmax = 0.0;
    struct Environment* env = ac->env;
    struct PowellCoeff* pc = NULL;

    for(i=0 ; i<nz ; i++) {
        r = ac->blade[0].elem[i].r;
        rmax = ac->blade[0].rmax;
        if (fabs(r/rmax)<1E-3) { /* @todo: need to set r_max, this is solved for a single blade, needs to change */
            for(j=0 ; j<nt ; j++) {
                ac->fr[ij(i,j,nt)] = 0.0;
                ac->ft[ij(i,j,nt)] = 0.0;
                ac->fz[ij(i,j,nt)] = 0.0;
            } 
        } else {
            ac->index = i;
            ac->passed_msg = msg;
            ac->passed_ierr = ierr;

            if (ac->pc==NULL) {
                success = powell_allocate(ac, nt, msg, ierr); CHECKERRK(FATAL_26,"");
            }
            pc = ac->pc;
            pc->info = __cminpack_func__(hybrd)(ac_func,
                                                ac,
                                                pc->n,
                                                ac->work->w,
                                                ac->work->residual,
                                                pc->xtol,
                                                pc->maxfev,
                                                pc->ml,
                                                pc->mu,
                                                pc->epsfcn,
                                                pc->diag,
                                                pc->mode,
                                                pc->factor,
                                                pc->nprint,
                                                &pc->nfev,
                                                pc->fjac,
                                                pc->ldfjac,
                                                pc->r,
                                                pc->lr,
                                                pc->qtf,
                                                pc->wa1,
                                                pc->wa2,
                                                pc->wa3,
                                                pc->wa4);
            
            ac_func(ac, pc->n, ac->work->w, ac->work->residual, -9);
            if(ac->dstall) {
                for (j=0 ; j<ac->nb ; j++) {
                    dstall_reset(ac->blade[j].elem[i].d);
                }
            }
            ac->passed_msg = NULL;
            ac->passed_ierr = NULL;
            const double chord = ac->blade[0].elem[i].chord; /* @todo: this is solved for a single blade, needs to change */
            const double delta = ac->blade[0].elem[i].delta; /* @todo: this is solved for a single blade, needs to change */
            for(j=0 ; j<nt ; j++) {                
                cn = ac->work->cn[j];
                ct = ac->work->ct[j];
                q = 0.5*(env->rho)*pow(ac->work->W[j],2);
                ac->fr[ij(i,j,nt)] = -cn*q*(chord);
                ac->ft[ij(i,j,nt)] = ct*q*(chord)/cos(delta);
                ac->fz[ij(i,j,nt)] = -cn*q*(chord)*tan(delta);
            }
        }
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


int
ac_func(void* ac_ptr, int n, const double* w, double* fvec, int iflag)
{
    struct ActuatorCylinder* ac = (struct ActuatorCylinder*)ac_ptr;
    struct Environment* env = ac->env;
    ERROR_CODE success = SAFE;
    const int index = ac->index;
    const int nt = ac->nt;
    int i = 0;
    char* msg = NULL;
    ERROR_CODE* ierr = NULL;
    double value = -999.9;
    
    msg = ac->passed_msg;
    ierr = ac->passed_ierr;
    if (iflag==0) {
        /* @todo: insert print statements here when nprint is positive */
        return 0;
    }
    success = ac_velocity(ac, &ac->blade[0].elem[index], index, env->V[index], w, msg, ierr); CHECKERRK(FATAL_26, "<ac_velocity>");
    if (ac->dstall && iflag==-9) { /* Beddoes dynamic stall model */
        success = dstall_get_coefficients(&ac->blade[0].elem[index], ac, msg, ierr); CHECKERRK(FATAL_26, "<dstall_get_coefficients>");
    } else { /* cd/cl look up table */
        for (i=0 ; i<nt ; i++) {
            success = af_get_cl(&value, ac->blade[0].elem[index].af, ac->work->alpha[i], ac->work->Re[i], msg, ierr); CHECKERRK(FATAL_26, "<af_get_cl>");
            ac->work->cl[i] = value;        
            success = af_get_cd(&value, ac->blade[0].elem[index].af, ac->work->alpha[i], ac->work->Re[i], msg, ierr); CHECKERRK(FATAL_26, "<af_get_cd>");
            ac->work->cd[i] = value;
            // printf("(DS=FALSE) %f  %f  %f\n",RAD2DEG*(ac->theta[i]), ac->work->cl[i], ac->work->cd[i]);
        }
    }   
    success = ac_error(ac, &ac->blade[0].elem[index], w, fvec, index, env->V[index], msg, ierr); CHECKERRK(FATAL_26, "ac_error");
    return 0;
    /* @todo: need to return proper errors */
CLEAN_UP:
    return -1;
}


ERROR_CODE 
ac_error(struct ActuatorCylinder* ac, struct BladeElement* elem, const double* w, double* fvec, const int index, const double Uinf, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;
    const int nt = ac->nt;
    double a = 0.0;
    double ka = 0.0;
    double c_thrust = 0.0;
    double* W = ac->work->W;
    
    double* integrand = ac->work->integrand;
    double temp = 0.0;
    const double chord = elem->chord;
    const double r = elem->r;
    const double delta = elem->delta;
    const double sigma = ac->nb*chord/(2*r);
    const double* theta = ac->theta;
    const double* cl = ac->work->cl;
    const double* cd = ac->work->cd;
    const int rows = ac->nt;
    double phi = 0.0;

    for (i=0 ; i<nt ; i++) {
        phi = ac->work->phi[i];
        ac->work->cn[i] = cl[i]*cos(phi) + cd[i]*sin(phi);
        ac->work->ct[i] = cl[i]*sin(phi) - cd[i]*cos(phi);
        ac->work->Pn[i] = sigma/(2*PI)*(ac->work->cn[i])*pow((W[i]/Uinf),2.0);
        integrand[i] = pow((W[i]/Uinf),2)*(ac->work->cn[i]*sin(theta[i]) - ac->work->ct[i]*cos(theta[i])/cos(delta));
    }
    c_thrust = (sigma/(2*PI))*ac_periodic_integrand(nt, integrand, theta, msg, ierr);
    if (2.0<c_thrust) {
        a = 0.5*(1.0 + sqrt(1.0 + c_thrust));
        ka = 1.0/(a-1.0);
    } else if (0.96<c_thrust) {
        a = (1.0/7.0)*(1 + 3*sqrt(7.0/2.0*c_thrust - 3.0));
        ka = 18.0*a / (7.0*a*a - 2.0*a + 4.0);
    } else {
        a = 0.5*(1.0 - sqrt(1.0 - c_thrust));
        ka = 1.0/(1.0-a);
    }
    for (i=0 ; i<2*nt ; i++) {
        temp = 0.0;
        for (j=0 ; j<nt ; j++) {
            temp += ac->A[ij(i,j,rows)]*ac->work->Pn[j];
        }
        fvec[i] = temp*ka- w[i];
    }
    return SAFE;
}


ERROR_CODE 
ac_velocity(struct ActuatorCylinder* ac, struct BladeElement* elem, const int index, const double Uinf, const double* w, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    const int nt = ac->nt;
    const double* theta = ac->theta;

    /* work variables */
    double* Vt = ac->work->Vt;
    double* Vn = ac->work->Vn;
    double* wx = ac->work->wx; /* should be an iterated input varaible */
    double* wy = ac->work->wy; /* should be an iterated input varaible */

    double* phi = ac->work->phi;
    double* alpha = ac->work->alpha; /* @todo: this should be moved to element */
    double* W = ac->work->W;
    double* Re = ac->work->Re;

    const double omega = ac->omega;
    const double twist = elem->twist;
    const double r = elem->r;
    struct Environment* env = ac->env;
    const double chord = elem->chord;
    const double mu = env->mu;
    const double rho = env->rho;   

    for (i=0 ; i<nt ; i++) { 
        wx[i] = w[i];
        wy[i] = w[i+nt];
    }

    for (i=0 ; i<nt ; i++) {
        Vt[i] = omega*r + Uinf*(1 + wx[i])*cos(theta[i]) + Uinf*wy[i]*sin(theta[i]);        
        Vn[i] = Uinf*(1+wx[i])*sin(theta[i]) - Uinf*wy[i]*cos(theta[i]);
        W[i] = sqrt(pow(Vt[i],2) + pow(Vn[i],2));
        phi[i] = atan2(Vn[i], Vt[i]);
        alpha[i] = phi[i] - twist;
        Re[i] = rho*chord*W[i]/mu; /* @todo? */
    }
    return SAFE;
}


double
ac_periodic_integrand(const int n, const double* fx, const double* theta, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    double out = 0.0;
    double d_theta = 0.0;
    
    for (i=1 ; i<n-1 ; i++) { 
        out += (theta[i]-theta[i-1])*(fx[i] + fx[i-1])/2.0;
    }
    d_theta = 2.0*theta[0];
    out += (d_theta*(fx[0] + fx[n-1])/2.0);
    return out;
}
