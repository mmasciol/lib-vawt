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


#include <stdlib.h>
#include <math.h>
#include "dstall.h"
#include "sys.h"
#include "af.h"
#include "blade.h"
#include "ac.h"
#include "numerics.h"
#include "./simclist/simclist.h"
#include "./bstring/bstrlib.h"


int
dstall_sign(const double x)
{
    return (x>0) - (x<0);
}


ERROR_CODE
dstall_initialize(struct Blade* blade, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    const int nz = blade->m;
    int i = 0;
    int na = -999;
    for (i=0 ; i<nz ; i++) {
        if (!blade->elem[i].d) {
            blade->elem[i].d = malloc(sizeof(struct DStall)); CHECK_MALLOC(blade->elem[i].d, FATAL_8, "<DStall>");
            na = blade->elem[i].af->na;
            blade->elem[i].d->na = na;
            blade->elem[i].d->fn = NULL;   
            blade->elem[i].d->ft = NULL;
            blade->elem[i].d->fn = malloc(na*sizeof(double)); CHECK_MALLOC(blade->elem[i].d->fn, FATAL_8, "DStall <fn>");
            blade->elem[i].d->ft = malloc(na*sizeof(double)); CHECK_MALLOC(blade->elem[i].d->ft, FATAL_8, "DStall <ft>");
            dstall_reset(blade->elem[i].d);
            success = dstall_initialize_fn(blade->elem[i].d, blade->elem[i].af, msg, ierr); CHECKERRK(FATAL_33, "");
            success = dstall_initialize_ft(blade->elem[i].d, blade->elem[i].af, msg, ierr); CHECKERRK(FATAL_33, "");
        }
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
dstall_reset(struct DStall* d)
{
    int i = 0;
    d->flow_shift = false;
    d->vortex = false;
    d->separation = false;
    d->super_sonic = false;
    d->first_entry = true;
    d->a1 = 0.3;  /* Leishman 2011 constant */
    d->a2 = 0.7;  /* Leishman 2011 constant */
    d->b1 = 0.14; /* Leishman 2011 constant */
    d->b2 = 0.53; /* Leishman 2011 constant */
    d->alpha_e = 0.0;
    d->dt = 0.0;
    d->tv = 6.0;
    d->tf = 3.0;
    d->c = 335.0;      /* speed of sound, m/sec */
    d->tp = 1.7;
    d->tv_l = 11.0;
    d->Cn_cp = 0.0;
    d->CLA = 0.0;
    d->CDA = 0.0;
    d->CMA = 0.0;
    d->alpha_n = 0.0;
    d->ds = 0.0;
    d->fk = 0.0;
    d->Cni_q = 0.0;
    d->fn_pp = 0.0;
    d->ft_pp = 0.0;   
    d->Cmi = 0.0;
    d->Cmq = 0.0;    
    d->ct = 0.0;
    d->cn = 0.0;      
    for (i=0 ; i<2 ; i++) {
        d->alpha[i] = 0.0;
        d->Cn_p[i] = 0.0;   
        d->Cn_pot[i] = 0.0;
        d->alpha_f[i] = 0.0;
        d->Xn[i] = 0.0;
        d->Yn[i] = 0.0;
        d->tau[i] = 0.0;
        d->adot[i] = 0.0;
        d->dn[i] = 0.0;
        d->dq[i] = 0.0;
        d->qx[i] = 0.0;
        d->dCn_p[i] = 0.0;
        d->Df_n[i] = 0.0;
        d->Dp[i] = 0.0;
        d->Cv[i] = 0.0;
        d->Cn_v[i] = 0.0;
        d->dqp[i] = 0.0;
        d->Df_t[i] = 0.0;
        d->ft_p[i] = 0.0;
        d->fn_p[i] = 0.0;
        d->Df_alpha_f_e[i] = 0.0;
    }
    return SAFE;
}


ERROR_CODE
dstall_initialize_fn(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    const int na = d->na;
    const double c_na = af->cna;
    double alpha = -999.9;
    double cl = -999.9;
    double cd = -999.9;
    int i = 0;    
    for (i=0 ; i<na ; i++){
        if (af->na-1<i) {
            CALL_CHECKERRK(FATAL_28, "Error setting DStall <fn>");
        }
        alpha = af->alpha[i];
        cl = af->cl[i];
        cd = af->cd[i];
        d->cn = cl*cos(alpha) + (cd - af->cd0)*sin(alpha);
        if (1.0E-6<fabs(c_na)) { /* @todo: check for machine epsilon */
            d->fn[i] = d->cn/c_na;
        } else {
            d->fn[i] = 1.0;
        }
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
dstall_initialize_ft(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    const int na = d->na;
    const double c_na = af->cna;
    double alpha = -999.9;
    double cl = -999.9;
    double cd = -999.9;
    int i = 0;
    for (i=0 ; i<na ; i++){
        if (af->na-1<i) {
            CALL_CHECKERRK(FATAL_28, "Error setting DStall <ft>");
        }
        
        alpha = af->alpha[i];
        cl = af->cl[i];
        cd = af->cd[i];
        d->ct = cl*sin(alpha) - (cd - af->cd0)*cos(alpha); 
        if (1.0E-6<fabs(c_na)) {
            d->ft[i] = d->ct/c_na;
        } else {
            d->ft[i] = 1.0;
        }
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
dstall_free(struct DStall* d)
{
    FREE_OBJ(d->fn);
    FREE_OBJ(d->ft);  
    return SAFE;
}


void
bed_update(struct DStall* d)
{
    d->alpha[1] = d->alpha[0];
    d->Df_alpha_f_e[1] = d->Df_alpha_f_e[0];
    d->adot[1] = d->adot[0];
    d->Xn[1] = d->Xn[0];
    d->Yn[1] = d->Yn[0];
    d->Cn_pot[1] = d->Cn_pot[0];
    d->Dp[1] = d->Dp[0];
    d->fn_p[1] = d->fn_p[0];
    d->ft_p[1] = d->ft_p[0];
    d->tau[1] = d->tau[0];
    d->Df_n[1] = d->Df_n[0];
    d->Df_t[1] = d->Df_t[0];
    d->dn[1] = d->dn[0];
    d->Cn_v[1] = d->Cn_v[0];
    d->Cv[1] = d->Cv[0];
    d->Cn_p[1] = d->Cn_p[0];
    d->dCn_p[1] = d->dCn_p[0];
    d->qx[1] = d->qx[0];
    d->dq[1] = d->dq[0];
    d->alpha_f[1] = d->alpha_f[0];
    d->dqp[1] = d->dqp[0];
    return;
}


ERROR_CODE
dstall_get_coefficients(struct BladeElement* elem, struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    int i = 0;
    const int nt = ac->nt;    
    const double* alpha = ac->work->alpha;
    const double* Re = ac->work->Re;
    double coeff;
    double cla = 0.0;
    double cda = 0.0;
    const double dtheta = (ac->theta[1] - ac->theta[0]);
    const double dt = dtheta/ac->omega;
    
    for (i=0 ; i<nt ; i++) {
        if (elem->d->first_entry==true) {
            elem->d->first_entry = false; /* @todo: this might need to be reset for performance analysis */
            success = dstall_beddoes_init(elem, dt, elem->alpha, msg, ierr); /* @todo: error */
            success = af_get_cl(&coeff, elem->af, alpha[i], Re[i], msg, ierr); CHECKERRK(FATAL_26, "Using cd/cl look-up table for the first dynamic stall calculation <af_get_cl>");
            ac->work->cl[i] = coeff;                   
            success = af_get_cd(&coeff, elem->af, alpha[i], Re[i], msg, ierr); CHECKERRK(FATAL_26, "Using cd/cl look-up table for the first dynamic stall calculation <af_get_cd>");
            ac->work->cd[i] = coeff;
            bed_update(elem->d);
        } else {        
            success = dstall_beddoes(elem, alpha[i], ac->work->W[i], &cla, &cda, msg, ierr);
            ac->work->cl[i] = cla;
            ac->work->cd[i] = cda;
            bed_update(elem->d);
        }        
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
dstall_beddoes(struct BladeElement* elem, const double ALPHA, const double U, double* cla, double* cda, char* msg, ERROR_CODE* ierr)
{    
    /*
     * uses the Beddoes dynamic stall model
     *  the routine is entered with an angle of attack
     *  and returns CL and CD.
     * This routine is used regardless of whether the element
     *  is actually in dynamic stall state.
     * 
     * VARIABLES:
     *  W2    = Relative velocity squared over blade element
     *  J     = Index which identifies the blade element
     *  ALPHA = Angle of attack in radians
     *  CLA   = Lift coeff. which is calculated by the routine
     *  CDA   = Drag coeff. which is calculated by the routine
     *  CMA   = Moment coeff. which is calculated by the routine
     */
    ERROR_CODE success = SAFE;
    struct DStall* d = elem->d;
    struct Airfoil* af = elem->af;
    const double cd0 = af->cd0;
     
    /* @todo: need an option here to interpolate for numtiple airfoil tables 
     *        and find CNA as a function of Re using 2-D interpolation
     * Check to see if element has multiple airfoil tables, then interpolate values
     * of constants based on the current location.
     */
     
    d->alpha_n = ALPHA;    
    success = dstall_attach(d, af, elem->chord, U, msg, ierr); CHECKERRK(FATAL_35,"Attached flow");
    success = dstall_separate(d, af, msg, ierr); CHECKERRK(FATAL_35,"Separated flow");    
    success = dstall_vortex(d, msg, ierr); CHECKERRK(FATAL_35,"Vortex shedding");

    *cla = d->cn*cos(d->alpha_n) + d->ct*sin(d->alpha_n);
    *cda = d->cn*sin(d->alpha_n) - d->ct*cos(d->alpha_n) + cd0;
    // *CMA = af->PMC; /* @todo: this does not exist yet */
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
dstall_beddoes_init(struct BladeElement* elem, const double dt, double ALPHA, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    struct DStall* d = elem->d;
    struct Airfoil* af = elem->af;
    const double aol = af->aol;
    const double CNA1 = af->cna;;
    double coeff = 0.0;
    int idx = 0;
    const int na = d->na;
    
    d->alpha[0] = ALPHA;
    d->alpha_f[0] = ALPHA;
    d->dt = dt;
    
    /* @todo: need an option here to interpolate for numtiple airfoil tables 
     *        and find CNA as a function of Re using 2-D interpolation
     */
    d->Cn_pot[0] = CNA1*(ALPHA-aol);
    d->Cn_p[0] = d->Cn_pot[0];

    // *ALPHA = MIN(MAX(*ALPHA , af->alpha[0]), af->alpha[na-1]);
    idx = get_index(ALPHA, af->alpha, na);
    if (idx<0) {
        b = bformat("Could not find an angle of attack withint the airfoil range AoA = <%f>", ALPHA);
        CALL_CHECKERRK(FATAL_34,(const char*)b->data);        
    }

    if (idx==0) {
        coeff = 0.0;
    } else if (idx==na-1) {
        idx = idx-1;
        coeff = 1.0;
    } else {
        coeff = (af->alpha[idx] - ALPHA)/(af->alpha[idx] - af->alpha[idx+1]);
    }

    /* @todo: need an option here to interpolate for numtiple airfoil tables 
     *        and find CNA as a function of Re using 2-D interpolation
     */
    success = dstall_set_fn_p(d, idx, aol, coeff, msg, ierr);    
    success = dstall_set_ft_p(d, idx, aol, coeff, msg, ierr);
    success = dstall_set_Cv(d, msg, ierr);
    return SAFE;
CLEAN_UP:
    bdestroy(b);
    return FATAL;
}



ERROR_CODE
dstall_set_ft_p(struct DStall* d, const int idx, const double aol, const double coeff, char* msg, ERROR_CODE* ierr)
{
    double val = 0.0;
    int sign = 0;
    
    d->ft_p[0] = d->ft[idx] - (d->ft[idx] - d->ft[idx+1])*coeff;
    if (fabs(d->alpha_f[0]-aol)<1.0E-10) { /* @todo: 1E-10 is enough? */
        d->ft_p[0] = 1.0;
    } else {
        if (fabs(d->alpha_f[0])<1.0E-10) { /* @todo: should this be machine epsilon? */
            d->ft_p[0] = 1.0;
        } else {
            val = d->ft_p[0]/((d->alpha_f[0]-aol)*d->alpha_f[0]);            
            sign = dstall_sign(val);
            d->ft_p[0] = val*val*sign;
            if (d->ft_p[0]>1.0) {
                d->ft_p[0] = 1.0;
            }
            if (d->ft_p[0]<-1.0) {
                d->ft_p[0] = -1.0;
            }
        }
    }
    return SAFE;
}


ERROR_CODE
dstall_set_fn_p(struct DStall* d, const int idx, const double aol, const double coeff, char* msg, ERROR_CODE* ierr)
{
    double val = 0.0;
    int sign = 0;
    
    d->fn_p[0] = d->fn[idx] - (d->fn[idx] - d->fn[idx+1])*coeff;
    if (fabs(d->alpha_f[0]-aol)<1.0E-10) { /* @todo: 1E-10 is enough? */
        d->fn_p[0] = 1.0;
    } else {
        val = 2*sqrt(fabs(d->fn_p[0]/(d->alpha_f[0]-aol))) - 1;
        sign = dstall_sign(val);
        d->fn_p[0] = val*val*sign;
        if (d->fn_p[0]>1.0) {
            d->fn_p[0] = 1.0;
        }
        if (d->fn_p[0]<-1.0) {
            d->fn_p[0] = -1.0;
        }
    }
    return SAFE;
}


ERROR_CODE
dstall_set_Cv(struct DStall* d, char* msg, ERROR_CODE* ierr)
{
    int sign = 0;
    double gamma;
    sign = dstall_sign(d->fn_p[0]);
    gamma = sqrt(fabs(d->fn_p[0]))*sign + 1.0;
    d->fk = 0.25*gamma*gamma;
    d->Cv[0] = d->Cn_pot[0]*(1-d->fk);
    return SAFE;
}



ERROR_CODE
dstall_attach(struct DStall* d, struct Airfoil* af, const double chord, const double U, char* msg, ERROR_CODE* ierr)    
{
    ERROR_CODE success = SAFE;
    double beta = 0.0;
    double BS = 0.0;
    double Cni = 0.0;
    double Cnq = 0.0;
    double CO = 0.0;
    double da = 0.0;
    double PRP = 0.0;
    double X = 0.0;
    double k = 0.0;
    double mach = 0.0;
    double val = 0.0;
    double sign = 0.0;
    const double c = d->c;
    const double aol = af->aol;
    const double cna = af->cna;
    const double a1 = d->a1;
    const double a2 = d->a2;
    const double b1 = d->b1;
    const double b2 = d->b2;
    const double dt = d->dt;
    
    if (fabs(d->alpha_n)<=PI/2) { /*PIBY2 = PI/2 ?*/
        d->alpha[0] = d->alpha_n;
    } else if (d->alpha_n>PI/2) { /*PIBY2 = PI/2 ?*/
        d->alpha[0] = PI - d->alpha_n;
    } else {
        d->alpha[0] = -PI - d->alpha_n;
    }
    mach = U/(d->c);   /* mach number */
 
    /* Check to see that the element is not supersonic */
    if (!d->super_sonic && 1.0<=mach) {
        mach = 0.7;
        d->super_sonic = true;
        CALL_CHECKERRK(WARNING_3,"");
    } else if (d->super_sonic==true && mach<1.0) { 
        d->super_sonic = false;
        CALL_CHECKERRK(WARNING_3,"Super sonic conditions have ended");
    }

    // d->DT = CurrentTime - OldTime; /* @todo: this was removed from FAST 8, probably better to set this outside this function */
    beta = 1.0 - pow(mach,2);                                                             /* beta^2, look at Eq. 3.5 */
    d->ds = 2.0*dt*U/chord;                                                 /* Eq. 1.5b in Damiani */
    BS = beta*(d->ds);
    k = 0.75/((1-mach ) + PI*beta*pow(mach,2)*(a1*b1+a2*b2));                                     /* a1*b1 + a2*b2 = 0.413, Eq. 3.15 */
    X = (dt*c)/(chord*k);  /* note that Ti (Eq. 1.11) = chord/AS */      /* exponent in Eq. 3.21 */
    CO = k*chord/(c*mach);
 
    da = d->alpha[0] - d->alpha[1];                                                  /* change in angle of attack? */
    d->adot[0] = da/dt;                                                      /* Eq. 3.19 (d alpha)/dt = alpha dot (a time derivative) */ 
    PRP = (d->adot[0]*chord)/U;                                                /* Eq. 1.8 in Damiani, delta_alpha*chord/Vrel */
    PRP = dstall_saturation(PRP, 0.03, 0.1);
    d->adot[0] = PRP*U/chord;

    d->dn[0] = (d->dn[1])*exp(-X) + (d->adot[0] - d->adot[1])*exp(-0.5*X);    /* Eq. 3.21 (?) typo apparently (or Eq. 1.19 in Damiani) */
    Cni = 4*CO*(d->adot[0] - d->dn[0]);
    d->Cmi = -0.25*Cni;
    
    d->qx[0] = (d->adot[0] - d->adot[1])*chord/(U*dt);                     /* chord/(VREL*dt) */
    d->dq[0] = (d->dq[1])*exp(-X) + (d->qx[0] - d->qx[1])*exp(-0.5*X);
    Cnq = -CO*(d->qx[0] - d->dq[0]);                                              /* CO = XKA*chord/(d->AS)/XM; = ka*Ti/M */
    d->dqp[0] = d->dqp[1]*exp(-X/k) + (d->qx[0] - d->qx[1])*exp(-0.5*X/k);  /* Bottom Eq. 1.26 with mod (km,q is different), X/XKA = ((d->DT)*(d->AS)/chord/XKA)/(XKA) = dt/(Ti*ka*ka) */ 
    
    d->Cmq = -0.25*Cnq - (k*CO/3)*(d->qx[0] - d->dqp[0]);                       /* This Eq. 1.28 in Damiani, works out to   -(1/4)*CNQ-(ka*ka*Ti/(3*XM))*(d->QX-d->DQP) */
    val = Cni + Cnq;
    sign = dstall_sign(val);
    d->Cni_q = MIN(fabs(Cni+Cnq),1.0)*sign;
    
    d->Xn[0] = d->Xn[1]*exp(-0.14*BS) + 0.3*da*exp(-0.07*BS);    /* Eq. 3.11 */
    d->Yn[0] = d->Yn[1]*exp(-0.53*BS) + 0.7*da*exp(-0.265*BS);   /* Eq. 3.12 */
    
    d->alpha_e = d->alpha[0] - d->Yn[0] - d->Xn[0];                         /* Eq. 3.10 -- Eq. 1.15 in Damiani. Note that aol is subtracted in the next line */
    d->Cn_cp = cna*(d->alpha_e-aol);                                       /* Part of Eq. 1.48, 1.36 in Damiani , alpha_e = AE-AOL1 */
    d->Cn_pot[0] = d->Cn_cp + d->Cni_q;                                /* Eq. 3.23 */
    d->ct = d->Cn_pot[0]*(d->alpha_e);
    return SAFE;
CLEAN_UP:
    return FATAL;
}



double
dstall_saturation(const double X, const double VAL, const double SLOPE )
{
    double sat = 0.0;

    if (fabs(X)<=VAL) {
        sat = X;
    } else if (X>VAL) {
        sat = SLOPE*X + VAL*(1-SLOPE);
    } else {
        sat = SLOPE*X - VAL*(1-SLOPE);
    }
    return sat;
}


//    SUBROUTINE VORTEX( J, IBlade, AE )
ERROR_CODE
dstall_vortex(struct DStall* d, char* msg, ERROR_CODE* ierr)
{
    // double CMV = 0.0;
    double t_shedding = 0.0;
    double tv = d->tv;;
    
    d->Cv[0] = d->Cn_cp*(1-d->fk);                                                                               /* Eq. 1.47 in Damiani */

    if (d->tau[0]<1) {
        d->vortex = true;
        if (d->flow_shift) {
            d->vortex = false;
        }
    } else {
        d->vortex = false;
    }

    if (d->vortex) {
        d->Cn_v[0] = d->Cn_v[1]*exp(-(d->ds)/tv) + (d->Cv[0] - d->Cv[1])*exp(-0.5*(d->ds)/tv);  /* Eq. 1.45 in Damiani */
    } else {
        /* conditions for this to be met:  */
        d->Cn_v[0] = d->Cn_v[1]*exp(-(d->ds)/(tv*0.5));                                                    /* Eq. 1.48 in Damiani */
    }
    
    d->cn += d->Cn_v[0];                                                                              /* is this Eq. 1.50 in Damiani? */
    d->ct += d->Cn_v[0]*(d->alpha_e)*(1 - d->tau[0]);                                                           /* Eq. 1.52 in Damiani */
    // CMV = -0.2*(1 - cos(PI*(d->TAU)))*(d->CNV);
    // af->PMC = af->PMC + CMV + d->Cmi + d->Cmq;

    t_shedding = 2.0*(1.0- d->fn_pp)/(0.19); /* St = 0.19 */                                                               /* Eq. 1.51 in Damiani */

    // IF ( TAU(J,IBLADE) .GT. 1. + TSH/TVL .AND. .NOT. SHIFT) THEN
    if (((d->tau[0])>=(1+t_shedding/d->tv_l)) && (!d->flow_shift)) {  /* shift = LESF in AeroDyn 15 */
        d->tau[0] = 0;
        d->separation = false;
    }

    if (d->tau[0]>1) {
        if (d->alpha[0]<0) {
            if (d->dCn_p[0]<=0 && d->dCn_p[1]>=0) {
                d->separation = false;
                d->tau[0] = 0.;
            }
            if (d->alpha[1]>0) {
                d->separation = false;
                d->tau[0] = 0.;
            }
       } else {
            if ((d->dCn_p[0]>=0) && (d->dCn_p[1]<=0)) {
                d->separation = false;
                d->tau[0] = 0.;
            }

            if (d->alpha[1]<0) {
                d->separation = false;
                d->tau[0] = 0;
            }
        }
    }
    return SAFE;
}


ERROR_CODE
dstall_separate(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    double AFF = 0.0;
    double coeff = 0.0;
    double angle = 0.0;
    const int na = af->na;
    const double aol = af->aol;
    const double cna = af->cna;
    double tf_e = 0.0;
    int idx = 0;

    tf_e  = d->tf;
    d->Dp[0] = d->Dp[1]*exp(-(d->ds)/(d->tp)) + (d->Cn_pot[0] - d->Cn_pot[1])*exp(-0.5*(d->ds)/(d->tp));   /* Eq. 1.34 in Damiani */
    d->Cn_p[0] = d->Cn_pot[0] - d->Dp[0];    /* Eq. 1.33 from Damiani */
    d->dCn_p[0] = d->Cn_p[0] - d->Cn_p[1];     /* delta C_n prime, used to dertermine if AoA is increasing/decreasing */

    /*
     * USE CNPD to determine if AOA is increasing or decreasing.
     * Vortex lift decays more rapidly for decreasing AOA.
     */

    if ((d->alpha[0])*(d->dCn_p[0])<0.0 ) {
        d->flow_shift = true;
    } else {
        d->flow_shift = false;
    }

    AFF = d->Cn_p[0]/cna + aol;    /* Eq. 1.32 in Damiani */

    if(fabs(d->alpha_n)<=PI/2) { /*PIBY2 = PI/2 ?*/
        d->alpha_f[0] =  AFF;
    } else if (d->alpha_n>PI/2) { /*PIBY2 = PI/2 ?*/
        d->alpha_f[0] = PI - AFF;
    } else {
        d->alpha_f[0] = -PI - AFF;
    }

    angle = d->alpha_f[0];
    angle = fmod(angle,(2.0*PI));
    if (PI<angle) {
        angle = angle - 2.0*PI;
    } else if (angle<-PI) {
        angle = angle + 2.0*PI;
    }
    d->alpha_f[0] = angle;
    
    d->alpha_f[0] = MIN(MAX(d->alpha_f[0], af->alpha[0]), af->alpha[na-1] ); /* @todo: NLIFT is af->len, i think */
    idx = get_index(d->alpha_f[0], af->alpha, na);
    
    if (idx==0) {
        idx = 0;
    } else if (idx==na-1) {
        idx = idx-1;
    }    
    coeff = (af->alpha[idx] - d->alpha_f[0])/(af->alpha[idx] - af->alpha[idx+1]);

    success = dstall_set_fn_p(d, idx, aol, coeff, msg, ierr);
    success = dstall_set_ft_p(d, idx, aol, coeff, msg, ierr); /* success = dstall_set_ftbc(d, i, aol, coeff, msg, ierr); */
    success = dstall_check_leading_edge_separation(d, af, msg, ierr);
     
    if (d->flow_shift) {
        tf_e = 1.5*tf_e;
    }

    success = dstall_normal_coefficient(d, tf_e, msg, ierr);
    success = dstall_tangent_coefficient(d, tf_e, msg, ierr);     
    success = dstall_pitch_moment(d, af, tf_e, msg, ierr); /* The pitch moment */
    return SAFE;
}


ERROR_CODE
dstall_check_leading_edge_separation(struct DStall* d, struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    const double cns = af->cns;
    const double cnsl = af->cnsl;
    const double ds = d->ds;
    const double tv_l = d->tv_l;
    if (d->Cn_p[0]>cns) {
        d->separation = true;
    }
    if (d->Cn_p[0]<cnsl) {
        d->separation = true;
    }

    if (d->separation) {
        d->tau[0] = d->tau[1] + ds/tv_l; /* ds = 2*dt*U/chord; */
    }
    return SAFE;
}



ERROR_CODE
dstall_normal_coefficient(struct DStall* d, const double tf_e, char* msg, ERROR_CODE* ierr)
{
    const double ds = d->ds;
    int sign = 0;
    double gamma = 0.0;

    d->Df_n[0] = d->Df_n[1]*exp(-ds/tf_e) + (d->fn_p[0] - d->fn_p[1])*exp(-0.5*ds/tf_e); /* Eq. 1.34 in Damiani */
    d->fn_pp = d->fn_p[0] - d->Df_n[0];                        /* Eq. 1.34 in Damiani, FSP = f_prime, FP = f_prime_prime */
    sign = dstall_sign(d->fn_pp);
    gamma = sqrt(fabs(d->fn_pp))*sign + 1;                    
    d->fk = 0.25*gamma*gamma;                              /* Eq. 108 in Moriarty (right hand side), also */
    d->cn = d->Cn_cp*d->fk + d->Cni_q;
    return SAFE;
}


ERROR_CODE
dstall_tangent_coefficient(struct DStall* d, const double tf_e, char* msg, ERROR_CODE* ierr)
{
    const double ds = d->ds;
    int sign = 0;
    double gamma = 0.0;
    d->Df_t[0] = d->Df_t[1]*exp(-ds/tf_e) + (d->ft_p[0] - d->ft_p[1])*exp(-0.5*ds/tf_e);   /* Eq. 1.34 in Daminani */
    d->ft_pp = d->ft_p[0] - d->Df_t[0];
    sign = dstall_sign(d->ft_pp);
    gamma = sqrt(fabs(d->ft_pp))*sign;
    d->ct = d->ct*gamma;                                 /* Eq. 130 in Damiani */
    return SAFE;
}


ERROR_CODE /* retain for later use */
dstall_pitch_moment(struct DStall* d, struct Airfoil* af, const double tf_e, char* msg, ERROR_CODE* ierr)
{
    // /* for the pitch moment */
    // d->Df_alpha_f_e[0] = d->Df_alpha_f_e[1]*exp(-d->ds/(0.1*TFE)) + (d->alpha_f[0] - d->alpha_f[1]) * exp(-0.5*d->ds/(0.1*TFE)); /* Eq. 1.41 in Damiani */
    // alpha_f_p = d->alpha_f[0] - d->Df_alpha_f_e[0];                                                                          /* Eq. 1.41 in Damiani */
    // 
    // alpha_f_p = MIN(MAX(alpha_f_p, af->alpha[0]), af->alpha[na-1] ); /* @todo: not sure if this is needed if alpha cover -180<x<180 */
    // i = get_index(alpha_f_p, af->alpha, na);
    // 
    // if (i==0) {
    //     coeff = 0.0;
    // } else if (i==na-1) {
    //     i = i-1;
    //     coeff = 1.0;
    // } else {
    //     coeff = (af->alpha[i] - alpha_f_p)/(af->alpha[i] - af->alpha[i+1]);
    // }
    // af->PMC = af->CM[I1] - ((af->CM[I1] - af->CM[I1P1])*P1);  /* linear interpolation times a constant */
    return SAFE;
}
