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


#ifndef _BLADE_H
#define _BLADE_H

#include <stdio.h>
#include "error.h"


struct ActuatorCylinder;
struct PowellCoeff;
struct DStall;
struct Airfoil;


struct BladeElement {
    struct Airfoil* af; /* reference to airfoil */
    struct DStall* d; /* because we need to preserve time history of fluid evolution over blade element chord */
    double chord;
    double r;
    double W;     /* relative velocity */   
    double phi;   /* wind direction in the inertial frame of length netheta */
    double alpha; /* angle of attack */
    double Re;   
    double twist;
    double swept_offset;
    double the; /* this is the theta value for the local element in the global frame, the = theta + azimuth + swept_offset */
    double delta;
};


struct Blade {
    struct BladeElement* elem;
    double azimuth;
    double* Rb; /* body frame */ //Rp_along_blade = ac->coeff->Rp_along_blade;
    double* Tb; /* body frame */ //Tp_along_blade = ac->coeff->Tp_along_blade;
    double* Zb; /* body frame */ //_along_blade = ac->coeff->Zp_along_blade;
    double* Ri; /* inertial frame */ //_rotated = ac->coeff->Rp_rotated;
    double* Ti; /* inertial frame */ //_rotated = ac->coeff->Tp_rotated;
    double* Q;
    double Qavg;
    double rmax;
    int n; /* @todo: what is this for?? */
    int m; /* number of blade elements */
    void** vtk;
};


ERROR_CODE blade_set_number(struct ActuatorCylinder* ac, const int num, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_initialize(struct Blade* blade, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_allocate_memory(struct Blade* blade, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_free(struct Blade* blade, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_set_azimuth(struct Blade* blade, const double azimuth, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_set_delta(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_set_rmax(struct ActuatorCylinder* ac, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_increment_azimuth(struct ActuatorCylinder* ac, const double dazimuth, char* msg, ERROR_CODE* ierr);
ERROR_CODE blade_increment_element_azimuth(struct ActuatorCylinder* ac, const double theta, char* msg, ERROR_CODE* ierr);


    
#endif /* _BLADE_H */
