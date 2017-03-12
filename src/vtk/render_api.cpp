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


#include <vector>
#include <stdio.h>

#include "render_engine.h"
#include "line_generator.h"
#include "render_api.h"


#define checkpoint() printf("Checkpoint: Line %d in file %s\n", __LINE__, __FILE__);



#ifdef _WIN32
#  include "stdbool.h"
#  define CALL __declspec( dllexport ) //define something for Windows (64-bit)
#elif _WIN64
#  include "stdbool.h"
#  define CALL __declspec( dllexport ) //define something for Windows (64-bit)
#else
#  ifdef __cplusplus
#    define CALL extern "C"  
#  else
#    define CALL 
#  endif
#  include <stdbool.h>
#endif

CALL void*
new_point()
{    
    Points* point = new Points();
    return point;
}


CALL void*
new_vtk()
{
    VawtVTK* vtk = new VawtVTK();
    return vtk;
}


CALL void
delete_point(void* p)
{
    Points* point = static_cast<Points*>(p);
    delete point;
}


CALL void
delete_vtk(void* engine)
{
    VawtVTK* vtk = static_cast<VawtVTK*>(engine);
    delete vtk;
}


CALL int
C_vtk_airfoil_hook(void* p, void* engine, const double chord)
{
    int error = 0;
    // const double cp = 0.25;
    VawtVTK* vtk = static_cast<VawtVTK*>(engine);    
    Points* point = static_cast<Points*>(p);
    
    error = point->read_file(chord);
    error = point->add_points();
    if (error!=0) {
        return 0;
    }
    point->set_polydata();
    point->set_mapper_actor();
    vtk->add_actor(point->get_actor());
    return 0;
}


CALL int
C_vtk_airfoil_displace(void* p, const double x, const double y, const double z)
{
    Points* point = static_cast<Points*>(p);
    point->translate(x,y,z);
    return 0;
}


CALL int
C_vtk_airfoil_rotate(void* p, const double phi, const double the, const double psi)
{
    Points* point = static_cast<Points*>(p);
    point->rotate(phi,the,psi);
    return 0;
}


CALL int
C_vtk_airfoil_transform(void* p)
{
    Points* point = static_cast<Points*>(p);
    point->transform();
    return 0;
}


CALL void
vtk_start(void* engine)
{
    VawtVTK* vtk = static_cast<VawtVTK*>(engine);
    vtk->set_background(0.0, 0.0, 0.0);
    vtk->start();
    return;
}


