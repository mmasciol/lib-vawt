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


CALL void* new_point();
CALL void* new_vtk();

CALL void delete_point(void* p);
CALL void delete_vtk(void* engine);

CALL int C_vtk_airfoil_hook(void* p, void* engine, const double chord);
CALL int C_vtk_airfoil_displace(void* p, const double x, const double y, const double z);
CALL int C_vtk_airfoil_rotate(void* p, const double phi, const double the, const double psi);
CALL int C_vtk_airfoil_transform(void* p);


CALL void vtk_start(void* engine);


