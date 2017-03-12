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


#ifndef ERROR_H
#define ERROR_H


typedef enum ERROR_CODE {
    SAFE     ,  // = 0 by default
    NOTE     ,  // = 1
    WARNING  ,  // = 2
    ERROR    ,  // = 3
    FATAL    ,  // = 4
  
    /* These are used internally to the program and are used to
     * map text to the specific error code */
    FATAL_4   , /* Airfoil file failed to open */
    FATAL_5   , /* Could not initialize model */
    FATAL_6   , /* Initialization operation terminated prematurely reading the blade aerodynamics file */
    FATAL_7   , /* Error parsing input file */
    FATAL_8   , /* Failed to allocate memory */
    FATAL_9   , /* Inconsistent value */
    FATAL_10  , /* Inconsistent number of blade foil entries */
    FATAL_11  , /* Could not blend airfoils */
    FATAL_12  , /* Failed to generate AC matrix */
    FATAL_13  , /* Failed to evaluate the actuator cylinder forces */
    FATAL_14  , /* Assignment error */
    FATAL_15  , /* Error in airfoil extrapolation */
    FATAL_16  , /* Drag (Cd) coefficient cannot be negative */
    FATAL_17  , /* AoA is outside the range for airfoil extrapolation */
    FATAL_18  , /* Failed to create extrapolated table */
    FATAL_19  , /* Failed to pack extrapolated airfoils */
    FATAL_20  , /* Failed to preallocate actuator cylinder memory */
    FATAL_21  , /* Failed to calculate average power */
    FATAL_22  , /* Failed to set inflow wind speed */  
    FATAL_23  , /* Memory is not allocated */
    FATAL_24  , /* Failed to set VAWT rotational speed */
    FATAL_25  , /* Interpolation failure */
    FATAL_26  , /* Failed during function evaluation in solver */
    FATAL_27  , /* Integration error */
    FATAL_28  , /* Array size is exceeded */
    FATAL_29  , /* Failed to set initial theta offset */
    FATAL_30  , /* Could not allocate blades */
    FATAL_31  , /* Failed to deallocate memory */
    FATAL_32  , /* Could not set blade azimuth */
    FATAL_33  , /* Could not initialize the dynamic stall model */
    FATAL_34  , /* Index is out of range */
    FATAL_35  , /* Dynamic stall error */
    FATAL_36  , /* Error assigning the dynamic stall coefficient */
    FATAL_37  , /* Exceeded iteration count */
    ERROR_1   , /*  */
    WARNING_1 , /* Maximum drag coefficient value in aerodynamic table is larger than maximum requested */
    WARNING_2 , /* Failed to deallocate memory */
    WARNING_3 , /* Mach number is above 1 */
    WARNING_4 , /* Inconsistent value */
    WARNING_5 , /* Library is compiled without the VTK engine. Cannot render the geometry */
    NOTE_1    , /* Estimating maximum drag coefficient from te aspect ratio */
    NOTE_2    , 
} ERROR_CODE ; 


void ierr_reset(char* msg, ERROR_CODE* ierr);
void ierr_set(char* msg, ERROR_CODE* ierr, const ERROR_CODE new_ierr, const char* file, int line);
void ierr_msg_set(char* msg, ERROR_CODE* ierr, const ERROR_CODE new_ierr, const char* file, int line, const char* in_string, ...);
void copy_string(char* target, unsigned char* source);
ERROR_CODE print_help_to_screen();


#endif /* ERROR_H */
