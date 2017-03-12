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


#include "error.h"
#include "sys.h"
#include "./bstring/bstrlib.h"


const char ERROR_STRING[][256] = {
    /* SAFE       */  "",
    /* NOTE       */  "",
    /* WARNING    */  "",
    /* ERROR      */  "",
    /* FATAL      */  "",
    /* FATAL_4    */  "Airfoil file failed to open", 
    /* FATAL_5    */  "Could not initialize model",
    /* FATAL_6    */  "Initialization operation terminated prematurely reading the blade aerodynamics file",
    /* FATAL_7    */  "Error parsing input file",
    /* FATAL_8    */  "Failed to allocate memory",
    /* FATAL_9    */  "Inconsistent value",
    /* FATAL_10   */  "Inconsistent number of blade foil entries",
    /* FATAL_11   */  "Could not blend airfoils",
    /* FATAL_12   */  "Failed to preallocate actuator cylinder matrix",
    /* FATAL_13   */  "Failed to evaluate the actuator cylinder forces",
    /* FATAL_14   */  "Assignment error",
    /* FATAL_15   */  "Error in airfoil extrapolation",
    /* FATAL_16   */  "Drag (Cd) coefficient cannot be negative",
    /* FATAL_17   */  "AoA is outside the range for airfoil extrapolation",
    /* FATAL_18   */  "Failed to create extrapolated table for one lift coefficient section",
    /* FATAL_19   */  "Failed to pack extrapolated airfoils",
    /* FATAL_20   */  "Failed to preallocate actuator cylinder memory",
    /* FATAL_21   */  "Failed to calculate average power",
    /* FATAL_22   */  "Failed to set inflow wind speed",
    /* FATAL_23   */  "Memory is not allocated",
    /* FATAL_24   */  "Failed to set VAWT rotational speed",
    /* FATAL_25   */  "Interpolation failure",
    /* FATAL_26   */  "Failed during function evaluation in solver",
    /* FATAL_27   */  "Integration error",
    /* FATAL_28   */  "Array size is exceeded",
    /* FATAL_29   */  "Failed to set initial theta offset",
    /* FATAL_30   */  "Could not allocate blades",
    /* FATAL_31   */  "Failed to deallocate memory",
    /* FATAL_32   */  "Could not set blade azimuth",
    /* FATAL_33   */  "Could not initialize the dynamic stall model",
    /* FATAL_34   */  "Index is out of range",
    /* FATAL_35   */  "Dynamic stall error",
    /* FATAL_36   */  "Error assigning the dynamic stall coefficient",
    /* FATAL_37   */  "2D iteration failure",
    /* ERROR_1    */  "",
    /* WARNING_1  */  "Maximum drag coefficient value in aerodynamic table is larger than maximum requested",
    /* WARNING_2  */  "Failed to deallocate memory",
    /* WARNING_3  */  "Mach number is above 1",
    /* WARNING_4  */  "Inconsistent value",
    /* WARNING_5  */  "Library is compiled without the VTK engine. Cannot render the geometry",
    /* NOTE_1     */  "Estimating maximum drag coefficient from the aspect ratio",
    /* NOTE_2     */  "",
};


void 
ierr_reset(char* msg, ERROR_CODE* ierr)
{
    *ierr = SAFE;
    msg[0] = '\0';
}


void 
ierr_set(char* msg, ERROR_CODE* ierr, const ERROR_CODE new_ierr, const char* file, int line)
{
    int error_number = 0;
    bstring note = NULL;
    bstring out_string = NULL;
    bstring message = bformat("%s", msg); /* first format msg to contain previously raised errors */
    
    if (new_ierr>=NOTE_1) { 
        error_number = new_ierr-NOTE_1 + 1;
        note = bformat("[INFO %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        //out_string = bformat("%s(%s:%d) %s.\n", note->data, file, line, ERROR_STRING[new_ierr]);
        out_string = bformat("%s(%s:%d) %s.\n", note->data, file, line, ERROR_STRING[new_ierr]);
        if (*ierr<=NOTE) {
            *ierr = NOTE;
        };
        
    } else if (new_ierr>=WARNING_1) { /* did not quite fail. Let users know what the error is */    
        error_number = new_ierr-WARNING_1+1;
        note = bformat("[WARN %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        out_string = bformat("%s(%s:%d) %s.\n", note->data, file, line, ERROR_STRING[new_ierr]);
        if (*ierr<=WARNING) {
            *ierr = WARNING;
        };
    } else if (new_ierr>=ERROR_1 ) { /* WT failed but recovered */    
        error_number = new_ierr-ERROR_1;    
        note = bformat("[ERROR %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        out_string = bformat("%s(%s:%d) %s.\n", note->data, file, line, ERROR_STRING[new_ierr]);
                if (*ierr<=ERROR) {
            *ierr = ERROR;
        };
    } else { /* WT failed and program must end prematurely */    
        error_number = new_ierr-1;
        note = bformat("[FATAL %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        out_string = bformat("%s(%s:%d) %s.\n", note->data, file, line, ERROR_STRING[new_ierr]);
                *ierr = FATAL;
    };
    
    bconcat(message, out_string);
    btrunc(message, ERROR_STRING_LENGTH-1);
    copy_string(msg, message->data);
    bdestroy(note);
    bdestroy(out_string);
    bdestroy(message);
}


void 
ierr_msg_set(char* msg, ERROR_CODE* ierr, const ERROR_CODE new_ierr, const char* file, int line, const char* in_string, ...)
{
    va_list arglist;
    bstring note = NULL;
    bstring out_string = NULL;
    bstring user_msg = NULL;
    bstring message = bformat("%s", msg); /* first format msg to contain previously raised errors */ 
    const int START_VSNBUFF = 16;    
    int error_number = 0;
    // int ret = 0;
    int r = 0;
    int n = 0;
    
    /* This is a re-implementation of the bstring library routines 'bformat(...)  
     * Take the variable argument list and create a string with it. This lets you
     * create a custom message to be rpinted to the terminal.
     */
    do { 
        n = (int)(2*strlen(in_string));
        if (n<START_VSNBUFF) {
            n = START_VSNBUFF;
        };
        user_msg = bfromcstralloc(n+2, "");
        if (!user_msg) {
            n = 1;
            user_msg = bfromcstralloc(n+2, "");
            if (!user_msg) {
                user_msg = NULL;
                break;
            };
        };
        while (1) {
            va_start(arglist, in_string);      
#           if !defined(_MSC_VER)
            r = vsnprintf((char*)user_msg->data, n+1, in_string, arglist); /* this is a copy of exvsnprintf in bstring library */
#           else
            r = vsnprintf_s((char*)user_msg->data, n, _TRUNCATE, in_string, arglist); /* windows way (or ISO C11 Annex K) way of doing things */
            /* This function works, but you need to specify the _CRT_SECURE_NO_WARNINGS compiler flag. Visual Studio hates this: 
             * r = vsnprintf((char*)user_msg->data, n + 1, in_string, arglist);
             */
#           endif
            va_end(arglist);
            user_msg->data[n] = (unsigned char)'\0';
            user_msg->slen = (int)strlen((char*)user_msg->data);
            if (user_msg->slen < n) {
                break;
            };
            if (r>n) {
                n = r;
            } else {
                n += n;
            };
            if (0!=balloc(user_msg, n+2)) {
                bdestroy(user_msg);
                break;
            };
        }; 
    } while (0);
  
    if (new_ierr>=NOTE_1) { 
        /* Notice for users */
        error_number = new_ierr - NOTE_1 + 1;
        note = bformat("[INFO %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        out_string = bformat("%s(%s:%d) %s. %s\n", note->data, file, line, ERROR_STRING[new_ierr], user_msg->data);
        if (*ierr<=NOTE) {
            *ierr = NOTE;
        };
    } else if (new_ierr>=WARNING_1) { 
        /* MAP did not quite fail. Let users know what the error is */    
        error_number = new_ierr - WARNING_1+1;
        note = bformat("[WARN %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        out_string = bformat("%s(%s:%d) %s. %s\n", note->data, file, line, ERROR_STRING[new_ierr], user_msg->data);
        if (*ierr<=WARNING) {
            *ierr = WARNING;
        };
    } else if (new_ierr>=ERROR_1 ) { 
        /* MAP failed but recovered */    
        error_number = new_ierr - ERROR_1;    
        note = bformat("[ERROR %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };        
        out_string = bformat("%s(%s:%d) %s. %s\n", note->data, file, line, ERROR_STRING[new_ierr], user_msg->data);
        if (*ierr<=ERROR) {
            *ierr = ERROR;
        };
    } else { 
        /* MAP failed and program must end prematurely */    
        error_number = new_ierr-1;
        note = bformat("[FATAL %d]", error_number);
        while (note->slen<12) {
            bconchar(note, ' ');
        };
        out_string = bformat("%s(%s:%d) %s. %s\n", note->data, file, line, ERROR_STRING[new_ierr], user_msg->data);
        *ierr = FATAL;
    };
    
    bconcat(message, out_string);
    btrunc(message, ERROR_STRING_LENGTH-1);
    copy_string(msg, message->data);
    bdestroy(note);
    bdestroy(out_string);
    bdestroy(message);
    bdestroy(user_msg);
}


void 
copy_string(char* target, unsigned char* source)
{
    while (*source) {
        *target = *source;
        source++;
        target++;
    };
    *target = '\0';
}


ERROR_CODE 
print_help_to_screen()
{
    print_machine_name_to_screen();
    printf("Input file section definitions:\n");
    printf("    \n");
    return SAFE;
}

