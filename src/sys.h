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


#ifndef _SYS_H
#define _SYS_H

#include <stdio.h>


#if defined(_WIN32) || defined(_WIN64)
#    include <Windows.h>
#    include <tchar.h>
#else
#    include <unistd.h>
#endif


#if defined(_MSC_VER)
typedef int bool;
#    define false 0
#    define true 1
#    define map_snprintf _snprintf
#    define map_strcat(a,b,c) strcat_s(a,b,c)
#    define MAP_STRCPY(a,b,c) strcpy_s(a,b,c)
#else
#    include <stdbool.h>
#    define map_snprintf snprintf
#    define map_strcat(a,b,c) strncat(a,c,b)
#    define MAP_STRCPY(a,b,c) strcpy(a,c)
#endif


#ifndef BUILD_DEFS_H
#    define BUILD_DEFS_H
#    define BUILD_YEAR_CH0 (__DATE__[ 7])
#    define BUILD_YEAR_CH1 (__DATE__[ 8])
#    define BUILD_YEAR_CH2 (__DATE__[ 9])
#    define BUILD_YEAR_CH3 (__DATE__[10])
#    define BUILD_MONTH_IS_JAN (__DATE__[0] == 'J' && __DATE__[1] == 'a' && __DATE__[2] == 'n')
#    define BUILD_MONTH_IS_FEB (__DATE__[0] == 'F')
#    define BUILD_MONTH_IS_MAR (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'r')
#    define BUILD_MONTH_IS_APR (__DATE__[0] == 'A' && __DATE__[1] == 'p')
#    define BUILD_MONTH_IS_MAY (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'y')
#    define BUILD_MONTH_IS_JUN (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'n')
#    define BUILD_MONTH_IS_JUL (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'l')
#    define BUILD_MONTH_IS_AUG (__DATE__[0] == 'A' && __DATE__[1] == 'u')
#    define BUILD_MONTH_IS_SEP (__DATE__[0] == 'S')
#    define BUILD_MONTH_IS_OCT (__DATE__[0] == 'O')
#    define BUILD_MONTH_IS_NOV (__DATE__[0] == 'N')
#    define BUILD_MONTH_IS_DEC (__DATE__[0] == 'D')
#    define BUILD_MONTH_CH0 (__DATE__[ 0])
#    define BUILD_MONTH_CH1 (__DATE__[ 1])
#    define BUILD_MONTH_CH2 (__DATE__[ 2])
#    define BUILD_DAY_CH0 ((__DATE__[4] >= '0') ? (__DATE__[4]) : '0')
#    define BUILD_DAY_CH1 (__DATE__[ 5])
#endif 


#ifdef DEBUG
#    define checkpoint() printf("Checkpoint: Line %d in file %s\n", __LINE__, __FILE__);
#else
#    define checkpoint() 
#endif 

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))    

#define MAX_INIT_VERSION_STRING_LENGTH 99
#define MAX_INIT_COMPILING_DATA_STRING_LENGTH 25
#define ERROR_STRING_LENGTH 1024
#define MACHINE_EPSILON 1e-16
#define TIME_BUFFER_SIZE 64
#define PROGNAME "diats"
#define PROGVERSION "0.10.00"
#define SPACE_LENGTH 12
#define DEG2RAD 0.017453292519943295
#define RAD2DEG 57.295779513082323
#define PI 3.14159265359
// #define LOG_ERR(code) do{ *ierr=code; } while (0);
#define CHECKERR(code) do{ if(*ierr!=SAFE) {ierr_set(msg, ierr, code, __FILE__, __LINE__);} } while(0);
#define CHECKERRQ(code, string) do{ ierr_msg_set(msg, ierr, code, __FILE__, __LINE__, string); } while(0);
#define CHECKERRK(code, string) do{if(success==FATAL) { ierr_msg_set(msg, ierr, code, __FILE__, __LINE__, string);  goto CLEAN_UP; } else if (success!=SAFE) { ierr_msg_set(msg, ierr, code, __FILE__, __LINE__, string); }} while(0);
#define CALL_CHECKERRK(code, string) do{success=FATAL; if(success==FATAL) {ierr_msg_set(msg, ierr, code, __FILE__, __LINE__, string);  goto CLEAN_UP; } else if (success!=SAFE) { ierr_msg_set(msg, ierr, code, __FILE__, __LINE__, string); }} while(0);
#define CHECK_MALLOC(var, code, string) do{ if(var==NULL) {ierr_msg_set(msg, ierr, code, __FILE__, __LINE__, string); goto CLEAN_UP;}} while(0);
#define FREE_OBJ(obj) do{ if(obj!=NULL) {free(obj); obj=NULL;};} while(0);

struct bstrList;
void print_machine_name_to_screen();
void split_words(struct bstrList* parsed, int* value, int* variable);
bool is_numeric(const char* str);
void copy_target_string(char* target, unsigned char* source);


#endif /* _SYS_H */
