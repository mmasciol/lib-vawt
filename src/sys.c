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


#include <stdlib.h>
#include "sys.h"
#include "./bstring/bstrlib.h"


void 
print_machine_name_to_screen() {
    printf( "%s Ver. %s ", PROGNAME, PROGVERSION); 
    printf( "%c",BUILD_MONTH_CH0 );// build month
    printf( "%c",BUILD_MONTH_CH1 );
    printf( "%c",BUILD_MONTH_CH2 );
    printf( "-" );
    printf( "%c",BUILD_DAY_CH0 );// build day
    printf( "%c",BUILD_DAY_CH1 );
    printf( "-" );
    printf( "%c",BUILD_YEAR_CH0 ); // build year 
    printf( "%c",BUILD_YEAR_CH1 );
    printf( "%c",BUILD_YEAR_CH2 );
    printf( "%c\n",BUILD_YEAR_CH3 );
}


void
split_words(struct bstrList* parsed, int* value, int* variable) 
{
    int n = 0;
    int count = 0;

    while (n<parsed->qty-1) { 
        if (parsed->entry[n]->slen && count<2) { 
            if (count!=0) {
                *variable = n;
                break;
            }
            count++;
        }
        n++;
    }
    return;
}


bool 
is_numeric(const char* str)
{
    char* p = NULL;
    if (str==NULL || *str=='\0' || isspace(*str)) {
        false;
    }
    strtod (str, &p);
    if (*p=='\0') {
        return true;
    } else {
        return false;
    }
}


void
copy_target_string(char* target, unsigned char* source)
{
    while (*source) {
        *target = *source;
        source++;
        target++;
    }
    *target = '\0';
    return;
}
