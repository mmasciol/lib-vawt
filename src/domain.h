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


#ifndef _DOMAIN_H
#define _DOMAIN_H

#include <stdio.h>
#include "error.h"
#include "./simclist/simclist.h"


struct Environment;


struct Domain {
    list_t ac;
    list_t af;
    list_t env;
};


ERROR_CODE domain_initialize(struct Domain* domain, char* msg, ERROR_CODE* ierr); 
ERROR_CODE domain_free(struct Domain* domain, char* msg, ERROR_CODE *ierr);


#endif /* _DOMAIN_H */
