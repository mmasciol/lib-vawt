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
#include "domain.h"
#include "ac.h"
#include "af.h"
#include "env.h"


ERROR_CODE
domain_initialize(struct Domain* domain, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    
    ac_initialize_list(domain, msg, ierr); CHECKERRK(FATAL_5,"List initialization <ActuatorCylinder>");
    af_initialize_list(domain, msg, ierr); CHECKERRK(FATAL_5,"List initialization <Airfoil>");
    env_initialize_list(domain, msg, ierr); CHECKERRK(FATAL_5,"List initialization <Environment>");
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
domain_free(struct Domain* domain, char* msg, ERROR_CODE *ierr)
{
    ERROR_CODE success = SAFE;
    struct ActuatorCylinder* ac_ = NULL;    
    struct Airfoil* af_ = NULL;    
    struct Environment* env_ = NULL;
    
    list_iterator_start(&domain->ac); 
    while (list_iterator_hasnext(&domain->ac)) { 
        ac_ = (struct ActuatorCylinder*)list_iterator_next(&domain->ac);
        success = ac_free(ac_, msg, ierr); CHECKERRK(WARNING_2, "<ActuatorCylinder>");
    }
    list_iterator_stop(&domain->ac);     

    list_iterator_start(&domain->af); 
    while (list_iterator_hasnext(&domain->af)) { 
        af_ = (struct Airfoil*)list_iterator_next(&domain->af);
        success = af_free(af_, msg, ierr); CHECKERRK(WARNING_2, "<Airfoil>");
    }
    list_iterator_stop(&domain->af);     

    list_iterator_start(&domain->env); 
    while (list_iterator_hasnext(&domain->env)) { 
        env_ = (struct Environment*)list_iterator_next(&domain->env);
        success = env_free(env_, msg, ierr); CHECKERRK(WARNING_2, "<Environment>");
    }
    list_iterator_stop(&domain->env);     
    
    list_destroy(&domain->ac);
    list_destroy(&domain->af);
    list_destroy(&domain->env);
    
    FREE_OBJ(domain);
    return SAFE;
CLEAN_UP:
    return FATAL;
}
