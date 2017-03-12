/****************************************************************
 * Original work by A. Ning (c) 2016                            *
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sys.h"
#include "af.h"
#include "ac.h"
#include "numerics.h"
#include "domain.h"
#include "./bstring/bstrlib.h"


/* need to change Airfoil to polar */
ERROR_CODE 
af_set_grid(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    return SAFE;
}


size_t
af_size_meter(const void *el) 
{
    return sizeof(struct Airfoil);
}


ERROR_CODE
af_set_beddoes_coefficients(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    success = af_get_cd0(af, msg, ierr); CHECKERRK(FATAL_36,"<Cd0>");                                                   
    success = af_get_aol(af, msg, ierr); CHECKERRK(FATAL_36,"<AoL>");                                                   
    success = af_get_cna(af, msg, ierr); CHECKERRK(FATAL_36,"<Cna>");                                                   
    success = af_get_cns(af, msg, ierr); CHECKERRK(FATAL_36,"<Cns>");                                                   
    success = af_get_cnsl(af, msg, ierr); CHECKERRK(FATAL_36,"<Cnsl>");                                                   
    return SAFE;
CLEAN_UP:
    return FATAL;
}

ERROR_CODE
af_get_aol(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    bstring b = NULL;
    const int na = af->na;
    int i = 0;
    double cn1 = -999.9;
    double cn2 = -999.9;
    double old_aol = af->aol;
    for (i=0 ; i<na-1 ; i++){
        cn1 = (af->cl[i])*cos(af->alpha[i]) + (af->cd[i] - af->cd0)*sin(af->alpha[i]);
        cn2 = (af->cl[i+1])*cos(af->alpha[i+1]) + (af->cd[i+1] - af->cd0)*sin(af->alpha[i+1]);
        if (cn1<=0.0 && 0.0<=cn2) { 
            af->aol = af->alpha[i];
        }            
    }
    if (1E-6<fabs(af->aol-old_aol)) {
        b = bformat("AoL(%s), Calculated=<%f> , airfoil Table=<%f>", af->type, af->aol, old_aol);     
        CHECKERRQ(WARNING_4,(const char*)b->data);
        bdestroy(b);
    }
    return SAFE;
}


ERROR_CODE
af_get_cd0(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    const int na = af->na;
    bstring b = NULL;
    int i = 0;
    double cd_min = af->cd[0];
    double old_cd0 = af->cd0;
    double old_aod = af->aod;
    for (i=0 ; i<na ; i++){
        if (af->cd[i]<cd_min) {
            cd_min = af->cd[i];
            af->aod = af->alpha[i];
        }
    }
    af->cd0 = cd_min;
    if (1E-6<fabs(af->cd0-old_cd0)) {
        b = bformat("Cd0(%s), Calculated=<%f> , airfoil Table=<%f>", af->cd0, old_cd0,af->type);     
        CHECKERRQ(WARNING_4,(const char*)b->data);
        bdestroy(b);
    }
    if (1E-6<fabs(af->aod-old_aod)) {
        b = bformat("AoD(%s), Calculated=<%f> , airfoil Table=<%f>", af->aod, old_aod,af->type);     
        CHECKERRQ(WARNING_4,(const char*)b->data);
        bdestroy(b);
    }
    return SAFE;
}


ERROR_CODE
af_get_cna(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    const int na = af->na;
    // const double c_na = af->cna;
    int i = 0;
    double cn1 = -999.9;
    double cn2 = -999.9;
    double delta = 0.0;
    bstring b = NULL;
    double old_cna = af->cna;
    for (i=0 ; i<na-1 ; i++){
        cn1 = (af->cl[i])*cos(af->alpha[i]) + (af->cd[i] - af->cd0)*sin(af->alpha[i]);
        cn2 = (af->cl[i+1])*cos(af->alpha[i+1]) + (af->cd[i+1] - af->cd0)*sin(af->alpha[i+1]);
        if(cn1<=0.0 && 0.0<=cn2) {
            delta = fabs(af->alpha[i]-af->alpha[i+1]);
            // printf("%f (%f)\n", (cn2-cn1)/delta, af->cna);
            af->cna = (cn2-cn1)/delta;
        }
    }    
    if (1E-6<fabs(af->cna-old_cna)) {
        b = bformat("Cna(%s), Calculated=<%f> , airfoil Table=<%f>", af->cna, old_cna,af->type);     
        CHECKERRQ(WARNING_4,(const char*)b->data);
        bdestroy(b);
    }
    return SAFE;
}


ERROR_CODE
af_get_cnsl(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    const int na = af->na;
    bstring b = NULL;
    // const double c_na = af->cna;
    int i = 0;
    double cn1 = -999.9;
    double old_cnsl = af->cnsl;
    for (i=0 ; i<na-1 ; i++){
        cn1 = (af->cl[i])*cos(af->alpha[i]) + (af->cd[i] - af->cd0)*sin(af->alpha[i]);
        if (-0.55<=af->alpha[i] && af->alpha[i]<=0.0) { /* positive static stall */
            if ((af->cl[i+1])>(af->cl[i])) { /* inflection point */
                af->cnsl = cn1;
                // printf("cnsl(%f) %f   <%f  %f>\n", af->alpha[i], cn1, af->cl[i+1], af->cl[i]);
                if (1E-6<fabs(af->cnsl-old_cnsl)) {
                    b = bformat("Cnsl(%s), Calculated=<%f> airfoil Table=<%f>", af->cnsl, old_cnsl,af->type);     
                    CHECKERRQ(WARNING_4,(const char*)b->data);
                    bdestroy(b);
                }
                return SAFE;
            }            
        }
    }    
    return SAFE;
}

ERROR_CODE
af_get_cns(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    const int na = af->na;
    // const double c_na = af->cna;
    int i = 0;
    bstring b = NULL;
    double cn1 = -999.9;
    double old_cns = af->cns;
    for (i=0 ; i<na-1 ; i++){
        cn1 = (af->cl[i])*cos(af->alpha[i]) + (af->cd[i] - af->cd0)*sin(af->alpha[i]);
        if (0.0<=af->alpha[i] && af->alpha[i]<=0.55) { /* positive static stall */
            if ((af->cl[i+1])<(af->cl[i])) { /* inflection point */
                // printf("cn(%f) %f   <%f  %f>\n", af->alpha[i], cn1, af->cl[i+1], af->cl[i]);
                af->cns = cn1;
                if (1E-6<fabs(af->cns-old_cns)) {
                    b = bformat("Cns(%s), Calculated=<%f> airfoil Table=<%f>", af->cns, old_cns,af->type);     
                    CHECKERRQ(WARNING_4,(const char*)b->data);
                    bdestroy(b);
                }
                return SAFE;
            }            
        }
    }    
    return SAFE;
}


ERROR_CODE 
af_initialize_list(struct Domain* domain, char* msg, ERROR_CODE* ierr)
{    
    list_init(&domain->af);  
    list_attributes_copy(&domain->af, af_size_meter, 1); 
    return SAFE;
}


ERROR_CODE 
af_read_file(list_t* restrict af, char file_name[255], char type_name[32], char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    FILE* fp = NULL;
    bstring line = NULL;
    //int warnings = 0;
    struct tagbstring tokens; 
    struct bstrList* parsed = NULL;
    int var = 0;
    int val = 0;
    struct Airfoil* af_ptr = NULL;

    af_ptr = (struct Airfoil*)list_get_at(af, list_size(af)-1);
    if (!af_ptr) {    
        b = bformat("Could not access airfoil <%d>", list_size(af)+1);     
        CHECKERRQ(FATAL_4,(const char*)b->data);
        bdestroy(b);
        return FATAL;
    }

    fp = fopen(file_name, "r");
    if (!fp) {    
        b = bformat("Trying to find <%s>", file_name);     
        CHECKERRQ(FATAL_4,(const char*)b->data);
        bdestroy(b);
        return FATAL;
    }
   
    copy_target_string(af_ptr->type, (unsigned char*)type_name);
 
    cstr2tbstr(tokens,"\t\n\r ");       
    for (success=SAFE ; NULL!=(line=bgets((bNgetc)fgetc, fp, (char)'\n')) ; bdestroy(line)) {
        parsed = bsplits(line, &tokens);
        split_words(parsed, &val, &var);         
        if (biseqcstrcaseless(parsed->entry[0],"RE")) {
            success = af_set_Re(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"AOL")) {
            success = af_set_aol(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CNA")) {
            success = af_set_cna(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CNS")) {
            success = af_set_cns(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CNSL")) {
            success = af_set_cnsl(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"AOD")) {
            success = af_set_aod(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CD0")) {
            success = af_set_cd0(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"AL")) {
            success = af_set_alpha(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CL")) {
            success = af_set_cl(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CD")) {
            success = af_set_cd(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        } else if (biseqcstrcaseless(parsed->entry[0],"CM")) {
            success = af_set_cm(af_ptr, parsed, msg, ierr); 
            CHECKERRK(FATAL_6, "");
        }
        success = bstrListDestroy(parsed);
    }
    fclose (fp); 
    return SAFE; 
CLEAN_UP:
    success = bstrListDestroy(parsed);
    bdestroy(line);
    fclose (fp);  
    return FATAL;
}


ERROR_CODE
af_set_Re(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr) 
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->Re = (double)atof((const char*)word->entry[i]->data);
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}


ERROR_CODE
af_set_aol(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr) /* Zero Cn angle of attack (word is in degrees) */
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->aol = (double)atof((const char*)word->entry[i]->data);
                af->aol *= DEG2RAD;
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}


ERROR_CODE
af_set_cna(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->cna = (double)atof((const char*)word->entry[i]->data);
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}

ERROR_CODE
af_set_cns(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->cns = (double)atof((const char*)word->entry[i]->data);
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}

ERROR_CODE
af_set_cnsl(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->cnsl = (double)atof((const char*)word->entry[i]->data);
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}


ERROR_CODE
af_set_aod(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr) 
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->aod = (double)atof((const char*)word->entry[i]->data);
                af->aod *= DEG2RAD;
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}


ERROR_CODE
af_set_cd0(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr) 
{    
    int i = 0;
    bstring b = NULL;
    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                af->cd0 = (double)atof((const char*)word->entry[i]->data);
                return SAFE;
            } else {
                b = bformat("Variable %s with value <%s>", word->entry[0]->data, word->entry[i]->data);     
                CHECKERRQ(FATAL_7, (const char*)b->data);                
                bdestroy(b);
                return FATAL;
            }
        }
    }
    return FATAL;
}


ERROR_CODE
af_set_alpha(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int count = 0;
    double* ptr = NULL;
    double val;
    bstring b = NULL;

    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                val = (double)atof((const char*)word->entry[i]->data);
                if (af->alpha==NULL) {
                    af->alpha = malloc(sizeof(double));
                } else {
                    ptr = realloc(af->alpha, (count+1)*sizeof(double));
                    if (ptr==NULL) {
                        b = bformat("<%s> in foil file.", word->entry[0]->data);
                        CHECKERRQ(FATAL_8, (const char*)b->data);                
                        bdestroy(b);
                        return FATAL;
                    } else {
                        af->alpha = ptr;
                    }
                }
                af->alpha[count] = val*DEG2RAD;
                count++;        
            } else {
                b = bformat("Value for %s is not a numeric entry <%s> in file <%s>.", word->entry[0]->data, word->entry[i]->data);
                CHECKERRQ(FATAL_9, (const char*)b->data);
                bdestroy(b);
                return FATAL;
            }
        }
    }
    af->na = count;
    return SAFE;
}


ERROR_CODE
af_set_cl(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int count = 0;
    // double* ptr = NULL;
    double val;
    bstring b = NULL;
    const int na = af->na;

    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                val = (double)atof((const char*)word->entry[i]->data);        
                if (af->cl==NULL) {
                    af->cl = malloc(na*sizeof(double));
                    if (!af->cl) {
                        b = bformat("<%s> in foil file.", word->entry[0]->data);
                        CHECKERRQ(FATAL_8, (const char*)b->data);                
                        bdestroy(b);
                        return FATAL;
                    }
                }
                af->cl[count] = val;
                count++;        
            } else {
                b = bformat("Value <%s> is not a numeric entry for %s.", word->entry[i]->data, word->entry[0]->data);
                CHECKERRQ(FATAL_9, (const char*)b->data);
                bdestroy(b);
                return FATAL;
            }
        }
    }
    
    if (na!=count) {    
        b = bformat("%s has %d entries in file <%s>. %d are required", word->entry[0]->data, count, na);
        CHECKERRQ(FATAL_10, (const char*)b->data);
        bdestroy(b);
        return FATAL;
    }
    return SAFE;
}


ERROR_CODE
af_set_cd(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int count = 0;
    // double* ptr = NULL;
    double val;
    bstring b = NULL;
    const int na = af->na;

    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                val = (double)atof((const char*)word->entry[i]->data);        
                if (af->cd==NULL) {
                    af->cd = malloc(na*sizeof(double));
                    if (!af->cd) {
                        b = bformat("<%s> in foil file.", word->entry[0]->data);
                        CHECKERRQ(FATAL_8, (const char*)b->data);                
                        bdestroy(b);
                        return FATAL;
                    }
                }
                af->cd[count] = val;
                count++;        
            } else {
                b = bformat("Value <%s> is not a numeric entry for %s.", word->entry[i]->data, word->entry[0]->data);
                CHECKERRQ(FATAL_9, (const char*)b->data);
                bdestroy(b);
                return FATAL;
            }
        }
    }
    
    if (na!=count) {    
        b = bformat("%s has %d entries in file <%s>. %d are required", word->entry[0]->data, count, na);
        CHECKERRQ(FATAL_10, (const char*)b->data);
        bdestroy(b);
        return FATAL;
    }
    return SAFE;
}


ERROR_CODE
af_set_cm(struct Airfoil* af, struct bstrList* word, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int count = 0;
    // double* ptr = NULL;
    double val;
    bstring b = NULL;
    const int na = af->na;

    for (i=1 ; i<word->qty-1 ; i++) { /* if the string length is not 0 */
        if (word->entry[i]->slen!=0) {
            if (is_numeric((const char*)word->entry[i]->data)) {        
                val = (double)atof((const char*)word->entry[i]->data);        
                if (af->cm==NULL) {
                    af->cm = malloc(na*sizeof(double));
                    if (!af->cm) {
                        b = bformat("<%s> in foil file.", word->entry[0]->data);
                        CHECKERRQ(FATAL_8, (const char*)b->data);                
                        bdestroy(b);
                        return FATAL;
                    }
                }
                af->cm[count] = val;
                count++;        
            } else {
                b = bformat("Value <%s> is not a numeric entry for %s.", word->entry[i]->data, word->entry[0]->data);
                CHECKERRQ(FATAL_9, (const char*)b->data);
                bdestroy(b);
                return FATAL;
            }
        }
    }
    
    if (na!=count) {    
        b = bformat("%s has %d entries in file <%s>. %d are required", word->entry[0]->data, count, na);
        CHECKERRQ(FATAL_10, (const char*)b->data);
        bdestroy(b);
        return FATAL;
    }
    return SAFE;
}


ERROR_CODE 
af_free(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    FREE_OBJ(af->alpha);
    FREE_OBJ(af->cl);
    FREE_OBJ(af->cd);
    FREE_OBJ(af->cm);
    return SAFE;
}


ERROR_CODE 
af_initialize(struct Airfoil* af, char* msg, ERROR_CODE* ierr)
{
    af->Re = -999.9;
    af->aol = -999.9;
    af->cna = -999.9;
    af->cns = -999.9;
    af->cnsl = -999.9;
    af->aod = -999.9;
    af->cd0 = -999.9;
    af->alpha = NULL;;
    af->cl = NULL;;
    af->cd = NULL;;
    af->cm = NULL;;
    af->na = 0;    
    af->type[0] = '\0';
    af->vtk_point = NULL;
    return SAFE;
}


ERROR_CODE
af_new_foil(list_t* restrict af, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    struct Airfoil new_af;
    struct Airfoil* af_ = NULL;

    success = af_initialize(&new_af, msg, ierr); CHECKERRK(FATAL_5,"Failed to initialize Airfoil");
    list_append(af, &new_af);   
    af_ = (struct Airfoil*)list_get_at(af, list_size(af)-1); CHECK_MALLOC(af_, FATAL_8, "<Airfoil>");
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
af_blend(list_t* restrict af, struct Airfoil* af1, struct Airfoil* af2, const double weight, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    double* alpha = NULL;
    double* cl1 = NULL;
    double* cl2 = NULL;
    double* cd1 = NULL;
    double* cd2 = NULL;
    double* cm1 = NULL;
    double* cm2 = NULL;
    double min_alpha = -999.9;
    double max_alpha = -999.9;
    int i = 0;
    int n = 0;
    struct Airfoil* blend = NULL;

    blend = (struct Airfoil*)list_get_at(af, list_size(af)-1);
    CHECK_MALLOC(blend, FATAL_8, "Blended foil <Airfoil>");

    alpha = af_merge_aoa(af1->alpha, af1->na, af2->alpha, af2->na, &n, msg, ierr); 
    if (!alpha) {
        CALL_CHECKERRK(FATAL_8, "<alpha>");
    }

    min_alpha = min(af1->alpha, af1->na)>min(af2->alpha, af2->na) ? min(af1->alpha, af1->na) : min(af2->alpha, af2->na);
    max_alpha = max(af1->alpha, af1->na)<max(af2->alpha, af2->na) ? max(af1->alpha, af1->na) : max(af2->alpha, af2->na);
    
    for(i=0 ; i<n ; i++) {        
        if (min_alpha<=alpha[i] && alpha[i]<=max_alpha) {
            blend->na++;
            af_increase_array(&blend->alpha, blend->na);
            blend->alpha[blend->na-1] = alpha[i];
        }
    }

    cl1 = linear_interpolation(blend->alpha, blend->na, af1->alpha, af1->cl, af1->na, msg, ierr); CHECK_MALLOC(cl1, FATAL_8, "CL1");
    cl2 = linear_interpolation(blend->alpha, blend->na, af2->alpha, af2->cl, af2->na, msg, ierr); CHECK_MALLOC(cl2, FATAL_8, "CL2");
    cd1 = linear_interpolation(blend->alpha, blend->na, af1->alpha, af1->cd, af1->na, msg, ierr); CHECK_MALLOC(cd1, FATAL_8, "CD1");
    cd2 = linear_interpolation(blend->alpha, blend->na, af2->alpha, af2->cd, af2->na, msg, ierr); CHECK_MALLOC(cd2, FATAL_8, "CD2");
    cm1 = linear_interpolation(blend->alpha, blend->na, af1->alpha, af1->cm, af1->na, msg, ierr); CHECK_MALLOC(cm1, FATAL_8, "CM1");
    cm2 = linear_interpolation(blend->alpha, blend->na, af2->alpha, af2->cm, af2->na, msg, ierr); CHECK_MALLOC(cm2, FATAL_8, "CM2");

    // # linearly blend
    blend->cl = malloc((blend->na)*sizeof(double)); CHECK_MALLOC(blend->cl, FATAL_8, "CL");
    blend->cd = malloc((blend->na)*sizeof(double)); CHECK_MALLOC(blend->cd, FATAL_8, "CD");
    blend->cm = malloc((blend->na)*sizeof(double)); CHECK_MALLOC(blend->cm, FATAL_8, "CM");

    blend->Re = af1->Re + weight*(af2->Re-af1->Re);    
    for (i=0 ; i<blend->na ; i++) {
        blend->cl[i] = cl1[i] + weight*(cl2[i]-cl1[i]);
        blend->cd[i] = cd1[i] + weight*(cd2[i]-cd1[i]);
        blend->cm[i] = cm1[i] + weight*(cm2[i]-cm1[i]);
    }
    blend->aol = (af1->aol) + weight*((af2->aol)-(af1->aol));  
    blend->cna = (af1->cna) + weight*((af2->cna)-(af1->cna));  
    blend->cns = (af1->cns) + weight*((af2->cns)-(af1->cns));  
    blend->cnsl = (af1->cnsl) + weight*((af2->cnsl)-(af1->cnsl)); 
    blend->aod = (af1->aod) + weight*((af2->aod)-(af1->aod));  
    blend->cd0 = (af1->cd0) + weight*((af2->cd0)-(af1->cd0));  
    
    FREE_OBJ(alpha);    
    FREE_OBJ(cl1); 
    FREE_OBJ(cl2);
    FREE_OBJ(cd1);
    FREE_OBJ(cd2);
    FREE_OBJ(cm1);
    FREE_OBJ(cm2);
    return SAFE;
CLEAN_UP:
    FREE_OBJ(alpha);
    FREE_OBJ(cl1); 
    FREE_OBJ(cl2);
    FREE_OBJ(cd1);
    FREE_OBJ(cd2);
    FREE_OBJ(cm1);
    FREE_OBJ(cm2);
    return FATAL;
}


double* 
af_merge_aoa(double* arr1, const int n1, double* arr2, const int n2, int* len, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    int j = 0;
    int k = 0;
    double* af = NULL;

    while ((i<n1) && (j<n2)) {
        if (arr1[i] < arr2[j]) {
            af_increase_array(&af, k);
            if (!af) {
                return NULL;
            }            
            af[k] = arr1[i];
            i++;
            k++;
        } else if (arr1[i] > arr2[j]) {
            af_increase_array(&af, k);
            if (!af) {
                return NULL;
            }            
            af[k] = arr2[j];
            j++;
            k++;
        } else {
            af_increase_array(&af, k);
            if (!af) {
                return NULL;
            }            
            af[k] = arr1[i];
            i++;
            j++;
            k++;
        }
    }
    
    if (i==n1) {
        while (j<n2) {
            af_increase_array(&af, k);
            if (!af) {
                return NULL;
            }            
            af[k] = arr2[j];
            j++;
            k++;
        }
    } else {
        while (i<n1) {
            af_increase_array(&af, k);
            if (!af) {
                return NULL;
            }            
            af[k] = arr1[i];
            i++;
            k++;
        }
    }
    *len = k;
    return af;
}


void 
af_increase_array(double** data, const int N)
{
    double* ptr = NULL;
    if (*data==NULL) {
        *data = malloc(sizeof(double));
    } else {
        ptr = realloc(*data, (N+1)*sizeof(double));
        if (ptr==NULL) {
            *data = NULL;
        } else {
            *data = ptr;  
        }
    }
    return;
}


ERROR_CODE
af_extrapolate(struct Airfoil* af, const double cd_max, const double aspect_ratio, const double cd_min, const int n_alpha, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    bstring b = NULL;
    bool logged = false;
    int i = 0;       
    double maximum_cd = af->cd[0];
    struct ExtrapolatedFoil nf;
    const int na = af->na;

    success = af_initialize_extrapolated_foil(&nf, msg, ierr);     
    nf.cl_adj = 0.7;
    nf.alpha_high = af->alpha[na-1]; 
    nf.cl_high = af->cl[na-1]; 
    nf.cd_high = af->cd[na-1]; 
    nf.cm_high = af->cm[na-1]; 
    nf.alpha_low = af->alpha[0]; 
    nf.cl_low = af->cl[0]; 
    nf.cd_low = af->cd[0]; 
    
    if (cd_min<0.0) {
        b = bformat("Cd <%f>", cd_min);        
        CALL_CHECKERRK(FATAL_16,(const char*)b->data);
    }
    
    if (0.0<aspect_ratio) {
        nf.cdmax = 1.11 + 0.018*(aspect_ratio);
        b = bformat("Overiding default max(Cd) of <%f> and using max(Cd) <%f>", cd_max, nf.cdmax);        
        CHECKERRQ(NOTE_1, (const char*)b->data);
        bdestroy(b);
    } 

    for (i=0 ; i<na ; i++) {
        if (maximum_cd<=af->cd[i]) {
            maximum_cd = af->cd[i];
        }
    }

    if (cd_max<maximum_cd) {
        nf.cdmax = maximum_cd;
        logged = true;
    } else {
        nf.cdmax = cd_max;
    }

    if (logged) {
        b = bformat("Using max(Cd) from aerodynamic table <%f>", nf.cdmax);        
        CHECKERRQ(WARNING_1, (const char*)b->data);
        bdestroy(b);
    }
    
    if (PI/2.0<nf.alpha_high) {
        b = bformat("max(AoA) <%f> deg", nf.alpha_high*RAD2DEG);        
        CALL_CHECKERRK(FATAL_17,(const char*)b->data);
    }
    
    if (nf.alpha_high<-PI/2.0) {
        b = bformat("min(AoA) <%f> deg", nf.alpha_low*RAD2DEG);        
        CALL_CHECKERRK(FATAL_17,(const char*)b->data);
    }
    
    nf.A = (nf.cl_high - nf.cdmax*sin(nf.alpha_high)*cos(nf.alpha_high))*sin(nf.alpha_high)/pow(cos(nf.alpha_high),2);
    nf.B = (nf.cd_high - nf.cdmax*sin(nf.alpha_high)*sin(nf.alpha_high))/cos(nf.alpha_high);
    success = af_extrapolate_foil_1(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <1>.");
    success = af_extrapolate_foil_2(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <2>."); 
    success = af_extrapolate_foil_3(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <3>.");
    success = af_extrapolate_foil_4(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <4>.");
    success = af_extrapolate_foil_5(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <5>.");
    success = af_extrapolate_foil_6(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <6>.");
    success = af_extrapolate_foil_7(&nf, n_alpha, msg, ierr); CHECKERRK(FATAL_18,"Section <7>.");
    nf.N += (af->na-1);
    // nf.N += af->na;
    success = af_extrapolate_pack(af, &nf, n_alpha, msg, ierr); CHECKERRK(FATAL_19,"");
    success = af_free_extrapolated_foil(&nf, msg, ierr);
    return SAFE;
CLEAN_UP:    
    success = af_free_extrapolated_foil(&nf, msg, ierr);
    bdestroy(b);
    return FATAL;
}


ERROR_CODE 
af_extrapolate_pack(struct Airfoil* af, struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    int i = 0;
    int m = 0;
    bstring b = NULL;
    nf->alpha = malloc((nf->N)*sizeof(double)); CHECK_MALLOC(nf->alpha, FATAL_8, "ALPHA");
    nf->cl = malloc((nf->N)*sizeof(double)); CHECK_MALLOC(nf->cl, FATAL_8, "CL");
    nf->cd = malloc((nf->N)*sizeof(double)); CHECK_MALLOC(nf->cd, FATAL_8, "CD");
     
    for (i=0 ; i<n ; i++) {
        nf->alpha[m] = nf->alpha7[i];
        nf->cl[m] = nf->cl7[i];
        nf->cd[m] = nf->cd7[i];
        m++;
    }

    for (i=0 ; i<n-1 ; i++) {
        nf->alpha[m] = nf->alpha6[i];
        nf->cl[m] = nf->cl6[i];
        nf->cd[m] = nf->cd6[i];
        m++;        
    }

    for (i=0 ; i<n-2 ; i++) {
        //for (i=0 ; i<n-1 ; i++) {
        nf->alpha[m] = nf->alpha5[i];
        nf->cl[m] = nf->cl5[i];
        nf->cd[m] = nf->cd5[i];
        m++;
    }

    if (nf->alpha4) {
        for (i=0 ; i<n-3 ; i++) {
            nf->alpha[m] = nf->alpha4[i];
            nf->cl[m] = nf->cl4[i];
            nf->cd[m] = nf->cd4[i];
            m++;
        }

    }

    for (i=0 ; i<af->na ; i++) {
        nf->alpha[m] = af->alpha[i];
        nf->cl[m] = af->cl[i];
        nf->cd[m] = af->cd[i];
        m++;        
    }

    // for (i=0 ; i<af->na ; i++) {
    //     nf->alpha[m] = af->alpha[i];
    //     printf(">>> %f\n",nf->alpha[m]);
    //     nf->cl[m] = af->cl[i];
    //     nf->cd[m] = af->cd[i];
    //     m++;        
    // }
    
    for (i=0 ; i<n-1 ; i++) {
        nf->alpha[m] = nf->alpha1[i];
        nf->cl[m] = nf->cl1[i];
        nf->cd[m] = nf->cd1[i];
        m++;
    }

    for (i=0 ; i<n-1 ; i++) {
        nf->alpha[m] = nf->alpha2[i];
        nf->cl[m] = nf->cl2[i];
        nf->cd[m] = nf->cd2[i];
        m++;
    }

    for (i=0 ; i<n-1 ; i++) {
        nf->alpha[m] = nf->alpha3[i];
        nf->cl[m] = nf->cl3[i];
        nf->cd[m] = nf->cd3[i];
        m++;
    }

    if (m!=nf->N) {
        checkpoint();
        b = bformat("Count of airfoil entries is not consistent: <%d> versus <%d>.", nf->N, m);        
        CALL_CHECKERRK(FATAL_10, (const char*)b->data)
    }

    af_increase_array(&af->alpha, (nf->N));
    af_increase_array(&af->cl, (nf->N));
    af_increase_array(&af->cd, (nf->N));
    af->na = (nf->N);
    for (i=0 ; i<(nf->N) ; i++) {
        af->alpha[i] = nf->alpha[i];
        af->cd[i] = nf->cd[i];
        af->cl[i] = nf->cl[i];
    }
    return SAFE;
CLEAN_UP:    
    bdestroy(b);
    return FATAL;
}

ERROR_CODE
af_extrapolate_foil_1(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    /* alpha_high <-> 90 */
    ERROR_CODE success = SAFE;
    int i = 0;
    const double delta = (PI/2.0-nf->alpha_high)/(n-1);
    const double step = nf->alpha_high;

    nf->N += (n-1);
    nf->alpha1 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->alpha1, FATAL_8, "NF->ALPHA1");
    nf->cd1 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cd1, FATAL_8, "cd1");
    nf->cl1 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cl1, FATAL_8, "cl1");
    for (i=1 ; i<n ; i++) { /* start at 1 or contanating with input foil data */
        nf->alpha1[i-1] = step + delta*i;
    }
    success = af_viterna_extrapolation(nf->alpha1, n, nf->cdmax, nf->cd1, nf->cl1, nf->A, nf->B,  1.0, msg, ierr); CHECKERRK(FATAL_25, "Viterna exprapolation for foil 1");
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_extrapolate_foil_2(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    /* 90 <-> 180-alpha_high */
    ERROR_CODE success = SAFE;
    int i = 0;
    const double delta = ((PI-nf->alpha_high)-(PI/2.0))/(n-1);
    const double step = (PI/2);

    nf->N += (n-1);
    nf->alpha2 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->alpha2, FATAL_8, "ALPHA2");
    nf->cd2 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cd2, FATAL_8, "cd2");
    nf->cl2 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cl2, FATAL_8, "cl2");
    for (i=1 ; i<n ; i++) { /* start at 1 or contanating with input foil data */
        nf->alpha2[i-1] = PI-(step + i*delta);
    }
    success = af_viterna_extrapolation(nf->alpha2, n, nf->cdmax, nf->cd2, nf->cl2, nf->A, nf->B, -nf->cl_adj, msg, ierr); CHECKERRK(FATAL_25, "Viterna exprapolation for foil 2");
    for (i=1 ; i<n ; i++) { /* post correction */ 
        nf->alpha2[i-1] = (step + i*delta);
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_extrapolate_foil_3(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    /* 180-alpha_high <-> 180 */
    ERROR_CODE success = SAFE;
    int i = 0;
    const double delta = (PI-(PI-nf->alpha_high))/(n-1);
    const double step = (PI-nf->alpha_high);

    nf->N += (n-1);
    nf->alpha3 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->alpha3, FATAL_8, "ALPHA3");
    nf->cd3 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cd3, FATAL_8, "cd3");
    nf->cl3 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cl3, FATAL_8, "cl3");
    for (i=1 ; i<n ; i++) { /* start at 1 or contanating with input foil data */
        nf->alpha3[i-1] = PI - (step + i*delta);
    }
    success = af_viterna_extrapolation(nf->alpha3, n, nf->cdmax, nf->cd3, nf->cl3, nf->A, nf->B,  1.0, msg, ierr); CHECKERRK(FATAL_25, "Viterna exprapolation for foil 3");
    for (i=1 ; i<n ; i++) { /* post correction */ 
        nf->alpha3[i-1] = (step + i*delta);
        nf->cl3[i-1] = (nf->alpha3[i-1]-PI)/nf->alpha_high*nf->cl_high*nf->cl_adj; /* override with linear variation */
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_extrapolate_foil_4(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    int i = 0;
    const double delta = (nf->alpha_low-(-nf->alpha_high))/(n-1);
    const double step = (-nf->alpha_high);

    nf->alpha4 = malloc((n-3)*sizeof(double)); CHECK_MALLOC(nf->alpha4, FATAL_8, "ALPHA4");
    nf->cd4 = malloc((n-3)*sizeof(double)); CHECK_MALLOC(nf->cd4, FATAL_8, "cd4");
    nf->cl4 = malloc((n-3)*sizeof(double)); CHECK_MALLOC(nf->cl4, FATAL_8, "cl4");
    if (nf->alpha_low<=-nf->alpha_high) {
        FREE_OBJ(nf->alpha4);
        FREE_OBJ(nf->cd4);
        FREE_OBJ(nf->cl4);        
        nf->alpha5_max = nf->alpha_low;
    } else {        
        nf->N += (n-3);
        for (i=1 ; i<n-2 ; i++) { /* start at 1 or contanating with input foil data */
            nf->alpha4[i-1] = step + (i*delta);
            nf->cl4[i-1] = -nf->cl_high*nf->cl_adj + (nf->alpha4[i-1]+nf->alpha_high)/(nf->alpha_low+nf->alpha_high)*(nf->cl_low+nf->cl_high*nf->cl_adj);
            nf->cd4[i-1] = nf->cd_low + (nf->alpha4[i-1]-nf->alpha_low)/(-nf->alpha_high-nf->alpha_low)*(nf->cd_high-nf->cd_low); 
        }
        nf->alpha5_max = -nf->alpha_high;
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_extrapolate_foil_5(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    /* -90 <-> -alpha_high */
    ERROR_CODE success = SAFE;
    int i = 0;
    const double delta = (nf->alpha5_max-(-PI/2.0))/(n-1);
    const double step = (-PI/2);

    nf->N += (n-1);
    nf->alpha5 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->alpha5, FATAL_8, "ALPHA5");
    nf->cd5 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cd5, FATAL_8, "cd5");    
    nf->cl5 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cl5, FATAL_8, "cl5");    
    for (i=1 ; i<n ; i++) { /* start at 1 or contanating with input foil data */
        nf->alpha5[i-1] = -(step + delta*i);
    }
    success = af_viterna_extrapolation(nf->alpha5, n, nf->cdmax, nf->cd5, nf->cl5, nf->A, nf->B, -nf->cl_adj, msg, ierr); CHECKERRK(FATAL_25, "Viterna exprapolation for foil 5");
    for (i=1 ; i<n ; i++) { /* post correction */ 
        nf->alpha5[i-1] = -nf->alpha5[i-1];
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
af_extrapolate_foil_6(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    /* -180+alpha_high <-> -90 */
    ERROR_CODE success = SAFE;
    int i = 0;
    const double delta = ((-PI/2.0) - (-PI+nf->alpha_high))/(n-1);
    const double step = (-PI+nf->alpha_high);
    
    nf->N += (n-1);
    nf->alpha6 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->alpha6, FATAL_8, "ALPHA6");
    nf->cd6 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cd6, FATAL_8, "cd6");
    nf->cl6 = malloc((n-1)*sizeof(double)); CHECK_MALLOC(nf->cl6, FATAL_8, "cl6");        
    for (i=1 ; i<n ; i++) { /* start at 1 or contanating with input foil data */
        nf->alpha6[i-1] = PI + (step + delta*i);
    }
    success = af_viterna_extrapolation(nf->alpha6, n, nf->cdmax, nf->cd6, nf->cl6, nf->A, nf->B, nf->cl_adj, msg, ierr); CHECKERRK(FATAL_25, "Viterna exprapolation for foil 6");
    for (i=1 ; i<n ; i++) { /* post correction */ 
        nf->alpha6[i-1] = (step + delta*i);
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_extrapolate_foil_7(struct ExtrapolatedFoil* nf, const int n, char* msg, ERROR_CODE* ierr)
{
    /* -180 <-> -180 + alpha_high */
    ERROR_CODE success = SAFE;
    int i = 0;
    const double delta = ((-PI+nf->alpha_high) - (-PI))/(n-1);
    const double step = (-PI);

    nf->N += n;
    nf->alpha7 = malloc((n)*sizeof(double)); CHECK_MALLOC(nf->alpha7, FATAL_8, "ALPHA7");
    nf->cd7 = malloc((n)*sizeof(double)); CHECK_MALLOC(nf->cd7, FATAL_8, "cd7");
    nf->cl7 = malloc((n)*sizeof(double)); CHECK_MALLOC(nf->cl7, FATAL_8, "cl7");    
    for (i=0 ; i<n ; i++) { /* start at 1 or contanating with input foil data */
        nf->alpha7[i] = PI+ (step + delta*i);
    }
    success = af_viterna_extrapolation(nf->alpha7, n+1, nf->cdmax, nf->cd7, nf->cl7, nf->A, nf->B, 1.0, msg, ierr); CHECKERRK(FATAL_25, "Viterna exprapolation for foil 7");
    for (i=0 ; i<n ; i++) { /* post correction */ 
        nf->alpha7[i] -= PI;
        nf->cl7[i] = (nf->alpha7[i]+PI)/nf->alpha_high*nf->cl_high*nf->cl_adj; /* override with linear variation */
    }
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_free_extrapolated_foil(struct ExtrapolatedFoil* nf, char* msg, ERROR_CODE* ierr)
{
    FREE_OBJ(nf->alpha1);
    FREE_OBJ(nf->alpha2);
    FREE_OBJ(nf->alpha3);
    FREE_OBJ(nf->alpha4);
    FREE_OBJ(nf->alpha5);
    FREE_OBJ(nf->alpha6);
    FREE_OBJ(nf->alpha7);
    FREE_OBJ(nf->cd1);
    FREE_OBJ(nf->cd2);
    FREE_OBJ(nf->cd3);
    FREE_OBJ(nf->cd4);
    FREE_OBJ(nf->cd5);
    FREE_OBJ(nf->cd6);
    FREE_OBJ(nf->cd7);
    FREE_OBJ(nf->cl1);
    FREE_OBJ(nf->cl2);
    FREE_OBJ(nf->cl3);
    FREE_OBJ(nf->cl4);
    FREE_OBJ(nf->cl5);
    FREE_OBJ(nf->cl6);
    FREE_OBJ(nf->cl7);
    FREE_OBJ(nf->alpha);
    FREE_OBJ(nf->cl);
    FREE_OBJ(nf->cd);
    return SAFE;
}

ERROR_CODE
af_initialize_extrapolated_foil(struct ExtrapolatedFoil* nf, char* msg, ERROR_CODE* ierr)
{
    nf->N = 0;
    nf->cl_adj = 0.7;
    nf->alpha_high = 0.0;
    nf->cl_high = 0.0;
    nf->cd_high = 0.0;
    nf->cm_high = 0.0;
    nf->alpha_low = 0.0;
    nf->cl_low = 0.0;
    nf->cd_low = 0.0;
    nf->cdmax = 0.0;
    nf->A = 0.0;
    nf->B = 0.0;
    nf->alpha5_max = 0.0;
    nf->alpha = NULL;
    nf->alpha1 = NULL;
    nf->alpha2 = NULL;
    nf->alpha3 = NULL;
    nf->alpha4 = NULL;
    nf->alpha5 = NULL;
    nf->alpha6 = NULL;
    nf->alpha7 = NULL;
    nf->cd = NULL;
    nf->cd1 = NULL;
    nf->cd2 = NULL;
    nf->cd3 = NULL;
    nf->cd4 = NULL;
    nf->cd5 = NULL;
    nf->cd6 = NULL;
    nf->cd7 = NULL;
    nf->cl = NULL;
    nf->cl1 = NULL;
    nf->cl2 = NULL;
    nf->cl3 = NULL;
    nf->cl4 = NULL;
    nf->cl5 = NULL;
    nf->cl6 = NULL;
    nf->cl7 = NULL;
    return SAFE;
}


ERROR_CODE 
af_get_cl(double* value, struct Airfoil* af, const double alpha, const double Re, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    int xi = 50;  /* index */ 
    int yi = 1;  /* index */
    int dim_n = 1;
    int dim_m = af->na;
    /* @todo: xi should be stored in the struct for quicker access at the array start point 
     * based on the previous iteration 
     */
    success = linear_bivariate_spline(value, alpha, Re, af->alpha, af->cl, dim_m, dim_n, &xi, &yi, msg, ierr); CHECKERRK(FATAL_25, "Linear bivariate spline");
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE 
af_get_cd(double* value, struct Airfoil* af, const double alpha, const double Re, char* msg, ERROR_CODE* ierr)
{
    ERROR_CODE success = SAFE;
    int xi = 50;  /* index */ 
    int yi = 1;  /* index */
    int dim_n = 1;
    int dim_m = af->na;
    /* @todo: xi should be stored in the struct for quicker access at the array start point 
     * based on the previous iteration 
     */
    success = linear_bivariate_spline(value, alpha, Re, af->alpha, af->cd, dim_m, dim_n, &xi, &yi, msg, ierr); CHECKERRK(FATAL_25, "Linear bivariate spline");
    return SAFE;
CLEAN_UP:
    return FATAL;
}


ERROR_CODE
af_viterna_extrapolation(double* alpha, const int n, const double cdmax, double* cd, double* cl, const double A, const double B, const double cl_adj, char* msg, ERROR_CODE* ierr)
{
    double* max = NULL; 
    int i = 0;

    max = malloc((n-1)*sizeof(double)); CHECK_MALLOC(max, FATAL_8, "max");
    for (i=0 ; i<n-1 ; i++) {        
        max[i] = alpha[i];
        if (max[i]<=0.0001) {
            max[i] = 0.0001;
        }
    }

    for (i=0 ; i<n-1 ; i++) {
        cl[i] = cdmax/2*sin(2.0*alpha[i]) + A*pow(cos(max[i]),2)/sin(max[i]);
        cl[i]*=cl_adj;
        cd[i] = cdmax*pow(sin(max[i]),2) + B*cos(max[i]);
    }    
    FREE_OBJ(max);
    return SAFE;
CLEAN_UP:
    FREE_OBJ(max);
    return FATAL;
}
