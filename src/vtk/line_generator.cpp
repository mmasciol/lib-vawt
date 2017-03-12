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


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkTransform.h"
#include "vtkLine.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#pragma GCC diagnostic pop

#include <stdio.h>
#include <vector>

#include "line_generator.h"
#include "../bstring/bstrlib.h"


bool 
Points::is_numeric(const char* str)
{
    char* p = NULL;
    if (str==NULL || *str=='\0' || isspace(*str)) {
        return false;
    }
    strtod (str, &p);
    if (*p=='\0') {
        return true;
    } else {
        return false;
    }
}

int
Points::translate(const double x, const double y, const double z)
{
    this->trans->Translate(x,y,z);
    return 0;
}


int
Points::rotate(const double x, const double y, const double z)
{    
    this->trans->RotateZ(z);
    this->trans->RotateY(y);
    this->trans->RotateX(x);
    return 0;
}


int
Points::transform()
{
    this->actor->SetUserTransform(this->trans);
    return 0;
}


int
Points::set_colors()
{
    unsigned int i = 0;
    unsigned char red[3] = { 255, 0, 0 };
    vtkUnsignedCharArray* colors = vtkUnsignedCharArray::New();
    colors->SetNumberOfComponents(3);
    for (i=0 ; i<x.size() ; i++) {
        colors->InsertNextTypedTuple(red);
    }
    linesPolyData->GetCellData()->SetScalars(colors);
    colors->Delete();
    return 0;
}


int
Points::set_mapper_actor()
{
    this->mapper->SetInputData(this->linesPolyData);
    this->actor->SetMapper(this->mapper);
    return 0;
}


int
Points::set_polydata()
{
    int i = 0;
    const int n = this->x.size()-1; 
    for(i=0 ; i<n ; i++) {
        this->line->GetPointIds()->SetId(0,i);
        this->line->GetPointIds()->SetId(1,i+1);
        this->lines->InsertNextCell(this->line);
    }          
    this->linesPolyData->SetPoints(this->points);
    this->linesPolyData->SetLines(this->lines);    
    return 0;
}


Points::Points() 
{
    this->trans = vtkTransform::New();
    this->linesPolyData = vtkPolyData::New();   
    this->mapper = vtkPolyDataMapper::New();
    this->actor = vtkActor::New();
    this->line = vtkLine::New();
    this->lines = vtkCellArray::New();
    this->points = vtkPoints::New();
}


Points::~Points() 
{
    this->trans->Delete();
    this->linesPolyData->Delete();
    this->mapper->Delete();
    this->actor->Delete();
    this->line->Delete();
    this->lines->Delete();
    this->points->Delete();
}


int
Points::read_file(const double chord)
{
    FILE *fp = fopen ("../foils/naca0015.pro", "r");
    struct tagbstring tokens; 
    struct bstrList* parsed = NULL;
    bstring text;
    int warnings = 0;
    unsigned char* x = NULL;
    unsigned char* y = NULL;
    const double cp_x = 0.25;
    int i = 0;

    cstr2tbstr(tokens,"\t\n\r "); /* token for splitting text into indivdual words is a tab and space */   
    for (warnings=0 ; NULL!=(text=bgets((bNgetc)fgetc, fp, (char)'\n')) ; bdestroy(text)) {
        parsed = bsplits(text, &tokens);
        x = NULL;
        y = NULL;
        for (i=0 ; i<parsed->qty ; i++) {
            if (x==NULL && Points::is_numeric((const char*)parsed->entry[i]->data) && parsed->entry[i]->slen) {
                x = parsed->entry[i]->data;                
            } else if (y==NULL && Points::is_numeric((const char*)parsed->entry[i]->data) && parsed->entry[i]->slen) {
                y = parsed->entry[i]->data;
            }
        }
        if (x && y) {
            this->x.push_back(chord*atof((const char*)x)-cp_x);
            this->y.push_back(chord*atof((const char*)y));
        } else if(x || y) {
            /* print warning, line is skipped; */
        }
        bstrListDestroy(parsed);
    }
    fclose (fp);
    return 0;
}


int
Points::add_points()
{
    int i = 0;
    const int n = this->x.size();
    double p0[3] = {0.0, 0.0, 0.0};
    for(i=0 ; i!=n ; i++) {
        p0[0] = this->x[i];
        p0[1] = this->y[i];
        this->points->InsertNextPoint(p0); 
    }
    return 0;
}
