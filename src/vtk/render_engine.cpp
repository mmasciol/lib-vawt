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
#include "vtkRenderWindow.h"
#include "vtkAxesActor.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkCaptionActor2D.h"
#include "vtkTextProperty.h"
#include "vtkTextActor.h"
#include "vtkLine.h"
#include "vtkPolyData.h"
#include "vtkProperty.h"
#pragma GCC diagnostic pop

#include "render_engine.h"


VawtVTK::VawtVTK()
{
    this->ren = vtkRenderer::New();
    this->renWin = vtkRenderWindow::New();
    this->axes = vtkAxesActor::New();
    this->iren = vtkRenderWindowInteractor::New();
    this->style = vtkInteractorStyleTrackballCamera::New();
    this->renWin->AddRenderer(this->ren);
    this->renWin->SetSize(600, 600);
    this->make_surface_mesh();
}


void
VawtVTK::make_surface_mesh()
{
    const int step = 5;
    const int m = 20;
    const int n = 20;
    int i = 0;
    double p[3] = {0.0, 0.0, -2.0};
    int pn = 0;
    this->surface_points = vtkPoints::New();
    this->surface_line = vtkLine::New();
    this->surface_lines = vtkCellArray::New();
    this->surface_linesPolyData = vtkPolyData::New();
    this->surface_mapper = vtkPolyDataMapper::New();
    this->surface_actor = vtkActor::New();

    for(i=-n ; i<n+1 ; i+=step) {
        p[0] = i;
        p[1] = -(m+step);
        pn++;
        this->surface_points->InsertNextPoint(p);
        
        p[0] = i;
        p[1] = (m+step);
        pn++;
        this->surface_points->InsertNextPoint(p);
        this->surface_line->GetPointIds()->SetId(0,pn-2);
        this->surface_line->GetPointIds()->SetId(1,pn-1);
        this->surface_lines->InsertNextCell(this->surface_line);
    }

    for(i=-n ; i<n+1 ; i+=step) {
        p[1] = i;
        p[0] = -(m+step);
        pn++;
        this->surface_points->InsertNextPoint(p);
        
        p[1] = i;
        p[0] = (m+step);
        pn++;
        this->surface_points->InsertNextPoint(p);
        this->surface_line->GetPointIds()->SetId(0,pn-2);
        this->surface_line->GetPointIds()->SetId(1,pn-1);
        this->surface_lines->InsertNextCell(this->surface_line);
    }

    this->surface_linesPolyData->SetPoints(this->surface_points);
    this->surface_linesPolyData->SetLines(this->surface_lines);    
    
    this->surface_mapper->SetInputData(this->surface_linesPolyData);
    this->surface_actor->SetMapper(this->surface_mapper);
    this->add_actor(this->surface_actor);
    this->surface_actor->GetProperty()->SetColor(0.0, 1.0, 1.0);
    this->surface_actor->GetProperty()->LightingOff();
    this->surface_actor->GetProperty()->SetOpacity(0.25);
}


VawtVTK::~VawtVTK()
{
    this->surface_points->Delete();
    this->surface_line->Delete();
    this->surface_lines->Delete();
    this->surface_linesPolyData->Delete();
    this->ren->Delete();
    this->renWin->Delete();
    this->axes->Delete();
    this->iren->Delete();
    this->style->Delete();
    this->surface_mapper->Delete();
    this->surface_actor->Delete();
    //this->planeSource->Delete();
    //this->plane_mapper->Delete();
    //this->plane_actor->Delete();
}


void
VawtVTK::add_actor(vtkActor* act)
{
    this->ren->AddActor(act);
}


void
VawtVTK::set_background(const double r, const double g, const double b)
{
    this->ren->SetBackground(r,g,b);
}


void
VawtVTK::start()
{    
    // this->ren->AddActor(plane_actor);
    this->axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(1,0,0);
    this->axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,1,0);
    this->axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(0,0,1);
    this->axes->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
    this->axes->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
    this->axes->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
    this->axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(20);
    axes->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(20);
    axes->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(20);
    
    this->axes->SetXAxisLabelText("");
    this->axes->SetYAxisLabelText("");
    this->axes->SetZAxisLabelText("");

    this->ren->AddActor(axes);

    this->iren->SetRenderWindow(renWin);

    this->iren->SetInteractorStyle(style);
    this->iren->Start();

}
