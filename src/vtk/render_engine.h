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


class vtkRenderer;
class vtkRenderWindow;
class vtkAxesActor;
class vtkRenderWindowInteractor;
class vtkInteractorStyleTrackballCamera;
class vtkActor;
class vtkPoints;
class vtkPolyData;
class vtkLine;
class vtkCellArray;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkActor;


class VawtVTK
{
private:    
    vtkRenderer* ren;
    vtkRenderWindow* renWin;
    vtkAxesActor* axes;
    vtkRenderWindowInteractor* iren;
    vtkInteractorStyleTrackballCamera* style;
    vtkPoints* surface_points;
    vtkCellArray* surface_lines;
    vtkLine* surface_line;
    vtkPolyData* surface_linesPolyData;
    vtkPolyDataMapper* surface_mapper;
    vtkActor* surface_actor;
public:
    VawtVTK();
    ~VawtVTK();
    void add_actor(vtkActor* act);
    void set_background(const double r, const double g, const double b);
    void start();
    void make_surface_mesh();

};
