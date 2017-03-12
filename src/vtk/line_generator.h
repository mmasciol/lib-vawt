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


#include <vector>


class vtkTransform;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkActor;
class vtkPoints;
class vtkCellArray;
class vtkLine;

class Points
{
private:
    vtkTransform* trans;
    vtkPolyData* linesPolyData;
    vtkPolyDataMapper* mapper;
    vtkActor* actor;
    vtkPoints* points;
    vtkCellArray* lines;
    vtkLine* line;
    std::vector<double> x;
    std::vector<double> y;
public:
    Points();
    ~Points();
    int read_file(const double chord);
    int add_points();
    int set_polydata();
    int set_mapper_actor();
    vtkActor* get_actor() {return this->actor;};
    int translate(const double x, const double y, const double z);
    int rotate(const double x, const double y, const double z);
    int transform();
    int set_colors();
    static bool is_numeric(const char* str);

};
