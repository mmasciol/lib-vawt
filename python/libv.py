'''
   Copyright (C) 2016 mdm                                    
                                                             
 Licensed to the Apache Software Foundation (ASF) under one  
 or more contributor license agreements.  See the NOTICE file
 distributed with this work for additional information       
 regarding copyright ownership.  The ASF licenses this file  
 to you under the Apache License, Version 2.0 (the           
 "License"); you may not use this file except in compliance  
 with the License.  You may obtain a copy of the License at  
                                                             
   http://www.apache.org/licenses/LICENSE-2.0                
                                                             
 Unless required by applicable law or agreed to in writing,  
 software distributed under the License is distributed on an 
 "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY      
 KIND, either express or implied.  See the License for the   
 specific language governing permissions and limitations           
 under the License.                                            
'''


import sys
from ctypes import *
import numpy as np
import os
from scipy.interpolate import UnivariateSpline, RectBivariateSpline
import math
np.set_printoptions(suppress=True)


_lib = cdll.LoadLibrary('../src/libv-0.10.00.so')                               
BSTRING_SIZE = c_size_t(0)


class Domain(Structure):
    _fields_ = []

    
class Environment(Structure):
    def __init__(self, obj, name=""):
        self.obj = obj
    _fields_ = [("V",POINTER(c_double)),
                ("ref_h",c_double),
                ("u_inf",c_double),
                ("rho",c_double),
                ("mu",c_double),
                ("shear",c_double)]


class ACWork(Structure):
    def __init__(self, obj, name=""):
        self.obj = obj
    _fields_ = [("w",POINTER(c_double)),
                ("wx",POINTER(c_double)),
                ("wy",POINTER(c_double)),
                ("Vt",POINTER(c_double)),
                ("Vn",POINTER(c_double)),
                ("cn",POINTER(c_double)),
                ("ct",POINTER(c_double)),
                ("Pn",POINTER(c_double)),
                ("integrand",POINTER(c_double)),
                ("cl",POINTER(c_double)),
                ("cd",POINTER(c_double)),
                ("residual",POINTER(c_double)),
                ("W",POINTER(c_double)),
                ("phi",POINTER(c_double)),
                ("alpha",POINTER(c_double)),
                ("Re",POINTER(c_double))]


class ACPerf(Structure):
    def __init__(self, obj, name=""):
        self.obj = obj
    _fields_ = [("fx",POINTER(c_double)),
                ("fy",POINTER(c_double)),
                ("fz",POINTER(c_double)),
                ("Q",POINTER(c_double)),
                ("Cp",c_double),
                ("Ct",c_double),
                ("Cx",c_double),
                ("Cy",c_double),
                ("Cz",c_double),
                ("n",c_int)]
                 

class ActuatorCylinder(Structure):
    def __init__(self, obj, name=""):
        self.obj = obj
    _fields_ = [("pc",c_void_p),
                ("work",POINTER(ACWork)),
                ("perf",POINTER(ACPerf)),
                ("blade",c_void_p),
                ("env",POINTER(Environment)),
                ("Am",POINTER(c_double)),
                ("theta",POINTER(c_double)),
                ("z",POINTER(c_double)),
                ("fr",POINTER(c_double)),
                ("ft ",POINTER(c_double)),
                ("fz",POINTER(c_double)),
                ("omega",c_double),
                ("the0",c_double),
                ("swept_area",c_double),
                ("nt",c_int),
                ("nz",c_int),
                ("nb",c_int),
                ("index",c_int),
                ("dstall", c_bool),
                ("msg", c_void_p),
                ("ierr", c_void_p)]
    

class FoilAerodynamics(Structure):
    def __init__(self, obj, name=""):
        self.obj = obj
    _fields_ = [("Re", c_double),
                ("aol", c_double),
                ("cna", c_double),
                ("cns", c_double),
                ("cnsl", c_double),
                ("aod", c_double),
                ("cd0", c_double),
                ("alpha", POINTER(c_double)),
                ("cl", POINTER(c_double)),
                ("cd", POINTER(c_double)),
                ("cm", POINTER(c_double)),
                ("len", c_int), 
                ("type", c_char*32)]


    
class Env(object):
    ierr = c_int(0)
    status = create_string_buffer(1024)
    env = None

    
    _lib.C_new_environment.argtype = [POINTER(Domain), c_char_p,  POINTER(c_int)]    
    _lib.C_new_environment.restype = POINTER(Environment)
    def __init__(self, domain):
        self.env = _lib.C_new_environment(domain, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        

    def __setattr__(self, attr, val):
        if attr=="env":
            super(Env, self).__setattr__(attr,val)
        else:
            if hasattr(self.env.contents, attr):
                setattr(self.env.contents, attr, val)

        
    # _lib.C_set_wind_velocity.argtypes = [POINTER(Environment), c_double, c_char_p, POINTER(c_int)]
    # def set_wind_velocity(self, u):
    #     _lib.C_set_wind_velocity(self.env, c_double(u), self.status, pointer(self.ierr))
    #     if self.ierr.value != 0 :
    #         print self.status.value        
    # 
    # _lib.C_set_ref_height.argtypes = [POINTER(Environment), c_double, c_char_p, POINTER(c_int)]
    # def set_ref_height(self, ref):
    #     _lib.C_set_ref_height(self.env, c_double(ref), self.status, pointer(self.ierr))
    #     if self.ierr.value != 0 :
    #         print self.status.value        
    #         
    #         
    # _lib.C_set_rho.argtypes = [POINTER(Environment), c_double, c_char_p, POINTER(c_int)]
    # def set_rho(self, rho):
    #     _lib.C_set_rho(self.env, c_double(rho), self.status, pointer(self.ierr))
    #     if self.ierr.value != 0 :
    #         print self.status.value        
    # 
    # 
    # _lib.C_set_mu.argtypes = [POINTER(Environment), c_double, c_char_p, POINTER(c_int)]
    # def set_mu(self, mu):
    #     _lib.C_set_mu(self.env, c_double(mu), self.status, pointer(self.ierr))
    #     if self.ierr.value != 0 :
    #         print self.status.value        
    # 
    # 
    # _lib.C_set_shear.argtypes = [POINTER(Environment), c_double, c_char_p, POINTER(c_int)]
    # def set_shear(self, shear):
    #     _lib.C_set_shear(self.env, c_double(shear), self.status, pointer(self.ierr))
    #     if self.ierr.value != 0 :
    #         print self.status.value        
    # 
    #         
    def env_get(self):
        return self.env


            
    
class Airfoil:    
    ierr = c_int(0)
    status = create_string_buffer(1024)
    af = None

    
    _lib.C_new_airfoil.argtype = [POINTER(Domain), c_char_p, c_char_p, c_char_p,  POINTER(c_int)]
    _lib.C_new_airfoil.restype = POINTER(FoilAerodynamics)
    def __init__(self, domain, open_this_file=None, type_name="None"):
        self.af = POINTER(FoilAerodynamics)
        if open_this_file!=None:
            foil_name = POINTER(c_char*255)
            passed_file_name = POINTER(c_char*255)
            passed_file_name = open_this_file+'\0'
            if type_name==None:
                foil_name = "9999" + "\0"
            else:
                foil_name = str(type_name) + "\0"
            self.af = _lib.C_new_airfoil(domain, passed_file_name, type_name, self.status, pointer(self.ierr))
            if self.ierr.value != 0 :
                print self.status.value        

    def __getattr__(self, attr):
        if hasattr(self.af.contents, attr):
            var = getattr(self.af.contents, attr)
            n = getattr(self.af.contents, "len")
        if type(var)==POINTER(c_double):
            return [var[i] for i in xrange(n)]
        else:
            return var


    _lib.C_extrapolate_airfoil.argtypes = [POINTER(FoilAerodynamics), c_double, c_double, c_double, c_int, c_char_p, POINTER(c_int)]
    def extrapolate(self, cd_max, aspect_ratio=-999.9, cd_min=0.001, n_alpha=15):
        #af = self.get_airfoil(index=index, name=name)
        _lib.C_extrapolate_airfoil(self.af, c_double(cd_max), c_double(aspect_ratio), c_double(cd_min), c_int(n_alpha), self.status, pointer(self.ierr)) 
        if self.ierr.value != 0 :
            print self.status.value        

            
    def af_get(self):
        return self.af


    _lib.C_set_blend_airfoils.argtype = [POINTER(Domain), POINTER(FoilAerodynamics), POINTER(FoilAerodynamics), c_double, c_char_p, c_char_p, POINTER(c_int)]
    _lib.C_set_blend_airfoils.restype = POINTER(FoilAerodynamics)
    def blend(self, domain, af1, af2, weight, name):
        self.af = _lib.C_set_blend_airfoils(domain, af1, af2, c_double(weight), c_char_p(name), self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        

        
    def get_data(self, name):
        var = getattr(self.af.contents, name)
        if type(var)==POINTER(c_double):
            n = getattr(self.af.contents, "len")
            return [var[i] for i in xrange(n)]
        else:
            return var


        
class Turbine(object):
    ierr = c_int(0)
    status = create_string_buffer(1024)
    ac = POINTER(ActuatorCylinder)
    perf = POINTER(ACPerf)
    work = POINTER(ACWork)

    _lib.C_new_actuator_cylinder.argtype = [POINTER(Domain), c_double, c_char_p,  POINTER(c_int)]    
    _lib.C_new_actuator_cylinder.restype = POINTER(ActuatorCylinder)
    def __init__(self, domain, nb, z, nt):
        arr = (c_double*len(z))(*z)
        self.ac = _lib.C_new_actuator_cylinder(domain, c_int(nb), c_int(nt), arr, c_int(len(z)), self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value

            
    def __setattr__(self, attr, val):
        if attr=="ac":
            super(Turbine, self).__setattr__(attr,val)
        elif attr=="perf":
            super(Turbine, self).__setattr__(attr,val)
        elif attr=="work":
            super(Turbine, self).__setattr__(attr,val)            
        else:
            if hasattr(self.ac.contents, attr):
                setattr(self.ac.contents, attr, val)
            elif hasattr(self.perf.contents, attr):
                setattr(self.perf.contents, attr, val)
            elif hasattr(self.work.contents, attr):
                setattr(self.work.contents, attr, val)


    def __getattr__(self, attr):
        self.work = getattr(self.ac.contents, "work")
        if hasattr(self.perf.contents, attr):
            var = getattr(self.perf.contents, attr)
            n = getattr(self.perf.contents, "n")
        elif hasattr(self.ac.contents, attr):
            var = getattr(self.ac.contents, attr)
            n = getattr(self.ac.contents, "nt")
        elif hasattr(self.work.contents, attr):
            var = getattr(self.work.contents, attr)
            n = getattr(self.ac.contents, "nt")
        if type(var)==POINTER(c_double):
            return [var[i] for i in xrange(n)]
        else:
            return var

            
    _lib.C_ac_commit.argtypes = [POINTER(ActuatorCylinder), POINTER(Environment), c_char_p, POINTER(c_int)]
    def commit(self, env):
        _lib.C_ac_commit(self.ac, env.env_get(), self.status, pointer(self.ierr))        
        if self.ierr.value != 0 :
            print self.status.value        


    _lib.C_set_dynamic_stall.argtypes = [POINTER(ActuatorCylinder), c_bool, c_char_p, POINTER(c_int)]
    def set_dstall(self, flag):        
        _lib.C_set_dynamic_stall(self.ac, c_bool(flag), self.status, pointer(self.ierr))        
        if self.ierr.value != 0 :
            print self.status.value        


    _lib.C_ac_solve.argtypes = [POINTER(ActuatorCylinder), c_char_p, POINTER(c_int)]
    def solve(self):
        _lib.C_ac_solve(self.ac, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        
        #var = getattr(self.ac.contents, "CP")
        #return var

    
    _lib.C_ac_instantaneous_force.argtypes = [POINTER(ActuatorCylinder), c_char_p, POINTER(c_int)]
    def force(self):
        _lib.C_ac_instantaneous_force(self.ac, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        
        #var = getattr(self.ac.contents, "CP")
        #return var
    

    _lib.C_ac_performance.argtypes = [POINTER(ActuatorCylinder), c_char_p, POINTER(c_int)]
    def performance(self):
        _lib.C_ac_performance(self.ac, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        
        self.perf = getattr(self.ac.contents, "perf")
        #var = getattr(perf.contents, "Cp")
        #return var

        
    # def get_data(self, name, key=None):
    #     var = getattr(self.perf.contents, name)
    #     if type(var)==POINTER(c_double):
    #         n = getattr(self.perf.contents, "n")
    #         return [var[i] for i in xrange(n)]
    #     else:
    #         return var


    # _lib.C_set_angular_velocity.argtypes = [POINTER(ActuatorCylinder), c_double, c_char_p, POINTER(c_int)]
    # def set_angular_velocity(self, omega):
    #     _lib.C_set_angular_velocity(self.ac, c_double(omega), self.status, pointer(self.ierr))
    #     if self.ierr.value != 0 :
    #         print self.status.value        

            
    _lib.C_set_initial_theta.argtypes = [POINTER(ActuatorCylinder), c_double, c_char_p, POINTER(c_int)]
    def set_initial_theta(self, the0):
        _lib.C_set_initial_theta(self.ac, c_double(the0), self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        
            

    _lib.C_set_airfoil_geometry.argtypes = [POINTER(ActuatorCylinder), POINTER(FoilAerodynamics), c_int, c_int, c_char_p, POINTER(c_int)]
    _lib.C_set_airfoil_chord.argtypes = [POINTER(ActuatorCylinder), c_double, c_int, c_int, c_char_p, POINTER(c_int)]
    _lib.C_set_airfoil_radius.argtypes = [POINTER(ActuatorCylinder), c_double, c_int, c_int, c_char_p, POINTER(c_int)]
    _lib.C_set_airfoil_twist.argtypes = [POINTER(ActuatorCylinder), c_double, c_int, c_int, c_char_p, POINTER(c_int)]
    _lib.C_set_airfoil_swept_offset.argtypes = [POINTER(ActuatorCylinder), c_double, c_int, c_int, c_char_p, POINTER(c_int)]
    def set_blade_foils(self, blade_num, af, chord, r, twist, swept):
        for i in xrange(len(af)):
            _lib.C_set_airfoil_geometry(self.ac, af[i], c_int(blade_num), c_int(i), self.status, pointer(self.ierr))            
            if self.ierr.value != 0 :
                print self.status.value        
            _lib.C_set_airfoil_chord(self.ac, c_double(chord[i]), c_int(blade_num), c_int(i), self.status, pointer(self.ierr))            
            if self.ierr.value != 0 :
                print self.status.value        
            _lib.C_set_airfoil_radius(self.ac, c_double(r[i]), c_int(blade_num), c_int(i), self.status, pointer(self.ierr))            
            if self.ierr.value != 0 :
                print self.status.value        
            _lib.C_set_airfoil_twist(self.ac, c_double(twist[i]), c_int(blade_num), c_int(i), self.status, pointer(self.ierr))            
            if self.ierr.value != 0 :
                print self.status.value        
            _lib.C_set_airfoil_swept_offset(self.ac, c_double(swept[i]), c_int(blade_num), c_int(i), self.status, pointer(self.ierr))            
            if self.ierr.value != 0 :
                print self.status.value        

                
    _lib.C_ac_get_stationary_performance.argtypes = [POINTER(ActuatorCylinder), c_double, c_char_p,  POINTER(c_int)]    
    _lib.C_ac_get_stationary_performance.restype = c_double
    def get_stationary_q(self, the):
        val = _lib.C_ac_get_stationary_performance(self.ac, c_double(the), self.status, pointer(self.ierr))            
        if self.ierr.value != 0 :
            print self.status.value        
        return val





class libv(object):
    ierr = c_int(0)
    status = create_string_buffer(1024)
    domain = POINTER(Domain)
    ac = []
    af = []
    env = []
    ren = None

    
    _lib.C_new_domain.argtype = [c_char_p,  POINTER(c_int)]
    _lib.C_new_domain.restype = POINTER(Domain)
    def __init__(self):
        self.domain = _lib.C_new_domain(self.status, pointer(self.ierr)) 
        if self.ierr.value != 0 :
            print self.status.value        


    def airfoil(self, fname, name=None):
        foil = Airfoil(self.domain, fname, name)
        self.af.append(foil)
        n = len(self.af)
        return self.af[n-1]

    
    def environment(self):
        env = Env(self.domain)
        self.env.append(env)
        n = len(self.env)
        return self.env[n-1]
    
        
    def vawt(self, nb, z, nt):
        actuator = Turbine(self.domain, nb, z, nt)
        self.ac.append(actuator)
        n = len(self.ac)
        return self.ac[n-1]

    
    _lib.C_set_airfoil_geometry.argtype = [c_void_p, POINTER(FoilAerodynamics), POINTER(c_char*255), c_char_p,  POINTER(c_int)]    
    def set_airfoil_geometry(self, geom_file, index=None, name=None):
        af = self.get_airfoil(index=index, name=name)
        f = POINTER(c_char*255)
        f = geom_file+'\0'
        _lib.C_set_airfoil_geometry(af, f, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        

            
    _lib.C_end_program.argtypes = [POINTER(Domain), c_char_p, POINTER(c_int)]
    def end(self):
        _lib.C_end_program(self.domain, self.status, pointer(self.ierr)) 
        if self.ierr.value != 0 :
            print self.status.value


    def get_airfoil(self, af=None, index=None, name=None):
        airfoil = POINTER(FoilAerodynamics)
        if index!=None:
            airfoil = _lib.C_get_airfoil(self.domain, c_int(index), self.status, pointer(self.ierr)) 
        elif name!=None:
            count = 0
            while 1:                
                airfoil = _lib.C_get_airfoil(self.domain, c_int(count), self.status, pointer(self.ierr))
                if airfoil.contents.type==name:
                    break
                count += 1
        else:
            return None
        return airfoil
    
            
    _lib.C_get_airfoil.argtype = [POINTER(Domain), c_int, c_char_p,  POINTER(c_int)]
    _lib.C_get_airfoil.restype = POINTER(FoilAerodynamics)
    def get_airfoil_data(self, var, index=None, name=None):
        airfoil = self.get_airfoil(index, name)
        N = airfoil.contents.len
        arr = getattr(airfoil.contents, var)
        return [arr[i] for i in xrange(N)]


    def airfoil_blend(self, weight, af1, af2, name="None"):
        foil = Airfoil(self.domain)
        foil.blend(self.domain, af1.af_get(), af2.af_get(), weight, name)
        self.af.append(foil)
        n = len(self.af)
        return self.af[n-1]

    
    _lib.C_set_renderer_engine.argtypes = [c_char_p, POINTER(c_int)]
    _lib.C_set_renderer_engine.restypes = c_void_p
    def render_engine(self):
        self.ren = c_void_p(_lib.C_set_renderer_engine(self.status, pointer(self.ierr)))
        if self.ierr.value != 0 :
            print self.status.value
            

    _lib.C_vawt_render.argtypes = [c_void_p, POINTER(Domain), c_char_p, POINTER(c_int)]
    def render(self):
        _lib.C_vawt_render(self.ren, self.domain, self.status, pointer(self.ierr))
        if self.ierr.value != 0 :
            print self.status.value        
