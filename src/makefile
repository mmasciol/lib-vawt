#   Copyright (C) 2015 mdm                                     
#                                                              
# Licensed to the Apache Software Foundation (ASF) under one   
# or more contributor license agreements.  See the NOTICE file 
# distributed with this work for additional information        
# regarding copyright ownership.  The ASF licenses this file   
# to you under the Apache License, Version 2.0 (the            
# "License"); you may not use this file except in compliance   
# with the License.  You may obtain a copy of the License at   
#                                                              
#   http://www.apache.org/licenses/LICENSE-2.0                 
#                                                              
# Unless required by applicable law or agreed to in writing,   
# software distributed under the License is distributed on an  
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       
# KIND, either express or implied.  See the License for the    
# specific language governing permissions and limitations            
# under the License.                                             

LIBV_VER:=0.10.00

BITS = -m64
PLATFORM = $(shell uname -s)
VPATH = cminpack:bstring:simclist:vtk


ifeq ($(OS),Windows_NT)
  DEL_CMD   = del
  LIB_EXT   = dll
  CFLAGS    = $(BITS) -g -std=c99 -DMAP_DLL_EXPORTS -DCMINPACK_NO_DLL  -DNDEBUG -D_WINDOWS -D_USRDLL -D_MINGW
  LDFLAGS   = $(BITS) -g -shared -Wl,--export-all-symbols
  LIB_FLAGS :=
else
  PLATFORM = $(shell uname -s)
  DEL_CMD   = rm -rf
  # GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)

  ifeq ($(PLATFORM),Darwin)
    LIB_EXT   = dylib
    CC_TOOLS  = clang
    CFLAGS    = -g -O3 -fno-omit-frame-pointer -fPIC -D DEBUG -Icminpack -Isimclist -Ivtk
    LDFLAGS   = -g -dynamiclib
    LIB_FLAGS := -lm
  else ifeq ($(PLATFORM),Linux)
    LIB_EXT   = so
    CC_TOOLS  = gcc
    LDFLAGS   = $(BITS) -g -O3 -shared -llapacke
    CFLAGS    = $(BITS) -g -O3 -fPIC -Wall -pedantic -std=c99 -D DEBUG -D WITH_LAPACK -D WITH_VTK -Icminpack
    LIB_FLAGS := -lm -llapacke /media/sf_Dropbox/offshore_wind/vawt/codes/v-lib/src/vtk/libvrender.so
    # CFLAGS    = $(BITS) -g -O3 -fPIC -Wall -pedantic -std=c99 -D DEBUG -D WITH_LAPACK -Icminpack
    # LIB_FLAGS := -lm -llapacke 
  endif
endif

DEBUG     = 
OBJ       = api.o \
	error.o \
	sys.o \
	domain.o \
	af.o \
	ac.o \
	env.o \
	blade.o \
	dstall.o \
	numerics.o \
	dpmpar.o \
	lmdif.o \
	hybrd.o \
	qform.o \
	r1updt.o \
	qrfac.o \
	fdjac2.o \
	dogleg.o \
	lmpar.o \
	qrsolv.o \
	enorm.o \
	fdjac1.o \
	r1mpyq.o \
	bstrlib.o \
	bstraux.o \
	simclist.o



all : $(OBJ)
	$(CC_TOOLS) $(LDFLAGS) -o libv-$(LIBV_VER).$(LIB_EXT) $(DEBUG) $(OBJ) $(LIB_FLAGS) 

.c.o :
	$(CC_TOOLS) -c $(CFLAGS) $<

clean:
	$(DEL_CMD) *.$(LIB_EXT) *.o *~

