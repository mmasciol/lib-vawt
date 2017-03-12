CPP     = g++
VTK_INC = /home/marco/vtk/VTK-Install/include/vtk-7.1/
VTK_LIB = /home/marco/vtk/VTK-Build/lib/

VPATH = ../bstring:../simclist
OPT     = -g -fPIC -O2 -Wall -Wno-long-long -pedantic
CFLAGS  = $(OPT) -I${VTK_INC}
LIBS    = -L${VTK_LIB} 
LFLAGS  = -lvtkCommonCore-7.1 \
	-lvtkRenderingCore-7.1 \
	-lvtkCommonDataModel-7.1 \
	-lvtkCommonTransforms-7.1 \
	-lvtkRenderingOpenGL2-7.1 \
	-lvtkInteractionStyle-7.1 \
	-lvtkRenderingFreeType-7.1 \
	-lvtkRenderingAnnotation-7.1 \

OBJ = render_api.o \
	render_engine.o \
	line_generator.o \
	bstrlib.o  

all : $(OBJ)
	$(CPP) $(LIBS) -shared -o libvrender.so $(OBJ) $(LFLAGS) 

.cpp.o :
	$(CPP) -c $(CFLAGS) $< 

clean:
	rm -rf *.o *~ *.so

memcheck:	
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-origins=yes ./test
