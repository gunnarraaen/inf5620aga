VPATH	= .


# target might be WINDOWS, OS_X or LINUX

TARGET = OS_X

# Home directory of the acw program

# point to the inclue directory
CINCLUDE =  CiRIOBase/CiRIOBase/source
SOURCEDIR = source

INCLUDES = -I/sw/include -I$(CINCLUDE)  -I$(SOURCEDIR)
# library dir
LIBDIR =  CiRIOBase/libs


# compiler specific flags
CFLAGS =  -O3 -D$(TARGET) -fopenmp

FFLAGS = -lglew  -framework GLUT -framework OpenGL, -framework Carbon -lm -lCiRIOBase -larmadillo -fopenmp


PROJECT = main

_obj 	=  main.o CWave.o CPoint.o CPoints.o WaveSolver.o


obj = $(patsubst %,$(SOURCEDIR)/%,$(_obj))

CC 	= g++ -m64

default: $(PROJECT)

$(PROJECT):  $(obj) 
	$(CC)  -L$(LIBDIR)  $(INCLUDES) -o $(PROJECT) $(obj) $(FFLAGS)  

%.o: %.cpp
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   

%.o: %.cu
	$(CC) -c -o $@ $^ $(INCLUDES) $(CFLAGS)   


clean:	
	rm *.o source/*.o work	


