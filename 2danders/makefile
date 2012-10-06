CXX = icpc
# Flags for optimized execs
CXXFLAGS = 
# Flags for debugging
# CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG

LIBS = -lm -larmadillo -fopenmp

default: wave

wave.o: wave.cpp
	${CXX} -c -O3 wave.cpp

wave: wave.o
	${CXX} -o wave wave.o
	
clean:
	rm wave.o
	rm wave
