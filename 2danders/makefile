CXX = g++
# Flags for optimized execs
CXXFLAGS = 
# Flags for debugging
# CXXFLAGS = -ftemplate-depth-30 -g -DBZ_DEBUG
OPTFLAG = -O3
LIBS = -lm -larmadillo -fopenmp

EXEC = wave
SRC = wave.cpp
OBJ = wave.o

$(EXEC): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OPTFLAG) -o $(EXEC) $(OBJ) $(LIBS)

$(OBJ):
	$(CXX) $(CXXFLAGS) $(OPTFLAG) -c $(SRC) #-fopenmp
	
clean:
	rm *.o
