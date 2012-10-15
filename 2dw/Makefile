## Basic Makefile for OpenMP
CC = g++
CFLAGS = -g -O3 -fopenmp
LD = g++
LDFLAGS = -g -fopenmp

PROGRAM = wave

all:	${PROGRAM}.out

${PROGRAM}.out: ${PROGRAM}.o
	${LD} ${LDFLAGS} $< -o ${PROGRAM}.out

${PROGRAM}.o: ${PROGRAM}.cpp
	${CC} ${CFLAGS} -c $< -o ${PROGRAM}.o

clean:
	rm -f ${PROGRAM}.o ${PROGRAM}.out
	rm -f *.o
	rm -f test.d0*
	rm -f _tmp*
