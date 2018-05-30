CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo


all: run

run: pair_fci.o jacobi.o
	${CC} ${CFLAGS} pair_fci.o jacobi.o -o run ${LIB}

pair_fci.o: pair_fci.cpp
	${CC} ${CFLAGS} -c pair_fci.cpp ${LIB}

jacobi.o: jacobi.cpp
	${CC} ${CFLAGS} -c jacobi.cpp ${LIB}

clean: 
	rm run *.o *.dat
