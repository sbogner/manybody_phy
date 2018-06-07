CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo


all: run

run: pairing_model.o jacobi.o
	${CC} ${CFLAGS} pairing_model.o jacobi.o -o run ${LIB}

pairing_model.o: pairing_model.cpp
	${CC} ${CFLAGS} -c pairing_model.cpp ${LIB}

jacobi.o: jacobi.cpp
	${CC} ${CFLAGS} -c jacobi.cpp ${LIB}

clean: 
	rm run *.o *.dat
