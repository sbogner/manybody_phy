CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo


all: run

run: pairing_model.o srg.o
	${CC} ${CFLAGS} pairing_model.o srg.o -o run ${LIB}

pairing_model.o: pairing_model.cpp
	${CC} ${CFLAGS} -c pairing_model.cpp ${LIB}

srg.o: srg.cpp
	${CC} ${CFLAGS} -c srg.cpp ${LIB}

jacobi.o: jacobi.cpp
	${CC} ${CFLAGS} -c jacobi.cpp ${LIB}

clean: 
	rm run *.o *.dat
