CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo


all: SRG MBPT2

SRG: pairing_model.o srg.o
	${CC} ${CFLAGS} pairing_model.o srg.o -o SRG ${LIB}

MBPT2: main.o system.o
	${CC} ${CFLAGS} main.o system.o -o MBPT2

main.o: main.cpp 
	${CC} ${CFLAGS} -c main.cpp

pairing_model.o: pairing_model.cpp
	${CC} ${CFLAGS} -c pairing_model.cpp ${LIB}

srg.o: srg.cpp
	${CC} ${CFLAGS} -c srg.cpp ${LIB}

jacobi.o: jacobi.cpp
	${CC} ${CFLAGS} -c jacobi.cpp ${LIB}

system.o: system.cpp
	${CC} ${CFLAGS} -c system.cpp

clean: 
	rm SRG MBPT2 
	rm *.o *.dat
