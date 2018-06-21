CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo


all: SRG MBPT2 TEST

SRG: pairing_model.o srg.o
	${CC} ${CFLAGS} pairing_model.o srg.o -o SRG ${LIB}

MBPT2: main.o system.o
	${CC} ${CFLAGS} main.o system.o -o MBPT2

TEST: test_main.o unit_tests.o system.o
	${CC} ${CFLAGS} test_main.o unit_tests.o system.o -o TEST 

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

test_main.o: test_main.cpp
	${CC} ${CFLAGS} -c test_main.cpp

unit_tests.o: unit_tests.cpp
	${CC} ${CFLAGS} -c unit_tests.cpp

clean: 
	rm SRG MBPT2 TEST
	rm *.o *.dat
