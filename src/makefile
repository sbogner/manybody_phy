CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo


all: RUN TEST

RUN: main.o system.o solver.o srg.o
	${CC} ${CFLAGS} main.o system.o solver.o srg.o -o RUN ${LIB}

TEST: test_main.o unit_tests.o system.o
	${CC} ${CFLAGS} test_main.o unit_tests.o system.o -o TEST 



main.o: main.cpp 
	${CC} ${CFLAGS} -c main.cpp

srg.o: srg.cpp
	${CC} ${CFLAGS} -c srg.cpp ${LIB}

solver.o: solver.cpp
	${CC} ${CFLAGS} -c solver.cpp

system.o: system.cpp
	${CC} ${CFLAGS} -c system.cpp

test_main.o: test/test_main.cpp
	${CC} ${CFLAGS} -c test/test_main.cpp

unit_tests.o: test/unit_tests.cpp
	${CC} ${CFLAGS} -c test/unit_tests.cpp



clean: 
	rm *.o
	rm RUN TEST