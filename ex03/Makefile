CFLAGS=-std=c++11
all: main

main: main.cpp functions
	g++ main.cpp $(CFLAGS) -o main.out functions.out

functions: functions.cpp
	g++ -c $(CFLAGS) functions.cpp -o functions.out
