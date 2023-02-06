all: a.out clean

main.o: main.cpp
	g++ -c main.cpp -pthread 

a.out: main.o
	g++  main.o -pthread 

clean: main.o
	rm main.o
