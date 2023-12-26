all: a.out clean

main.o: main.cpp
	g++ -c main.cpp

a.out: main.o
	g++  main.o 

clean: main.o
	rm main.o
