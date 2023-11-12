all: a.out clean

chebyshev.o: chebyshev.cpp
	g++ -c chebyshev.cpp

help.o: help.cpp
	g++ -c help.cpp

main.o: main.cpp
	g++ -c main.cpp

a.out: main.o chebyshev.o help.o
	g++  main.o chebyshev.o help.o 

clean: main.o chebyshev.o help.o
	rm main.o chebyshev.o help.o