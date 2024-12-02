all: a.out clean

main.o: main.cpp
	mpic++ -c main.cpp

a.out: main.o
	mpic++  main.o 

clean: main.o
	rm main.o
