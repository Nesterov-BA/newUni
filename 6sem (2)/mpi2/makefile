all: a.out clean

mpi2.o: mpi2.cpp
	mpic++ -c mpi2.cpp

a.out: mpi2.o
	mpic++  mpi2.o 

clean: mpi2.o
	rm mpi2.o
