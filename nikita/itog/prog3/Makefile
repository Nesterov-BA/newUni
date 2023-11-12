all: main.o help.o task.o
	gcc main.o help.o task.o -lm -lpthread 

main.o: main.c help.h task.h
	gcc main.c -c -o main.o 

matrix.o: help.c help.h
	gcc help.c -c -o help.o 

task.o: task.c task.h
	gcc task.c -c -o task.o 
clean:
	rm -f *.o 
