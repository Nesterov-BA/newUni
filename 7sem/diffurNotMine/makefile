all: a.out

p1.o: p1.c
	gcc -c p1.c

a.out: p1.o
	gcc p1.o -lm

clean:
	rm p1.o a.out