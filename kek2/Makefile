prog: main.o fun.o
	gcc main.o fun.o -o prog -lm

main.o: main.c fun.h	
	gcc -g -Wall -Wextra -Werror -pedantic -std=c99 -c main.c -lm

fun.o: fun.c fun.h
	gcc -g -Wall -Wextra -Werror -pedantic -std=c99 -c fun.c -lm

clean:
	rm -f *.o prog