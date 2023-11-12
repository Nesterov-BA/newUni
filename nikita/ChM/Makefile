prog: main.o fun.o
		gcc main.o fun.o -o prog

main.o: main.c fun.h
		gcc -g -Wall -Wextra -Werror -pedantic -std=c99 -c main.c

fun.o: fun.c fun.h
		gcc -g -Wall -Wextra -Werror -pedantic -std=c99 -c fun.c

clean:
		rm -f *.o prog