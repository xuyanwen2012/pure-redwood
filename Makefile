CC = gcc
CFLAGS = -O3 -Wall -Wextra -g 

main: main.o
	$(CC) $(CFLAGS) main.o -lm

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f a.out main.o
