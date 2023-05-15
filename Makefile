CC = gcc
CFLAGS = -Wall -Wextra -g

main: main.o
	$(CC) $(CFLAGS) main.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f a.out main.o
