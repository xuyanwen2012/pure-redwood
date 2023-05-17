CC = gcc
CFLAGS = -Wall -Wextra 
DEBUG_FLAG = -g -DREDWOOD_DEBUG

all: main

debug: main.c
	$(CC) $(CFLAGS) $(DEBUG_FLAG) main.c -lm -o debug.out

main: main.c
	$(CC) $(CFLAGS) -O3 main.c -lm

clean:
	rm -f a.out main.o
