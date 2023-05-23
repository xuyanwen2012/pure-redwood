CC = gcc
CFLAGS = -Wall -Wextra 
DEBUG_FLAG = -g -DREDWOOD_DEBUG

all: main

debug: main.c
	$(CC) $(CFLAGS) $(DEBUG_FLAG) main.c -lm -o debug.out

barnes: barnes.c
	$(CC) $(CFLAGS) $(DEBUG_FLAG) -O3 barnes.c -lm

nn: nn.c
	$(CC) $(CFLAGS) -O3 nn.c -lm

clean:
	rm -f a.out main.o
