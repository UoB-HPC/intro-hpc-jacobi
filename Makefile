CC = cc
CFLAGS = -std=c99 -Wall
LDFLAGS = -lm

jacobi: jacobi.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
