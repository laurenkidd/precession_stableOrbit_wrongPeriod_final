CC = gcc
CFLAGS = -g -lm -Wall -std=c99
OBJECTS = prec_euler_zeros.o
all : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o prec_euler_zeros


prec_euler_zeros: .c
	$(CC) $(CFLAGS) -c prec_euler_zeros.c
