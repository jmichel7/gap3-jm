# This is the Makefile for Guava binaries

# Set this to the location for Guava binaries
BINDIR = ../bin

CC = gcc
CFLAGS = -O2 -Wall
FILES = ../bin/leonconv

$(BINDIR)/leonconv: leonconv.c
	$(CC) $(CFLAGS) -o $(BINDIR)/leonconv leonconv.c

clean :
	rm -f $(FILES)
