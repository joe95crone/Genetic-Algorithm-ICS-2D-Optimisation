# Makefile for SPEA2

# Compiler
CC = g++

# Compiler options
CFLAGS = -g -Wall

# all object files
SEL_OBJECTS = spea2_functions.o spea2_io.o spea2.o

spea2 : $(SEL_OBJECTS)
	$(CC) $(CFLAGS) -lm $(SEL_OBJECTS) -o spea2

spea2_functions.o : spea2_functions.cpp spea2.hpp
	$(CC) $(CFLAGS) -c spea2_functions.cpp

spea2_io.o : spea2_io.cpp spea2.hpp
	$(CC) $(CFLAGS) -c spea2_io.cpp

spea2.o : spea2.cpp spea2.hpp
	$(CC) $(CFLAGS) -c spea2.cpp

clean :
	rm -f *~ *.o

install : spea2
	cp spea2 $(HOME)/bin/