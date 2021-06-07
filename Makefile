# Copyright (C) 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
# 2003, 2004, 2005, 2006, 2007, 2008, 2009  Free Software Foundation,
# Inc. This Makefile is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

CC =  gcc
CXX = g++
CFLAGS = -O2 -m64 -fopenmp
INCDIR = -I.
LDFLAGS = $(CFLAGS) -lm
TARGET  = ra_omp
OBJECTS = RandomAccess_omp.o verification.o

.SUFFIXES: .c

.c.o:
	$(CC) $(INCDIR) $(CFLAGS) -c $<

all:  RAOMP

RAOMP: $(OBJECTS)
	$(CC) $(INCDIR) $(CFLAGS) $(OBJECTS) -o $(TARGET) $(LDFLAGS)

clean:
	/bin/rm -f $(OBJECTS) core

cleanall:
	/bin/rm -f $(OBJECTS) $(TARGET) core

