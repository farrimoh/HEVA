CC = g++
DEBUG = #-g 
GSLINC=-I/include`
GSL=-L/lib
CFLAGS = -Wall $(DEBUG) -O3 -ffast-math -funroll-loops
LFLAGS = -Wall $(DEBUG) $(GSL) -lgsl -lgslcblas -lm


assembly:  tri_tri_intersect.cpp tri_tri_intersect.hpp geometry.cpp geometry.hpp montecarlo.cpp montecarlo.hpp run_assembly.cpp
	
	$(CC) $(CFLAGS) -c tri_tri_intersect.cpp
	$(CC) $(CFLAGS) -c geometry.cpp
	$(CC) $(CFLAGS) -c montecarlo.cpp
	$(CC) $(CFLAGS) -c run_assembly.cpp
	$(CC) $(CFLAGS) tri_tri_intersect.o geometry.o montecarlo.o run_assembly.o -o assemble $(LFLAGS)



clean:
	rm -fr *.o
