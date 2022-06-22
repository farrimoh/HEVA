CC = g++
DEBUG = #-g 
GSLINC=-I/include
GSL=-L/lib
CFLAGS = -Wall $(DEBUG) -O3 -ffast-math -funroll-loops
LFLAGS = -Wall $(DEBUG) $(GSL) -lgsl -lgslcblas -lm


drug:  tri_tri_intersect.cpp tri_tri_intersect.h geometry.cpp geometry.h montecarlo.cpp montecarlo.h run_assembly_drug.cpp
	
	$(CC) $(CFLAGS) -c tri_tri_intersect.cpp
	$(CC) $(CFLAGS) -c geometry.cpp
	$(CC) $(CFLAGS) -c montecarlo.cpp
	$(CC) $(CFLAGS) -c run_assembly.cpp
	$(CC) $(CFLAGS) tri_tri_intersect.o geometry.o montecarlo.o run_assembly.o -o drug $(LFLAGS)

clean:
	rm -fr *.o
