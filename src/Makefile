CC = C:\w64devkit\bin\g++.exe
DEBUG = -g -v
GSLINC=-I C:\msys64\mingw64\include
GSL=-L C:\msys64\mingw64\lib
CFLAGS = -Wall $(DEBUG) -O3 -ffast-math -funroll-loops -mconsole
LFLAGS = -Wall $(DEBUG) $(GSL) -lgsl -lgslcblas -lm -v -mconsole


assembly:  tri_tri_intersect.cpp tri_tri_intersect.hpp geometry.cpp geometry.hpp montecarlo.cpp montecarlo.hpp run_assembly.cpp

	$(CC) $(CFLAGS) $(GSLINC)-c tri_tri_intersect.cpp
	$(CC) $(CFLAGS) -c geometry.cpp
	$(CC) $(CFLAGS) -c montecarlo.cpp
	$(CC) $(CFLAGS) -c run_assembly.cpp
	$(CC) -mconsole tri_tri_intersect.o geometry.o montecarlo.o run_assembly.o $(LFLAGS) -o assemble.exe



clean:
	del /F /Q *.o assemble.exe
