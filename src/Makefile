CXX ?= g++
DEBUG ?= -g
CPPFLAGS ?= -I/usr/include
CXXFLAGS ?= -Wall $(DEBUG) -O3 -ffast-math -funroll-loops
LDFLAGS ?= -L/usr/lib/x86_64-linux-gnu
LDLIBS ?= -lgsl -lgslcblas -lm

TARGET := assemble
SOURCES := tri_tri_intersect.cpp geometry.cpp montecarlo.cpp seed_initialization.cpp run_assembly.cpp simulation_config.cpp simulation_setup.cpp simulation_runner.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: assembly

assembly: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
