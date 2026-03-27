CXX ?= g++
DEBUG ?= -g
GSL_CFLAGS ?= $(shell pkg-config --cflags gsl 2>/dev/null)
GSL_LIBS ?= $(shell pkg-config --libs gsl 2>/dev/null)
CPPFLAGS ?= $(GSL_CFLAGS)
CXXFLAGS ?= -Wall $(DEBUG) -O3 -ffast-math -funroll-loops
LDFLAGS ?=
LDLIBS ?= $(if $(strip $(GSL_LIBS)),$(GSL_LIBS),-lgsl -lgslcblas -lm)
DEPFLAGS ?= -MMD -MP

ASSEMBLE_TARGET := assemble
PREP_TARGET := prepare_initial_frame

ASSEMBLE_SOURCES := tri_tri_intersect.cpp geometry.cpp montecarlo.cpp run_assembly.cpp simulation_config.cpp simulation_setup.cpp simulation_runner.cpp
ASSEMBLE_OBJECTS := $(ASSEMBLE_SOURCES:.cpp=.o)
PREP_OBJECTS := tri_tri_intersect.o geometry.o prepare_initial_frame.o
OBJECTS := $(ASSEMBLE_OBJECTS) $(PREP_OBJECTS)
DEPS := $(OBJECTS:.o=.d)

all: assembly

assembly: $(ASSEMBLE_TARGET) $(PREP_TARGET)

$(ASSEMBLE_TARGET): $(ASSEMBLE_OBJECTS)
	$(CXX) $(CXXFLAGS) $(ASSEMBLE_OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

$(PREP_TARGET): $(PREP_OBJECTS)
	$(CXX) $(CXXFLAGS) $(PREP_OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(DEPS) $(ASSEMBLE_TARGET) $(PREP_TARGET)

-include $(DEPS)
