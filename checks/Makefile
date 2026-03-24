CXX ?= g++
DEBUG ?= -g
CPPFLAGS ?= -I../src
CXXFLAGS ?= -Wall $(DEBUG) -O0
LDFLAGS ?=
LDLIBS ?=

TARGET := test_simulation_config
SOURCES := test_simulation_config.cpp ../src/simulation_config.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

.PHONY: smoke_run

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

smoke_run:
	$(MAKE) -C ../src
	bash ./smoke_restart_run.sh ../src/assemble ./fixtures/restart_lammps.dat

clean:
	rm -f $(OBJECTS) $(TARGET)
