CXX ?= g++
DEBUG ?= -g
CPPFLAGS ?= -I../src
CXXFLAGS ?= -Wall $(DEBUG) -O0
LDFLAGS ?=
LDLIBS ?=

TARGET := test_simulation_config
SOURCES := test_simulation_config.cpp ../src/simulation_config.cpp
OBJECTS := $(SOURCES:.cpp=.o)
INVARIANT_TARGET := test_geometry_invariants

all: $(TARGET)

.PHONY: smoke_run regression_run invariants_run capacity_policy_run open_shell_run test_all

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(INVARIANT_TARGET): test_geometry_invariants.cpp geometry_invariants.cpp ../src/tri_tri_intersect.cpp ../src/geometry.cpp ../src/montecarlo.cpp ../src/simulation_config.cpp ../src/simulation_setup.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -o $@

smoke_run:
	$(MAKE) -C ../src
	bash ./smoke_restart_run.sh ../src/assemble ./fixtures/restart_lammps.dat

invariants_run: $(INVARIANT_TARGET)
	./$(INVARIANT_TARGET)

regression_run:
	$(MAKE) -C ../src
	bash ./regression_reference.sh ../src/assemble ./fixtures/restart_lammps.dat

capacity_policy_run:
	$(MAKE) -C ../src
	bash ./index_capacity_policy.sh ../src/assemble ./fixtures/restart_lammps.dat

open_shell_run:
	$(MAKE) -C ../src
	bash ./open_shell_max_sweeps.sh ../src/assemble ./fixtures/open_shell_restart_lammps.dat

test_all: $(TARGET) $(INVARIANT_TARGET)
	./$(TARGET)
	./$(INVARIANT_TARGET)
	$(MAKE) capacity_policy_run
	$(MAKE) open_shell_run
	$(MAKE) smoke_run
	$(MAKE) regression_run

clean:
	rm -f $(OBJECTS) $(TARGET) $(INVARIANT_TARGET)
