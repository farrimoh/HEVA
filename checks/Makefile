CXX ?= g++
DEBUG ?= -g
GSL_CFLAGS ?= $(shell pkg-config --cflags gsl 2>/dev/null)
GSL_LIBS ?= $(shell pkg-config --libs gsl 2>/dev/null)
CPPFLAGS ?= -I../src
CXXFLAGS ?= -Wall $(DEBUG) -O0
LDFLAGS ?=
LDLIBS ?= $(if $(strip $(GSL_LIBS)),$(GSL_LIBS),-lgsl -lgslcblas -lm)
DEPFLAGS ?= -MMD -MP

TARGET := test_simulation_config
SOURCES := test_simulation_config.cpp ../src/simulation_config.cpp
OBJECTS := $(SOURCES:.cpp=.o)
DEPS := $(OBJECTS:.o=.d)
INVARIANT_TARGET := test_geometry_invariants

all: $(TARGET)

.PHONY: smoke_run regression_run invariants_run capacity_policy_run open_shell_run relaxation_run seed_init_run test_all

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(INVARIANT_TARGET): test_geometry_invariants.cpp geometry_invariants.cpp ../src/tri_tri_intersect.cpp ../src/geometry.cpp ../src/montecarlo.cpp ../src/simulation_config.cpp ../src/simulation_setup.cpp
	$(CXX) $(CPPFLAGS) $(GSL_CFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

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

relaxation_run:
	$(MAKE) -C ../src
	bash ./relaxation_workflow_run.sh ../src/assemble ./fixtures/open_shell_restart_lammps.dat

seed_init_run:
	$(MAKE) -C ../src
	bash ./seed_init_modes.sh ../src/assemble

test_all: $(TARGET) $(INVARIANT_TARGET)
	./$(TARGET)
	./$(INVARIANT_TARGET)
	$(MAKE) capacity_policy_run
	$(MAKE) open_shell_run
	$(MAKE) relaxation_run
	$(MAKE) seed_init_run
	$(MAKE) smoke_run
	$(MAKE) regression_run

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET) $(INVARIANT_TARGET)

-include $(DEPS)
