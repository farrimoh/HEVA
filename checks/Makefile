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

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) $(LDLIBS) -o $@

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
