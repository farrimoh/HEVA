#ifndef SIMULATION_SETUP_HPP
#define SIMULATION_SETUP_HPP

#include "simulation_config.hpp"

class geometry;

void apply_simulation_config(const SimulationConfig &config, geometry &g);

#endif
