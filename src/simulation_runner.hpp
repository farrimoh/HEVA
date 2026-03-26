#ifndef SIMULATION_RUNNER_HPP
#define SIMULATION_RUNNER_HPP

#include "simulation_config.hpp"

#include <gsl/gsl_rng.h>

#include <cstdio>
#include <ctime>
#include <string>

class geometry;

struct SimulationRunStats
{
    int frame;
    int monomeradded;
    int dimeradded;
    int monomerremoved;
    int dimerremoved;
    int drugadded;
    int drugremoved;
    int typechanged;
    int fusion;
    int fission;
    int wedgefusion;
    int wedgefission;
    int boundtri;
    int binding;
    int unbinding;
    int deletednorate;
    int lastNhe;
    int lastNheGrowth;
    int npace;
    double avgpace;
    int avgAddInterval;

    SimulationRunStats();
};

struct SimulationLoopSettings
{
    int minHEUpdateNeigh;
    int minhe_fission;
    int freq_vis;
    int freq_log;
    int freq_out;
    int initial_equilibration_steps;
    int final_equilibration_steps;
    bool relaxationSweepsAreRelative;
    unsigned long maxSweeps;
    unsigned long cgSampleStartSweep;
    unsigned long cgSampleEvery;
    std::string cgSampleOutputPath;

    SimulationLoopSettings();
};

enum SimulationStopReason
{
    SIMULATION_STOP_CLOSED = 0,
    SIMULATION_STOP_MAX_SWEEPS = 1,
    SIMULATION_STOP_OVERLAP_ERROR = 2,
    SIMULATION_STOP_MIXED_MORPH = 3,
    SIMULATION_STOP_STALLED_GROWTH = 4,
    SIMULATION_STOP_TOO_LONG = 5,
    SIMULATION_STOP_TOO_LARGE = 6
};

SimulationLoopSettings make_simulation_loop_settings(const SimulationConfig &config);
void initialize_from_restart(geometry &g, gsl_rng *rng, const char *filename, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings);
void initialize_from_initial_frame_compat(geometry &g, const char *filename, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings);
void initialize_from_seed(geometry &g, gsl_rng *rng, const char *seed_config, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings);
SimulationStopReason run_relaxation_loop(geometry &g, gsl_rng *rng, FILE *ofile, unsigned long seed, time_t start_time, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings);
SimulationStopReason run_simulation_loop(geometry &g, gsl_rng *rng, FILE *ofile, unsigned long seed, time_t start_time, unsigned long &sweep, double ks0, SimulationRunStats &stats, const SimulationLoopSettings &settings);
void finalize_simulation(geometry &g, gsl_rng *rng, FILE *ofile, unsigned long seed, time_t start_time, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings, SimulationStopReason stop_reason);

#endif
