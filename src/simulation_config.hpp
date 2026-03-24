#ifndef SIMULATION_CONFIG_HPP
#define SIMULATION_CONFIG_HPP

#include <string>

struct SimulationConfig
{
    unsigned long seed;
    unsigned long indexCapacity;
    unsigned long maxSweeps;
    std::string restartPath;
    std::string outputDir;
    double epsilon0;
    double kappa0;
    double kappaPhi0;
    double theta0;
    double theta1;
    double gb0;
    double muCd;
    double ks0;
    double dmu;
    double dg;
    double mudrug;
    double gdrug0;
    double kd0;
    double dg12;
    double dg01;
    double dg20;
    double dg33;
    double dg00;
    double dgother;
};

std::string assemble_usage();
bool parse_simulation_config(int argc, char **argv, SimulationConfig &config, std::string &error_message);

#endif
