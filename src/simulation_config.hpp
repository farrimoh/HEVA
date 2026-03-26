#ifndef SIMULATION_CONFIG_HPP
#define SIMULATION_CONFIG_HPP

#include <string>

struct CapsidGeometryConfig
{
    double epsilon0;
    double kappa0;
    double kappaPhi0;
    double theta0;
    double theta1;
    double gb0;
    double dg12;
    double dg01;
    double dg20;
    double dg33;
    double dg00;
    double dgother;
};

struct SimulationParameterConfig
{
    double muCd;
    double ks0;
    double dmu;
    double dg;
};

struct DrugParameterConfig
{
    double mudrug;
    double gdrug0;
    double kd0;
};

struct InitializationConfig
{
    std::string mode;
    std::string seedConfig;
    std::string restartPath;
};

struct RuntimeConfig
{
    unsigned long seed;
    unsigned long maxSweeps;
    std::string outputDir;
    std::string workflow;
};

struct EngineConfig
{
    unsigned long indexCapacity;
    std::string runMode;
};

struct SimulationConfig
{
    CapsidGeometryConfig capsidGeometry;
    SimulationParameterConfig simulation;
    DrugParameterConfig drug;
    InitializationConfig initialization;
    RuntimeConfig runtime;
    EngineConfig engine;
};

std::string assemble_usage();
bool parse_simulation_config(int argc, char **argv, SimulationConfig &config, std::string &error_message);
std::string render_simulation_config(const SimulationConfig &config);

#endif
