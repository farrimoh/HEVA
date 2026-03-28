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
    double theta2;
    double theta3;
    double l0;
    double l1;
    double phi33;
    double phi12;
    double phi01;
    double phi20;
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
    double dmud;
};

struct CoreParameterConfig
{
    bool enabled;
    unsigned long maxBonds;
    double epsilonLJ;
    double sigmaLJ;
};

struct InitializationConfig
{
    std::string mode;
    std::string seedConfig;
    std::string path;
};

struct RuntimeConfig
{
    unsigned long seed;
    unsigned long maxSweeps;
    std::string outputDir;
    std::string workflow;
    bool resume;
};

struct EngineConfig
{
    unsigned long indexCapacity;
    std::string runMode;
};

struct CGParamOptConfig
{
    unsigned long sampleStartSweep;
    unsigned long sampleEvery;
    std::string sampleOutputPath;
};

struct SimulationConfig
{
    CapsidGeometryConfig capsidGeometry;
    SimulationParameterConfig simulation;
    DrugParameterConfig drug;
    CoreParameterConfig core;
    InitializationConfig initialization;
    RuntimeConfig runtime;
    EngineConfig engine;
    CGParamOptConfig cgParamOpt;
};

std::string assemble_usage();
bool parse_simulation_config(int argc, char **argv, SimulationConfig &config, std::string &error_message);
std::string render_simulation_config(const SimulationConfig &config);

#endif
