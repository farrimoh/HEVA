#include "simulation_config.hpp"

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>

namespace
{
struct FieldSpec
{
    int index;
    const char *name;
    double *value;
};

bool parse_unsigned_long(const char *text, unsigned long &value)
{
    if (text == nullptr || *text == '\0')
    {
        return false;
    }

    char *end = nullptr;
    errno = 0;
    unsigned long parsed = std::strtoul(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0')
    {
        return false;
    }

    value = parsed;
    return true;
}

bool parse_double_value(const char *text, double &value)
{
    if (text == nullptr || *text == '\0')
    {
        return false;
    }

    char *end = nullptr;
    errno = 0;
    double parsed = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !std::isfinite(parsed))
    {
        return false;
    }

    value = parsed;
    return true;
}
}

std::string assemble_usage()
{
    return "usage: ./assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 dmu dummydg mudrug gdrug kd0 dg12 dg01 dg20 dg33 dg00 dgother [--index-capacity N] [--max-sweeps N] [--restart PATH] [--output-dir PATH]";
}

bool parse_simulation_config(int argc, char **argv, SimulationConfig &config, std::string &error_message)
{
    if (argc < 21)
    {
        std::ostringstream message;
        message << "Expected at least 20 arguments after the executable name, got " << (argc - 1) << ".\n"
                << assemble_usage();
        error_message = message.str();
        return false;
    }

    config.indexCapacity = 1000000UL;
    config.maxSweeps = 0UL;
    config.restartPath = "restart_lammps.dat";
    config.outputDir = ".";

    if (!parse_unsigned_long(argv[1], config.seed))
    {
        error_message = "Invalid value for seed.\n" + assemble_usage();
        return false;
    }

    FieldSpec fields[] = {
        {2, "epsilon0", &config.epsilon0},
        {3, "kappa0", &config.kappa0},
        {4, "kappaPhi0", &config.kappaPhi0},
        {5, "theta0", &config.theta0},
        {6, "theta1", &config.theta1},
        {7, "LnK", &config.gb0},
        {8, "muCd", &config.muCd},
        {9, "ks0", &config.ks0},
        {10, "dmu", &config.dmu},
        {11, "dummydg", &config.dg},
        {12, "mudrug", &config.mudrug},
        {13, "gdrug", &config.gdrug0},
        {14, "kd0", &config.kd0},
        {15, "dg12", &config.dg12},
        {16, "dg01", &config.dg01},
        {17, "dg20", &config.dg20},
        {18, "dg33", &config.dg33},
        {19, "dg00", &config.dg00},
        {20, "dgother", &config.dgother},
    };

    const int field_count = sizeof(fields) / sizeof(fields[0]);
    for (int i = 0; i < field_count; ++i)
    {
        if (!parse_double_value(argv[fields[i].index], *fields[i].value))
        {
            error_message = std::string("Invalid value for ") + fields[i].name + ".\n" + assemble_usage();
            return false;
        }
    }

    for (int i = 21; i < argc; )
    {
        const std::string option = argv[i];
        if (option == "--index-capacity")
        {
            if (i + 1 >= argc || !parse_unsigned_long(argv[i + 1], config.indexCapacity) || config.indexCapacity == 0 || config.indexCapacity > static_cast<unsigned long>(INT_MAX))
            {
                error_message = "Invalid value for --index-capacity.\n" + assemble_usage();
                return false;
            }
            i += 2;
        }
        else if (option == "--max-sweeps")
        {
            if (i + 1 >= argc || !parse_unsigned_long(argv[i + 1], config.maxSweeps))
            {
                error_message = "Invalid value for --max-sweeps.\n" + assemble_usage();
                return false;
            }
            i += 2;
        }
        else if (option == "--restart")
        {
            if (i + 1 >= argc || argv[i + 1] == nullptr || *argv[i + 1] == '\0')
            {
                error_message = "Invalid value for --restart.\n" + assemble_usage();
                return false;
            }
            config.restartPath = argv[i + 1];
            i += 2;
        }
        else if (option == "--output-dir")
        {
            if (i + 1 >= argc || argv[i + 1] == nullptr || *argv[i + 1] == '\0')
            {
                error_message = "Invalid value for --output-dir.\n" + assemble_usage();
                return false;
            }
            config.outputDir = argv[i + 1];
            i += 2;
        }
        else
        {
            error_message = std::string("Unknown option: ") + option + ".\n" + assemble_usage();
            return false;
        }
    }

    error_message.clear();
    return true;
}
