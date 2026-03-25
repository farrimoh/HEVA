#include "simulation_config.hpp"

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct FieldSpec
{
    int index;
    const char *name;
    double *value;
};

struct ConfigFieldSpec
{
    const char *name;
    double *value;
    bool *was_set;
};

std::string trim(const std::string &text)
{
    std::size_t start = text.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
    {
        return "";
    }

    std::size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(start, end - start + 1);
}

bool is_absolute_path(const std::string &path)
{
    return !path.empty() && path[0] == '/';
}

std::string join_paths(const std::string &left, const std::string &right)
{
    if (left.empty() || left == ".")
    {
        return right;
    }
    if (left == "/")
    {
        return "/" + right;
    }
    if (left[left.size() - 1] == '/')
    {
        return left + right;
    }
    return left + "/" + right;
}

std::string parent_directory(const std::string &path)
{
    std::size_t slash = path.find_last_of('/');
    if (slash == std::string::npos)
    {
        return ".";
    }
    if (slash == 0)
    {
        return "/";
    }
    return path.substr(0, slash);
}

std::string resolve_relative_to(const std::string &base_dir, const std::string &path)
{
    if (path.empty() || is_absolute_path(path))
    {
        return path;
    }
    return join_paths(base_dir, path);
}

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

bool parse_unsigned_long_string(const std::string &text, unsigned long &value)
{
    return parse_unsigned_long(text.c_str(), value);
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

bool parse_double_string(const std::string &text, double &value)
{
    return parse_double_value(text.c_str(), value);
}

void set_default_config(SimulationConfig &config)
{
    config.indexCapacity = 1000000UL;
    config.maxSweeps = 0UL;
    config.initMode = "restart";
    config.seedConfig = "triangle";
    config.restartPath = "restart_lammps.dat";
    config.outputDir = ".";
}

bool is_valid_seed_config(const std::string &value)
{
    return value == "triangle" || value == "pentamer" || value == "hexamer";
}

bool apply_optional_argument(const std::string &option, const std::string &value, SimulationConfig &config, std::string &error_message)
{
    if (option == "--index-capacity")
    {
        if (!parse_unsigned_long_string(value, config.indexCapacity) || config.indexCapacity == 0 || config.indexCapacity > static_cast<unsigned long>(INT_MAX))
        {
            error_message = "Invalid value for --index-capacity.\n" + assemble_usage();
            return false;
        }
        return true;
    }

    if (option == "--max-sweeps")
    {
        if (!parse_unsigned_long_string(value, config.maxSweeps))
        {
            error_message = "Invalid value for --max-sweeps.\n" + assemble_usage();
            return false;
        }
        return true;
    }

    if (option == "--restart")
    {
        if (value.empty())
        {
            error_message = "Invalid value for --restart.\n" + assemble_usage();
            return false;
        }
        config.restartPath = value;
        return true;
    }

    if (option == "--init")
    {
        if (value == "restart")
        {
            config.initMode = value;
            return true;
        }
        if (value == "seed")
        {
            config.initMode = value;
            return true;
        }
        if (is_valid_seed_config(value))
        {
            config.initMode = "seed";
            config.seedConfig = value;
            return true;
        }
        error_message = "Invalid value for --init. Expected restart, seed, triangle, pentamer, or hexamer.\n" + assemble_usage();
        return false;
    }

    if (option == "--seed-config")
    {
        if (!is_valid_seed_config(value))
        {
            error_message = "Invalid value for --seed-config. Expected triangle, pentamer, or hexamer.\n" + assemble_usage();
            return false;
        }
        config.seedConfig = value;
        return true;
    }

    if (option == "--output-dir")
    {
        if (value.empty())
        {
            error_message = "Invalid value for --output-dir.\n" + assemble_usage();
            return false;
        }
        config.outputDir = value;
        return true;
    }

    error_message = std::string("Unknown option: ") + option + ".\n" + assemble_usage();
    return false;
}

bool parse_optional_arguments(int argc, char **argv, int start_index, SimulationConfig &config, std::string &error_message)
{
    for (int i = start_index; i < argc; )
    {
        const std::string option = argv[i];
        if (i + 1 >= argc || argv[i + 1] == nullptr)
        {
            error_message = std::string("Missing value for ") + option + ".\n" + assemble_usage();
            return false;
        }

        if (!apply_optional_argument(option, argv[i + 1], config, error_message))
        {
            return false;
        }
        i += 2;
    }
    return true;
}

bool load_simulation_config_file(const std::string &config_path, SimulationConfig &config, std::string &error_message)
{
    std::ifstream input(config_path.c_str());
    if (!input)
    {
        error_message = std::string("Failed to open config file: ") + config_path;
        return false;
    }

    bool seed_set = false;
    bool epsilon0_set = false;
    bool kappa0_set = false;
    bool kappaPhi0_set = false;
    bool theta0_set = false;
    bool theta1_set = false;
    bool gb0_set = false;
    bool muCd_set = false;
    bool ks0_set = false;
    bool dmu_set = false;
    bool dg_set = false;
    bool mudrug_set = false;
    bool gdrug0_set = false;
    bool kd0_set = false;
    bool dg12_set = false;
    bool dg01_set = false;
    bool dg20_set = false;
    bool dg33_set = false;
    bool dg00_set = false;
    bool dgother_set = false;

    ConfigFieldSpec numeric_fields[] = {
        {"epsilon0", &config.epsilon0, &epsilon0_set},
        {"kappa0", &config.kappa0, &kappa0_set},
        {"kappaPhi0", &config.kappaPhi0, &kappaPhi0_set},
        {"theta0", &config.theta0, &theta0_set},
        {"theta1", &config.theta1, &theta1_set},
        {"gb0", &config.gb0, &gb0_set},
        {"LnK", &config.gb0, &gb0_set},
        {"muCd", &config.muCd, &muCd_set},
        {"ks0", &config.ks0, &ks0_set},
        {"dmu", &config.dmu, &dmu_set},
        {"dg", &config.dg, &dg_set},
        {"mudrug", &config.mudrug, &mudrug_set},
        {"gdrug0", &config.gdrug0, &gdrug0_set},
        {"kd0", &config.kd0, &kd0_set},
        {"dg12", &config.dg12, &dg12_set},
        {"dg01", &config.dg01, &dg01_set},
        {"dg20", &config.dg20, &dg20_set},
        {"dg33", &config.dg33, &dg33_set},
        {"dg00", &config.dg00, &dg00_set},
        {"dgother", &config.dgother, &dgother_set},
    };

    const std::string config_dir = parent_directory(config_path);

    std::string raw_line;
    int line_number = 0;
    while (std::getline(input, raw_line))
    {
        ++line_number;
        std::size_t comment_pos = raw_line.find('#');
        std::string line = trim(raw_line.substr(0, comment_pos));
        if (line.empty())
        {
            continue;
        }

        std::size_t equals_pos = line.find('=');
        if (equals_pos == std::string::npos)
        {
            std::ostringstream message;
            message << "Invalid config line " << line_number << ": expected key=value.";
            error_message = message.str();
            return false;
        }

        const std::string key = trim(line.substr(0, equals_pos));
        const std::string value = trim(line.substr(equals_pos + 1));
        if (key.empty() || value.empty())
        {
            std::ostringstream message;
            message << "Invalid config line " << line_number << ": expected non-empty key and value.";
            error_message = message.str();
            return false;
        }

        if (key == "seed")
        {
            if (!parse_unsigned_long_string(value, config.seed))
            {
                error_message = "Invalid value for seed in config file.";
                return false;
            }
            seed_set = true;
            continue;
        }

        if (key == "indexCapacity" || key == "index_capacity")
        {
            if (!parse_unsigned_long_string(value, config.indexCapacity) || config.indexCapacity == 0 || config.indexCapacity > static_cast<unsigned long>(INT_MAX))
            {
                error_message = "Invalid value for indexCapacity in config file.";
                return false;
            }
            continue;
        }

        if (key == "maxSweeps" || key == "max_sweeps")
        {
            if (!parse_unsigned_long_string(value, config.maxSweeps))
            {
                error_message = "Invalid value for maxSweeps in config file.";
                return false;
            }
            continue;
        }

        if (key == "restart" || key == "restartPath" || key == "restart_path")
        {
            config.restartPath = resolve_relative_to(config_dir, value);
            continue;
        }

        if (key == "init" || key == "initMode" || key == "init_mode")
        {
            if (value == "restart" || value == "seed")
            {
                config.initMode = value;
                continue;
            }
            if (!is_valid_seed_config(value))
            {
                error_message = "Invalid value for init in config file.";
                return false;
            }
            config.initMode = "seed";
            config.seedConfig = value;
            continue;
        }

        if (key == "seedConfig" || key == "seed_config")
        {
            if (!is_valid_seed_config(value))
            {
                error_message = "Invalid value for seedConfig in config file.";
                return false;
            }
            config.seedConfig = value;
            continue;
        }

        if (key == "outputDir" || key == "output_dir")
        {
            config.outputDir = resolve_relative_to(config_dir, value);
            continue;
        }

        bool matched_numeric = false;
        const int field_count = sizeof(numeric_fields) / sizeof(numeric_fields[0]);
        for (int i = 0; i < field_count; ++i)
        {
            if (key == numeric_fields[i].name)
            {
                if (!parse_double_string(value, *numeric_fields[i].value))
                {
                    error_message = std::string("Invalid value for ") + key + " in config file.";
                    return false;
                }
                *numeric_fields[i].was_set = true;
                matched_numeric = true;
                break;
            }
        }

        if (!matched_numeric)
        {
            error_message = std::string("Unknown config key: ") + key;
            return false;
        }
    }

    if (!seed_set || !epsilon0_set || !kappa0_set || !kappaPhi0_set || !theta0_set || !theta1_set || !gb0_set ||
        !muCd_set || !ks0_set || !dmu_set || !dg_set || !mudrug_set || !gdrug0_set || !kd0_set || !dg12_set ||
        !dg01_set || !dg20_set || !dg33_set || !dg00_set || !dgother_set)
    {
        error_message = "Config file is missing one or more required fields.";
        return false;
    }

    return true;
}
}

std::string assemble_usage()
{
    return "usage: ./assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 dmu dummydg mudrug gdrug kd0 dg12 dg01 dg20 dg33 dg00 dgother [--index-capacity N] [--max-sweeps N] [--init MODE] [--seed-config NAME] [--restart PATH] [--output-dir PATH]\n"
           "   or: ./assemble --config PATH [--index-capacity N] [--max-sweeps N] [--init MODE] [--seed-config NAME] [--restart PATH] [--output-dir PATH]";
}

bool parse_simulation_config(int argc, char **argv, SimulationConfig &config, std::string &error_message)
{
    if (argc >= 2 && std::string(argv[1]) == "--config")
    {
        if (argc < 3 || argv[2] == nullptr || *argv[2] == '\0')
        {
            error_message = "Missing value for --config.\n" + assemble_usage();
            return false;
        }
        set_default_config(config);
        if (!load_simulation_config_file(argv[2], config, error_message))
        {
            error_message += "\n" + assemble_usage();
            return false;
        }
        if (!parse_optional_arguments(argc, argv, 3, config, error_message))
        {
            return false;
        }
        error_message.clear();
        return true;
    }

    if (argc < 21)
    {
        std::ostringstream message;
        message << "Expected at least 20 arguments after the executable name, got " << (argc - 1) << ".\n"
                << assemble_usage();
        error_message = message.str();
        return false;
    }

    set_default_config(config);

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

    if (!parse_optional_arguments(argc, argv, 21, config, error_message))
    {
        return false;
    }

    error_message.clear();
    return true;
}
