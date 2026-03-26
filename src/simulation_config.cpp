#include "simulation_config.hpp"

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

namespace
{
const unsigned long kMinimumIndexCapacity = 1000000UL;
const unsigned long kTestIndexCapacity = 1000000UL;
const unsigned long kRunIndexCapacity = 2000000UL;

struct NumericFieldSpec
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
    config.engine.indexCapacity = 0UL;
    config.engine.runMode.clear();
    config.runtime.maxSweeps = 0UL;
    config.runtime.outputDir = ".";
    config.runtime.workflow = "assembly";
    config.initialization.mode = "restart";
    config.initialization.seedConfig = "triangle";
    config.initialization.restartPath = "restart_lammps.dat";
}

bool is_valid_seed_config(const std::string &value)
{
    return value == "triangle" || value == "pentamer" || value == "hexamer";
}

bool is_valid_run_mode(const std::string &value)
{
    return value == "test" || value == "run" || value == "extended";
}

bool is_valid_runtime_workflow(const std::string &value)
{
    return value == "assembly" || value == "relaxation";
}

bool validate_index_capacity_value(unsigned long value)
{
    return value >= kMinimumIndexCapacity && value <= static_cast<unsigned long>(INT_MAX);
}

bool finalize_runtime_workflow(SimulationConfig &config, std::string &error_message)
{
    if (config.runtime.workflow.empty())
    {
        config.runtime.workflow = "assembly";
    }

    if (!is_valid_runtime_workflow(config.runtime.workflow))
    {
        error_message = "Invalid value for runtime.workflow. Expected assembly or relaxation.\n" + assemble_usage();
        return false;
    }

    if (config.runtime.workflow == "relaxation" && config.runtime.maxSweeps == 0UL)
    {
        error_message = "Relaxation workflow requires runtime.maxSweeps (or --max-sweeps) to be greater than zero.\n" + assemble_usage();
        return false;
    }

    return true;
}

bool finalize_run_mode(SimulationConfig &config, std::string &error_message)
{
    std::string mode = config.engine.runMode;
    if (mode.empty())
    {
        mode = config.engine.indexCapacity == 0UL ? "run" : "extended";
    }

    if (!is_valid_run_mode(mode))
    {
        error_message = "Invalid value for --run-mode. Expected test, run, or extended.\n" + assemble_usage();
        return false;
    }

    if (mode == "test")
    {
        if (config.engine.indexCapacity != 0UL)
        {
            error_message = "The test run mode uses a fixed preset and cannot be combined with --index-capacity or engine.indexCapacity.\n" + assemble_usage();
            return false;
        }
        config.engine.indexCapacity = kTestIndexCapacity;
    }
    else if (mode == "run")
    {
        if (config.engine.indexCapacity != 0UL)
        {
            error_message = "The run mode uses a fixed preset and cannot be combined with --index-capacity or engine.indexCapacity.\n" + assemble_usage();
            return false;
        }
        config.engine.indexCapacity = kRunIndexCapacity;
    }
    else
    {
        if (config.engine.indexCapacity == 0UL)
        {
            error_message = "Extended run mode requires --index-capacity (or engine.indexCapacity in the config file).\n" + assemble_usage();
            return false;
        }
        if (!validate_index_capacity_value(config.engine.indexCapacity))
        {
            std::ostringstream message;
            message << "Invalid value for index capacity. Extended mode requires a value between "
                    << kMinimumIndexCapacity << " and " << INT_MAX << ".\n"
                    << assemble_usage();
            error_message = message.str();
            return false;
        }
    }

    config.engine.runMode = mode;
    return true;
}

bool apply_optional_argument(const std::string &option, const std::string &value, SimulationConfig &config, std::string &error_message)
{
    if (option == "--run-mode")
    {
        if (!is_valid_run_mode(value))
        {
            error_message = "Invalid value for --run-mode. Expected test, run, or extended.\n" + assemble_usage();
            return false;
        }
        config.engine.runMode = value;
        return true;
    }

    if (option == "--index-capacity")
    {
        if (!parse_unsigned_long_string(value, config.engine.indexCapacity) || config.engine.indexCapacity == 0UL || config.engine.indexCapacity > static_cast<unsigned long>(INT_MAX))
        {
            error_message = "Invalid value for --index-capacity.\n" + assemble_usage();
            return false;
        }
        return true;
    }

    if (option == "--max-sweeps")
    {
        if (!parse_unsigned_long_string(value, config.runtime.maxSweeps))
        {
            error_message = "Invalid value for --max-sweeps.\n" + assemble_usage();
            return false;
        }
        return true;
    }

    if (option == "--workflow")
    {
        if (!is_valid_runtime_workflow(value))
        {
            error_message = "Invalid value for --workflow. Expected assembly or relaxation.\n" + assemble_usage();
            return false;
        }
        config.runtime.workflow = value;
        return true;
    }

    if (option == "--restart")
    {
        if (value.empty())
        {
            error_message = "Invalid value for --restart.\n" + assemble_usage();
            return false;
        }
        config.initialization.restartPath = value;
        return true;
    }

    if (option == "--init")
    {
        if (value == "restart" || value == "seed")
        {
            config.initialization.mode = value;
            return true;
        }
        if (is_valid_seed_config(value))
        {
            config.initialization.mode = "seed";
            config.initialization.seedConfig = value;
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
        config.initialization.seedConfig = value;
        return true;
    }

    if (option == "--output-dir")
    {
        if (value.empty())
        {
            error_message = "Invalid value for --output-dir.\n" + assemble_usage();
            return false;
        }
        config.runtime.outputDir = value;
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

bool handle_section_numeric_field(const std::string &section, const std::string &key, const std::string &value, NumericFieldSpec *fields, int field_count, std::string &error_message)
{
    if (section.empty())
    {
        error_message = "Config keys must appear inside a named section.";
        return false;
    }

    for (int i = 0; i < field_count; ++i)
    {
        if (key == fields[i].name)
        {
            if (!parse_double_string(value, *fields[i].value))
            {
                error_message = std::string("Invalid value for ") + section + "." + key + ".";
                return false;
            }
            *fields[i].was_set = true;
            return true;
        }
    }

    return false;
}

bool load_simulation_config_file(const std::string &config_path, SimulationConfig &config, std::string &error_message)
{
    std::ifstream input(config_path.c_str());
    if (!input)
    {
        error_message = std::string("Failed to open config file: ") + config_path;
        return false;
    }

    bool runtime_seed_set = false;

    bool epsilon0_set = false;
    bool kappa0_set = false;
    bool kappaPhi0_set = false;
    bool theta0_set = false;
    bool theta1_set = false;
    bool gb0_set = false;
    bool dg12_set = false;
    bool dg01_set = false;
    bool dg20_set = false;
    bool dg33_set = false;
    bool dg00_set = false;
    bool dgother_set = false;

    bool muCd_set = false;
    bool ks0_set = false;
    bool dmu_set = false;
    bool dg_set = false;

    bool mudrug_set = false;
    bool gdrug0_set = false;
    bool kd0_set = false;

    NumericFieldSpec capsid_fields[] = {
        {"epsilon0", &config.capsidGeometry.epsilon0, &epsilon0_set},
        {"kappa0", &config.capsidGeometry.kappa0, &kappa0_set},
        {"kappaPhi0", &config.capsidGeometry.kappaPhi0, &kappaPhi0_set},
        {"theta0", &config.capsidGeometry.theta0, &theta0_set},
        {"theta1", &config.capsidGeometry.theta1, &theta1_set},
        {"gb0", &config.capsidGeometry.gb0, &gb0_set},
        {"dg12", &config.capsidGeometry.dg12, &dg12_set},
        {"dg01", &config.capsidGeometry.dg01, &dg01_set},
        {"dg20", &config.capsidGeometry.dg20, &dg20_set},
        {"dg33", &config.capsidGeometry.dg33, &dg33_set},
        {"dg00", &config.capsidGeometry.dg00, &dg00_set},
        {"dgother", &config.capsidGeometry.dgother, &dgother_set},
    };

    NumericFieldSpec simulation_fields[] = {
        {"muCd", &config.simulation.muCd, &muCd_set},
        {"ks0", &config.simulation.ks0, &ks0_set},
        {"dmu", &config.simulation.dmu, &dmu_set},
        {"dg", &config.simulation.dg, &dg_set},
    };

    NumericFieldSpec drug_fields[] = {
        {"mudrug", &config.drug.mudrug, &mudrug_set},
        {"gdrug0", &config.drug.gdrug0, &gdrug0_set},
        {"kd0", &config.drug.kd0, &kd0_set},
    };

    const std::string config_dir = parent_directory(config_path);
    std::string current_section;

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

        if (line[0] == '[')
        {
            if (line[line.size() - 1] != ']')
            {
                std::ostringstream message;
                message << "Invalid config line " << line_number << ": malformed section header.";
                error_message = message.str();
                return false;
            }

            current_section = trim(line.substr(1, line.size() - 2));
            if (current_section != "capsid_geometry" &&
                current_section != "simulation" &&
                current_section != "drug" &&
                current_section != "init" &&
                current_section != "runtime" &&
                current_section != "engine")
            {
                error_message = std::string("Unknown config section: ") + current_section;
                return false;
            }
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

        if (current_section.empty())
        {
            error_message = "Config keys must appear inside a named section.";
            return false;
        }

        if (current_section == "capsid_geometry")
        {
            if (!handle_section_numeric_field(current_section, key, value, capsid_fields, sizeof(capsid_fields) / sizeof(capsid_fields[0]), error_message))
            {
                error_message = std::string("Unknown config key in [capsid_geometry]: ") + key;
                return false;
            }
            continue;
        }

        if (current_section == "simulation")
        {
            if (!handle_section_numeric_field(current_section, key, value, simulation_fields, sizeof(simulation_fields) / sizeof(simulation_fields[0]), error_message))
            {
                error_message = std::string("Unknown config key in [simulation]: ") + key;
                return false;
            }
            continue;
        }

        if (current_section == "drug")
        {
            if (!handle_section_numeric_field(current_section, key, value, drug_fields, sizeof(drug_fields) / sizeof(drug_fields[0]), error_message))
            {
                error_message = std::string("Unknown config key in [drug]: ") + key;
                return false;
            }
            continue;
        }

        if (current_section == "init")
        {
            if (key == "mode")
            {
                if (value == "restart" || value == "seed")
                {
                    config.initialization.mode = value;
                    continue;
                }
                if (is_valid_seed_config(value))
                {
                    config.initialization.mode = "seed";
                    config.initialization.seedConfig = value;
                    continue;
                }
                error_message = "Invalid value for init.mode.";
                return false;
            }

            if (key == "seedShape")
            {
                if (!is_valid_seed_config(value))
                {
                    error_message = "Invalid value for init.seedShape.";
                    return false;
                }
                config.initialization.seedConfig = value;
                continue;
            }

            if (key == "restartPath")
            {
                config.initialization.restartPath = resolve_relative_to(config_dir, value);
                continue;
            }

            error_message = std::string("Unknown config key in [init]: ") + key;
            return false;
        }

        if (current_section == "runtime")
        {
            if (key == "seed")
            {
                if (!parse_unsigned_long_string(value, config.runtime.seed))
                {
                    error_message = "Invalid value for runtime.seed.";
                    return false;
                }
                runtime_seed_set = true;
                continue;
            }

            if (key == "maxSweeps")
            {
                if (!parse_unsigned_long_string(value, config.runtime.maxSweeps))
                {
                    error_message = "Invalid value for runtime.maxSweeps.";
                    return false;
                }
                continue;
            }

            if (key == "outputDir")
            {
                config.runtime.outputDir = resolve_relative_to(config_dir, value);
                continue;
            }

            if (key == "workflow")
            {
                if (!is_valid_runtime_workflow(value))
                {
                    error_message = "Invalid value for runtime.workflow.";
                    return false;
                }
                config.runtime.workflow = value;
                continue;
            }

            error_message = std::string("Unknown config key in [runtime]: ") + key;
            return false;
        }

        if (current_section == "engine")
        {
            if (key == "profile")
            {
                if (!is_valid_run_mode(value))
                {
                    error_message = "Invalid value for engine.profile.";
                    return false;
                }
                config.engine.runMode = value;
                continue;
            }

            if (key == "indexCapacity")
            {
                if (!parse_unsigned_long_string(value, config.engine.indexCapacity) || config.engine.indexCapacity == 0UL || config.engine.indexCapacity > static_cast<unsigned long>(INT_MAX))
                {
                    error_message = "Invalid value for engine.indexCapacity.";
                    return false;
                }
                continue;
            }

            error_message = std::string("Unknown config key in [engine]: ") + key;
            return false;
        }
    }

    if (!runtime_seed_set ||
        !epsilon0_set || !kappa0_set || !kappaPhi0_set || !theta0_set || !theta1_set || !gb0_set ||
        !dg12_set || !dg01_set || !dg20_set || !dg33_set || !dg00_set || !dgother_set ||
        !muCd_set || !ks0_set || !dmu_set || !dg_set ||
        !mudrug_set || !gdrug0_set || !kd0_set)
    {
        error_message = "Config file is missing one or more required fields.";
        return false;
    }

    return true;
}
}

std::string assemble_usage()
{
    std::ostringstream usage;
    usage << "usage: ./assemble --config PATH "
          << "[--run-mode MODE] [--index-capacity N] [--max-sweeps N] [--workflow MODE] [--init MODE] [--seed-config NAME] [--restart PATH] [--output-dir PATH]\n"
          << "config sections: [capsid_geometry], [simulation], [drug], [init], [runtime], [engine]\n"
          << "runtime workflows: assembly (default), relaxation (move_vertex only)\n"
          << "engine profiles: test=" << kTestIndexCapacity
          << ", run=" << kRunIndexCapacity
          << ", extended unlocks --index-capacity >= " << kMinimumIndexCapacity;
    return usage.str();
}

std::string render_simulation_config(const SimulationConfig &config)
{
    std::ostringstream output;
    output.setf(std::ios::fixed);
    output.precision(6);

    output << "[capsid_geometry]\n";
    output << "epsilon0 = " << config.capsidGeometry.epsilon0 << "\n";
    output << "kappa0 = " << config.capsidGeometry.kappa0 << "\n";
    output << "kappaPhi0 = " << config.capsidGeometry.kappaPhi0 << "\n";
    output << "theta0 = " << config.capsidGeometry.theta0 << "\n";
    output << "theta1 = " << config.capsidGeometry.theta1 << "\n";
    output << "gb0 = " << config.capsidGeometry.gb0 << "\n";
    output << "dg12 = " << config.capsidGeometry.dg12 << "\n";
    output << "dg01 = " << config.capsidGeometry.dg01 << "\n";
    output << "dg20 = " << config.capsidGeometry.dg20 << "\n";
    output << "dg33 = " << config.capsidGeometry.dg33 << "\n";
    output << "dg00 = " << config.capsidGeometry.dg00 << "\n";
    output << "dgother = " << config.capsidGeometry.dgother << "\n\n";

    output << "[simulation]\n";
    output << "muCd = " << config.simulation.muCd << "\n";
    output << "ks0 = " << config.simulation.ks0 << "\n";
    output << "dmu = " << config.simulation.dmu << "\n";
    output << "dg = " << config.simulation.dg << "\n\n";

    output << "[drug]\n";
    output << "mudrug = " << config.drug.mudrug << "\n";
    output << "gdrug0 = " << config.drug.gdrug0 << "\n";
    output << "kd0 = " << config.drug.kd0 << "\n\n";

    output << "[init]\n";
    output << "mode = " << config.initialization.mode << "\n";
    output << "seedShape = " << config.initialization.seedConfig << "\n";
    output << "restartPath = " << config.initialization.restartPath << "\n\n";

    output << "[runtime]\n";
    output << "seed = " << config.runtime.seed << "\n";
    output << "maxSweeps = " << config.runtime.maxSweeps << "\n";
    output << "outputDir = " << config.runtime.outputDir << "\n";
    output << "workflow = " << config.runtime.workflow << "\n\n";

    output << "[engine]\n";
    output << "profile = " << config.engine.runMode << "\n";
    if (config.engine.runMode == "extended")
    {
        output << "indexCapacity = " << config.engine.indexCapacity << "\n";
    }

    return output.str();
}

bool parse_simulation_config(int argc, char **argv, SimulationConfig &config, std::string &error_message)
{
    if (argc < 2 || std::string(argv[1]) != "--config")
    {
        error_message = "HEVA now requires --config PATH using the sectioned config format.\n" + assemble_usage();
        return false;
    }

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

    if (!finalize_run_mode(config, error_message))
    {
        return false;
    }

    if (!finalize_runtime_workflow(config, error_message))
    {
        return false;
    }

    error_message.clear();
    return true;
}
