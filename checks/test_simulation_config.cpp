#include "simulation_config.hpp"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

namespace
{
void write_file(const std::string &path, const std::string &contents)
{
    std::ofstream output(path.c_str());
    assert(output.good());
    output << contents;
}
}

int main()
{
    char arg0[] = "assemble";
    SimulationConfig config;
    std::string error_message;

    char *missing_config_argv[] = {arg0};
    bool parsed = parse_simulation_config(1, missing_config_argv, config, error_message);
    assert(!parsed);
    assert(error_message.find("requires --config PATH") != std::string::npos);

    const std::string good_config_path = "simulation_config_test.in";
    write_file(
        good_config_path,
        "[capsid_geometry]\n"
        "epsilon0 = 4200\n"
        "kappa0 = 40\n"
        "kappaPhi0 = 800\n"
        "theta0 = 0.240\n"
        "theta1 = 0.480\n"
        "gb0 = -9.8\n"
        "dg12 = 0.3\n"
        "dg01 = 0.1\n"
        "dg20 = -0.1\n"
        "dg33 = 0.0\n"
        "dg00 = -0.6\n"
        "dgother = -0.95\n"
        "\n"
        "[simulation]\n"
        "muCd = -11.5\n"
        "ks0 = 0.02\n"
        "dmu = -4.5\n"
        "dg = 0.0\n"
        "\n"
        "[drug]\n"
        "mudrug = -100.0\n"
        "gdrug0 = 10.0\n"
        "kd0 = 0.0\n"
        "\n"
        "[init]\n"
        "mode = restart\n"
        "restartPath = restart_lammps.dat\n"
        "\n"
        "[runtime]\n"
        "seed = 7\n"
        "outputDir = out\n"
        "\n"
        "[engine]\n"
        "profile = run\n");

    char cfg_opt[] = "--config";
    char cfg_path[] = "simulation_config_test.in";
    char *config_argv[] = {arg0, cfg_opt, cfg_path};
    parsed = parse_simulation_config(3, config_argv, config, error_message);
    assert(parsed);
    assert(error_message.empty());
    assert(config.runtime.seed == 7UL);
    assert(config.engine.indexCapacity == 2000000UL);
    assert(config.engine.runMode == "run");
    assert(config.runtime.maxSweeps == 0UL);
    assert(config.initialization.mode == "restart");
    assert(config.initialization.seedConfig == "triangle");
    assert(config.initialization.restartPath == "restart_lammps.dat");
    assert(config.runtime.outputDir == "out");
    assert(config.capsidGeometry.epsilon0 == 4200.0);
    assert(config.capsidGeometry.kappaPhi0 == 800.0);
    assert(config.capsidGeometry.dg00 == -0.6);
    assert(config.simulation.muCd == -11.5);
    assert(config.simulation.ks0 == 0.02);
    assert(config.drug.mudrug == -100.0);
    assert(config.drug.gdrug0 == 10.0);

    char override_opt1[] = "--max-sweeps";
    char override_val1[] = "25";
    char override_opt2[] = "--restart";
    char override_val2[] = "fixtures/restart_lammps.dat";
    char override_opt3[] = "--init";
    char override_val3[] = "seed";
    char override_opt4[] = "--seed-config";
    char override_val4[] = "hexamer";
    char override_opt5[] = "--output-dir";
    char override_val5[] = "runs/smoke";
    char override_opt6[] = "--run-mode";
    char override_val6[] = "extended";
    char override_opt7[] = "--index-capacity";
    char override_val7[] = "2000000";
    char *with_overrides_argv[] = {
        arg0, cfg_opt, cfg_path,
        override_opt1, override_val1,
        override_opt2, override_val2,
        override_opt3, override_val3,
        override_opt4, override_val4,
        override_opt5, override_val5,
        override_opt6, override_val6,
        override_opt7, override_val7,
    };
    parsed = parse_simulation_config(17, with_overrides_argv, config, error_message);
    assert(parsed);
    assert(config.runtime.maxSweeps == 25UL);
    assert(config.initialization.mode == "seed");
    assert(config.initialization.seedConfig == "hexamer");
    assert(config.initialization.restartPath == "fixtures/restart_lammps.dat");
    assert(config.runtime.outputDir == "runs/smoke");
    assert(config.engine.runMode == "extended");
    assert(config.engine.indexCapacity == 2000000UL);

    char fixture_path[] = "fixtures/smoke_config.in";
    char cfg_override_opt[] = "--max-sweeps";
    char cfg_override_val[] = "5";
    char *fixture_argv[] = {arg0, cfg_opt, fixture_path, cfg_override_opt, cfg_override_val};
    parsed = parse_simulation_config(5, fixture_argv, config, error_message);
    assert(parsed);
    assert(config.runtime.seed == 521759UL);
    assert(config.runtime.maxSweeps == 5UL);
    assert(config.engine.indexCapacity == 1000000UL);
    assert(config.engine.runMode == "extended");
    assert(config.initialization.mode == "restart");
    assert(config.initialization.seedConfig == "triangle");
    assert(config.initialization.restartPath == "fixtures/restart_lammps.dat");
    assert(config.runtime.outputDir == "fixtures/generated-output");

    const std::string bad_config_path = "simulation_config_bad.in";
    write_file(
        bad_config_path,
        "[capsid_geometry]\n"
        "epsilon0 = 4200\n");

    char bad_cfg_path[] = "simulation_config_bad.in";
    char *bad_cfg_argv[] = {arg0, cfg_opt, bad_cfg_path};
    parsed = parse_simulation_config(3, bad_cfg_argv, config, error_message);
    assert(!parsed);
    assert(error_message.find("missing one or more required fields") != std::string::npos);

    const std::string bad_extended_path = "simulation_config_extended_bad.in";
    write_file(
        bad_extended_path,
        "[capsid_geometry]\n"
        "epsilon0 = 4200\n"
        "kappa0 = 40\n"
        "kappaPhi0 = 800\n"
        "theta0 = 0.240\n"
        "theta1 = 0.480\n"
        "gb0 = -9.8\n"
        "dg12 = 0.3\n"
        "dg01 = 0.1\n"
        "dg20 = -0.1\n"
        "dg33 = 0.0\n"
        "dg00 = -0.6\n"
        "dgother = -0.95\n"
        "\n"
        "[simulation]\n"
        "muCd = -11.5\n"
        "ks0 = 0.02\n"
        "dmu = -4.5\n"
        "dg = 0.0\n"
        "\n"
        "[drug]\n"
        "mudrug = 0.0\n"
        "gdrug0 = 0.0\n"
        "kd0 = 0.0\n"
        "\n"
        "[runtime]\n"
        "seed = 7\n"
        "\n"
        "[engine]\n"
        "profile = extended\n");

    char bad_extended_cfg[] = "simulation_config_extended_bad.in";
    char *bad_extended_argv[] = {arg0, cfg_opt, bad_extended_cfg};
    parsed = parse_simulation_config(3, bad_extended_argv, config, error_message);
    assert(!parsed);
    assert(error_message.find("Extended run mode requires") != std::string::npos);

    assert(assemble_usage().find("--config PATH") != std::string::npos);
    assert(assemble_usage().find("[capsid_geometry]") != std::string::npos);
    assert(assemble_usage().find("[simulation]") != std::string::npos);
    assert(assemble_usage().find("[drug]") != std::string::npos);
    assert(assemble_usage().find("[runtime]") != std::string::npos);
    assert(assemble_usage().find("[engine]") != std::string::npos);

    assert(render_simulation_config(config).find("[capsid_geometry]") != std::string::npos);
    assert(render_simulation_config(config).find("profile = extended") != std::string::npos);

    std::remove(good_config_path.c_str());
    std::remove(bad_config_path.c_str());
    std::remove(bad_extended_path.c_str());

    std::cout << "simulation_config checks passed" << std::endl;
    return 0;
}
