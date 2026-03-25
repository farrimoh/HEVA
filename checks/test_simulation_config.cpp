#include "simulation_config.hpp"

#include <cassert>
#include <iostream>
#include <string>

int main()
{
    char arg0[] = "assemble";
    char arg1[] = "7";
    char arg2[] = "4200";
    char arg3[] = "40";
    char arg4[] = "800";
    char arg5[] = "0.240";
    char arg6[] = "0.480";
    char arg7[] = "-9.8";
    char arg8[] = "-11.5";
    char arg9[] = "0.02";
    char arg10[] = "-4.5";
    char arg11[] = "0.0";
    char arg12[] = "-100.0";
    char arg13[] = "10.0";
    char arg14[] = "0.0";
    char arg15[] = "0.3";
    char arg16[] = "0.1";
    char arg17[] = "-0.1";
    char arg18[] = "0.0";
    char arg19[] = "-0.6";
    char arg20[] = "-0.95";

    char *argv[] = {
        arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
    };

    SimulationConfig config;
    std::string error_message;
    bool parsed = parse_simulation_config(21, argv, config, error_message);
    assert(parsed);
    assert(error_message.empty());
    assert(config.seed == 7UL);
    assert(config.indexCapacity == 2000000UL);
    assert(config.runMode == "run");
    assert(config.maxSweeps == 0UL);
    assert(config.initMode == "restart");
    assert(config.seedConfig == "triangle");
    assert(config.restartPath == "restart_lammps.dat");
    assert(config.outputDir == ".");
    assert(config.epsilon0 == 4200.0);
    assert(config.kappaPhi0 == 800.0);
    assert(config.dg00 == -0.6);
    assert(config.muCd == -11.5);
    assert(config.ks0 == 0.02);

    char bad_arg3[] = "not-a-number";
    char *bad_argv[] = {
        arg0, arg1, arg2, bad_arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
    };

    SimulationConfig bad_config;
    parsed = parse_simulation_config(21, bad_argv, bad_config, error_message);
    assert(!parsed);
    assert(error_message.find("kappa0") != std::string::npos);

    parsed = parse_simulation_config(20, argv, bad_config, error_message);
    assert(!parsed);
    assert(error_message.find("Expected at least 20 arguments") != std::string::npos);

    char opt1[] = "--index-capacity";
    char val1[] = "2000000";
    char opt2[] = "--max-sweeps";
    char val2[] = "25";
    char opt3[] = "--restart";
    char val3[] = "fixtures/restart_lammps.dat";
    char opt4[] = "--init";
    char val4[] = "seed";
    char opt5[] = "--seed-config";
    char val5[] = "hexamer";
    char opt6[] = "--output-dir";
    char val6[] = "runs/smoke";
    char *with_options_argv[] = {
        arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
        opt1, val1, opt2, val2, opt3, val3, opt4, val4, opt5, val5, opt6, val6,
    };
    parsed = parse_simulation_config(33, with_options_argv, config, error_message);
    assert(parsed);
    assert(config.indexCapacity == 2000000UL);
    assert(config.runMode == "extended");
    assert(config.maxSweeps == 25UL);
    assert(config.initMode == "seed");
    assert(config.seedConfig == "hexamer");
    assert(config.restartPath == "fixtures/restart_lammps.dat");
    assert(config.outputDir == "runs/smoke");

    char cfg_opt[] = "--config";
    char cfg_path[] = "fixtures/smoke_config.in";
    char cfg_override_opt[] = "--max-sweeps";
    char cfg_override_val[] = "5";
    char *config_argv[] = {arg0, cfg_opt, cfg_path, cfg_override_opt, cfg_override_val};
    parsed = parse_simulation_config(5, config_argv, config, error_message);
    assert(parsed);
    assert(config.seed == 521759UL);
    assert(config.maxSweeps == 5UL);
    assert(config.indexCapacity == 1000000UL);
    assert(config.runMode == "extended");
    assert(config.initMode == "restart");
    assert(config.seedConfig == "triangle");

    const std::string expected_restart = "fixtures/restart_lammps.dat";
    const std::string expected_output_dir = "fixtures/generated-output";
    assert(config.restartPath == expected_restart);
    assert(config.outputDir == expected_output_dir);

    assert(assemble_usage().find("epsilon0") != std::string::npos);
    assert(assemble_usage().find("--config") != std::string::npos);
    assert(assemble_usage().find("--run-mode") != std::string::npos);
    assert(assemble_usage().find("--index-capacity") != std::string::npos);
    assert(assemble_usage().find("--max-sweeps") != std::string::npos);
    assert(assemble_usage().find("--init") != std::string::npos);
    assert(assemble_usage().find("--seed-config") != std::string::npos);
    assert(assemble_usage().find("--restart") != std::string::npos);
    assert(assemble_usage().find("--output-dir") != std::string::npos);

    char mode_opt[] = "--run-mode";
    char mode_val[] = "test";
    char *test_mode_argv[] = {
        arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
        mode_opt, mode_val,
    };
    parsed = parse_simulation_config(23, test_mode_argv, config, error_message);
    assert(parsed);
    assert(config.runMode == "test");
    assert(config.indexCapacity == 1000000UL);

    char extended_mode_val[] = "extended";
    char *bad_extended_mode_argv[] = {
        arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
        mode_opt, extended_mode_val,
    };
    parsed = parse_simulation_config(23, bad_extended_mode_argv, config, error_message);
    assert(!parsed);
    assert(error_message.find("Extended run mode requires") != std::string::npos);

    char run_mode_val[] = "run";
    char tiny_capacity_val[] = "10";
    char *bad_run_capacity_argv[] = {
        arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
        mode_opt, run_mode_val, opt1, tiny_capacity_val,
    };
    parsed = parse_simulation_config(25, bad_run_capacity_argv, config, error_message);
    assert(!parsed);
    assert(error_message.find("run mode uses a fixed preset") != std::string::npos);

    char *bad_small_capacity_argv[] = {
        arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10,
        arg11, arg12, arg13, arg14, arg15, arg16, arg17, arg18, arg19, arg20,
        opt1, tiny_capacity_val,
    };
    parsed = parse_simulation_config(23, bad_small_capacity_argv, config, error_message);
    assert(!parsed);
    assert(error_message.find("Extended mode requires a value between 1000000") != std::string::npos);

    std::cout << "simulation_config checks passed" << std::endl;
    return 0;
}
