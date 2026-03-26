/*
 * run_assembly.cpp
 *
 * Created on: Apr 25, 2019
 * Author: Farri Mohajerani
 *
 * MIT License
 *
 * Copyright (c) 2022 Farri Mohajerani
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "geometry.hpp"
#include "simulation_config.hpp"
#include "simulation_runner.hpp"
#include "simulation_setup.hpp"

#include <gsl/gsl_rng.h>

#include <cerrno>
#include <climits>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <time.h>
#include <unistd.h>

using namespace std;

namespace
{
bool is_absolute_path(const std::string &path)
{
    return !path.empty() && path[0] == '/';
}

std::string join_paths(const std::string &left, const std::string &right)
{
    if (left.empty())
    {
        return right;
    }
    if (left == "/")
    {
        return "/" + right;
    }
    if (!left.empty() && left[left.size() - 1] == '/')
    {
        return left + right;
    }
    return left + "/" + right;
}

std::string resolve_path_from(const std::string &base_dir, const std::string &path)
{
    if (path.empty() || path == ".")
    {
        return base_dir;
    }
    if (is_absolute_path(path))
    {
        return path;
    }
    return join_paths(base_dir, path);
}

bool ensure_directory_exists(const std::string &path, std::string &error_message)
{
    if (path.empty())
    {
        error_message = "Output directory path cannot be empty.";
        return false;
    }

    std::string current = is_absolute_path(path) ? "/" : "";
    std::size_t start = is_absolute_path(path) ? 1 : 0;

    while (start <= path.size())
    {
        std::size_t end = path.find('/', start);
        std::string part = path.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (!part.empty())
        {
            current = current.empty() ? part : join_paths(current, part);
            if (mkdir(current.c_str(), 0777) != 0 && errno != EEXIST)
            {
                error_message = std::string("Failed to create output directory '") + path + "': " + std::strerror(errno);
                return false;
            }
        }

        if (end == std::string::npos)
        {
            break;
        }
        start = end + 1;
    }

    return true;
}

std::string relative_to_output_dir_or_absolute(const std::string &output_dir, const std::string &path)
{
    const std::string prefix = output_dir == "/" ? "/" : output_dir + "/";
    if (path == output_dir)
    {
        return ".";
    }
    if (path.find(prefix) == 0)
    {
        return path.substr(prefix.size());
    }
    return path;
}

SimulationConfig make_persisted_run_config(const SimulationConfig &config, const std::string &restart_path, const std::string &output_dir)
{
    SimulationConfig persisted = config;
    persisted.runtime.outputDir = ".";
    if (persisted.initialization.mode == "restart" || persisted.initialization.mode == "legacy_lammps" || persisted.initialization.mode == "initial_frame")
    {
        persisted.initialization.restartPath = relative_to_output_dir_or_absolute(output_dir, restart_path);
    }
    else
    {
        persisted.initialization.restartPath = "restart_lammps.dat";
    }
    return persisted;
}

void write_invocation(FILE *stream, const SimulationConfig &config, unsigned long sweep)
{
    const std::string rendered = render_simulation_config(config);
    fprintf(stream, "# normalized run configuration\n%s", rendered.c_str());
    fprintf(stream, "\n# effective engine.profile=%s engine.indexCapacity=%lu sweep=%lu\n",
            config.engine.runMode.c_str(), config.engine.indexCapacity, sweep);
}
}

int main(int argc, char **argv)
{
    time_t timer1;
    time(&timer1);

    SimulationConfig config;
    std::string config_error;
    if (!parse_simulation_config(argc, argv, config, config_error))
    {
        std::cerr << config_error << std::endl;
        return -1;
    }

    char cwd_buffer[PATH_MAX];
    if (getcwd(cwd_buffer, sizeof(cwd_buffer)) == nullptr)
    {
        std::cerr << "Failed to resolve the current working directory: " << std::strerror(errno) << std::endl;
        return -1;
    }

    const std::string launch_dir = cwd_buffer;
    const std::string restart_path = resolve_path_from(launch_dir, config.initialization.restartPath);
    const std::string output_dir = resolve_path_from(launch_dir, config.runtime.outputDir);

    if ((config.initialization.mode == "restart" || config.initialization.mode == "legacy_lammps" || config.initialization.mode == "initial_frame") && access(restart_path.c_str(), R_OK) != 0)
    {
        std::cerr << "Initialization file is not readable: " << restart_path << std::endl;
        return -1;
    }

    std::string path_error;
    if (!ensure_directory_exists(output_dir, path_error))
    {
        std::cerr << path_error << std::endl;
        return -1;
    }

    const gsl_rng_type *t = gsl_rng_taus2;
    gsl_rng *rng = gsl_rng_alloc(t);
    srand((unsigned)config.runtime.seed);
    gsl_rng_set(rng, config.runtime.seed);
    cout << "HERE " << endl;

    geometry g;
    set_output_directory(output_dir);
    apply_simulation_config(config, g);

    cout << "l_thermal_sigma is " << g.l_thermal_sigma << endl;
    cout << "l_thermal_kappa is " << g.l_thermal_kappa << endl;
    cout << "theta_thermal_kappa is " << g.theta_thermal_kappa << endl;
    cout << "gaussian sigma " << g.gaussian_sigma << endl;

    FILE *ofile;
    FILE *fi;
    const std::string energy_path = join_paths(output_dir, "energy.dat");
    const std::string run_config_path = join_paths(output_dir, "run_config.out");
    const std::string parameters_path = join_paths(output_dir, "parameters_run.out");
    const SimulationConfig persisted_config = make_persisted_run_config(config, restart_path, output_dir);
    g.dump_parameters();
    ofile = fopen(energy_path.c_str(), "a");
    fprintf(stderr, " log files exist\n");
    cout << " file openned" << endl;

    fi = fopen(run_config_path.c_str(), "w");
    const std::string rendered_config = render_simulation_config(persisted_config);
    fprintf(fi, "%s", rendered_config.c_str());
    fclose(fi);

    fi = fopen(parameters_path.c_str(), "a");
    fprintf(fi, "%s\n", assemble_usage().c_str());
    write_invocation(fi, persisted_config, 0UL);

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {
            if (g.gb[i][j] != 0)
            {
                fprintf(fi, "half_edge %d -> %d : %.3f kT \n", i, j, g.gb[i][j]);
                fprintf(stderr, "half_edge %d -> %d : %.3f kT \n", i, j, g.gb[i][j]);
            }
        }
    }

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {
            if (g.gdrug[i][j] != 0)
            {
                fprintf(fi, "half_edge %d -> %d : %.3f kT \n", i, j, g.gdrug[i][j]);
                fprintf(stderr, "half_edge %d -> %d : %.3f kT \n", i, j, g.gdrug[i][j]);
            }
        }
    }

    fprintf(fi, "theta_thermal_kappa %.5f\n", g.theta_thermal_kappa);
    fprintf(fi, "l_thermal_kappa %.5f\n", g.l_thermal_kappa);
    fprintf(fi, "l_thermal_epsilon %.5f\n", g.l_thermal_sigma);
    fprintf(fi, "gaussian sigma %.5f\n", g.gaussian_sigma);
    fflush(fi);
    fclose(fi);

    unsigned long sweep = 0;
    fprintf(stderr, "%s\n", assemble_usage().c_str());
    write_invocation(stderr, persisted_config, sweep);
    cout << "# WORKFLOW " << config.runtime.workflow << endl;

    for (int j = 0; j < 4; j++)
    {
        fprintf(stderr, "%d  g.epsilon %f g.kappa %f g.kappaPhi %f g.l0 %f  g.theta0 %f g.phi0 %f\n", j, g.epsilon[j], g.kappa[j], g.kappaPhi[j], g.l0[j], g.theta0[j], g.phi0[j]);
    }

    SimulationRunStats stats;
    SimulationLoopSettings settings = make_simulation_loop_settings(config);
    if (config.initialization.mode == "restart")
    {
        initialize_from_restart(g, rng, restart_path.c_str(), sweep, stats, settings);
    }
    else if (config.initialization.mode == "legacy_lammps" || config.initialization.mode == "initial_frame")
    {
        initialize_from_initial_frame_compat(g, restart_path.c_str(), sweep, stats, settings);
    }
    else
    {
        initialize_from_seed(g, rng, config.initialization.seedConfig.c_str(), sweep, stats, settings);
    }
    SimulationStopReason stop_reason = SIMULATION_STOP_MAX_SWEEPS;
    if (config.runtime.workflow == "relaxation")
    {
        stop_reason = run_relaxation_loop(g, rng, ofile, config.runtime.seed, timer1, sweep, stats, settings);
    }
    else
    {
        stop_reason = run_simulation_loop(g, rng, ofile, config.runtime.seed, timer1, sweep, config.simulation.ks0, stats, settings);
    }
    finalize_simulation(g, rng, ofile, config.runtime.seed, timer1, sweep, stats, settings, stop_reason);

    fclose(ofile);
    gsl_rng_free(rng);
    return 0;
}
