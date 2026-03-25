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

gsl_rng *r;

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
    const std::string restart_path = resolve_path_from(launch_dir, config.restartPath);
    const std::string output_dir = resolve_path_from(launch_dir, config.outputDir);

    if (config.initMode == "restart" && access(restart_path.c_str(), R_OK) != 0)
    {
        std::cerr << "Restart file is not readable: " << restart_path << std::endl;
        return -1;
    }

    std::string path_error;
    if (!ensure_directory_exists(output_dir, path_error))
    {
        std::cerr << path_error << std::endl;
        return -1;
    }

    if (chdir(output_dir.c_str()) != 0)
    {
        std::cerr << "Failed to enter output directory '" << output_dir << "': " << std::strerror(errno) << std::endl;
        return -1;
    }

    const gsl_rng_type *t = gsl_rng_taus2;
    r = gsl_rng_alloc(t);
    srand((unsigned)config.seed);
    gsl_rng_set(r, config.seed);
    cout << "HERE " << endl;

    geometry g;
    apply_simulation_config(config, g);

    cout << "l_thermal_sigma is " << g.l_thermal_sigma << endl;
    cout << "l_thermal_kappa is " << g.l_thermal_kappa << endl;
    cout << "theta_thermal_kappa is " << g.theta_thermal_kappa << endl;
    cout << "gaussian sigma " << g.gaussian_sigma << endl;

    FILE *ofile;
    FILE *fi;
    FILE *paramfile;
    g.dump_parameters();
    ofile = fopen("energy.dat", "a");
    if (access("energy.dat", F_OK) != -1)
    {
        fprintf(stderr, " log files exist\n");
        cout << " file openned" << endl;
    }
    else
    {
        paramfile = fopen("allparam.dat", "a");
        fprintf(paramfile, "# seed, g.epsilon[0], g.kappa[0], g.theta0[0],g.theta0[1], g.gb0, g.mu[0], g.mu[1], dgother,  dg12, dg01,dg20,dg33,dg00 , ks0,g.theta0[2],kd0,gdrug0, g.mudrug \n");
        fprintf(paramfile, "%lu %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f  %.6f %.3f %.6f %.3f %.3f",
                config.seed, g.epsilon[0], g.kappa[0], g.theta0[0], g.theta0[1], g.gb0, g.mu[0], g.mu[1], config.dgother, config.dg12, config.dg01, config.dg20, config.dg33, config.dg00, config.ks0, g.theta0[2], config.kd0, config.gdrug0, g.mudrug);
        fclose(paramfile);
    }

    fi = fopen("parameters_run.out", "a");
    fprintf(fi, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 dmu dummydg mudrug gdrug kd0 dg12 dg01 dg20 dg33 dg00 dgother [--index-capacity N] [--max-sweeps N] [--init MODE] [--seed-config NAME] [--restart PATH] [--output-dir PATH]\n");
    fprintf(fi, "./source/assemble %lu %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.3f %.3f --index-capacity %lu --max-sweeps %lu --init %s --seed-config %s --restart %s --output-dir %s\n",
            config.seed, g.epsilon[0], g.kappa[0], config.kappaPhi0, g.theta0[0], g.theta0[1], g.gb0, g.mu[0], config.ks0, config.dmu, g.dg, g.mudrug, config.gdrug0, config.kd0, config.dg12, config.dg01, config.dg20, config.dg33, config.dg00, config.dgother, config.indexCapacity, config.maxSweeps, config.initMode.c_str(), config.seedConfig.c_str(), restart_path.c_str(), output_dir.c_str());

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
    fprintf(stderr, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK LnZ ks0 muAB mudrugdrugProb kd0 sweep [--index-capacity N] [--max-sweeps N] [--init MODE] [--seed-config NAME] [--restart PATH] [--output-dir PATH]\n");
    fprintf(stderr, "./source/assemble %lu %f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %lu --index-capacity %lu --max-sweeps %lu --init %s --seed-config %s --restart %s --output-dir %s\n",
            config.seed, g.epsilon[0], g.kappa[0], config.kappaPhi0, g.theta0[0], g.theta0[1], g.gb0, g.mu[0], config.ks0, g.mu[1], g.dg, g.mudrug, config.gdrug0, config.kd0, sweep, config.indexCapacity, config.maxSweeps, config.initMode.c_str(), config.seedConfig.c_str(), restart_path.c_str(), output_dir.c_str());

    for (int j = 0; j < 4; j++)
    {
        fprintf(stderr, "%d  g.epsilon %f g.kappa %f g.kappaPhi %f g.l0 %f  g.theta0 %f g.phi0 %f\n", j, g.epsilon[j], g.kappa[j], g.kappaPhi[j], g.l0[j], g.theta0[j], g.phi0[j]);
    }

    SimulationRunStats stats;
    SimulationLoopSettings settings = make_simulation_loop_settings(config);
    if (config.initMode == "restart")
    {
        initialize_from_restart(g, r, restart_path.c_str(), sweep, stats, settings);
    }
    else
    {
        initialize_from_seed(g, r, config.seedConfig.c_str(), sweep, stats, settings);
    }
    SimulationStopReason stop_reason = run_simulation_loop(g, r, ofile, config.seed, timer1, sweep, config.ks0, stats, settings);
    finalize_simulation(g, r, ofile, config.seed, timer1, sweep, stats, settings, stop_reason);

    fclose(ofile);
    gsl_rng_free(r);
    return 0;
}
