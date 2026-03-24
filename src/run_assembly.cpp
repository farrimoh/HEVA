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

#include <iostream>
#include <string>
#include <time.h>
#include <unistd.h>

using namespace std;

gsl_rng *r;

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
    fprintf(fi, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 dmu dummydg mudrug gdrug kd0 dg12 dg01 dg20 dg33 dg00 dgother [indexCapacity] [maxSweeps]\n");
    fprintf(fi, "./source/assemble %lu %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.3f %.3f %lu %lu\n",
            config.seed, g.epsilon[0], g.kappa[0], config.kappaPhi0, g.theta0[0], g.theta0[1], g.gb0, g.mu[0], config.ks0, config.dmu, g.dg, g.mudrug, config.gdrug0, config.kd0, config.dg12, config.dg01, config.dg20, config.dg33, config.dg00, config.dgother, config.indexCapacity, config.maxSweeps);

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
    fprintf(stderr, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK LnZ ks0 muAB mudrugdrugProb kd0 sweep [indexCapacity] [maxSweeps]\n");
    fprintf(stderr, "./source/assemble %lu %f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %lu %lu %lu\n",
            config.seed, g.epsilon[0], g.kappa[0], config.kappaPhi0, g.theta0[0], g.theta0[1], g.gb0, g.mu[0], config.ks0, g.mu[1], g.dg, g.mudrug, config.gdrug0, config.kd0, sweep, config.indexCapacity, config.maxSweeps);

    for (int j = 0; j < 4; j++)
    {
        fprintf(stderr, "%d  g.epsilon %f g.kappa %f g.kappaPhi %f g.l0 %f  g.theta0 %f g.phi0 %f\n", j, g.epsilon[j], g.kappa[j], g.kappaPhi[j], g.l0[j], g.theta0[j], g.phi0[j]);
    }

    SimulationRunStats stats;
    SimulationLoopSettings settings = make_simulation_loop_settings(config);
    const char filename[] = "restart_lammps.dat";

    initialize_from_restart(g, r, filename, sweep, stats, settings);
    SimulationStopReason stop_reason = run_simulation_loop(g, r, ofile, config.seed, timer1, sweep, config.ks0, stats, settings);
    finalize_simulation(g, r, ofile, config.seed, timer1, sweep, stats, settings, stop_reason);

    fclose(ofile);
    gsl_rng_free(r);
    return 0;
}
