#include "simulation_setup.hpp"

#include "geometry.hpp"

#include <cmath>

void apply_simulation_config(const SimulationConfig &config, geometry &g)
{
    g.initialize(4);
    g.all_neigh = 0;

    g.epsilon[0] = config.epsilon0;
    g.epsilon[1] = g.epsilon[0];
    g.epsilon[2] = g.epsilon[0];
    g.epsilon[3] = g.epsilon[0];

    g.kappa[0] = config.kappa0;
    g.kappa[1] = g.kappa[0];
    g.kappa[2] = g.kappa[0];
    g.kappa[3] = g.kappa[0];

    g.kappaPhi[0] = config.kappaPhi0;
    g.kappaPhi[1] = config.kappaPhi0;
    g.kappaPhi[2] = config.kappaPhi0;
    g.kappaPhi[3] = config.kappaPhi0;

    g.theta0[0] = config.theta0;
    g.theta0[1] = config.theta1;
    g.theta0[2] = g.theta0[0];
    g.theta0[3] = g.theta0[0];

    g.gb0 = config.gb0;
    g.mu[0] = config.muCd;
    g.mu[3] = g.mu[0];
    g.mu[1] = g.mu[0] + config.dmu;
    g.mu[2] = g.mu[1];

    g.dg = config.dg;
    g.mudrug = config.mudrug;

    double alp = 1;
    g.l_thermal_kappa = std::sqrt((3.0 * g.l0[0] * g.l0[0] * alp * g.T / (2.0 * g.kappa[0])));
    g.theta_thermal_kappa = std::sqrt(2.0 * (alp * g.T / (g.kappa[0])));
    g.l_thermal_sigma = std::sqrt(2.0 * (alp * g.T / g.epsilon[0]));
    g.gaussian_sigma = 0.5 * g.l_thermal_kappa;

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {
            if (i == 1 && j == 2)
                g.gb[i][j] = (1 + config.dg12) * g.gb0;
            else if (i == 0 && j == 1)
                g.gb[i][j] = (1 + config.dg01) * g.gb0;
            else if (i == 3 && j == 1)
                g.gb[i][j] = (1 + config.dg01) * g.gb0;
            else if (i == 2 && j == 0)
                g.gb[i][j] = (1 + config.dg20) * g.gb0;
            else if (i == 2 && j == 3)
                g.gb[i][j] = (1 + config.dg20) * g.gb0;
            else if (i == 3 && j == 3)
                g.gb[i][j] = (1 + config.dg33) * g.gb0;
            else if (i == 0 && j == 0)
                g.gb[i][j] = (1 + config.dg00) * g.gb0;
            else if (i == 0 && j == 3)
                g.gb[i][j] = (1 + config.dg00) * g.gb0;
            else if (i == 3 && j == 0)
                g.gb[i][j] = (1 + config.dg00) * g.gb0;
            else
                g.gb[i][j] = (1 + config.dgother) * g.gb0;
        }
    }

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {
            if (j == 3 || j == 0)
                g.gdrug[i][j] = config.gdrug0;
            else
                g.gdrug[i][j] = 0;
        }
    }

    g.l0[0] = 1.05;
    g.l0[1] = .95;
    g.l0[2] = .95;
    g.l0[3] = 1.05;

    g.phi0[0] = 1.05;
    g.phi0[1] = 1.17;
    g.phi0[2] = .98;
    g.phi0[3] = 1.05;

    g.xi = .5;
    g.T = 1;
    g.Nd = 0;
}
