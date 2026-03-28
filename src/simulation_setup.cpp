#include "simulation_setup.hpp"

#include "geometry.hpp"

#include <cmath>

void apply_simulation_config(const SimulationConfig &config, geometry &g)
{
    g.initialize(4, config.engine.indexCapacity);
    g.all_neigh = 0;
    g.T = 1;
    g.xi = .5;
    g.Nd = 0;

    g.l0[0] = config.capsidGeometry.l0;
    g.l0[1] = config.capsidGeometry.l1;
    g.l0[2] = .95;
    g.l0[3] = 1.05;

    g.phi0[0] = config.capsidGeometry.phi33;
    g.phi0[1] = config.capsidGeometry.phi12;
    g.phi0[2] = config.capsidGeometry.phi01;
    g.phi0[3] = config.capsidGeometry.phi20;

    g.epsilon[0] = config.capsidGeometry.epsilon0;
    g.epsilon[1] = g.epsilon[0];
    g.epsilon[2] = g.epsilon[0];
    g.epsilon[3] = g.epsilon[0];

    g.kappa[0] = config.capsidGeometry.kappa0;
    g.kappa[1] = g.kappa[0];
    g.kappa[2] = g.kappa[0];
    g.kappa[3] = g.kappa[0];

    g.kappaPhi[0] = config.capsidGeometry.kappaPhi0;
    g.kappaPhi[1] = config.capsidGeometry.kappaPhi0;
    g.kappaPhi[2] = config.capsidGeometry.kappaPhi0;
    g.kappaPhi[3] = config.capsidGeometry.kappaPhi0;

    g.theta0[0] = config.capsidGeometry.theta0;
    g.theta0[1] = config.capsidGeometry.theta1;
    g.theta0[2] = config.capsidGeometry.theta2;
    g.theta0[3] = config.capsidGeometry.theta3;

    g.gb0 = config.capsidGeometry.gb0;
    g.mu[0] = config.simulation.muCd;
    g.mu[3] = g.mu[0];
    g.mu[1] = g.mu[0] + config.simulation.dmu;
    g.mu[2] = g.mu[1];

    g.dg = config.simulation.dg;
    g.mudrug = config.drug.mudrug;
    g.dmud = config.drug.dmud;
    g.core_enabled = config.core.enabled;
    g.maxbondsRNA = static_cast<int>(config.core.maxBonds);
    g.epsilon_lj = config.core.epsilonLJ;
    g.sigma_lj = config.core.sigmaLJ;

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
                g.gb[i][j] = (1 + config.capsidGeometry.dg12) * g.gb0;
            else if (i == 0 && j == 1)
                g.gb[i][j] = (1 + config.capsidGeometry.dg01) * g.gb0;
            else if (i == 3 && j == 1)
                g.gb[i][j] = (1 + config.capsidGeometry.dg01) * g.gb0;
            else if (i == 2 && j == 0)
                g.gb[i][j] = (1 + config.capsidGeometry.dg20) * g.gb0;
            else if (i == 2 && j == 3)
                g.gb[i][j] = (1 + config.capsidGeometry.dg20) * g.gb0;
            else if (i == 3 && j == 3)
                g.gb[i][j] = (1 + config.capsidGeometry.dg33) * g.gb0;
            else if (i == 0 && j == 0)
                g.gb[i][j] = (1 + config.capsidGeometry.dg00) * g.gb0;
            else if (i == 0 && j == 3)
                g.gb[i][j] = (1 + config.capsidGeometry.dg00) * g.gb0;
            else if (i == 3 && j == 0)
                g.gb[i][j] = (1 + config.capsidGeometry.dg00) * g.gb0;
            else
                g.gb[i][j] = (1 + config.capsidGeometry.dgother) * g.gb0;
        }
    }

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {
            if (j == 3 || j == 0)
                g.gdrug[i][j] = config.drug.gdrug0;
            else
                g.gdrug[i][j] = 0;
        }
    }
}
