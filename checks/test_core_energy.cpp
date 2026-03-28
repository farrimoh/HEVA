#include "geometry.hpp"
#include "simulation_config.hpp"
#include "simulation_setup.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

namespace
{
SimulationConfig make_base_config()
{
    SimulationConfig config;
    config.capsidGeometry.epsilon0 = 4200.0;
    config.capsidGeometry.kappa0 = 40.0;
    config.capsidGeometry.kappaPhi0 = 800.0;
    config.capsidGeometry.theta0 = 0.24;
    config.capsidGeometry.theta1 = 0.48;
    config.capsidGeometry.theta2 = 0.24;
    config.capsidGeometry.theta3 = 0.24;
    config.capsidGeometry.l0 = 1.05;
    config.capsidGeometry.l1 = 0.95;
    config.capsidGeometry.phi33 = 1.05;
    config.capsidGeometry.phi12 = 1.17;
    config.capsidGeometry.phi01 = 0.98;
    config.capsidGeometry.phi20 = 1.05;
    config.capsidGeometry.gb0 = -9.48;
    config.capsidGeometry.dg12 = 0.3;
    config.capsidGeometry.dg01 = 0.1;
    config.capsidGeometry.dg20 = -0.1;
    config.capsidGeometry.dg33 = 0.0;
    config.capsidGeometry.dg00 = -0.8;
    config.capsidGeometry.dgother = -0.95;

    config.simulation.muCd = -10.8;
    config.simulation.ks0 = 0.02;
    config.simulation.dmu = -3.5;
    config.simulation.dg = 0.1;

    config.drug.mudrug = 0.0;
    config.drug.gdrug0 = 0.0;
    config.drug.kd0 = 0.0;

    config.core.enabled = false;
    config.core.maxBonds = 1000UL;
    config.core.epsilonLJ = 1.5;
    config.core.sigmaLJ = 1.2;
    config.engine.indexCapacity = 1000000UL;
    config.engine.runMode = "extended";
    config.runtime.seed = 7UL;
    config.runtime.maxSweeps = 0UL;
    config.runtime.outputDir = ".";
    config.runtime.workflow = "assembly";
    config.runtime.resume = false;
    config.initialization.mode = "seed";
    config.initialization.seedConfig = "triangle";
    config.initialization.path.clear();
    config.cgParamOpt.sampleStartSweep = 0UL;
    config.cgParamOpt.sampleEvery = 0UL;
    config.cgParamOpt.sampleOutputPath = "data.dat";
    return config;
}

double total_core_binding(geometry &g)
{
    double total = 0.0;
    for (std::vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
    {
        g.update_half_edge(it->id);
        total += g.RNA_bind_energy(g.heidtoindex[it->id]);
    }
    return total;
}
}

int main()
{
    SimulationConfig baseline_config = make_base_config();
    geometry baseline;
    apply_simulation_config(baseline_config, baseline);
    make_initial_triangle(baseline);
    baseline.update_boundary();
    baseline.update_normals();
    const double baseline_energy = baseline.compute_energy();

    SimulationConfig core_config = make_base_config();
    core_config.core.enabled = true;
    geometry with_core;
    apply_simulation_config(core_config, with_core);
    make_initial_triangle(with_core);
    with_core.update_boundary();
    with_core.update_normals();

    const double binding_total = total_core_binding(with_core);
    const double core_energy = with_core.compute_energy();

    assert(std::fabs(binding_total) > 1e-9);
    assert(std::fabs((core_energy - baseline_energy) - binding_total) < 1e-9);

    std::cout << "core energy checks passed" << std::endl;
    return 0;
}
