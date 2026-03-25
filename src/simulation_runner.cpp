#include "simulation_runner.hpp"

#include "geometry.hpp"
#include "montecarlo.hpp"
#include "seed_initialization.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

SimulationRunStats::SimulationRunStats()
    : frame(0),
      monomeradded(0),
      dimeradded(0),
      monomerremoved(0),
      dimerremoved(0),
      drugadded(0),
      drugremoved(0),
      typechanged(0),
      fusion(0),
      fission(0),
      wedgefusion(0),
      wedgefission(0),
      boundtri(0),
      binding(0),
      unbinding(0),
      deletednorate(0),
      lastNhe(0),
      lastNheGrowth(0),
      npace(0),
      avgpace(0.0),
      avgAddInterval(10000)
{
}

SimulationLoopSettings::SimulationLoopSettings()
    : minHEUpdateNeigh(150),
      minhe_fission(50),
      freq_vis(100),
      freq_log(100),
      freq_out(1000),
      initial_equilibration_steps(1000),
      final_equilibration_steps(1000),
      maxSweeps(0UL)
{
}

namespace
{
int elapsed_seconds(time_t start_time)
{
    time_t timer2;
    time(&timer2);
    return static_cast<int>(difftime(timer2, start_time));
}

void print_energy_state(const geometry &g, double energy, int frame, const char *label)
{
    cout << "# ENERGY " << label << " " << energy << " # ENERGY PER DIMER " << 2 * energy / g.Nhe << " # FRAME " << frame << " #" << endl;
    cout << "# NHE " << g.Nhe << " # NHESURF " << g.boundary.size() << " # NVSURF " << g.boundaryv.size() << " # NV_BONDSURF " << g.boundaryvbond.size() << " #" << endl;
}

double write_log_snapshot(geometry &g, FILE *ofile, unsigned long sweep, unsigned long seed, time_t start_time)
{
    if (g.Nhe == 6)
    {
        recenter(g);
    }

    g.update_boundary();
    g.check_odd_neigh();

    double energy = g.compute_energy();
    int seconds = elapsed_seconds(start_time);

    dump_analysis(g, ofile, sweep, seed, seconds);
    dump_lammps_data_file(g, 22222222);
    dump_lammps_data_dimers(g, 11111111);
    dump_restart_lammps_data_file(g, sweep);
    return energy;
}

void print_periodic_summary(const geometry &g, double energy, unsigned long sweep, time_t start_time, const SimulationRunStats &stats)
{
    int seconds = elapsed_seconds(start_time);
    cout << "###################################################################" << endl;
    cout << " ################ RUN TIME " << seconds << " SECONDS ###############" << endl;
    cout << " ############# SWEEP " << sweep << "##############" << endl;
    cout << "######### FRAME " << stats.frame << " ##############" << endl;
    cout << "######### ENERGY " << energy << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * energy / g.Nhe << " ##############" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << g.boundary.size() << " ##############" << endl;
    cout << "#########  NVSURF " << g.boundaryv.size() << " ##############" << endl;
    cout << "#########  NV_BONDSURF " << g.boundaryvbond.size() << " ##############" << endl;
    cout << "#########  NV5 " << g.Nv5 << " ##############" << endl;
    cout << "#########  MONOMER ADDED " << stats.monomeradded << " ##############" << endl;
    cout << "#########  MONOMER REMOVED " << stats.monomerremoved << " ##############" << endl;
    cout << "#########  DIMER ADDED " << stats.dimeradded << " ##############" << endl;
    cout << "#########  DIMER REMOVED " << stats.dimerremoved << " ##############" << endl;
    cout << "#########  no rate REMOVED " << stats.deletednorate << " ##############" << endl;
    cout << "#########  Surface bound " << stats.binding << " ##############" << endl;
    cout << "#########  Surface Unbound " << stats.unbinding << " ##############" << endl;
    cout << "#########  DrugAdded " << stats.drugadded << " ##############" << endl;
    cout << "#########  DrugRemoved " << stats.drugremoved << " ##############" << endl;
    cout << "#########  ND " << g.Nd << " ##############" << endl;
    cout << "#########  TYPE CHANGED " << stats.typechanged << " ##############" << endl;
    cout << "#########  WEDGE FUSION " << stats.wedgefusion << " ##############" << endl;
    cout << "#########  WEDGE FISSION " << stats.wedgefission << " ##############" << endl;
    cout << "#########  FUSION " << stats.fusion << " ##############" << endl;
    cout << "#########  FISSION " << stats.fission << " ##############" << endl;
    cout << "#########  FUSION HALFEDGES " << g.fusionhe.size() << " ##############" << endl;
    cout << "#########  WEDGE FUSION HALFEDGES " << g.fusionwedgehe.size() << " ##############" << endl;
    cout << "#########  ALL NEIGH " << g.all_neigh << " ##############" << endl;
    cout << "#########  Nboundary " << g.Nboundary << " ##############" << endl;
    cout << "#########  Bound Triangle " << stats.boundtri << " ##############" << endl;
    cout << "#########  Acceptance vmove " << (1.0 * g.accepted_vmove) / (1.0 * (g.accepted_vmove + g.rejected_vmove)) << "#################" << endl;
    cout << "#########  Nvlast " << g.Nvlast << " Nhelast " << g.Nhelast << " ################" << endl;
    cout << "#########  T4 " << g.NCD_T4_in << " T3 " << g.NCD_T3_in << " ################" << endl;
    cout << "#########  avgAddInterval " << stats.avgAddInterval << " ###############" << endl;
}

void print_final_summary(const geometry &g, double energy, unsigned long sweep, time_t start_time, const SimulationRunStats &stats, const char *stop_label)
{
    int seconds = elapsed_seconds(start_time);
    cout << " ################  FULL RUN TIME " << seconds << " SECONDS ###############" << endl;
    cout << " ################  STOP REASON " << stop_label << " ###############" << endl;
    cout << " ############# SWEEP " << sweep << "##############" << endl;
    cout << "######### ENERGY " << energy << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * energy / g.Nhe << " ##############" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << g.boundary.size() << " ##############" << endl;
    cout << "#########  NVSURF " << g.boundaryv.size() << " ##############" << endl;
    cout << "#########  NV_BONDSURF " << g.boundaryvbond.size() << " ##############" << endl;
    cout << "#########  NV5 " << g.Nv5 << " ##############" << endl;
    cout << "#########  MONOMER ADDED " << stats.monomeradded << " ##############" << endl;
    cout << "#########  MONOMER REMOVED " << stats.monomerremoved << " ##############" << endl;
    cout << "#########  DIMER ADDED " << stats.dimeradded << " ##############" << endl;
    cout << "#########  DIMER REMOVED " << stats.dimerremoved << " ##############" << endl;
    cout << "#########  no rate REMOVED " << stats.deletednorate << " ##############" << endl;
    cout << "#########  Surface bound " << stats.binding << " ##############" << endl;
    cout << "#########  Surface Unbound " << stats.unbinding << " ##############" << endl;
    cout << "#########  DrugAdded " << stats.drugadded << " ##############" << endl;
    cout << "#########  DrugRemoved " << stats.drugremoved << " ##############" << endl;
    cout << "#########  ND " << g.Nd << " ##############" << endl;
    cout << "#########  TYPE CHANGED " << stats.typechanged << " ##############" << endl;
    cout << "#########  WEDGE FUSION " << stats.wedgefusion << " ##############" << endl;
    cout << "#########  WEDGE FISSION " << stats.wedgefission << " ##############" << endl;
    cout << "#########  FUSION " << stats.fusion << " ##############" << endl;
    cout << "#########  FISSION " << stats.fission << " ##############" << endl;
    cout << "#########  FUSION HALFEDGES " << g.fusionhe.size() << " ##############" << endl;
    cout << "#########  WEDGE FUSION HALFEDGES " << g.fusionwedgehe.size() << " ##############" << endl;
    cout << "#########  ALL NEIGH " << g.all_neigh << " ##############" << endl;
    cout << "#########  Nboundary " << g.Nboundary << " ##############" << endl;
    cout << "#########  Bound Triangle " << stats.boundtri << " ##############" << endl;
}

const char *stop_reason_label(SimulationStopReason stop_reason)
{
    switch (stop_reason)
    {
    case SIMULATION_STOP_CLOSED:
        return "closed";
    case SIMULATION_STOP_MAX_SWEEPS:
        return "max_sweeps";
    case SIMULATION_STOP_OVERLAP_ERROR:
        return "overlap_error";
    case SIMULATION_STOP_MIXED_MORPH:
        return "mixed_morph";
    case SIMULATION_STOP_STALLED_GROWTH:
        return "stalled_growth";
    case SIMULATION_STOP_TOO_LONG:
        return "too_long";
    case SIMULATION_STOP_TOO_LARGE:
        return "too_large";
    default:
        return "unknown";
    }
}

void write_stop_snapshot(geometry &g, FILE *ofile, unsigned long sweep, unsigned long seed, time_t start_time)
{
    int seconds = elapsed_seconds(start_time);
    dump_analysis(g, ofile, sweep, seed, seconds);
    dump_restart_lammps_data_file(g, sweep);
}
}

SimulationLoopSettings make_simulation_loop_settings(const SimulationConfig &config)
{
    SimulationLoopSettings settings;
    settings.maxSweeps = config.runtime.maxSweeps;
    settings.final_equilibration_steps = 10 * settings.freq_log;
    if (config.initialization.mode != "restart")
    {
        settings.initial_equilibration_steps = 0;
    }
    return settings;
}

void initialize_from_restart(geometry &g, gsl_rng *rng, const char *filename, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings)
{
    sweep = read_restart_lammps_data_file(g, const_cast<char *>(filename));

    g.update_boundary();
    g.update_neigh();

    double energy = g.compute_energy();
    print_energy_state(g, energy, stats.frame, "before equilibration");

    for (int rstep = 0; rstep < settings.initial_equilibration_steps; rstep++)
    {
        move_vertex(g, rng);
        g.update_boundary();
    }

    energy = g.compute_energy();
    print_energy_state(g, energy, stats.frame, "after equilibration");
    fprintf(stderr, "Graph initialized.\n");

    if (g.Nboundary == 1)
    {
        dump_restart_lammps_data_file(g, sweep);
    }
}

void initialize_from_seed(geometry &g, gsl_rng *rng, const char *seed_config, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings)
{
    sweep = 0;

    const std::string config = seed_config == nullptr ? "" : seed_config;
    if (!make_seed(g, rng, config))
    {
        std::cerr << "Failed to initialize seed configuration: " << config << std::endl;
        std::exit(-1);
    }

    g.update_boundary();
    g.update_neigh();

    double energy = g.compute_energy();
    print_energy_state(g, energy, stats.frame, "before equilibration");

    for (int rstep = 0; rstep < settings.initial_equilibration_steps; rstep++)
    {
        move_vertex(g, rng);
        g.update_boundary();
    }

    energy = g.compute_energy();
    print_energy_state(g, energy, stats.frame, "after equilibration");
    fprintf(stderr, "Graph initialized.\n");
    dump_restart_lammps_data_file(g, sweep);
}

SimulationStopReason run_simulation_loop(geometry &g, gsl_rng *rng, FILE *ofile, unsigned long seed, time_t start_time, unsigned long &sweep, double ks0, SimulationRunStats &stats, const SimulationLoopSettings &settings)
{
    while (g.Nsurf > 0)
    {
        if (settings.maxSweeps > 0 && sweep >= settings.maxSweeps)
        {
            return SIMULATION_STOP_MAX_SWEEPS;
        }

        int ind = 0;
        int e = 0;
        int ss = 0;
        int ssadd = 0;
        double ps_attempt = 0.0;

        if (g.Nhe == 6)
        {
            ind = gsl_rng_uniform_int(rng, g.boundary.size());
            int hh = g.boundary[ind];
            int x = attempt_change_edge_type_tri(g, hh, rng);
            if (x >= 0)
            {
                stats.typechanged++;
            }
        }
        else
        {
            for (int nc = 0; nc < g.Nhe / 2; nc++)
            {
                int ind1 = gsl_rng_uniform_int(rng, g.Nhe);
                int e1 = g.he[ind1].id;
                int x = attempt_change_edge_type(g, e1, rng);
                if (x >= 0)
                {
                    stats.typechanged++;
                }
            }
        }

        ps_attempt = ks0 * g.Nsurf;
        if (gsl_rng_uniform(rng) < ps_attempt)
        {
            ind = gsl_rng_uniform_int(rng, g.boundary.size());
            e = g.boundary[ind];
            if (g.check_inside_overlap(e) > 0)
            {
                ssadd = attempt_add_monomer_dimer(g, e, rng);
                if (ssadd > 1)
                {
                    stats.dimeradded++;
                }
                else if (ssadd > 0)
                {
                    stats.monomeradded++;
                }
                g.update_boundary();
            }
        }

        move_vertex(g, rng);

        if (g.Nhe > 6)
        {
            ps_attempt = ks0 * g.Nsurf;
            if (gsl_rng_uniform(rng) < ps_attempt)
            {
                if (sweep > 0 && g.Nhe > 6)
                {
                    ind = gsl_rng_uniform_int(rng, g.boundary.size());
                    e = g.boundary[ind];
                    if (g.no_bond_boundary(e) > 0)
                    {
                        ss = attempt_remove_monomer_dimer(g, e, rng);
                        if (ss > 1)
                        {
                            stats.dimerremoved++;
                        }
                        else if (ss > 0)
                        {
                            stats.monomerremoved++;
                        }
                        g.update_boundary();
                    }
                }
            }

            if (g.Nhe > 15)
            {
                int cc = check_bind_triangle(g);
                if (cc > 0)
                {
                    cout << "bound triangle" << endl;
                    g.update_boundary();
                    stats.boundtri += cc;
                }
                else
                {
                    ind = gsl_rng_uniform_int(rng, g.boundary.size());
                    int hh = g.boundary[ind];
                    int tt = -1;
                    if (g.no_bond_boundary(hh) > 0)
                    {
                        tt = attempt_bind_wedge_dimer(g, hh, rng);
                        if (tt > 0)
                        {
                            stats.binding++;
                        }
                    }
                    g.update_boundary();

                    if (g.boundaryvbond.size() > 0)
                    {
                        ind = gsl_rng_uniform_int(rng, g.boundary.size());
                        hh = g.boundary[ind];
                        if ((g.is_bond_in_boundary(hh) > 0) || (g.is_bond_out_boundary(hh) > 0))
                        {
                            tt = attempt_unbind_wedge_dimer(g, hh, rng);
                            if (tt > 0)
                            {
                                stats.unbinding++;
                            }
                            g.update_boundary();
                        }
                    }

                    if (g.Nhe > settings.minhe_fission && g.Nsurf > 3)
                    {
                        int movetype = gsl_rng_uniform_int(rng, 4);
                        switch (movetype)
                        {
                        case 0:
                            if (g.all_neigh > 0)
                            {
                                int ff = attempt_wedge_fusion(g, rng);
                                g.update_boundary();
                                if (ff > 0)
                                {
                                    stats.wedgefusion++;
                                }
                            }
                            break;

                        case 1:
                        {
                            int ff = attempt_wedge_fission(g, rng);
                            g.update_boundary();
                            if (ff > 0)
                            {
                                stats.wedgefission++;
                            }
                        }
                        break;

                        case 2:
                            if (g.all_neigh > 0)
                            {
                                int ff = attempt_fusion(g, rng);
                                g.update_boundary();
                                if (ff > 0)
                                {
                                    stats.fusion++;
                                }
                            }
                            break;

                        case 3:
                        {
                            int ff = attempt_fission(g, rng);
                            g.update_boundary();
                            if (ff > 0)
                            {
                                stats.fission++;
                            }
                        }
                        break;
                        }
                    }
                }
            }
        }

        double energy = 0.0;
        if (sweep % settings.freq_log == 0)
        {
            energy = write_log_snapshot(g, ofile, sweep, seed, start_time);
        }

        if (sweep % settings.freq_out == 0)
        {
            print_periodic_summary(g, energy, sweep, start_time, stats);
        }

        int cc = check_bind_triangle(g);
        if (cc > 0)
        {
            cout << "bound triangle" << endl;
            g.update_boundary();
            stats.boundtri += cc;
        }

        g.check_odd_neigh();

        if (g.Nhe > 70 && g.Nhe < settings.minHEUpdateNeigh && sweep % 10000 == 0)
        {
            double thispace = float(g.Nhe - stats.lastNhe) / 10000.0;
            cout << "thispace " << thispace << endl;

            stats.avgpace = (stats.npace * stats.avgpace + thispace) / (stats.npace + 1.0);
            stats.avgAddInterval = int(pow(10, (-1 * int(floor(log10(stats.avgpace))))));
            cout << "avgpace" << stats.avgpace << " avgAddInterval " << stats.avgAddInterval << endl;

            stats.npace++;
            stats.lastNhe = g.Nhe;
        }

        if (g.Nhe > settings.minHEUpdateNeigh && sweep % 10000 == 0)
        {
            g.update_neigh();
            if (g.find_overlap_all() < 0)
            {
                cout << "error overlap" << endl;
                dump_lammps_data_dimers(g, 5555555);
                return SIMULATION_STOP_OVERLAP_ERROR;
            }
        }

        if (g.Nhe > settings.minHEUpdateNeigh && sweep % (10 * stats.avgAddInterval) == 0)
        {
            update_geometry_parameters(g);
            if (g.Nhe - stats.lastNhe <= 2)
            {
                if ((g.Nhe >= 220 && g.NCD_T4_in >= 26 && g.NCD_T3_in >= 3 && g.Nsurf > 10) ||
                    (g.Nhe >= 160 && g.NCD_T4_in >= 3 && g.NCD_T3_in >= 16 && g.Nsurf > 10) ||
                    (g.Nhe >= 200 && g.NCD_T4_in >= 5 && g.NCD_T3_in >= 5 && g.Nsurf > 10))
                {
                    cout << "STOP for now - mixed morph" << endl;
                    g.update_boundary();
                    dump_lammps_data_dimers(g, 44444444);
                    dump_lammps_data_dimers(g, 11111111);
                    write_stop_snapshot(g, ofile, sweep, seed, start_time);
                    return SIMULATION_STOP_MIXED_MORPH;
                }
            }

            if (sweep % (100 * stats.avgAddInterval) == 0)
            {
                if (g.NCD_T4_in > 0 && g.NCD_T3_in > 0 && abs(g.Nhe - stats.lastNheGrowth) <= 4)
                {
                    fprintf(stderr, "STOP for now - not growing\n");
                    g.update_boundary();
                    dump_lammps_data_dimers(g, 333333333);
                    dump_lammps_data_dimers(g, 11111111);
                    write_stop_snapshot(g, ofile, sweep, seed, start_time);
                    return SIMULATION_STOP_STALLED_GROWTH;
                }
                stats.lastNheGrowth = g.Nhe;
            }

            stats.lastNhe = g.Nhe;
        }

        if (sweep == 200000000)
        {
            fprintf(stderr, "STOP for now - too long\n");
            g.update_boundary();
            dump_lammps_data_dimers(g, 77777777);
            write_stop_snapshot(g, ofile, sweep, seed, start_time);
            return SIMULATION_STOP_TOO_LONG;
        }

        if (g.Nhe >= 310 || g.Nv >= 65)
        {
            fprintf(stderr, "STOP for now - too large\n");
            g.update_boundary();
            dump_lammps_data_dimers(g, 88888888);
            dump_lammps_data_dimers(g, 11111111);
            write_stop_snapshot(g, ofile, sweep, seed, start_time);
            return SIMULATION_STOP_TOO_LARGE;
        }

        sweep++;
    }

    return SIMULATION_STOP_CLOSED;
}

void finalize_simulation(geometry &g, gsl_rng *rng, FILE *ofile, unsigned long seed, time_t start_time, unsigned long &sweep, SimulationRunStats &stats, const SimulationLoopSettings &settings, SimulationStopReason stop_reason)
{
    if (stop_reason == SIMULATION_STOP_CLOSED)
    {
        dump_restart_lammps_data_file(g, sweep);

        for (int rstep = 0; rstep < settings.final_equilibration_steps; rstep++)
        {
            move_vertex(g, rng);
            g.update_boundary();

            int ind1 = gsl_rng_uniform_int(rng, g.Nhe);
            int e1 = g.he[ind1].id;
            int x = attempt_change_edge_type(g, e1, rng);
            if (x >= 0)
            {
                stats.typechanged++;
            }

            if (sweep % settings.freq_out == 0)
            {
                int seconds = elapsed_seconds(start_time);
                cout << "###################################################################" << endl;
                cout << " ################  RUN TIME " << seconds << " SECONDS ###############" << endl;
                cout << " ############# SWEEP " << sweep << "##############" << endl;

                double energy = g.compute_energy();
                cout << "######### ENERGY " << energy << " ##############" << endl;
                cout << "######### ENERGY PER DIMER " << 2 * energy / g.Nhe << " ##############" << endl;
            }

            if (sweep % settings.freq_log == 0)
            {
                int seconds = elapsed_seconds(start_time);
                dump_analysis(g, ofile, sweep, seed, seconds);
            }

            sweep++;
        }
    }
    else
    {
        fprintf(stderr, "STOP for now - %s\n", stop_reason_label(stop_reason));
        g.update_boundary();
    }

    dump_lammps_data_file(g, 22222222);
    dump_lammps_data_dimers(g, 11111111);
    dump_restart_lammps_data_file(g, sweep);

    double energy = g.compute_energy();
    print_final_summary(g, energy, sweep, start_time, stats, stop_reason_label(stop_reason));

    int seconds = elapsed_seconds(start_time);
    dump_analysis(g, ofile, sweep, seed, seconds);

    FILE *finalfile = fopen(output_file_path("last.dat").c_str(), "w");
    dump_analysis(g, finalfile, sweep, seed, seconds);
    fclose(finalfile);
}
