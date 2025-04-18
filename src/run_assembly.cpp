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
#include "montecarlo.hpp"
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#define PI 3.14159265
using namespace std;

gsl_rng *r;
int main(int argc, char **argv)
{
    time_t timer1, timer2;
    int seconds;
    time(&timer1);

    if (argc != 21)
    {
        fprintf(stderr, "%d", argc);
        fprintf(stderr, "usage: ./assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 dmu dummydg mudrug gdrug kd0 dg12 dg01 dg20 dg33 dg00 dgother\n");
        exit(-1);
    }

    int minHE_update_neigh = 150;
    int lastNhe = 0;
    int lastNheGrowth=0;
    int npace=0;
    double avgpace=0;
    int avgAddInterval=10000;

    // initialize the rng
    const gsl_rng_type *t;
    t = gsl_rng_taus2;
    r = gsl_rng_alloc(t);
    long unsigned int seed = atoi(argv[1]);
    srand((unsigned)seed);
    gsl_rng_set(r, seed);
    cout << "HERE " << endl;
    //const double pi = 4*atan(1);
    geometry g;
    g.initialize(4);
    g.all_neigh = 0;

    

    g.epsilon[0] = atof(argv[2]);
    g.epsilon[1] = g.epsilon[0];
    g.epsilon[2] = g.epsilon[0];
    g.epsilon[3] = g.epsilon[0];
    g.kappa[0] = atof(argv[3]);
    g.kappa[1] = g.kappa[0];
    g.kappa[2] = g.kappa[0];
    g.kappa[3] = g.kappa[0];
    double Phikappa = atof(argv[4]);
    g.kappaPhi[0] = Phikappa; // CD-CD DC-DC and all other
    g.kappaPhi[1] = Phikappa; //BA-AB
    g.kappaPhi[2] = Phikappa; //AB-CD and AB-DC
    g.kappaPhi[3] = Phikappa; //with drug
    g.theta0[0] = atof(argv[5]);
    g.theta0[1] = atof(argv[6]);
    g.theta0[2] = g.theta0[0]; //0.1;  // CD-CD
    g.theta0[3] = g.theta0[0]; //0.1 ; //349 ;

    g.gb0 = atof(argv[7]);

    /* MU parameters */
    double dmu = atof(argv[10]);
    g.mu[0] = atof(argv[8]);
    g.mu[3] = g.mu[0];
    g.mu[1] = g.mu[0] + dmu;
    g.mu[2] = g.mu[1];

    double ks0 = atof(argv[9]);

    g.dg = atof(argv[11]);

    /* Drug parameteres */
    g.mudrug = atof(argv[12]);
    double gdrug0 = atof(argv[13]);
    double kd0 = atof(argv[14]);


    // gaussian
    double alp=1;
    
    g.l_thermal_kappa = sqrt((3.0*g.l0[0]*g.l0[0]*alp*g.T/(2.0*g.kappa[0])));
    g.theta_thermal_kappa = sqrt(2.0*(alp*g.T/(g.kappa[0])));
    g.l_thermal_sigma = sqrt(2.0*(alp*g.T/g.epsilon[0]));
    //g.l_thermal_sigma = g.l_thermal_kappa;//sqrt(2.0*(alp*g.T/g.epsilon[0]));
    g.gaussian_sigma = 0.5 * g.l_thermal_kappa; 

    cout << "l_thermal_sigma is " << g.l_thermal_sigma<<endl;
    cout << "l_thermal_kappa is " << g.l_thermal_kappa<<endl;
    cout << "theta_thermal_kappa is " << g.theta_thermal_kappa<<endl;
    cout << "gaussian sigma " << g.gaussian_sigma <<endl;
    //exit(-1);

    /* GB parameteres */
    double dg12 = atof(argv[15]); //BA-AB
    double dg01 = atof(argv[16]); //CD-BA
    double dg20 = atof(argv[17]); //AB-CD
    double dg33 = atof(argv[18]);
    double dg00 = atof(argv[19]);
    double dgother = atof(argv[20]);

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {

            if (i == 1 && j == 2) // BA-AB
                g.gb[i][j] = (1 + dg12) * g.gb0;
            else if (i == 0 && j == 1) // CD-AB (T3)
                g.gb[i][j] = (1 + dg01) * g.gb0;
            else if (i == 3 && j == 1) // DC-AB
                g.gb[i][j] = (1 + dg01) * g.gb0;
            else if (i == 2 && j == 0) // AB-DC
                g.gb[i][j] = (1 + dg20) * g.gb0;
            else if (i == 2 && j == 3) //
                g.gb[i][j] = (1 + dg20) * g.gb0;
            else if (i == 3 && j == 3) //
                g.gb[i][j] = (1 + dg33) * g.gb0;
            else if (i == 0 && j == 0) //
                g.gb[i][j] = (1 + dg00) * g.gb0;
            else if (i == 0 && j == 3)
                g.gb[i][j] = (1 + dg00) * g.gb0;
            else if (i == 3 && j == 0)
                g.gb[i][j] = (1 + dg00) * g.gb0;
            else
                g.gb[i][j] = (1 + dgother) * g.gb0;
        }
    }

    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {

            if (j == 3 || j == 0)
                g.gdrug[i][j] = gdrug0;
            else
                g.gdrug[i][j] = 0;
        }
    }

    long unsigned int sweep = 0;

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

    int ind, e;

    // set up an output file
    //ofstream *efile, *finalfile, *fi, *paramfile;
    FILE *ofile, *finalfile, *fi, *paramfile;
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
        // XXX Update this
        fprintf(paramfile, "# seed, g.epsilon[0], g.kappa[0], g.theta0[0],g.theta0[1], g.gb0, g.mu[0], g.mu[1], dgother,  dg12, dg01,dg20,dg33,dg00 , ks0,g.theta0[2],kd0,gdrug0, g.mudrug \n");
        fprintf(paramfile, "%lu %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f  %.6f %.3f %.6f %.3f %.3f",
                seed, g.epsilon[0], g.kappa[0], g.theta0[0], g.theta0[1], g.gb0, g.mu[0], g.mu[1], dgother, dg12, dg01, dg20, dg33, dg00, ks0, g.theta0[2], kd0, gdrug0, g.mudrug);
        fclose(paramfile);
    }
    fi = fopen("parameters_run.out", "a");
    fprintf(fi, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 dmu dummydg mudrug gdrug kd0 dg12 dg01 dg20 dg33 dg00 dgother\n");
    //fprintf(fi, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCD ks0 dmu dummydg mudrug gdrug kd0                                 \n");
    fprintf(fi, "./source/assemble %lu %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.3f %.3f\n",
            seed, g.epsilon[0], g.kappa[0], Phikappa, g.theta0[0], g.theta0[1], g.gb0, g.mu[0], ks0, dmu, g.dg, g.mudrug, gdrug0, kd0, dg12, dg01, dg20, dg33, dg00, dgother);

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

    double ee = 0;
    fprintf(stderr, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK LnZ ks0 muAB mudrugdrugProb kd0 sweep\n");
    fprintf(stderr, "./source/assemble %lu %f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %lu\n", seed, g.epsilon[0], g.kappa[0], Phikappa, g.theta0[0], g.theta0[1], g.gb0, g.mu[0], ks0, g.mu[1], g.dg, g.mudrug, gdrug0, kd0, sweep);

    for (int j = 0; j < 4; j++)
    {
        fprintf(stderr, "%d  g.epsilon %f g.kappa %f g.kappaPhi %f g.l0 %f  g.theta0 %f g.phi0 %f\n", j, g.epsilon[j], g.kappa[j], g.kappaPhi[j], g.l0[j], g.theta0[j], g.phi0[j]);
    }
    int frame = 0;
    int monomeradded = 0;
    int dimeradded = 0;
    int monomerremoved = 0;
    int dimerremoved = 0;
    int drugadded = 0;
    int drugremoved = 0;
    int typechanged = 0;
    int fusion = 0;
    int fission = 0;
    int wedgefusion = 0;
    int wedgefission = 0;
    int boundtri = 0;

    char filename[30] = "restart_lammps.dat";

    sweep = read_restart_lammps_data_file(g, filename);

    g.update_boundary();
    g.update_neigh();
    //make_seed(g,r);
    ee = g.compute_energy();
    cout << "# ENERGY before equilibration " << ee << " # ENERGY PER DIMER " << 2 * ee / g.Nhe << " # FRAME " << frame << " #" << endl;
    cout << "# NHE " << g.Nhe << " # NHESURF " << g.boundary.size() << " # NVSURF " << g.boundaryv.size() << " # NV_BONDSURF " << g.boundaryvbond.size() << " #" << endl;

    for (int rstep = 0; rstep < 1000; rstep++)
    { //relaxing the shell
        move_vertex(g, r);
        g.update_boundary();
    }

    ee = g.compute_energy();
    cout << "# ENERGY after equilibration " << ee << " # ENERGY PER DIMER " << 2 * ee / g.Nhe << " # FRAME " << frame << " #" << endl;
    cout << "# NHE " << g.Nhe << " # NHESURF " << g.boundary.size() << " # NVSURF " << g.boundaryv.size() << " # NV_BONDSURF " << g.boundaryvbond.size() << " #" << endl;
    //dump_lammps_data_file(g, 0);
    fprintf(stderr, "Graph initialized.\n");

    if (g.Nboundary == 1)
        dump_restart_lammps_data_file(g, sweep);
    //exit(-1);
    int binding = 0;
    int unbinding = 0;

    int deletednorate = 0;

    int ssadd = 0;
    //int ssremove = 0;

    int runhpc = 1;
    int test = 0;

    int freq_vis = 1000;
    int freq_log = 1000;
    int freq_out = 10000; // shoud be >=freq_log
    if (runhpc == 1)
    {
        freq_vis = 100;
        freq_log = 100;
        freq_out = 1000;
    }
    if (test == 1)
    {
        freq_vis = 1;
        freq_log = 1;
        freq_out = 1;
    }
    int minhe_fission = 50;
    //int kb0=1;

    double ps_attempt = 0;
    int ss = 0;
    while (g.Nsurf > 0)

    {

        //move all vertices

       //move_vertex(g, r);

        //g.check_odd_neigh();

        //change edge type

        if (g.Nhe==6){
            ind = gsl_rng_uniform_int(r, g.boundary.size());
                //cout <<"ind id" << ind<<endl;
            int hh = g.boundary[ind];
            int x = attempt_change_edge_type_tri(g, hh, r);
            if (x >= 0)
                typechanged++; 
        }
        else{
            for (int nc = 0; nc < g.Nhe/2; nc++)
            {
                int ind1 = gsl_rng_uniform_int(r, g.Nhe);
                int e1 = g.he[ind1].id;
                int x = -1;
                x = attempt_change_edge_type(g, e1, r);
                if (x >= 0)
                    typechanged++;
            }
        }

    
       
        // attempt add monomer_dimer
        
        //if (g.Nhe>6)
        //{
            ps_attempt = ks0 * g.Nsurf;
            if (gsl_rng_uniform(r) < ps_attempt)
            {
                ind = gsl_rng_uniform_int(r, g.boundary.size());
                e = g.boundary[ind];
                if (g.check_inside_overlap(e) > 0)
                {
                    ssadd = attempt_add_monomer_dimer(g, e, r);
                    if (ssadd > 1)
                        dimeradded++;
                    else if (ssadd > 0)
                        monomeradded++;
                    ssadd = -1;
                    g.update_boundary();
                }
            }
        /*}
        else {
                ind = gsl_rng_uniform_int(r, g.boundary.size());
                e = g.boundary[ind];
                ssadd = attempt_add_monomer_dimer(g, e, r);
                if (ssadd > 1)
                    dimeradded++;
                else if (ssadd > 0)
                    monomeradded++;
                ssadd = -1;
                g.update_boundary();
        }*/
        

        //double ps_move = ks0 * g.Nv ;
        //if (gsl_rng_uniform(r) < ps_move)
        //{
        move_vertex(g, r);
        //}


        if (g.Nhe>6) {


            


            //remove monomer dimer
            ps_attempt = ks0 * g.Nsurf;
            if (gsl_rng_uniform(r) < ps_attempt)
            {
                if (sweep > 0 && g.Nhe > 6)
                {
                    ind = gsl_rng_uniform_int(r, g.boundary.size());
                    e = g.boundary[ind];
                    if (g.no_bond_boundary(e) > 0)
                    {
                        ss = attempt_remove_monomer_dimer(g, e, r);
                        if (ss > 1)
                            dimerremoved++;
                        else if (ss > 0)
                            monomerremoved++;
                        ss = 0;
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
                    boundtri += cc;
                }
                else 
                {
                    //bind wedge
                    //if (gsl_rng_uniform(r) < pb_attempt){
                    ind = gsl_rng_uniform_int(r, g.boundary.size());
                    int hh = g.boundary[ind];
                    int tt = -1;
                    if (g.no_bond_boundary(hh) > 0)
                    {
                        tt = attempt_bind_wedge_dimer(g, hh, r);
                        if (tt > 0)
                            binding++;
                        tt = -1;
                    }
                    g.update_boundary();


                    //unbind wedge
                    if (g.boundaryvbond.size() > 0) // &&  gsl_rng_uniform(r) < pb_attempt)
                    {
                        ind = gsl_rng_uniform_int(r, g.boundary.size());
                        int hh = g.boundary[ind];
                        if ((g.is_bond_in_boundary(hh) > 0) || (g.is_bond_out_boundary(hh) > 0))
                        {
                            //cout <<" trying unbind vv is " << vv <<endl;
                            int tt = attempt_unbind_wedge_dimer(g, hh, r);
                            if (tt > 0)
                                unbinding++;
                            tt = -1;
                            g.update_boundary();
                        }
                    }


                    if (g.Nhe > minhe_fission && g.Nsurf > 3){ // dont try if only last triangle is open

                        int movetype = gsl_rng_uniform_int(r, 4);
                        switch (movetype)
                        {
                
                        case 0:
                            //wedge fusion
                            if (g.all_neigh > 0 )
                            {
                            //double kf0=1;
                            //double pf_attempt=kf0*g.Nsurf;
                            //if (gsl_rng_uniform(r) < pf_attempt)
                            //{
                            int ff = attempt_wedge_fusion(g, r);
                            g.update_boundary();
                            if (ff > 0)
                                wedgefusion++;
                            ff = 0;
                            //}
                            }
                            break;

                        case 1:
                            //wedge fission
                            {
                                //pf_attempt=kf0; //g.Nsurf
                                //if (gsl_rng_uniform(r) < pf_attempt)
                                //{
                                int ff = attempt_wedge_fission(g, r);
                                g.update_boundary();
                                if (ff > 0)
                                    wedgefission++;
                                ff = 0;
                                //}
                            }
                            break;

                        case 2:
                            //fusion
                            if (g.all_neigh > 0 )
                            {
                                //double kf0=1;
                                //double pf_attempt=kf0;//g.Nsurf;
                                //if (gsl_rng_uniform(r) < pf_attempt)
                                //{
                                int ff = attempt_fusion(g, r);
                                g.update_boundary();
                                if (ff > 0)
                                    fusion++;
                                ff = 0;
                            }
                            break;

                        case 3:
                            //fission
                            
                            {
                            //pf_attempt=kf0*.25; //g.Nsurf
                            //if (gsl_rng_uniform(r) < pf_attempt)
                            //{
                            int ff = attempt_fission(g, r);
                            g.update_boundary();
                            if (ff > 0)
                                fission++;
                            ff = 0;
                            //}
                            }
                            break;

                        }
                    }


                }

            }
        }
        // attempt add monomer_dimer
        //case 3:
            /*
            if (g.Nhe>6)
            {
                ps_attempt = ks0 * g.Nsurf;
            if (gsl_rng_uniform(r) < ps_attempt)
            {
                ind = gsl_rng_uniform_int(r, g.boundary.size());
                e = g.boundary[ind];
                if (g.check_inside_overlap(e) > 0)
                {
                    ssadd = attempt_add_monomer_dimer(g, e, r);
                    if (ssadd > 1)
                        dimeradded++;
                    else if (ssadd > 0)
                        monomeradded++;
                    ssadd = -1;
                    g.update_boundary();
                }
            }
            }
            else {
                ind = gsl_rng_uniform_int(r, g.boundary.size());
                e = g.boundary[ind];
                ssadd = attempt_add_monomer_dimer(g, e, r);
                if (ssadd > 1)
                    dimeradded++;
                else if (ssadd > 0)
                    monomeradded++;
                ssadd = -1;
                g.update_boundary();
            }
        //    break;

        double ps_move = ks0 * g.Nhe ;
        if (gsl_rng_uniform(r) < ps_move)
        {
            move_vertex(g, r);
        }
        */
        //if (g.Nhe>6){

       /*for (int ii=0; ii<10;ii++)
        {
        int movetype = gsl_rng_uniform_int(r, 10);
        //cout <<"try moves"<<endl;
        //if (g.boundary.size() <= 3 && sweep > 10000 && g.Nhe > 10)
        //    movetype = -1;
        switch (movetype)
        {
        //attempt bind wedge
        case 0:
            if (g.Nhe > 15) // &&  gsl_rng_uniform(r) < pb_attempt)
            {
                ind = gsl_rng_uniform_int(r, g.boundary.size());
                //cout <<"ind id" << ind<<endl;
                int hh = g.boundary[ind];
                int tt = -1;
                if (g.no_bond_boundary(hh) > 0)
                {
                    tt = attempt_bind_wedge_dimer(g, hh, r);
                    if (tt > 0)
                        binding++;
                    tt = -1;
                }
                g.update_boundary();
            }
            break;

        //attempt unbind wedge
        case 1:
            if (g.boundaryvbond.size() > 0) // &&  gsl_rng_uniform(r) < pb_attempt)
            {
                ind = gsl_rng_uniform_int(r, g.boundary.size());
                int hh = g.boundary[ind];
                if ((g.is_bond_in_boundary(hh) > 0) || (g.is_bond_out_boundary(hh) > 0))
                {
                    //cout <<" trying unbind vv is " << vv <<endl;
                    int tt = attempt_unbind_wedge_dimer(g, hh, r);
                    if (tt > 0)
                        unbinding++;
                    tt = -1;
                    g.update_boundary();
                }
            }
            break;

        // attempt remove monomer_dimer
        case 2:

            ps_attempt = ks0 * g.Nsurf ;
            if (gsl_rng_uniform(r) < ps_attempt)
            {
                if (sweep > 0 && g.Nhe > 6)
                {
                    ind = gsl_rng_uniform_int(r, g.boundary.size());
                    e = g.boundary[ind];
                    if (g.no_bond_boundary(e) > 0)
                    {
                        ss = attempt_remove_monomer_dimer(g, e, r);
                        if (ss > 1)
                            dimerremoved++;
                        else if (ss > 0)
                            monomerremoved++;
                        ss = 0;
                        g.update_boundary();
                    }
                }
            }
            break;

        

        //attempt wedge fusion
        case 3:
            if (g.all_neigh > 0 && g.Nhe > minhe_fission && g.Nsurf > 3)
            {
                //double kf0=1;
                //double pf_attempt=kf0;//g.Nsurf;
                //if (gsl_rng_uniform(r) < pf_attempt)
                //{
                int ff = attempt_wedge_fusion(g, r);
                g.update_boundary();
                if (ff > 0)
                    wedgefusion++;
                ff = 0;
                //}
            }
            break;

        // attempt_wedge fission
        case 4:
            if (g.all_neigh > 0 && g.Nhe > minhe_fission && g.Nsurf > 3)
            {
                //pf_attempt=kf0; //g.Nsurf
                //if (gsl_rng_uniform(r) < pf_attempt)
                //{
                int ff = attempt_wedge_fission(g, r);
                g.update_boundary();
                if (ff > 0)
                    wedgefission++;
                ff = 0;
                //}
            }
            break;

        //attempt fusion  with kg0=1
        case 5:
            if (g.all_neigh > 0 && g.Nhe > minhe_fission && g.Nsurf > 3)
            {
                //double kf0=1;
                //double pf_attempt=kf0;//g.Nsurf;
                //if (gsl_rng_uniform(r) < pf_attempt)
                //{
                int ff = attempt_fusion(g, r);
                g.update_boundary();
                if (ff > 0)
                    fusion++;
                ff = 0;
            }
            break;

        //attempt fission
        case 6:
            if (g.all_neigh > 0 && g.Nhe > minhe_fission && g.Nsurf > 3)
            {
                //pf_attempt=kf0*.25; //g.Nsurf
                //if (gsl_rng_uniform(r) < pf_attempt)
                //{
                int ff = attempt_fission(g, r);
                g.update_boundary();
                if (ff > 0)
                    fission++;
                ff = 0;
                //}
            }
            break;
        

        case 7:
            move_vertex(g, r);
            break;


        case 8:
            {
            ps_attempt = ks0 * g.Nsurf;
            if (gsl_rng_uniform(r) < ps_attempt)
                {
                    ind = gsl_rng_uniform_int(r, g.boundary.size());
                    e = g.boundary[ind];
                    if (g.check_inside_overlap(e) > 0)
                    {
                        ssadd = attempt_add_monomer_dimer(g, e, r);
                        if (ssadd > 1)
                            dimeradded++;
                        else if (ssadd > 0)
                            monomeradded++;
                        ssadd = -1;
                        g.update_boundary();
                    }
                }
            }
            break;

        case 9:
        {
            for (int nc = 0; nc < g.Nhe; nc++)
            {
                int ind1 = gsl_rng_uniform_int(r, g.Nhe);
                int e1 = g.he[ind1].id;
                int x = -1;
                x = attempt_change_edge_type(g, e1, r);
                if (x >= 0)
                    typechanged++;
            }
        }
        break;


        }
        }*/
        /*case 8:
            move_vertex(g, r);
            break;

            //g.check_odd_neigh();

            //change edge type
        case 9:
        {
            for (int nc = 0; nc < g.Nhe; nc++)
            {
                int ind1 = gsl_rng_uniform_int(r, g.Nhe);
                int e1 = g.he[ind1].id;
                int x = -1;
                x = attempt_change_edge_type(g, e1, r);
                if (x >= 0)
                    typechanged++;
            }
        }
        break;*/

            /*case 8:
            if (g.Nsurf > 0){
                int p_attempt = kb0 * g.Nsurf;
                if (gsl_rng_uniform(r) < p_attempt)
                {
                    ind = gsl_rng_uniform_int(r, g.boundary.size());
                    e = g.boundary[ind];
                    int ssbind = attempt_bind_triangle(g, e, r);
                    if (ssbind>0){
                        cout <<"bound triangle"<<endl;
                    }
                    g.update_boundary();
                }
            }
            break;

        case 9:
            
            int p_attempt = kb0 * g.Nhe;
            if ( gsl_rng_uniform(r) < p_attempt)
            {
                ind = gsl_rng_uniform_int(r, g.Nhe);
                e = g.he[ind].id;
                int ssunbind = attempt_unbind_triangle(g, e, r);
                if (ssunbind>0){
                    cout << "unbound triangle"<<endl;
                }
                g.update_boundary();
            }
            break;*/

            /****************** DRUG ****************************/
            /* uncomment this part for drug 
                
                case 8:
                    double d_attempt = kd0 * g.Nhe/2;

                    if (gsl_rng_uniform(r) < d_attempt)
                    {

                        ind = gsl_rng_uniform_int(r, g.Nhe);
                        //fprintf(stderr, "Attempt add drug, Nd=%d index=%d\n", g.Nd, ind );
                        e = g.he[ind].id;
                        //if ((g.he[ind].type == 0 || g.he[ind].type == 3))
                        //{ //no_bond_boundary(e)>0) {
                            //fprintf(stderr, "Attempt delete monomer\n" );
                            ss = attempt_add_drug(g, e, r);
                            if (ss > 0)
                            { //cout << "drug added " << endl;
                                drugadded++;
                            }

                            //g.update_index();

                            ss = 0;
                        //}
                    }
                    //g.update_boundary();
                    break;

                    //remove drug

                case 9:    
                    if (gsl_rng_uniform(r) < d_attempt && g.Nd > 0)
                    {

                        ind = gsl_rng_uniform_int(r, g.Nhe);
                        //fprintf(stderr, "Attempt delete drug, Nd=%d index=%d\n", g.Nd, ind );
                        e = g.he[ind].id;
                        //if (g.he[ind].type == 0 || g.he[ind].type == 3)
                        //{ //no_bond_boundary(e)>0) {
                            //fprintf(stderr, "Attempt delete monomer\n" );
                            ss = attempt_remove_drug(g, e, r);
                            if (ss > 0)
                            { //cout << "drug removed " << endl;
                                drugremoved++;
                            }
                            
                            ss = 0;
                        //}
                    }
                    
                    g.update_boundary();
                    */
        //}
        //}
        /****************** OUTPUT ****************************/

        double ee = 0;
        if (sweep % (freq_log) == 0 ) 
        {
            
            if (g.Nhe==6) recenter(g);
            g.update_boundary();
            g.check_odd_neigh();
            ee = g.compute_energy();

            ////dump_lammps_traj(g, int(sweep));
            //dump_lammps_traj_dimers(g, sweep);

            time(&timer2);
            seconds = difftime(timer2, timer1);

            dump_analysis(g, ofile, sweep, seed, seconds);
            dump_lammps_data_file(g, 22222222);
            //dump_lammps_traj_restart(g, sweep);
            dump_lammps_data_dimers(g, 11111111);

            dump_restart_lammps_data_file(g, sweep);
        }

        if (sweep % freq_out == 0 )
        {
            time(&timer2);
            seconds = difftime(timer2, timer1);
            cout << "###################################################################" << endl;
            cout << " ################ RUN TIME " << seconds << " SECONDS ###############" << endl;
            cout << " ############# SWEEP " << sweep << "##############" << endl;
            cout << "######### FRAME " << frame << " ##############" << endl;
            cout << "######### ENERGY " << ee << " ##############" << endl;
            cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
            cout << "#########  NHE " << g.Nhe << " #######################" << endl;
            cout << "#########  NHESURF " << g.boundary.size() << " ##############" << endl;
            cout << "#########  NVSURF " << g.boundaryv.size() << " ##############" << endl;
            cout << "#########  NV_BONDSURF " << g.boundaryvbond.size() << " ##############" << endl;
            cout << "#########  NV5 " << g.Nv5 << " ##############" << endl;
            cout << "#########  MONOMER ADDED " << monomeradded << " ##############" << endl;
            cout << "#########  MONOMER REMOVED " << monomerremoved << " ##############" << endl;
            cout << "#########  DIMER ADDED " << dimeradded << " ##############" << endl;
            cout << "#########  DIMER REMOVED " << dimerremoved << " ##############" << endl;
            cout << "#########  no rate REMOVED " << deletednorate << " ##############" << endl;
            cout << "#########  Surface bound " << binding << " ##############" << endl;
            cout << "#########  Surface Unbound " << unbinding << " ##############" << endl;
            cout << "#########  DrugAdded " << drugadded << " ##############" << endl;
            cout << "#########  DrugRemoved " << drugremoved << " ##############" << endl;
            cout << "#########  ND " << g.Nd << " ##############" << endl;
            cout << "#########  TYPE CHANGED " << typechanged << " ##############" << endl;
            cout << "#########  WEDGE FUSION " << wedgefusion << " ##############" << endl;
            cout << "#########  WEDGE FISSION " << wedgefission << " ##############" << endl;
            cout << "#########  FUSION " << fusion << " ##############" << endl;
            cout << "#########  FISSION " << fission << " ##############" << endl;
            cout << "#########  FUSION HALFEDGES " << g.fusionhe.size() << " ##############" << endl;
            cout << "#########  WEDGE FUSION HALFEDGES " << g.fusionwedgehe.size() << " ##############" << endl;
            cout << "#########  ALL NEIGH " << g.all_neigh << " ##############" << endl;
            cout << "#########  Nboundary " << g.Nboundary << " ##############" << endl;
            cout << "#########  Bound Triangle " << boundtri << " ##############" << endl;
            cout << "#########  Acceptance vmove " << (1.0 * g.accepted_vmove) / (1.0 * (g.accepted_vmove + g.rejected_vmove)) << "#################" << endl;
            cout << "#########  Nvlast " << g.Nvlast << " Nhelast " << g.Nhelast << " ################" << endl;
            cout << "#########  T4 " << g.NCD_T4_in << " T3 " << g.NCD_T3_in << " ################" << endl;
            cout << "#########  avgAddInterval "<<avgAddInterval<<" ###############"<<endl;
        }
        /*int boundaryvsize = g.boundaryv.size();
            if ((boundaryvsize > 10) && (boundaryvsize - surfclosev(g) <= 3))
            {

                fprintf(stderr, "Almost closed\n");
                g.update_boundary();
                dump_lammps_data_file(g, 99999999);
                exit(-1);
            }*/

        int cc = check_bind_triangle(g);
        if (cc > 0)
        {
            cout << "bound triangle" << endl;
            g.update_boundary();
            boundtri += cc;
        }



        /*if (sweep % (freq_vis)==0) {
            
            //dump_lammps_traj_dimers(g, int(sweep)); 
            
        }*/
        g.check_odd_neigh();
        //g.update_neigh();

        
        if (g.Nhe>70 && g.Nhe < minHE_update_neigh && sweep % 10000 == 0)
        {
            double thispace=float(g.Nhe-lastNhe)/10000.0; //pace of adding Nhe per sweep
            cout << "thispace "<< thispace <<endl; 
            
            avgpace=(npace*avgpace+thispace)/(npace+1.0); //average pace of adding 

            avgAddInterval=int(pow(10,(-1 * int(floor(log10(avgpace))) ) ) );
            cout << "avgpace" << avgpace << " avgAddInterval " <<avgAddInterval << endl;
            
            npace++;
            lastNhe=g.Nhe;
        }

        if (g.Nhe > minHE_update_neigh && sweep % 10000 == 0)
        {
            //g.check_odd_neigh();
            g.update_neigh();
            //g.find_overlap_all();
            if (g.find_overlap_all() < 0)
            {
                cout << "error overlap" << endl;
                dump_lammps_data_dimers(g, 5555555);
                exit(-1);
            }
        }
        //see if capsid is growing or it is stalled in mixed morphology
        if (g.Nhe > minHE_update_neigh && sweep % (10*avgAddInterval) == 0)
        {    
            update_geometry_parameters(g);
            if ( g.Nhe - lastNhe<=2 ){
                if ((g.Nhe >= 220 && g.NCD_T4_in >=26 && g.NCD_T3_in >= 3  && g.Nsurf > 10 ) || 
                    (g.Nhe >= 160 && g.NCD_T4_in >= 3 && g.NCD_T3_in >=16  && g.Nsurf > 10) || 
                    (g.Nhe >= 200 && g.NCD_T4_in >= 5 && g.NCD_T3_in >=5  && g.Nsurf > 10) )
                    {
                        cout << "STOP for now - mixed morph" << endl;
                        g.update_boundary();
                        //dump_lammps_traj_dimers(g, int(sweep));
                        dump_lammps_data_dimers(g, 44444444);
                        dump_lammps_data_dimers(g, 11111111);
                        time(&timer2);
                        seconds = difftime(timer2, timer1);
                        dump_analysis(g, ofile, sweep, seed, seconds);
                        exit(-1);
                    }
            }
            
            
            if (sweep % (100*avgAddInterval) == 0)
            {
                
                if (g.NCD_T4_in>0 && g.NCD_T3_in>0 && abs( g.Nhe - lastNheGrowth)<=4 ){
                    fprintf(stderr, "STOP for now - not growing\n");
                    g.update_boundary();
                    //dump_lammps_traj_dimers(g, int(sweep));
                    dump_lammps_data_dimers(g, 333333333);
                    dump_lammps_data_dimers(g, 11111111);
                    dump_restart_lammps_data_file(g, sweep);
                    time(&timer2);
                    seconds = difftime(timer2, timer1);
                    dump_analysis(g, ofile, sweep, seed, seconds);
                    exit(-1);
                }
                lastNheGrowth = g.Nhe;

            }

            lastNhe = g.Nhe;
            
        }

        if (sweep == 200000000)
        {

            fprintf(stderr, "STOP for now - too long\n");
            g.update_boundary();
            //dump_lammps_traj_dimers(g, int(sweep));
            //dump_lammps_traj_restart(g, int(sweep));
            dump_lammps_data_dimers(g, 77777777);
            dump_restart_lammps_data_file(g, sweep);
            time(&timer2);
            seconds = difftime(timer2, timer1);
            dump_analysis(g, ofile, sweep, seed, seconds);
            exit(-1);
        }

        if (g.Nhe >= 310 || g.Nv >= 65)
        {

            fprintf(stderr, "STOP for now - too large\n");
            g.update_boundary();
            //dump_lammps_traj_dimers(g, int(sweep));
            dump_lammps_data_dimers(g, 88888888);
            dump_lammps_data_dimers(g, 11111111);
            dump_restart_lammps_data_file(g, sweep);
            time(&timer2);
            seconds = difftime(timer2, timer1);
            dump_analysis(g, ofile, sweep, seed, seconds);
            exit(-1);
        }
        sweep++;
    }
    //dump_lammps_traj_restart(g, int(sweep));
    //dump_lammps_traj_dimers(g, int(sweep));

    //if (g.Nboundary == 1)
    dump_restart_lammps_data_file(g, sweep);
    //equilibrating final structure
    for (int rstep = 0; rstep < (10 * freq_log); rstep++)
    {
        move_vertex(g, r);
        g.update_boundary();
        int ind1 = gsl_rng_uniform_int(r, g.Nhe);
        int e1 = g.he[ind1].id;
        int x = -1;

        //change edge type
        x = attempt_change_edge_type(g, e1, r);
        if (x >= 0)
            typechanged++;

        if (sweep % freq_vis == 0)
        {

            //dump_lammps_traj_dimers(g, int(sweep));
        }
        if (sweep % freq_out == 0)
        {
            time(&timer2);
            seconds = difftime(timer2, timer1);
            cout << "###################################################################" << endl;
            cout << " ################  RUN TIME " << seconds << " SECONDS ###############" << endl;
            cout << " ############# SWEEP " << sweep << "##############" << endl;

            double ee = g.compute_energy();
            cout << "######### ENERGY " << ee << " ##############" << endl;
            cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
        }
        if (sweep % freq_log == 0)
        {

            time(&timer2);
            seconds = difftime(timer2, timer1);

            dump_analysis(g, ofile, sweep, seed, seconds);
        }
        sweep++;
    }

    //dump_lammps_traj_dimers(g, frame++);
    dump_lammps_data_file(g, 22222222);
    dump_lammps_data_dimers(g, 11111111);
    dump_restart_lammps_data_file(g, sweep);

    time(&timer2);
    seconds = difftime(timer2, timer1);
    cout << " ################  FULL RUN TIME " << seconds << " SECONDS ###############" << endl;
    cout << " ############# SWEEP " << sweep << "##############" << endl;

    ee = g.compute_energy();
    cout << "######### ENERGY " << ee << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << g.boundary.size() << " ##############" << endl;
    cout << "#########  NVSURF " << g.boundaryv.size() << " ##############" << endl;
    cout << "#########  NV_BONDSURF " << g.boundaryvbond.size() << " ##############" << endl;
    cout << "#########  NV5 " << g.Nv5 << " ##############" << endl;
    cout << "#########  MONOMER ADDED " << monomeradded << " ##############" << endl;
    cout << "#########  MONOMER REMOVED " << monomerremoved << " ##############" << endl;
    cout << "#########  DIMER ADDED " << dimeradded << " ##############" << endl;
    cout << "#########  DIMER REMOVED " << dimerremoved << " ##############" << endl;
    cout << "#########  no rate REMOVED " << deletednorate << " ##############" << endl;
    cout << "#########  Surface bound " << binding << " ##############" << endl;
    cout << "#########  Surface Unbound " << unbinding << " ##############" << endl;
    cout << "#########  DrugAdded " << drugadded << " ##############" << endl;
    cout << "#########  DrugRemoved " << drugremoved << " ##############" << endl;
    cout << "#########  ND " << g.Nd << " ##############" << endl;
    cout << "#########  TYPE CHANGED " << typechanged << " ##############" << endl;
    cout << "#########  WEDGE FUSION " << wedgefusion << " ##############" << endl;
    cout << "#########  WEDGE FISSION " << wedgefission << " ##############" << endl;
    cout << "#########  FUSION " << fusion << " ##############" << endl;
    cout << "#########  FISSION " << fission << " ##############" << endl;
    cout << "#########  FUSION HALFEDGES " << g.fusionhe.size() << " ##############" << endl;
    cout << "#########  WEDGE FUSION HALFEDGES " << g.fusionwedgehe.size() << " ##############" << endl;
    cout << "#########  ALL NEIGH " << g.all_neigh << " ##############" << endl;
    cout << "#########  Nboundary " << g.Nboundary << " ##############" << endl;
    cout << "#########  Bound Triangle " << boundtri << " ##############" << endl;

    time(&timer2);
    seconds = difftime(timer2, timer1);

    dump_analysis(g, ofile, sweep, seed, seconds);
    finalfile = fopen("last.dat", "w");
    dump_analysis(g, finalfile, sweep, seed, seconds);

    fclose(ofile);
    fclose(finalfile);
    //fclose(outfile);
    return 0;
}
