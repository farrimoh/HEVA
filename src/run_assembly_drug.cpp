#include "geometry.h"
#include "montecarlo.h"
#include <iostream>
#include <vector>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265
using namespace std;

gsl_rng *r;
int main(int argc, char **argv)
{
    time_t timer1, timer2;
    double seconds;
    time(&timer1);
    // parse inputs
    /*if (argc!=12) {
	fprintf(stderr,"usage: ./assemble seed epsilon0 epsilon1 kappa0 kappa1 theta0 theta1 LnK LnZ ks0 hexpent\n");
	exit(-1);
	}*/
    if (argc != 15)
    {
        fprintf(stderr, "usage: ./assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCd ks0 muAB dg mudrug gdrug kd0\n");
        exit(-1);
    }
    
    
    int minHE_update_neigh = 20;


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
    g.all_neigh=0;
    g.epsilon[0] = atof(argv[2]);
    g.epsilon[1] = g.epsilon[0];
    g.epsilon[2] = g.epsilon[0];
    g.epsilon[3] = g.epsilon[0];
    g.kappa[0] = atof(argv[3]);
    g.kappa[1] = g.kappa[0];
    g.kappa[2] = g.kappa[0];
    g.kappa[3] = g.kappa[0];
    double Phikappa = atof(argv[4]);
    g.kappaPhi[0] = Phikappa;     // CD-CD DC-DC and all other
    g.kappaPhi[1] = Phikappa;     //BA-AB
    g.kappaPhi[2] = Phikappa;     //AB-CD and AB-DC
    g.kappaPhi[3] = Phikappa ; //with drug
    g.theta0[0] = atof(argv[5]);
    g.theta0[1] = atof(argv[6]);
    g.theta0[2] = g.theta0[0]; //0.1;  // CD-CD
    g.theta0[3] = g.theta0[0]; //0.1 ; //349 ;

    double gb0 = atof(argv[7]);

    g.dg = atof(argv[11]);
    
    for (int i = 0; i < g.Ntype; i++)
    {
        for (int j = 0; j < g.Ntype; j++)
        {

            if (i == 1 && j == 2)
                g.gb[i][j] = gb0 + (1 * g.dg * gb0);
            else if (i == 0 && j == 1)
                g.gb[i][j] = gb0 + (1 * g.dg * gb0);
            else if (i == 3 && j == 1)
                g.gb[i][j] = gb0 + (1 * g.dg * gb0);
            else if (i == 2 && j == 0)
                g.gb[i][j] = gb0 + (-1 * g.dg * gb0);
            else if (i == 2 && j == 3)
                g.gb[i][j] = gb0 + (-1 * g.dg * gb0);
            else if (i == 3 && j == 3)
                g.gb[i][j] = gb0 + (0 * g.dg * gb0);
            else if (i == 0 && j == 0)
                g.gb[i][j] = .3*gb0;//+ (-3 * g.dg * gb0);
            else if (i == 0 && j == 3)
                g.gb[i][j] = .3* gb0;// + (-3 * g.dg * gb0);
            else if (i == 3 && j == 0)
                g.gb[i][j] = .3 *gb0;// + (-3 * g.dg * gb0);
            else
                g.gb[i][j] = .3 *gb0;
        }
    }

    g.mudrug= atof(argv[12]); 
    double gdrug0 = atof(argv[13]); 


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

    g.mu[0] = atof(argv[8]);
    g.mu[3] = g.mu[0];
    g.mu[1] = atof(argv[10]);
    g.mu[2] = g.mu[1];
    
    double ks0 = atof(argv[9]);
    double kd0 = atof(argv[14]);
    

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
    FILE *efile, *finalfile, *fi; // , *eefile; ,
    g.dump_parameters();
    efile = fopen("energy.dat", "w");
    //outfile = fopen("log.dat","w");
    //finalfile = fopen("last.dat","w");
    fprintf(efile, "#sweep, seed, g.epsilon[0], g.kappa[0], g.theta0[0], gb0, g.mu[0], g.mu[1], g.dg, ks0,g.theta0[2],kd0,g.mudrug , g.compute_energy(), g.Nsurf, g.Nd, g.Nv, g.Nhe / 2\n");

    fi = fopen("parameters_run.out", "w");
    fprintf(fi, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK muCD ks0 muAB mudrug gdrug kd0\n");
    fprintf(fi, "./source/assemble %lu %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.6f\n", seed, g.epsilon[0], g.kappa[0], Phikappa, g.theta0[0], g.theta0[1], gb0, g.mu[0], ks0, g.mu[1], g.dg, g.mudrug , gdrug0, kd0);
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
    fflush(fi);

    double ee = 0;
    fprintf(stderr, "./source/assemble seed epsilon0 kappa0 kappaPhi0 theta0 theta1 LnK LnZ ks0 muAB mudrugdrugProb kd0\n");
    fprintf(stderr, "./source/assemble %lu %f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f\n", seed, g.epsilon[0], g.kappa[0], Phikappa, g.theta0[0], g.theta0[1], gb0, g.mu[0], ks0, g.mu[1], g.dg, g.mudrug , gdrug0, kd0);

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
    int fusion=0;
    int fission=0;




    //make_initial_pentamer(g);
    make_initial_triangle(g);
    g.update_neigh();
    //make_seed(g,r);
    ee = g.compute_energy();
    cout << "######### ENERGY after equilibration " << ee << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
    cout << " ############# FRAME " << frame << "##############" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << g.surfheid.size() << " ##############" << endl;
    //dump_lammps_data_file(g, frame++);
    for (int rstep = 0; rstep < 1000; rstep++)
    { //relaxing the shell
        move_vertex(g, r);
        g.update_surface();
        
    }
    
    ee = g.compute_energy();
    cout << "######### ENERGY after equilibration " << ee << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
    cout << " ############# FRAME " << frame << "##############" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << g.surfheid.size() << " ##############" << endl;
    cout << "#########  NVSURF " << g.surfv.size() << " ##############" << endl;
    cout << "#########  NV_BONDSURF " << g.surfvbond.size() << " ##############" << endl;
    dump_lammps_data_file(g, frame++);
    fprintf(stderr, "Graph initialized.\n");
    

    int binding = 0;
    int unbinding = 0;

    int deletednorate = 0;
    int sweep = 0;
    int ssadd = 0;
    //int ssremove = 0;
    
    int runhpc=1;

    int freq_vis=10000;
    int freq_log=1000;
    int freq_out=10000; // shoud be >=freq_log

    if (runhpc==1){ 
        freq_vis=freq_vis*1000;
    }
          

    while (g.Nsurf > 0)
    //for (int sw=0;sw<1001; sw++)
    {
        
        if (g.Nsurf == 3 && sweep > 1000 && g.Nhe > 10) //LAST TRINGLE
        { //complete capsid
            int hindex0 = g.heidtoindex[g.surfheid[0]];
            int hindex1 = g.heidtoindex[g.surfheid[1]];
            int hindex2 = g.heidtoindex[g.surfheid[2]];
            if (g.he[hindex0].vout == g.he[hindex1].vin)
            {
                g.set_prev_next(g.surfheid[0], g.surfheid[2], g.surfheid[1]);
                g.set_prev_next(g.surfheid[1], g.surfheid[0], g.surfheid[2]);
                g.set_prev_next(g.surfheid[2], g.surfheid[1], g.surfheid[0]);
            }
            else if (g.he[hindex0].vout == g.he[hindex2].vin)
            {
                g.set_prev_next(g.surfheid[0], g.surfheid[1], g.surfheid[2]);
                g.set_prev_next(g.surfheid[1], g.surfheid[2], g.surfheid[0]);
                g.set_prev_next(g.surfheid[2], g.surfheid[0], g.surfheid[1]);
            }

            g.update_surface();
            dump_lammps_data_file(g, frame++);
            time(&timer2);
            seconds = difftime(timer2, timer1);         
            fprintf(efile, "%d %lu %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d\n", sweep, seed,seconds, g.epsilon[0], g.kappa[0], g.theta0[0], gb0, g.mu[0], g.mu[1], g.dg, ks0,g.theta0[2],kd0,g.mudrug , g.compute_energy(), g.Nd,g.Nsurf, g.Nv, g.Nhe / 2);
            fflush(efile);
            break;
        }

        //time(&t1);
        //cout<< "move_vertex"<<endl;
        //if (g.Nhe > 1)
        //{
            //for (int rstep = 0; rstep < g.Nhe ; rstep++)
            //{
        move_vertex(g, r);

                //g.update_surface();
            //}
            /*for ( unsigned int rstep = 0; rstep < g.surfv.size() *5; rstep++)
            {
                ind = gsl_rng_uniform_int(r, g.surfv.size());
                e = g.surfv[ind];
                move_one_vertex(g, e, r);
                //g.update_surface();
            }*/
        //}
        g.check_odd_neigh();
        //g.update_surface();
        //g.update_neigh()
        // CHANGE TYPE ATTEMPT

        //double pc_attempt = 1;//.01 * (g.Nhe / 2); //ks0*g.Nhe;

       // if (gsl_rng_uniform(r) < pc_attempt)
        //{
        //for (int step = 0; step <  g.Nhe; step++)
        //{
            int ind1 = gsl_rng_uniform_int(r, g.Nhe);
            int e1 = g.he[ind1].id;
            int x = -1;
            //int et=g.he[ind1].type;
            //cout<< "try changetype ind1 " << ind1 << " e "<< e1 <<endl;
            x = attempt_change_edge_type(g, e1, r);
            if (x >= 0)
            {
                typechanged++;

                //fprintf(stderr,"  %d edge type %d CHANGED!!!!!!!!!!!!!!!!! to %d \n",ind,et,x);
            }
            //ssadd=0;
            //ssremove=0;
            //g.update_surface();
        //}
        //}
        //g.update_surface();
        //for (int step = 0; step < 2*g.Nsurf; step++)
        /*for (int step = 0; step < 2; step++)
        {
            int ind2 = gsl_rng_uniform_int(r, g.Nsurf);
            int e2 = g.surfheid[ind2];
            int x2 = -1;
            //int et=g.he[ind1].type;
            //cout<< "try changetype ind1 " << ind1 << " e "<< e1 <<endl;
            x2 = attempt_change_edge_type(g, e2, r);
            if (x2 >= 0)
            {
                typechanged++;

                //fprintf(stderr,"  %d edge type %d CHANGED!!!!!!!!!!!!!!!!! to %d \n",ind,et,x);
            }
           
        }*/
        //g.update_surface();
        //}

        //BIND -UNBIND ATTEMPT

        //if (ssadd > 0 || ssremove > 0 || (g.Nhe > 50 && sweep % 100 == 0))
        //double kb0=1;
        //double pb_attempt=kb0*g.surfv.size();
        if (g.Nhe>15)// && gsl_rng_uniform(r) < pb_attempt)
        {

            //cout <<" trying bind "<<endl;//g.surfv.size()
            ind = gsl_rng_uniform_int(r, g.surfv.size());
            //cout <<"ind id" << ind<<endl;
            int vv = g.surfv[ind];
            int tt = -1;
            if(g.is_bond_vsurface(vv)<0){
                //cout <<" trying bind vid is " << vv <<endl;
                tt = attempt_bind_wedge_dimer(g, vv, r);
               
                if (tt > 0)
                {
                    binding++;
                    tt=-1;
                } //cout << "Bound " << e <<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!frame " << frame<<endl;  g.update_surface(); dump_lammps_data_file(g, frame++); }// exit(-1);}
                //cout << "try binding"<<endl;
                
            }
            g.update_surface();
            //dump_lammps_data_file(g,frame++);
        }
        if (g.surfvbond.size()>0)// &&  gsl_rng_uniform(r) < pb_attempt)
        {

            ind = gsl_rng_uniform_int(r, g.surfv.size());
            //cout <<"ind id" << ind<<endl;
            int vv = g.surfv[ind];

            if (g.is_bond_vsurface(vv)>0) { 
            
                //cout <<" trying unbind vv is " << vv <<endl;
                int tt = attempt_unbind_wedge_dimer(g, vv, r);

                if (tt > 0)
                {
                    unbinding++;
                } //cout << "UNNNNBound "<<endl;

                            
                g.update_surface();
            }
        }
            

        int ss = 0;

        // REMOVE DIMER ATTEMPT

        double ps_attempt = ks0 * g.Nsurf;

        if (gsl_rng_uniform(r) < ps_attempt)
        {
            if (sweep > 0 && g.Nhe > 6)
            {
                ind = gsl_rng_uniform_int(r, g.surfheid.size());
                //fprintf(stderr, "Attempt delete monomer dimer, Nsurf=%d index=%d\n", g.Nsurf, ind );
                e = g.surfheid[ind];
                if (g.no_bond_surface(e) > 0)
                {
                    //fprintf(stderr, "Attempt delete monomer\n" );
                    ss = attempt_remove_monomer_dimer(g, e, r);
                    if (ss > 1)
                    {
                        dimerremoved++;
                        //if (g.Nhe>minHE_update_neigh) { g.update_neigh();}
                    }
                    else if (ss > 0)
                    {
                        monomerremoved++;
                        //if (g.Nhe>minHE_update_neigh) { g.update_neigh();}
                    }
                    
                    ss = 0;
                    g.update_surface();
                    //g.update_neigh();
                }
                

            }
            
        }
        
       // ADD DIMER ATTEMPT
        
        ps_attempt = ks0 * g.Nsurf;
        if (gsl_rng_uniform(r) < ps_attempt)
        {

            ind = gsl_rng_uniform_int(r, g.surfheid.size());
            e = g.surfheid[ind];
            
            //if g.check_overlap_g(g.heidtoindex())
                //fprintf(stderr, "Now Run insert monomer_dimer for ind %d edge %d\n" ,ind ,e);
            ssadd = attempt_add_monomer_dimer(g, e, r);

            if (ssadd > 1)
            {
                dimeradded++;
                //if (g.Nhe>minHE_update_neigh) { g.update_neigh();}
            }
            else if (ssadd > 0)
            {
                monomeradded++;
            }
                
            g.update_surface();
            //g.update_neigh();    
        }
       
        

        //ATTEMPT FUSION_FISSION
        
        if (g.all_neigh>0 && g.Nhe>20 && g.Nsurf>3) {
            //double kf0=1;
            //double pf_attempt=kf0;//*g.Nsurf;
            //if (gsl_rng_uniform(r) < pf_attempt)
            //{    
                int ff =attempt_vertex_fusion(g, r);
                g.update_surface(); 
                
                if (ff > 0) { 
                    fusion++;                    
                } 
            //}
        }
        if (g.all_neigh>0 && g.Nhe>20 && g.Nsurf>3) {
            //pf_attempt=kf0*.25; //*g.Nsurf
            //if (gsl_rng_uniform(r) < pf_attempt)
            //{    

                int ff = attempt_vertex_fission(g, r);
                
                g.update_surface();   
                
                
                if (ff > 0){    
                    fission++; 
                    
                }
            //}
            
        }
        
        
        //ADD- REMOVE DRUG

        double d_attempt = kd0 * g.Nhe/2;

        if (gsl_rng_uniform(r) < d_attempt)
        {

            ind = gsl_rng_uniform_int(r, g.Nhe);
            //fprintf(stderr, "Attempt add drug, Nd=%d index=%d\n", g.Nd, ind );
            e = g.he[ind].id;
            //if ((g.he[ind].type == 0 || g.he[ind].type == 3))
            //{ //no_bond_surface(e)>0) {
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
        //g.update_surface();
        
        //REMOVE DRUG

        if (gsl_rng_uniform(r) < d_attempt && g.Nd > 0)
        {

            ind = gsl_rng_uniform_int(r, g.Nhe);
            //fprintf(stderr, "Attempt delete drug, Nd=%d index=%d\n", g.Nd, ind );
            e = g.he[ind].id;
            //if (g.he[ind].type == 0 || g.he[ind].type == 3)
            //{ //no_bond_surface(e)>0) {
                //fprintf(stderr, "Attempt delete monomer\n" );
                ss = attempt_remove_drug(g, e, r);
                if (ss > 0)
                { //cout << "drug removed " << endl;
                    drugremoved++;
                }
                
                ss = 0;
            //}
        }
        
        //g.update_surface();

        
        /****************** OUNTPUT ****************************/
        

            
        double ee=0;
        if (sweep % (freq_log)==0) 
        {
            g.update_surface();
            g.check_odd_neigh();
            ee = g.compute_energy();
            time(&timer2);
            seconds = difftime(timer2, timer1);         
            fprintf(efile, "%d %lu %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d\n", sweep, seed,seconds, g.epsilon[0], g.kappa[0], g.theta0[0], gb0, g.mu[0], g.mu[1], g.dg, ks0,g.theta0[2],kd0,g.mudrug , g.compute_energy(), g.Nd,g.Nsurf, g.Nv, g.Nhe / 2);
            fflush(efile);
            dump_lammps_data_file(g, 22222222);
        }

        if (sweep % freq_out == 0){
                        
            cout << "###################################################################" << endl;
            cout << " ################ RUN TIME " << seconds << " SECONDS ###############" << endl;
            cout << " ############# SWEEP " << sweep << "##############" << endl;
            cout << "######### FRAME " << frame << " ##############" << endl;
            cout << "######### ENERGY " << ee << " ##############" << endl;
            cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
            cout << "#########  NHE " << g.Nhe << " #######################" << endl;
            cout << "#########  NHESURF " << g.surfheid.size() << " ##############" << endl;
            cout << "#########  NVSURF " << g.surfv.size() << " ##############" << endl;
            cout << "#########  NV_BONDSURF " << g.surfvbond.size() << " ##############" << endl;
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
            cout << "#########  FUSION " << fusion << " ##############" << endl;
            cout << "#########  FISSION " << fission << " ##############" << endl;
           
                     
            int surfvsize=g.surfv.size();
            if ((surfvsize>10 ) &&  (surfvsize-surfclosev(g)<=3))
            {

                fprintf(stderr, "Almost closed\n");
                g.update_surface();
                dump_lammps_data_file(g, 999999);
                exit(-1);
            }            

        }

        if (sweep % (freq_vis)==0) {
            recenter(g);
            dump_lammps_data_file(g, frame++); 
        }

        if (g.Nhe>minHE_update_neigh && sweep % 10000 == 0)
        {
            g.check_odd_neigh();
        }


        if (sweep == 100000000)
        {

            fprintf(stderr, "STOP for now - too long\n");
            g.update_surface();
            dump_lammps_data_file(g, frame++);
            exit(-1);
        }

        sweep++;
        
    }

    //equlibrating final structure
    for (int rstep = 0; rstep < (10*freq_out); rstep++)
    {
        move_vertex(g, r);
        g.update_surface();
        if (sweep % freq_vis == 0)
        {

            dump_lammps_data_file(g, frame++);
        }
        if (sweep % freq_out == 0){
            time(&timer2);
            seconds = difftime(timer2, timer1);
            cout << "###################################################################" << endl;
            cout << " ################  RUN TIME " << seconds << " SECONDS ###############" << endl;
            cout << " ############# SWEEP " << sweep << "##############" << endl;

            double ee = g.compute_energy();
            cout << "######### ENERGY " << ee << " ##############" << endl;
            cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
        }
        if (sweep % freq_log == 0){

            time(&timer2);
            seconds = difftime(timer2, timer1);         
            fprintf(efile, "%d %lu %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d\n", sweep, seed,seconds, g.epsilon[0], g.kappa[0], g.theta0[0], gb0, g.mu[0], g.mu[1], g.dg, ks0,g.theta0[2],kd0,g.mudrug , g.compute_energy(), g.Nd,g.Nsurf, g.Nv, g.Nhe / 2);
            fflush(efile);
        }
        sweep++;
    }

    dump_lammps_data_file(g, frame++);
    dump_lammps_data_file(g, 22222222);
    time(&timer2);
    seconds = difftime(timer2, timer1);
    cout << " ################  FULL RUN TIME " << seconds << " SECONDS ###############" << endl;
    cout << " ############# SWEEP " << sweep << "##############" << endl;

    ee = g.compute_energy();
    cout << "######### ENERGY " << ee << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << g.surfheid.size() << " ##############" << endl;
    cout << "#########  NVSURF " << g.surfv.size() << " ##############" << endl;
    cout << "#########  NV_BONDSURF " << g.surfvbond.size() << " ##############" << endl;
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
    cout << "#########  FUSION " << fusion << " ##############" << endl;
    cout << "#########  FISSION " << fission << " ##############" << endl;

    time(&timer2);
    seconds = difftime(timer2, timer1);         
    fprintf(efile, "%d %lu %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d\n", sweep, seed,seconds, g.epsilon[0], g.kappa[0], g.theta0[0], gb0, g.mu[0], g.mu[1], g.dg, ks0,g.theta0[2],kd0,g.mudrug , g.compute_energy(), g.Nd,g.Nsurf, g.Nv, g.Nhe / 2);
    finalfile = fopen("last.dat", "w");
    

    fprintf(finalfile, "%d %lu %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %d %d\n", sweep, seed,seconds, g.epsilon[0], g.kappa[0], g.theta0[0], gb0, g.mu[0], g.mu[1], g.dg, ks0,g.theta0[2],kd0,g.mudrug , g.compute_energy(), g.Nd,g.Nsurf, g.Nv, g.Nhe / 2);
    dump_lammps_data_file(g, 11111111);
    //dump_lammps_data_file(g, frame + sweep + 100);

    fclose(efile);
    fclose(finalfile);
    //fclose(outfile);
    return 0;
}
