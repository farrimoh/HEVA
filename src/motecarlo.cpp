 /*
 * MonteCarlo.cpp
 *
 *  Created on: May 3, 2019
 *      Author: farri
 */
#include "geometry.hpp"
#include "montecarlo.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#define PI 3.14159265
void move_vertex(geometry &g, gsl_rng *r)
{
    double *newv = new double[3];
    double *oldv = new double[3];

    int overlapflag = -1;
    
    //cout << "in move_vertex"<<endl;
    for (unsigned int i = 0; i < g.v.size(); i++)
    { 
        overlapflag = -1;
        
        double e1 = 0, e2 = 0, de = 0;
        int ind = gsl_rng_uniform_int(r, g.v.size());
        /* current neigh 
        vector<int> vecupdate;      
        for (vector<int>::iterator it =g.v[ind].vneigh.begin(); it != g.v[ind].vneigh.end(); ++it)
        {	
            vecupdate.push_back(*it);
        }*/
        
        /*TEMP*/
        
        e1+=g.vertex_energy(g.v[ind].vid);
        
        oldv[0] = g.v[ind].co[0];
        oldv[1] = g.v[ind].co[1];
        oldv[2] = g.v[ind].co[2];
        g.move_p(oldv, newv, r);
        
        g.v[ind].co[0] = newv[0];
        g.v[ind].co[1] = newv[1];
        g.v[ind].co[2] = newv[2];
        
        for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ++ithe)
        { //update the geometry
            int heindex = g.heidtoindex[*ithe];
           
            g.update_half_edge(*ithe);
            g.update_half_edge(g.he[heindex].opid);
            
            if (g.he[heindex].previd != -1)
            {
                g.update_half_edge(g.he[heindex].previd);
                g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].previd]].opid);
            }

            if (g.he[heindex].nextid != -1)
            {
                g.update_half_edge(g.he[heindex].nextid);
                g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].nextid]].opid);
            }
            
            
        }
        
        g.update_normals();
        e2+=g.vertex_energy(g.v[ind].vid); 
         

        /* first update old neighbors 
        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
        {	
            g.update_neigh_vertex(*it);
        }
        vecupdate.clear();*/
        /* now update ne */
        /* save vecupdate */
        /*g.update_neigh_vertex(g.v[ind].vid);
        if (g.v[ind].vneigh.size() > 0)
        {
            for (vector<int>::iterator it =g. v[ind].vneigh.begin(); it != g.v[ind].vneigh.end(); ++it)
            {	
                g.update_neigh_vertex(*it);
                vecupdate.push_back(*it);
            }
        }*/
        if ( g.Nhe>20) {
            if (g.check_overlap_g(g.v[ind].vid) < 0) { overlapflag = 1; } 
            
        }

            

        de = e2 - e1;
       
        double crit = exp((-de) / g.T);
       
        
        if (gsl_rng_uniform(r) > crit || overlapflag ==1)
        {
            /* not accepted putting things back */
            g.v[ind].co[0] = oldv[0];
            g.v[ind].co[1] = oldv[1];
            g.v[ind].co[2] = oldv[2];
            for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ithe++)
            { 

                int heindex = g.heidtoindex[*ithe];
                
                g.update_half_edge(*ithe);
                g.update_half_edge(g.he[heindex].opid);
                if (g.he[heindex].previd != -1)
                {
                    g.update_half_edge(g.he[heindex].previd);
                    g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].previd]].opid);
                }

                if (g.he[heindex].nextid != -1)
                {
                    g.update_half_edge(g.he[heindex].nextid);
                    g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].nextid]].opid);
                }


            }
            /* first update vecupdate 
            for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
            {	
                g.update_neigh_vertex(*it);
            }
            vecupdate.clear();*/

            /* now update to old neighborhood 
            g.update_neigh_vertex(g.v[ind].vid); 
            if (g.v[ind].vneigh.size() > 0)
            {
                for (vector<int>::iterator it =g. v[ind].vneigh.begin(); it != g.v[ind].vneigh.end(); ++it)
                {	
                    g.update_neigh_vertex(*it);
                }
            }*/
          

        }
        
        g.update_normals();
        //cout << "move vertex  index" << ind <<endl;
        //g.check_odd_neigh();
        
    }
    //after update neigh // XXX this soulb be updated to just update the surface and surf-neighbors
    //if (g.Nhe>100) g.update_neigh();
    //else {
    for (vector<int>::iterator it = g.surfv.begin(); it != g.surfv.end(); ++it)
        {
            g.update_neigh_vertex(*it);
            for (vector<int>::iterator vit =g. v[g.vidtoindex[*it]].vneigh.begin(); vit != g.v[g.vidtoindex[*it]].vneigh.end(); ++vit)
                {	
                    g.update_neigh_vertex(*vit);
                }
        }
    
    
    //}
    
    /*for (unsigned int i = 0; i < g.v.size(); i++)
    {    
        if (g.check_overlap_g(g.v[i].vid) < 0) { 
             
            dump_lammps_data_file(g, 5555555);
            cout << "after move vertex" <<endl;
            cout <<"OVERLAP vid " << g.v[i].vid << "vindex is " << g.vidtoindex[g.v[i].vid] << endl;
            //g.find_overlap_g(g.v[ind].vid);
            exit(-1); 
        } 
    }  */  
    delete[] newv;
    delete[] oldv;
}

/*int move_one_vertex(geometry &g, int vid0, gsl_rng *r)
{
    int ind = g.vidtoindex[vid0];
    double *newv = new double[3];
    double *oldv = new double[3];
    double e1 = 0, e2 = 0, de = 0;
    int overlapflag = -1;
    for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ++ithe)
    {
        //cout << "here"<<endl;
        int heindex = g.heidtoindex[*ithe];
        //cout << "heindex0 " << heindex<<endl;
        int opindex = g.heidtoindex[g.he[heindex].opid];
        //cout << "opindex0 " << opindex<<endl;

        //cout << "previndex0 " << previndex<<endl;
        e1 += g.stretch_energy(heindex);
        //if (
        e1 += g.bend_energy(heindex);
        e1 += g.dimer_bend_energy(opindex);
        e1 += g.dimer_bend_energy(heindex);

        if (g.he[heindex].previd != -1)
        {

            int previndex = g.heidtoindex[g.he[heindex].previd];
            e1 += g.bend_energy(previndex) + g.dimer_bend_energy(previndex);
        }
    }
    // cout << "e1 is" <<e1;
    //cout << " e1 is " << e1 <<endl;
    oldv[0] = g.v[ind].co[0];
    oldv[1] = g.v[ind].co[1];
    oldv[2] = g.v[ind].co[2];
    g.move_v(oldv, newv, r);

    if (g.check_overlap_centerv(newv) < 0)
    {
        return -1;
    }
    g.v[ind].co[0] = newv[0];
    g.v[ind].co[1] = newv[1];
    g.v[ind].co[2] = newv[2];

    for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ++ithe)
    { //update the geometry
        //if *ithe
        int heindex = g.heidtoindex[*ithe];
        //cout << " updating edge after move" << *ithe << endl;
        g.update_half_edge(*ithe);
        //cout << " updating op edge " << g.he[heindex].opid << endl;
        g.update_half_edge(g.he[heindex].opid);
        if (g.he[heindex].previd != -1)
        {

            g.update_half_edge(g.he[heindex].previd);
            g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].previd]].opid);
        }
    }
    //g.update_surface();
    g.update_normals();
    
    for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ++ithe)
    { //
        //cout << " calculate energy  for vertex "<< ind << " *ithe " << *ithe <<endl;
        int heindex = g.heidtoindex[*ithe];
        //cout << " heidex is g.heidtoindex[*ithe]; "<< heindex << endl;
        int opindex = g.heidtoindex[g.he[heindex].opid];
        //cout << " opid is " << " g.he[heindex].opid; "<< opindex << endl;
        e2 += g.stretch_energy(heindex);
        e2 += g.bend_energy(heindex);
        e2 += g.dimer_bend_energy(heindex) + g.dimer_bend_energy(opindex);
        if (g.he[heindex].previd != -1)
        {

            int previndex = g.heidtoindex[g.he[heindex].previd];
            e2 += g.bend_energy(previndex) + g.dimer_bend_energy(previndex);
        }

        if (g.check_overlap_hesurf(*ithe) < 0)
        {
            overlapflag = 1;
        }
    }

    de = e2 - e1;
    //cout << " de is " << de << endl;
    //  cout << " T is " << g.T << endl;
    double crit = exp((-de) / g.T);
    //cout << "crit move vertex is " << crit <<endl;
    if (crit < 0)
    {
        cout << "ERRRRORRRR!!!!" << endl;
        exit(-1);
    }
    if (gsl_rng_uniform(r) > crit || overlapflag == 1)
    {
        g.v[ind].co[0] = oldv[0];
        g.v[ind].co[1] = oldv[1];
        g.v[ind].co[2] = oldv[2];
        for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ithe++)
        { 

            int heindex = g.heidtoindex[*ithe];
            //cout << " updating edge " << *ithe << endl;
            g.update_half_edge(*ithe);
            //cout << " updating edge " << *ithe << endl;
            g.update_half_edge(g.he[heindex].opid);
            if (g.he[heindex].previd != -1)
            {
                g.update_half_edge(g.he[heindex].previd);
                g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].previd]].opid);
            }
        }
    }
    //g.update_surface();
    g.update_normals();

    delete[] newv;
    delete[] oldv;
    return 0;
}*/

int attempt_add_monomer_dimer(geometry &g, int heid0, gsl_rng *r) //!!! Should update with pre_oipen wedge
{
    //cout << " in attempt_add_monomer-dimer" <<endl;
    if (g.is_surface(heid0) < 0)
    {
        cout << " cannot add not on the surface !" << endl;
        exit(-1);
    }
    //double gbb=0;

    int etypenew1 = -1;
    int etypenew2 = -1;
   

    int heindex0 = g.heidtoindex[heid0];
    int etypeheid0 = g.he[heindex0].type;

    //monomer with next

    //if (g.he[heindex0].nextid!=-1 && g.he[heindex0].previd!=-1) {cout <<" wrong not on surface" <<endl; exit(-1); }
    int xid = -1;
    
    int vid0 = g.is_bond_out_surface(heid0);
    int vin0 = g.is_bond_in_surface(heid0);

    //adding monomer
    g.update_fusion_pairs();

    if (vid0 >= 0 && g.v[g.vidtoindex[g.he[heindex0].vin]].hein.size() < 6 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[heindex0].nextid]].vin]].hein.size() < 6)
    {

        xid = g.he[heindex0].nextid;
        int xidindex = g.heidtoindex[xid];
        if ( (g.v[g.vidtoindex[g.he[heindex0].vin]].fusion_vid!=-1) ||  (g.v[g.vidtoindex[g.he[xidindex].vout]].fusion_vid!=-1)) { return -1;}

        if (g.is_bond_in_surface(xid) != vid0 || g.is_bond_out_surface(xid) != -1)
        {
            cout << "wrong bound on surface" << endl;
            exit(-1);
        }
        else
        {
            //cout << "try add monomer" <<endl;
            int etypenew = -1;

            int etypexid = g.he[g.heidtoindex[xid]].type;
            // SELECTED TYPES

            if ((etypeheid0 == 2) && ((etypexid == 0) || (etypexid == 3)))
            {
                etypenew = 1;
            }
            else if (((etypeheid0 == 3) || (etypeheid0 == 0)) && (etypexid == 1))
            {
                etypenew = 2;
            }
            else if (((etypeheid0 == 3) && (etypexid == 3)))
            {
                etypenew = 3;
            }
            else if (((etypeheid0 == 0) && (etypexid == 0)) || ((etypeheid0 == 1) && (etypexid == 2)))
            {
                etypenew = 0;
            }
            else if ((etypeheid0 == 1) && (etypexid == 2))
            {
                if (gsl_rng_uniform(r) < 0.5)
                {
                    etypenew = 0;
                }
                else
                {
                    etypenew = 3;
                }
            }
            else
            {
                etypenew = gsl_rng_uniform_int(r, 4);
            }
            // TYPES BASED ON Concentration
            /*if (gsl_rng_uniform(r) < g.cdProb) {
                if (gsl_rng_uniform(r) < 0.5) {etypenew=0;}
                else {etypenew=3; }
            } 
            else{
                if (gsl_rng_uniform(r) < 0.5) {etypenew=1;}
                else {etypenew=2; } 

            }  */
            // etypenew=gsl_rng_uniform_int(r,4);

            bool drug0 = 0;
            //bool drug1=0;
            /*if ((etypenew==3) ||(etypenew==0) ) {
            if (gsl_rng_uniform(r) < g.drugProb) {drug0=1; }
            if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;  }
            }*/
            //cout << "g.he[heindex0].din" << g.he[heindex0].din<<endl;
            double gbb = g.find_dg(etypenew, etypeheid0, g.he[heindex0].din);
            gbb += g.find_dg(etypexid, etypenew, drug0);
            //cout << "add_monomer with prev" <<endl;
            

            


            //double e1=g.bend_energy(heindex0) + g.bend_energy(xidindex);
            int x = g.add_monomer(heid0, xid, etypenew);
            int Nhelast2index = g.heidtoindex[g.Nhelast - 2];
            if (x > 0)
            {
                //double de=g.monomer_energy(g.Nhelast-1);
                int voutid=g.he[heindex0].vout;
                int vinid=g.he[xidindex].vin;
                vector<int> vecupdate;
                vecupdate.push_back(vinid);
                vecupdate.push_back(voutid);
                g.update_half_edge(heid0);
                g.update_half_edge(g.he[heindex0].opid);
                g.update_half_edge(xid);
                g.update_half_edge(g.he[xidindex].opid);
                g.update_half_edge(g.Nhelast - 1);
                //g.he[g.heidtoindex(g.Nhelast-1)].din=drug1;
                g.update_half_edge(g.Nhelast - 2);
                //g.he[g.heidtoindex[g.Nhelast-2]].din=drug0;
                g.update_surface();
                double de = g.stretch_energy(Nhelast2index);
                de += g.dimer_bend_energy(Nhelast2index) + g.dimer_bend_energy(xidindex);
                de += g.bend_energy(heindex0) + g.bend_energy(xidindex);//-e1; // g.dimer_bend_energy(heindex0);
                //de+=  ; //g.monomer_energy(heid0);
                //cout << " de is " << de <<endl;
                //double e2=g.compute_energy();
                //de += g.gb * 3 - g.mu;
                //gbb=gb0next+gb0prev;
                de += gbb - g.mu[etypenew];
                /*int vinid=g.he[g.heidtoindex[g.Nhelast - 1]].vin;
                int vindex0=g.vidtoindex[vinid];
                g.update_neigh_vertex(vinid);
                if (g.v[vindex0].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                    {	
                        g.update_neigh_vertex(*it);
                    }
                }
                int voutid=g.he[g.heidtoindex[g.Nhelast - 1]].vout;
                vindex0=g.vidtoindex[voutid];
                if (g.v[vindex0].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                    {	
                        g.update_neigh_vertex(*it);
                    }
                }*/

                //double crit = 2.*g.z*g.K*g.K*g.K*exp((-de)/g.T);
                double crit = exp((-de) / g.T) / 2;
                int overlapflag=-1;
                
                if (gsl_rng_uniform(r) < crit ) 
                {
                    for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                    {	
                        g.update_neigh_vertex(*it);
                    }
                    
                    for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                    {	
                        if (g.check_overlap_g(*it)<0) overlapflag=1;

                    }
                    
                    if (overlapflag==-1) {
                        //cout << "monomer added no overlap" <<endl;
                        vecupdate.clear();
                        return 1;
                    }
                }
                else
                {
                    int nextid0 = g.he[Nhelast2index].nextid;
                    int previd0 = g.he[Nhelast2index].previd;
                    //int xidindex=g.heidtoindex[xid];
                    if (nextid0 == -1 || previd0 == -1)
                    {
                        cout << "not accepted in add monomer!" << endl;
                        exit(-1);
                    }
                    if (g.delete_edge(g.Nhelast - 1) > 0)
                    {
                        
                        g.he[heindex0].previd = -1;
                        g.he[xidindex].nextid = -1;
                        
                        g.update_half_edge(heid0);
                        g.update_half_edge(g.he[heindex0].opid);
                        g.update_half_edge(xid);
                        g.update_half_edge(g.he[xidindex].opid);
                        g.update_index();
                        //cout << " MONOMER REMOVED AFTER addition xid1" <<endl;

                        
                        /*vindex0=g.vidtoindex[vinid];
                        
                        if (g.v[vindex0].vneigh.size() > 0)
                        {
                            for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                            {	
                                g.update_neigh_vertex(*it);
                            }
                        }
                        vindex0=g.vidtoindex[voutid];
                        if (g.v[vindex0].vneigh.size() > 0)
                        {
                            for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                            {	
                                g.update_neigh_vertex(*it);
                            }
                        }*/
                        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                        vecupdate.clear();
                        return -1;
                    }
                    else
                    {
                        cout << " could not delete_monomer HERE 4444";
                        exit(-1);
                    }
                }
            }
            else
            {
                
                //cout << "could not add monomer" <<endl;
                return -1;
            }
        }
    }
    
    //adding monomer
    else if (vin0 >= 0 && g.v[g.vidtoindex[g.he[heindex0].vout]].hein.size() < 6 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[heindex0].previd]].vin]].hein.size() < 6)
    {

        xid = g.he[heindex0].previd;
        int xidindex=g.heidtoindex[xid];
        if ( (g.v[g.vidtoindex[g.he[heindex0].vout]].fusion_vid!=-1) ||  (g.v[g.vidtoindex[g.he[xidindex].vin]].fusion_vid!=-1)) { return -1;}

        if (g.is_bond_out_surface(xid) != vin0 || g.is_bond_in_surface(xid) != -1)
        {
            cout << "wrong bound on surface" << endl;
            exit(-1);
        }
        else
        {
           
            int etypexid = g.he[g.heidtoindex[xid]].type;
            int etypenew = -1;
            //cout  << "add_monomer with next" <<endl;
            //BASED on types
            if ((etypexid == 2) && ((etypeheid0 == 0) || (etypeheid0 == 3)))
            {
                etypenew = 1;
            }
            else if (((etypexid == 3) || (etypexid == 0)) && (etypeheid0 == 1))
            {
                etypenew = 2;
            }
            else if (((etypexid == 3) && (etypeheid0 == 3)))
            {
                etypenew = 3;
            }
            else if (((etypexid == 0) && (etypeheid0 == 0)) || ((etypexid == 1) && (etypeheid0 == 2)))
            {
                etypenew = 0;
            }
            else if ((etypexid == 1) && (etypeheid0 == 2))
            {
                if (gsl_rng_uniform(r) < .5)
                {
                    etypenew = 0;
                }
                else
                {
                    etypenew = 3;
                }
            }
            else
            {
                etypenew = gsl_rng_uniform_int(r, 4);
            }
            /* if (gsl_rng_uniform(r) < g.cdProb) {
                if (gsl_rng_uniform(r) < 0.5) {etypenew=0;}
                else {etypenew=3; }
            } 
            else{
                if (gsl_rng_uniform(r) < 0.5) {etypenew=1;}
                else {etypenew=2; } 
            }*/
            //CHOOSING RANDOM
            //etypenew=gsl_rng_uniform_int(r,4);

            //else { cout << }
            bool drug0 = 0;
            //bool drug1=0;
            /*if ((etypenew==3) ||(etypenew==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug0=1;}
                if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;}
            }*/
            double gbb = g.find_dg(etypenew, etypexid, g.he[g.heidtoindex[xid]].din);
            gbb += g.find_dg(etypeheid0, etypenew, drug0);
            
            //double e1=g.bend_energy(heindex0) + g.bend_energy(xidindex);
            int x = g.add_monomer(xid, heid0, etypenew);
            if (x > 0)
            {
                int vinid=g.he[heindex0].vin;
                int voutid=g.he[xidindex].vout;
                vector<int> vecupdate;
                vecupdate.push_back(vinid);
                vecupdate.push_back(voutid);
                //double de = (g.stretch_energy(g.heidtoindex[g.Nhelast-2]) + g.bend_energy(g.heidtoindex[g.Nhelast-2]));
                g.update_half_edge(heid0);
                g.update_half_edge(g.he[g.heidtoindex[heid0]].opid);
                g.update_half_edge(xid);
                g.update_half_edge(g.he[g.heidtoindex[xid]].opid);
                g.update_half_edge(g.Nhelast - 1);
                //g.he[g.heidtoindex(g.Nhelast-1)].din=drug1;
                g.update_half_edge(g.Nhelast - 2);
                //g.he[g.heidtoindex[g.Nhelast-2]].din=drug0;
                g.update_surface();
                double de = g.stretch_energy(g.heidtoindex[g.Nhelast - 2]);
                de += g.dimer_bend_energy(g.heidtoindex[g.Nhelast - 2]) + g.dimer_bend_energy(heindex0);
                de += g.bend_energy(heindex0) + g.bend_energy(xidindex);//-e1;

                
                //cout << " crit is " << crit << endl;
                de += gbb - g.mu[etypenew];

                

                /*int vinid=g.he[g.heidtoindex[g.Nhelast - 1]].vin;
                int vindex0=g.vidtoindex[vinid];
                
                if (g.v[vindex0].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                    {	
                        g.update_neigh_vertex(*it);
                    }
                }
                int voutid=g.he[g.heidtoindex[g.Nhelast - 1]].vout;
                vindex0=g.vidtoindex[voutid];
                if (g.v[vindex0].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                    {	
                        g.update_neigh_vertex(*it);
                    }
                }*/

                double crit = exp((-de) / g.T) / 2;
                int overlapflag=-1;
                if (gsl_rng_uniform(r) < crit ) 
                {
                    for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                    {	
                        g.update_neigh_vertex(*it);
                    }
                    
                    for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                    {	
                        if (g.check_overlap_g(*it)<0) overlapflag=1;

                    }
                    
                    if (overlapflag==-1) {
                        //cout << "monomer added no overlap" <<endl;
                        vecupdate.clear();
                        return 1;
                    }
                }
                else
                {

                    int nextid0 = g.he[g.heidtoindex[g.Nhelast - 2]].nextid;
                    int previd0 = g.he[g.heidtoindex[g.Nhelast - 2]].previd;
                    if (nextid0 == -1 || previd0 == -1)
                    {
                        cout << "not accepted in add monomer!" << endl;
                        exit(-1);
                    }
                    //cout << "!!!!!!!" <<endl;
                    int xidindex = g.heidtoindex[xid];
                    if (g.delete_edge(g.he[g.Nhe - 1].id) > 0)
                    {
                        
                        g.he[heindex0].nextid = -1;
                        g.he[xidindex].previd = -1;
                       
                        g.update_half_edge(heid0);
                        g.update_half_edge(g.he[heindex0].opid);
                        g.update_half_edge(xid);
                        g.update_half_edge(g.he[xidindex].opid);
                        g.update_index();
                        //cout << " MONOMER REMOVED AFTER addition xid2" <<endl;


                        /*vindex0=g.vidtoindex[vinid];
                        
                        if (g.v[vindex0].vneigh.size() > 0)
                        {
                            for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                            {	
                                g.update_neigh_vertex(*it);
                            }
                        }
                        vindex0=g.vidtoindex[voutid];
                        if (g.v[vindex0].vneigh.size() > 0)
                        {
                            for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                            {	
                                g.update_neigh_vertex(*it);
                            }
                        }*/
                        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                        vecupdate.clear();
                        return -1;
                    }
                    else
                    {
                        cout << " could not delete_monomer HERE 555";
                        exit(-1);
                    }
                }
            }
            else
            {
                //cout << "could not add monomer" << endl;
                return -1;
            }
        }
    }
    else if (g.v[g.vidtoindex[g.he[heindex0].vin]].hein.size() < 6 && (g.v[g.vidtoindex[g.he[heindex0].vout]].hein.size()) < 6) //adding dimer
    {
        if (vin0>=0 || vid0>=0) {return(-1);}
       

        if (etypeheid0 == 0)
        {
            if (gsl_rng_uniform(r) < .5)
            {
                etypenew1 = 1;
                etypenew2 = 2;
            } //ba
            else
            {
                etypenew1 = 0;
                etypenew2 = 0;
            }
        }

        else if (etypeheid0 == 3)
        {
            if (gsl_rng_uniform(r) < .5)
            {
                etypenew1 = 3;
                etypenew2 = 3;
            } //ba
            else
            {
                etypenew1 = 1;
                etypenew2 = 2;
            }
        }

        else if (etypeheid0 == 1)
        {
            if (gsl_rng_uniform(r) < .5)
            {
                etypenew1 = 2;
                etypenew2 = 0;
            } //ba
            else
            {
                etypenew1 = 2;
                etypenew2 = 3;
            }
        }

        else if (etypeheid0 == 2)
        {
            if (gsl_rng_uniform(r) < .5)
            {
                etypenew1 = 0;
                etypenew2 = 1;
            }
            else
            {
                etypenew1 = 3;
                etypenew2 = 1;
            }
        }
        //cout << "adding dimer"<<endl;
        // Based on concentration
        /*    if (gsl_rng_uniform(r) < g.cdProb) {
                if (gsl_rng_uniform(r) < 0.5) {etypenew1=0;}
                else {etypenew1=3; }
            } 
            else{
                if (gsl_rng_uniform(r) < 0.5) {etypenew1=1;}
                else {etypenew1=2; } 

            } 

            if (gsl_rng_uniform(r) < g.cdProb) {
                if (gsl_rng_uniform(r) < 0.5) {etypenew2=0;}
                else {etypenew2=3; }
            } 
            else{
                if (gsl_rng_uniform(r) < 0.5) {etypenew2=1;}
                else {etypenew2=2; } 

            }*/
        //cout << "new types are etypenew1= " <<etypenew1 << " etypenew2 "  << etypenew2<<endl;
        //etypenew1=gsl_rng_uniform_int(r,4);
        //etypenew2=gsl_rng_uniform_int(r,4);
        bool drug1 = 0;
        bool drug2 = 0;
        //bool drug10=0;
        //bool drug20=0;
        /*if ((etypenew1==3) ||(etypenew1==0) ) {
            if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;}
            if (gsl_rng_uniform(r) < g.drugProb) {drug10=1;}
        }
         if ((etypenew2==3) ||(etypenew2==0) ) {
            if (gsl_rng_uniform(r) < g.drugProb) {drug2=1;}
            if (gsl_rng_uniform(r) < g.drugProb) {drug20=1;}
        }*/

        //cout << "types " << etype1 <<" " << etype2<<endl;
        //double gbb=g.find_gbb(etypeheid0,etypenew1,etypenew2);
        double gbb = g.find_dg(etypeheid0, etypenew1, drug1);
        gbb += g.find_dg(etypenew1, etypenew2, drug2);
        gbb += g.find_dg(etypenew2, etypeheid0, g.he[heindex0].din);
        //if
        //cout << "gbb is  " << gbb <<endl;
        int x = g.add_dimer(heid0, r, etypenew1, etypenew2);
    
        if (x < 0)
        {
            //cout << "could not add dimer " << endl;
            return (-1);
        }
        ///else {
        //cout << " dimer added , now check if it stays"  <<endl;
        // }

        g.update_half_edge(g.Nhelast - 2);
        //g.he[g.heidtoindex[g.Nhelast-2]].din=drug2; //double check in add_dimer
        g.update_half_edge(g.Nhelast - 1);
        //g.he[g.get_heindex(g.Nhelast - 1)].din=drug20;
        g.update_half_edge(g.Nhelast - 4);
        //g.he[g.get_heindex(g.Nhelast - 4)].din=drug1;
        g.update_half_edge(g.Nhelast - 3);
        //g.he[g.get_heindex(g.Nhelast - 3)].din=drug10;
        g.update_surface();
        

        ////cout<< "op of edge" << g.Nhelast-2 << " is " << g.he[g.heidtoindex[g.Nhelast-2]].opid<<endl;
        int index2 = g.heidtoindex[g.Nhelast - 2];
        int index4 = g.heidtoindex[g.Nhelast - 4];
        double de = g.stretch_energy(index2) + g.dimer_bend_energy(index2);
        de += g.stretch_energy(index4) + g.dimer_bend_energy(index4);
        de += g.bend_energy(heindex0) + g.dimer_bend_energy(heindex0);
       
        double vp = 4 / 3. * 4 * atan(1.) * g.xi * g.xi * g.xi;
       
        de += gbb - (g.mu[etypenew1] + g.mu[etypenew2]);
        
        double crit = 2 * vp * exp(-de / g.T);
        //cout << " crit is " << crit << endl;
        int overlapflag=-1;
        //g.update_index();
        g.update_neigh_vertex(g.Nvlast-1);
        //cout<<"updating neigh of added vetex"<<endl;
        /* save vecupdate */ 
        int vindex0=g.vidtoindex[g.Nvlast-1];
        vector<int> vecupdate;      
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {	
                vecupdate.push_back(*it);
                g.update_neigh_vertex(*it);
            }
        }

        if (g.check_overlap_g(g.Nvlast-1)<0) { overlapflag=1;}
        
        if (gsl_rng_uniform(r) < crit && overlapflag==-1 )
        {
            
            vecupdate.clear();
            return 2;
                
        }
        else
        { // removing the asded dimer

            if (g.delete_edge(g.Nhelast - 1) < 0)
            {
                cout << " could not delete_added_dimer HERE 777";
                exit(-1);
            }
            
            g.update_index();
            if (g.delete_edge(g.Nhelast - 3) < 0)
            {
                cout << " could not delete_dimer HERE 888";
                exit(-1);
            }
            
            g.update_index();
            

            if (g.delete_vertex(g.Nvlast - 1) > 0)
            {
                g.set_prev_next(heid0, -1, -1);
                         
                g.update_index();
                /* update neigh of vecupdate*/
                for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                {	
                    g.update_neigh_vertex(*it);
                }
                vecupdate.clear();

                //cout << " added dimer edges and vertices are deleted" <<endl;
                heindex0=g.heidtoindex[heid0];         
                // MAYBE NOT NEEDED XXX
                

                //g.update_surface();
                return 0;
            }
            else
            {
                cout << "HERE 1111";
                exit(-1);
            }
            //vecupdate.clear();
            //return -1;
        }
    }
    return -1;
}

int attempt_remove_monomer_dimer(geometry &g, int heid0, gsl_rng *r)
{
    //cout << "attempt_remove_monomer_dimer" << endl;
    //cout << "Nd is " <<g.Nd <<endl;

    if (!(g.is_surface(heid0) > 0))
    {
        cout << "not on surface cannot remove" << endl;
        exit(-1);
    }

    g.update_half_edge(heid0);
    int heindex0 = g.heidtoindex[heid0]; // index of this edge
    //int heid0type=g.he[heindex0].type;
    if ((g.he[heindex0].nextid != -1) || (g.he[heindex0].previd != -1))
    {
        //cout << " has next or previous cannot remove" << endl;
        return -1;
    }

    g.update_half_edge(g.he[heindex0].opid);

    int heopindex0 = g.heidtoindex[g.he[heindex0].opid]; // indexd of opposite edge
    int optype = g.he[heopindex0].type;
    int nextopid0 = g.he[heopindex0].nextid; // id of next of opposite edge
    int prevopid0 = g.he[heopindex0].previd;
    //cout << " nextopid0 " << nextopid0 << " prevopid0 " <<prevopid0 <<endl;
    if ((nextopid0 == -1) || (prevopid0 == -1))
    {
        cout << "wrong geometry" << endl;
        exit(-1);
    }
    int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
    int prevopindex0 = g.heidtoindex[prevopid0];
    int opnexttype = g.he[nextopindex0].type;
    int opprevtype = g.he[prevopindex0].type;

    if ((g.he[heindex0].din == 1) || g.he[heopindex0].din == 1)
    {
        //cout << " din=1 "<<endl;
        return (-1);
    }                                    // WITH DRUG NO REMOVAL
    int heid1 = g.he[nextopindex0].opid; // now back to this side
    int heid2 = g.he[prevopindex0].opid; // afetr vertex

    double gbb = 0;
    
    //DLETING MONOMER
    if (g.is_surface(heid1) < 0 && g.is_surface(heid2) < 0 && g.is_vsurface(g.he[nextopindex0].vout) < 0)
    { //delete monomer
        //cout << " deleting monomer" <<endl;
        //if (g.is_vsurface(g.he[g.heidtoindex[nextopid0]].vout) > 0) { cout << "v on surface wrong geometry!"<<endl; exit(-1);}
        if ((g.he[nextopindex0].din == 1) || g.he[prevopindex0].din == 1)
        {
            return (-1);
        }


        double de = -(g.stretch_energy(heindex0));
        //if (g.he[heindex0].previd!=-1) { de-=g.dimer_bend_energy(g.get_heindex(g.he[heindex0].previd)); }
        //if (g.he[heindex0].nextid!=-1) {de-=g.dimer_bend_energy(heindex0); }
        de -= (g.dimer_bend_energy(heopindex0) + g.dimer_bend_energy(prevopindex0));
        //de -= g.bend_energy(nextopindex0) + g.bend_energy(prevopindex0); //g.monomer_energy(heid0);
        
        gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
        gbb += g.find_dg(optype, opnexttype, g.he[nextopindex0].din);
        

        de -= (gbb - g.mu[g.he[heindex0].type]);
        double crit = 2 * exp(-de / g.T); ///(2.0*g.z*g.K*g.K*g.K);
        if (gsl_rng_uniform(r) < crit)
        {
            //g.Nd-=g.he[heopindex0].din;
            //g.Nd-=g.he[heindex0].din;
            //cout << "REMOVING MONOMER two drugs removed " << g.he[heopindex0].din + g.he[heindex0].din <<endl;

            //cout << "de is" <<de <<endl;
            //cout << "crit is  " << crit <<endl;
            //if ( g.he[heindex0].previd!=-1) {g.he[g.get_heindex(g.he[heindex0].previd)].nextid=-1;}
            //if (g.he[heindex0].nextid!=-1) {g.he[g.get_heindex(g.he[heindex0].nextid)].previd=-1; }
            //cout << "g.he[heopindex0].previ " << g.he[heopindex0].previd <<endl;
            //cout << "g.he[prevopindex0].nextid " << g.he[prevopindex0].nextid <<endl;
            //cout << "g.he[heopindex0].nextid " << g.he[heopindex0].nextid <<endl;
            //cout << "g.he[nextopindex0].previd" << g.he[nextopindex0].previd <<endl;
            if (g.he[heopindex0].previd != -1)
            {
                g.he[prevopindex0].nextid = -1;
            }
            if (g.he[heopindex0].nextid != -1)
            {
                g.he[nextopindex0].previd = -1;
            }
            int vidin=g.he[heindex0].vin;
            int vidout=g.he[heindex0].vout;
            int x = g.delete_edge(heid0);
            //g.update_edge();_

            if (x < 0)
            {
                cout << "could not delete " << endl;
                exit(-1);
            }
            //

            //g.set_prev_next(nextopid0, -1, prevopid0);
            //g.set_prev_next(prevopid0, nextopid0, -1);
            //cout << "monomer removed" <<endl;
            //cout << " 004 g.Nd is " <<g.Nd<<endl;
            g.update_index();
            g.update_neigh_vertex(vidin);
            g.update_neigh_vertex(vidout);
            
            //cout << " 005 g.Nd is " <<g.Nd<<endl;
            return 1;
        }
        else
        {
            //cout << " remove monomer not accepted  " << endl;
            //g.add_monomer(nextopid0, prevopid0,optype);

            return -1;
        }
    }
    else
    {
        int heindex1 = g.heidtoindex[heid1];
        int heindex2 = g.heidtoindex[heid2];
        //double *vco=new double[3];
        //DLETING DIMER
        if (g.is_surface(heid1) > 0 && g.he[heindex1].previd == -1) //delete dimer this and next
        {
            //cout <<"deleting dimer with next"<<endl;
            if (g.he[heindex1].din == 1 || g.he[nextopindex0].din == 1)
            {
                return (-1);
            } //NOREMOVAL DRUG

            // int yid=g.he[heindex0].nextid;
            int vi = g.he[heopindex0].vout;

            //gbb+= g.find_gbb(optype,opnexttype,opprevtype); //whole triangle
            gbb += g.find_dg(optype, opnexttype, g.he[g.heidtoindex[nextopid0]].din);
            gbb += g.find_dg(opnexttype, opprevtype, g.he[prevopindex0].din);
            gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
            double de = -(g.stretch_energy(heindex0) + g.stretch_energy(g.heidtoindex[heid1]));
            de -= g.bend_energy(g.heidtoindex[heid2]);

            de -= (g.dimer_bend_energy(heopindex0));
            de -= (g.dimer_bend_energy(nextopindex0) + g.dimer_bend_energy(prevopindex0));
            //if (yid!=-1) {  de-=g.dimer_bend_energy(heindex0);}

            //if (xid!=-1 && g.he[xindex].nextid==heid1) {
            //    gbb+=g.find_dg( g.he[xindex].type,g.he[g.heidtoindex[heid1]].type);
            //     de-=g.dimer_bend_energy(xindex);
            //xidconnect=1;
            //}

            de -= (gbb - (g.mu[g.he[heindex0].type] + g.mu[g.he[heindex1].type]));
            double vp = 4 / 3. * 4 * atan(1.) * g.xi * g.xi * g.xi;
            double crit = exp(-de / g.T) / (2 * vp); ///vp; //*g.z*g.K*g.K);
                //cout << " crit is " << crit << endl;
            if (gsl_rng_uniform(r) < crit)
            {
                vector<int> vecupdate;
                if (g.v[g.vidtoindex[vi]].vneigh.size()>0){
                    for (vector<int>::iterator it =g.v[g.vidtoindex[vi]].vneigh.begin(); it != g.v[g.vidtoindex[vi]].vneigh.end(); ++it)
                    {	
                        vecupdate.push_back(*it);
                    }
                }
                //cout << " 006 g.Nd is " <<g.Nd<<endl;

                //WITHDRUG No removal
                //g.Nd-=g.he[heindex0].din;
                //g.Nd-=g.he[heopindex0].din;
                //g.Nd-=g.he[heindex1].din;
                //g.Nd-=g.he[nextopindex0].din;

                //cout << "REMOVING DIMER four drugs removed " << g.he[heopindex0].din + g.he[heindex0].din +  g.he[heindex1].din + g.he[nextopindex0].din<<endl;
                //cout << "removing heid0 and heid1 " << heid0 <<" " << heid1 <<endl;
                int success = g.delete_edge(heid0);
                if (success < 0)
                {
                    cout << "!!!" << endl;
                    exit(-1);
                }
                g.update_index();
                success = g.delete_edge(heid1); // next of op
                if (success < 0)
                {
                    cout << "!!!" << endl;
                    exit(-1);
                }
                g.update_index();
                int x = g.delete_vertex(vi);
                if (x < 0)
                {
                    exit(-1);
                }
                g.set_prev_next(prevopid0, -1, -1); //prev of op
                g.update_index();
                for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                {	
                    g.update_neigh_vertex(*it);
                }
                vecupdate.clear();
                return 2;
            }
        }
        else if (g.is_surface(heid2) > 0 && g.he[heindex2].nextid == -1)
        { //delete dimer this and prev

            //cout << "deleting dimer with prev of opid " <<endl;
            //if (g.he[g.heidtoindex[heid2]].din==1) return -1;
            //if (g.he[heindex2].din==1 || g.he[g.heidtoindex[prevopid0]].din==1 ) { return(-1);} //WITHDRUG no removal
            if (g.he[heindex2].din == 1 || g.he[prevopindex0].din == 1)
            {
                return (-1);
            } //NOREMOVAL DRUG
            int vi = g.he[heopindex0].vin;

            //gbb+= g.find_gbb(optype,opnexttype,opprevtype); //whole triangle
            gbb += g.find_dg(optype, opnexttype, g.he[g.heidtoindex[nextopid0]].din);
            gbb += g.find_dg(opnexttype, opprevtype, g.he[g.heidtoindex[prevopid0]].din);
            gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
            double de = -(g.stretch_energy(heindex0) + g.stretch_energy(g.heidtoindex[heid2]));
            de -= g.bend_energy(g.heidtoindex[heid1]);

            de -= (g.dimer_bend_energy(heopindex0));
            de -= (g.dimer_bend_energy(g.heidtoindex[nextopid0]) + g.dimer_bend_energy(g.heidtoindex[prevopid0]));

            de -= (gbb - (g.mu[g.he[heindex0].type] + g.mu[g.he[heindex2].type]));
            double vp = 4 / 3. * 4 * atan(1.) * g.xi * g.xi * g.xi;
            double crit = exp(-de / g.T) / (2 * vp); ///vp; //*g.z*g.K*g.K);
                //cout << " crit is " << crit << endl;
            if (gsl_rng_uniform(r) < crit)
            {
                //cout << " 008 g.Nd is " <<g.Nd<<endl;

                //g.Nd-=g.he[heindex0].din;
                //g.Nd-=g.he[heopindex0].din;
                //g.Nd-=g.he[heindex2].din;
                //g.Nd-=g.he[prevopindex0].din;
                //cout << "REMOVING DIMER four drugs removed " << g.he[heopindex0].din + g.he[heindex0].din +  g.he[heindex2].din + g.he[prevopindex0].din<<endl;
                vector<int> vecupdate;
                if (g.v[g.vidtoindex[vi]].vneigh.size()>0){
                    for (vector<int>::iterator it =g.v[g.vidtoindex[vi]].vneigh.begin(); it != g.v[g.vidtoindex[vi]].vneigh.end(); ++it)
                    {	
                        vecupdate.push_back(*it);
                    }
                }
                int success = g.delete_edge(heid0);
                if (success < 0)
                {
                    cout << "!!!heid2-0" << endl;
                    exit(-1);
                }
                g.update_index();
                success = g.delete_edge(heid2); // next of op
                if (success < 0)
                {
                    cout << "!!!heid2-2" << endl;
                    exit(-1);
                }
                g.update_index();
                int x = g.delete_vertex(vi);
                if (x < 0)
                {
                    exit(-1);
                }
                g.set_prev_next(nextopid0, -1, -1); //prev of op
                //g.update_half_edge(nextopid0);
                //g.update_half_edge(heid1);

                //multvec(g.v[g.vidtoindex(vi)].co,1,vco);

                //double de= g.compute_energy()-e1;
                //gbb=gb0next+gbnextprev+gb0prev;

                //delete[] vco;
                g.update_index();
                for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
                {	
                    g.update_neigh_vertex(*it);
                }
                vecupdate.clear();
                return 2;
            }
        }
        //delete[] vco;
    }

    return -1;
}




int attempt_vertex_fusion(geometry &g, gsl_rng *r)
{
    g.update_index();
    //cout << " Starting Fusion g.Nv " << g.Nv << endl;
    g.update_fusion_pairs();
   
    if (g.fusionv.size() == 0){  return -1;}
    int ind = gsl_rng_uniform_int(r, g.fusionv.size());
    int vidi = g.fusionv[ind];
    int vindexi = g.vidtoindex[vidi];
    int vidj = g.v[vindexi].fusion_vid;
    int vindexj = g.vidtoindex[vidj];
    //cout << " vidi " <<vidi <<"vidj " << vidj <<endl;

    if (g.v[vindexj].fusion_vid != vidi)
    {
        cout << " fusion pair wrong " << endl;
        exit(-1);
    }

    int heidnext = -1;
    int heidm = -1;
    int heidprev = -1;
    int vidtemp = -1;

    //find edges heidprev heidm heinext
    //first case !
    //cout << "start with vidi " << vidi << endl;
    heidnext = g.v[vindexi].hesurfinid;
    //cout << "heidnext is g.v[vindexi].hesurfinid " << g.v[vindexi].hesurfinid <<endl;
    vidtemp = g.he[g.heidtoindex[heidnext]].vin;
    //cout  << " first inside vertex is g.he[g.heidtoindex[heidnext]].vin" << g.he[g.heidtoindex[heidnext]].vin <<endl;
    heidm = g.v[g.vidtoindex[vidtemp]].hesurfinid;
    //cout << "heidm is g.v[g.vidtoindex[vidtemp]].hesurfinid " << g.v[g.vidtoindex[vidtemp]].hesurfinid <<endl;
    vidtemp = g.he[g.heidtoindex[heidm]].vin;
    //cout  << " second inside vertex is g.he[g.heidtoindex[heidm]].vin" << g.he[g.heidtoindex[heidm]].vin <<endl;
    heidprev = g.v[g.vidtoindex[vidtemp]].hesurfinid;
    //cout << "heidprev is g.v[g.vidtoindex[vidtemp]].hesurfinid " << g.v[g.vidtoindex[vidtemp]].hesurfinid <<endl;
    vidtemp = g.he[g.heidtoindex[heidprev]].vin;
    //cout  << " third vertex is g.he[g.heidtoindex[heidprev]].vin" << g.he[g.heidtoindex[heidprev]].vin <<endl;

    if (vidtemp != vidj)
    {
        //cout << "  vidtemp " <<  vidtemp <<endl;
        //cout << "now start with vidj " << vidj << endl;
        heidnext = g.v[vindexj].hesurfinid;
        //cout << "heidnext is g.v[vindexi].hesurfinid " << g.v[vindexi].hesurfinid <<endl;
        vidtemp = g.he[g.heidtoindex[heidnext]].vin;
        //cout  << " first inside vertex is g.he[g.heidtoindex[heidnext]].vin" << g.he[g.heidtoindex[heidnext]].vin <<endl;
        heidm = g.v[g.vidtoindex[vidtemp]].hesurfinid;
        //cout << "heidm is g.v[g.vidtoindex[vidtemp]].hesurfinid " << g.v[g.vidtoindex[vidtemp]].hesurfinid <<endl;
        vidtemp = g.he[g.heidtoindex[heidm]].vin;
        //cout  << " second inside vertex is g.he[g.heidtoindex[heidm]].vin" << g.he[g.heidtoindex[heidm]].vin <<endl;
        heidprev = g.v[g.vidtoindex[vidtemp]].hesurfinid;
        //cout << "heidprev is g.v[g.vidtoindex[vidtemp]].hesurfinid " << g.v[g.vidtoindex[vidtemp]].hesurfinid <<endl;
        vidtemp = g.he[g.heidtoindex[heidprev]].vin;
        //cout  << " third vertex is g.he[g.heidtoindex[heidprev]].vin" << g.he[g.heidtoindex[heidprev]].vin <<endl;
        if (vidtemp != vidi)
        {
            cout << "  vidtemp " << vidtemp << endl;
            cout << "wrong pairs vidin vidout" << endl;
            exit(-1);
        }
    }

    int bondmnext = g.is_bond_in_surface(heidnext);
    int bondprevm = g.is_bond_in_surface(heidm);
    if (bondmnext==-1 && bondprevm==-1) { //cout << "bondmnex=-1 no fusion"<<endl; 
        return(-1);}


    int vindtemp = vindexi;
    if (vindexi > vindexj)
    {
        vindexi = vindexj;
        vindexj = vindtemp;
        vidi = g.v[vindexi].vid;
        vidj = g.v[vindexj].vid;
    }

    if (g.v[vindexj].fusion_vid != vidi)
    {
        cout << " fusion pair wrong " << endl;
        exit(-1);
    }

    // check g.v[vindexi].hein.size() + g.v[vindexj].hein.size()

    /* FUSION pair exist */
    double e1 = 0;

    e1 += g.vertex_energy(vidi);
    e1 += g.vertex_energy(vidj);
/* save vecupdate old neighbors except vidi , vidj */
    
    vector<int> vecupdate;      
    for (vector<int>::iterator it =g.v[vindexj].vneigh.begin(); it != g.v[vindexj].vneigh.end(); ++it)
    {	
        if (*it!=vidi) vecupdate.push_back(*it);
    }
    
    for (vector<int>::iterator it =g.v[vindexi].vneigh.begin(); it != g.v[vindexi].vneigh.end(); ++it)
    {	
        if (*it!=vidj) vecupdate.push_back(*it);
    }

    /* save old vertex i */
    VTX *vtxi;
    vtxi = new VTX;
    g.save_vtx(vidi, vtxi);

    /* save old vertex j */
    VTX *vtxj;
    vtxj = new VTX;
    g.save_vtx(vidj, vtxj);

    double *tempv = new double[3];
    double *newv = new double[3];
    centvec(g.v[vindexi].co, g.v[vindexj].co, newv); //tempv if random
    //no random position
    //g.move_p(tempv, newv, r);
    
    /* add_new vertex  vid is newvid */
    g.add_vertex(newv);
    int newvid = g.Nvlast - 1;
    /* update index */
    g.update_index();
    //cout << " newvid is " << newvid << " vidtoindex[newvid] is " << g.vidtoindex[newvid] << " Nv is " << g.Nv << endl;
    delete[] tempv;
    delete[] newv;

    //cout << "new vertex added" << endl;
    
    /* merging vertices  */
    //cout<< " After addition g.Nv " << g.Nv <<endl;
    int newvindex = g.vidtoindex[newvid];
 
    for (vector<int>::iterator ithe = g.v[vindexi].hein.begin(); ithe != g.v[vindexi].hein.end(); ++ithe)
    {

        g.v[newvindex].hein.push_back(*ithe);
        int heindex0 = g.heidtoindex[*ithe];
        g.he[heindex0].vout = newvid;
        g.he[g.heidtoindex[g.he[heindex0].opid]].vin = newvid;
        g.update_half_edge(g.he[heindex0].id);
        g.update_half_edge(g.he[heindex0].opid);
    }
    

    for (vector<int>::iterator ithe = g.v[vindexj].hein.begin(); ithe != g.v[vindexj].hein.end(); ++ithe)
    {

        g.v[newvindex].hein.push_back(*ithe);
        int heindex0 = g.heidtoindex[*ithe];
        g.he[heindex0].vout = newvid;
        g.he[g.heidtoindex[g.he[heindex0].opid]].vin = newvid;
        g.update_half_edge(g.he[heindex0].id);
        g.update_half_edge(g.he[heindex0].opid);
    }
    
    
    

    if (vindexi > vindexj)
    {
        g.delete_vertex(vidi);
        g.delete_vertex(vidj);
    }
    else
    {
        g.delete_vertex(vidj);
        g.delete_vertex(vidi);
    }

    


    //cout<< " After delete g.Nv " << g.Nv <<endl;
    g.update_index();
    
    
    //can do only for those not connected

    g.set_prev_next(heidm, heidprev, heidnext);
    g.set_prev_next(heidprev, heidnext, heidm);
    g.set_prev_next(heidnext, heidm, heidprev);
    //cout << " in fusion next prev is set " << endl;

    g.update_normals();
    // cout << " in fusion normald updated " << endl;
    //final configuration  // binding energies!
    double e2 = g.find_dg(g.he[g.heidtoindex[heidnext]].type, g.he[g.heidtoindex[heidprev]].type, g.he[g.heidtoindex[heidprev]].din);

    if (bondmnext < 0)
    {
        e2 += g.find_dg(g.he[g.heidtoindex[heidm]].type, g.he[g.heidtoindex[heidnext]].type, g.he[g.heidtoindex[heidnext]].din);
    }
    if (bondprevm < 0)
    {
        e2 += g.find_dg(g.he[g.heidtoindex[heidprev]].type, g.he[g.heidtoindex[heidm]].type, g.he[g.heidtoindex[heidm]].din);
    }
    //cout << " now new structure " << endl;
    //int overlap = -1;
    
    e2 += g.vertex_energy(newvid);


    double vp = g.xi * g.xi * g.xi;
    double de = e2 - e1;
    //cout <<"de is " <<de <<endl;
    double crit = exp(-de / g.T) / (vp * 2);
    if (g.Test_assembly==1) { crit=1;}
    //cout << "crit is "<<crit <<endl;
    //cout << " for now neglecting overlap " <<endl;
    int overlapflag=-1;
    //
    /* update neighbors of "neighbors of vtxi and vtxj"  */
    for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
    {	
        g.update_neigh_vertex(*it);
        
    }
    
    g.update_neigh_vertex(newvid); 
    int newvindex0=g.vidtoindex[newvid];
    if (g.v[newvindex0].vneigh.size() > 0)
    {
        for (vector<int>::iterator it =g. v[newvindex0].vneigh.begin(); it != g.v[newvindex0].vneigh.end(); ++it)
        {	
            g.update_neigh_vertex(*it);
            //vecupdate.push_back(*it);
        }
    }
    
    if (g.check_overlap_g(newvid)<0) { overlapflag=1;}
    else{
        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
        {	
            
            if (g.check_overlap_g(*it)<0) { overlapflag=1;}
        }
    }
    if (gsl_rng_uniform(r)<crit && overlapflag==-1 )
    {
        cout <<"fusion crit is met"<<endl;
          
        //g.check_odd_neigh();
        
        vecupdate.clear();
        delete vtxi;
        delete vtxj;
        return 1;
        
    }
    else
    {
        //cout << "fusion not accepted  now move things back to what it was " << endl; //

        g.he[g.heidtoindex[heidnext]].nextid = -1;
        g.he[g.heidtoindex[heidprev]].previd = -1;
        if (bondmnext < 0)
        { 
            g.he[g.heidtoindex[heidm]].nextid = -1;
            g.he[g.heidtoindex[heidnext]].previd = -1;
        }
        if (bondprevm < 0)
        {
            g.he[g.heidtoindex[heidprev]].nextid = -1;
            g.he[g.heidtoindex[heidm]].previd = -1;
        }

    
        //cout << "fusion not accepted deleting  newvid "<< newvid << endl;
        //vector<int> vecupdate;      
        int vindextemp=g.vidtoindex[newvid];
        //cout << "newvid  vindex is "<< vindextemp << endl;
        //cout << " now delete its neighbors " << endl;
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it =g. v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {	
                vecupdate.push_back(*it);
            }
        }

        
        //cout << "g.Nv " << g.Nv <<endl;
        //cout << "fusion not accepted adding vtxi back" << endl;
        g.add_vertex(vtxi->co);
       
        int lastindex = g.vidtoindex[g.Nvlast - 1];
        for (vector<int>::iterator ithe = vtxi->hein.begin(); ithe != vtxi->hein.end(); ++ithe)
        {

            int heindex = g.heidtoindex[*ithe];
            g.v[lastindex].hein.push_back(*ithe);
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;
            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);

        }
        //cout << "g.Nv " << g.Nv <<endl;
        //cout << "fusion not accepted adding vtxj back" << endl;
        g.add_vertex(vtxj->co);
        
        lastindex = g.vidtoindex[g.Nvlast - 1];
        for (vector<int>::iterator ithe = vtxj->hein.begin(); ithe != vtxj->hein.end(); ++ithe)
        {
            int heindex = g.heidtoindex[*ithe];
            g.v[lastindex].hein.push_back(*ithe);
            ////cout << "updating vout of he " << g.he[heindex].id << "index " << heindex << " *ithe " << *ithe <<endl;
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;
            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
        }

        //cout <<" deleting ne vertex" <<endl;
        g.delete_vertex(newvid);
        
        //cout <<" deleted , updating index" <<endl;
        //g.update_index();
        //cout << "g.Nv " << g.Nv <<endl;
       // cout << "fusion not accepted updating index here !!!!!!!!!!!! "<< endl;
        g.update_index();

        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
        {	
            
            //cout <<" *it"<< *it <<endl;
            if (*it !=newvid) { g.update_neigh_vertex(*it);}
        }
        g.update_neigh_vertex(g.Nvlast-1);
        int vindex0=g.vidtoindex[g.Nvlast-1];
        if (g.v[vindex0].vneigh.size() > 0)
	    {
		    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {	
                g.update_neigh_vertex(*it);
            }
        }	
        g.update_neigh_vertex(g.Nvlast-2);
        vindex0=g.vidtoindex[g.Nvlast-2];
	    if (g.v[vindex0].vneigh.size() > 0)
	    {
		    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {	
                g.update_neigh_vertex(*it);
            }
        }
        //cout<<"attempt_vertex_fusion not accepted" <<endl;
        //g.check_odd_neigh();
        delete vtxi;
        delete vtxj;
        vecupdate.clear();
        return 0;

    }
    return -1;
}



int attempt_vertex_fission(geometry &g, gsl_rng *r)
{
    int ind = gsl_rng_uniform_int(r, g.surfv.size());
    int vid0 = g.surfv[ind];

    int vindex0 = g.vidtoindex[vid0];

    //cout << "fission of vertex  " << vid0 <<endl;
    //cout << " in fission Nv is " << g.Nv << " Nvlast is " << g.Nvlast <<endl;

    if (g.v[vindex0].hein.size() < 4 || g.is_bond_vsurface(vid0) > 0)
    {
        return -1;
    }

    // choose where to break ->xid
    int indhe = gsl_rng_uniform_int(r, g.v[vindex0].hein.size());
    int xid = g.v[vindex0].hein[indhe];

    int xindex = g.heidtoindex[xid];
    int xidopid = g.he[xindex].opid;
    int xidnextid = g.he[xindex].nextid;
    int xidprevid = g.he[xindex].previd;

    //cout << " in fission xid is " << xid << "xidopid is" << xidopid  <<endl;
    if (g.is_surface(xid) > 0 || g.is_surface(xidopid) > 0)
    { 
        //cout << "fission return -1 111"<<endl;
        return -1;
    }

    if (xidnextid == -1)
    {
        cout << " why xidnextid  is -1" << endl;
        exit(-1);
    }
    int xidopidnext = g.he[g.heidtoindex[xidnextid]].opid;
    if (g.is_surface(xidopidnext) > 0)
    { 
        
        return -1;
    }

    if (g.is_vsurface(g.he[xindex].vin) > 0 || g.is_vsurface(g.he[g.heidtoindex[xidnextid]].vout) > 0)
    {
        
        return -1;
    }

    /* starting fission */

    //cout << "fission of vertex  " << vid0 <<endl;
    
    // saving neighbors of this
    vector<int> vecupdate;      
    int vindextemp=g.vidtoindex[vid0];
    for (vector<int>::iterator it =g. v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
    {	
        vecupdate.push_back(*it);
    }
    
    double e1 = 0;
    //energy of initial configuration
    e1 += g.vertex_energy(vid0);

    e1 += g.find_dg(g.he[g.heidtoindex[xid]].type, g.he[g.heidtoindex[xidnextid]].type, g.he[g.heidtoindex[xidnextid]].din);
    //int breakprevxid=-1;
    //int breaknextprev=-1;
    if (gsl_rng_uniform(r) < .5)
    { //BREAK prev- this
       // breakprevxid=g.he[g.heidtoindex[xidprevid]].vout;
        e1 += g.find_dg(g.he[g.heidtoindex[xidprevid]].type, g.he[g.heidtoindex[xid]].type, g.he[g.heidtoindex[xid]].din);
        g.he[g.heidtoindex[xidprevid]].nextid = -1;
        g.he[g.heidtoindex[xid]].previd = -1;
        //add a flag here
    }
    else
    {
        //breaknextprev=g.he[g.heidtoindex[xidnextid]].vout;
        //BREAK next-prev
        e1 += g.find_dg(g.he[g.heidtoindex[xidnextid]].type, g.he[g.heidtoindex[xidprevid]].type, g.he[g.heidtoindex[xidprevid]].din);
        g.he[g.heidtoindex[xidnextid]].nextid = -1;
        g.he[g.heidtoindex[xidprevid]].previd = -1;
    }
    
    
    // BREAK xid-next
    g.he[g.heidtoindex[xid]].nextid = -1;
    g.he[g.heidtoindex[xidnextid]].previd = -1;

    int newvid1 = -1;
    int newvid2 = -1;

    double *newv1 = new double[3];
    double *newv2 = new double[3];
    double *tempv = new double[3];
    g.move_p(g.v[vindex0].co, tempv, r);  // move one in one direction //
    multvec(tempv,.5,newv1);
    //g.move_p(g.v[vindex0].co, newv2, r);
    // move the other one in opposit direction
    subvec(newv1,g.v[vindex0].co,tempv);
    addvec(g.v[vindex0].co,tempv,newv2);

    int xidnextopindex = g.heidtoindex[g.he[g.heidtoindex[xidopid]].nextid];
    double len1 = veclen(g.he[xidnextopindex].hecent, newv1);

    double len2 = veclen(g.he[xidnextopindex].hecent, newv2);

    //cout << "len1 is "<< len1 << " len2 is " << len2 <<endl;

    /* store old vtx */
    VTX *vtxi;
    vtxi = new VTX;
    g.save_vtx(vid0, vtxi);

    
    //cout << "vtxi updated hein vtxi->hein.size()" <<vtxi->hein.size() <<endl;

    g.add_vertex(newv1);
    g.add_vertex(newv2);

    //cout << " in fission Nv is " << g.Nv << " Nvlast is " << g.Nvlast <<endl;
    if (len1 <= len2)
    {
        newvid1 = g.Nvlast - 2; //xid->yid.vout
        newvid2 = g.Nvlast - 1; //xidnextid->zid.vin
    }
    else
    {
        newvid1 = g.Nvlast - 1;
        newvid2 = g.Nvlast - 2;
    }

    delete[] newv1;
    delete[] newv2;
    delete[] tempv;
    int yid = xid;
    int yidopid = g.he[g.heidtoindex[yid]].opid;
    while (true)
    {
        if (g.is_surface(yidopid) > 0)
        { //cout << "yid  shoul be -1" << yid <<endl;
            break;
        };
        g.v[g.vidtoindex[newvid1]].hein.push_back(yid);
        yidopid = g.he[g.heidtoindex[yid]].opid;

        g.he[g.heidtoindex[yid]].vout = newvid1;
        g.he[g.heidtoindex[yidopid]].vin = newvid1;
        g.update_half_edge(yid);
        g.update_half_edge(yidopid);

        
        yid = g.he[g.heidtoindex[yidopid]].previd;
    }
    //cout << " newvid1 updated" <<endl;

    int zid = xidnextid;
    int zidopid = -1; //g.he[g.heidtoindex[yid]].opid;
    while (true)
    {

        zidopid = g.he[g.heidtoindex[zid]].opid;
        g.v[g.vidtoindex[newvid1]].hein.push_back(zidopid);

        g.he[g.heidtoindex[zid]].vin = newvid2;
        g.he[g.heidtoindex[zidopid]].vout = newvid2;
        g.update_half_edge(zid);
        g.update_half_edge(zidopid);

        zid = g.he[g.heidtoindex[zidopid]].nextid;
        if (zid == -1)
        { //cout << "zid  is -1" << zid <<endl;
            break;
        }
    }

    //cout << "now delete the vertex" <<endl;
    
    

    
  


    g.delete_vertex(vid0);
    g.update_surface();
    // calculate energy



    //g.update_normals();

    double e2 = 0;
    //energy of final
    e2 += g.vertex_energy(newvid1);
    e2 += g.vertex_energy(newvid2);

    double vp = g.xi * g.xi * g.xi;
    double de = e2 - e1;
    double crit = exp(-de / g.T) * (vp * 2);
    if (g.Test_assembly==1) { crit=1;}
    
    int overlapflag=-1;

    g.update_index();

    // update old neighbors
    for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
    {	
        
        if (*it!=vid0) g.update_neigh_vertex(*it);
    }
    
    //update new and neighbors
    g.update_neigh_vertex(g.Nvlast-2);
    vindextemp=g.vidtoindex[g.Nvlast-2];
    if (g.v[vindextemp].vneigh.size() > 0)
    {
        for (vector<int>::iterator it =g. v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
        {	
            g.update_neigh_vertex(*it);
        }
    }
    g.update_neigh_vertex(g.Nvlast-1);	
    vindextemp=g.vidtoindex[g.Nvlast-1];
    if (g.v[vindextemp].vneigh.size() > 0)
    {
        for (vector<int>::iterator it =g. v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
        {	
            g.update_neigh_vertex(*it);
            //vecupdate.push_back(*it);
        }
    }

    if (g.check_overlap_g(g.Nvlast-1)<0) { overlapflag=1;}
    else if (g.check_overlap_g(g.Nvlast-2)<0) { overlapflag=1;}
    else{
        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
        {	
            
            if (g.check_overlap_g(*it)<0) { overlapflag=1;}
        }
    }

    if (gsl_rng_uniform(r) < crit && overlapflag==-1)
    {
        cout << "fission accepted added vid " << newvid1 << " " << newvid2 << endl;
        //cout << "#############" << endl;
        
         
        delete vtxi;
        //g.check_odd_neigh();
        return 1;
            
    }
    else
    {

        //vector<int> vecupdate;      
        vindextemp=g.vidtoindex[g.Nvlast-1];
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it =g. v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {	
                vecupdate.push_back(*it);
            }
        }
        vindextemp=g.vidtoindex[g.Nvlast-2];
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it =g. v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {	
                vecupdate.push_back(*it);
            }
        }
        //cout << "fission not accepted " << endl; // now move things back to what it was

        g.add_vertex(vtxi->co);
        //update added neighbors
        
        
        //cout << "vtxi is added back -- updating indices g.Nvlast " << g.Nvlast  << " g.Nv "<< g.Nv <<endl;

        int lastindex = g.vidtoindex[g.Nvlast - 1];

        for (vector<int>::iterator ithe = vtxi->hein.begin(); ithe != vtxi->hein.end(); ++ithe)
        {
            g.v[lastindex].hein.push_back(*ithe);
            int heindex = g.heidtoindex[*ithe];
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;

            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
            //cout << " heid is " << g.he[heindex].id <<endl;
            //cout << " heid opid is " << g.he[heindex].opid <<endl;
        }

        g.set_prev_next(xid, xidprevid, xidnextid);
        g.set_prev_next(xidprevid, xidnextid, xid);
        g.set_prev_next(xidnextid, xid, xidprevid);

        int newvindex1 = g.vidtoindex[newvid1];
        int newvindex2 = g.vidtoindex[newvid2];
        //cout << " try deleting BEFORE!!!!! putting back the xtzi" <<endl;
        if (newvindex1 > newvindex2)
        {
            g.delete_vertex(g.v[newvindex1].vid);
            g.delete_vertex(g.v[newvindex2].vid);
        }
        else
        {
            g.delete_vertex(g.v[newvindex2].vid);
            g.delete_vertex(g.v[newvindex1].vid);
        }
        ///cout <<"no fission updating neighbor list"<<endl;
        //UPDATE INDICES
        
        //cout << "0000 fission" <<endl;
        g.update_index();
        for (vector<int>::iterator it =vecupdate.begin(); it != vecupdate.end(); ++it)
        {	
            
            //cout <<" *it"<< *it <<endl;
            if ((*it !=newvid1) && (*it !=newvid2)) { g.update_neigh_vertex(*it);}
        }
        
        //cout << "2222 fission" <<endl;
        
        g.update_neigh_vertex(g.Nvlast-1);
        int vindex0=g.vidtoindex[g.Nvlast-1];
        //cout << "updating fission vtx_neigh 3333 fission" <<endl;
	    if (g.v[vindex0].vneigh.size() > 0)
	    {
		    for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {	
                g.update_neigh_vertex(*it);
            }
        }	

        /*
        g.check_odd_neigh();
        
        for (unsigned int i = 0; i < g.v.size(); i++)
        {    
            if (g.check_overlap_g(g.v[i].vid) < 0) { 
                
                dump_lammps_data_file(g, 5555555);
                cout << "after move vertex" <<endl;
                cout <<"OVERLAP vid " << g.v[i].vid << "vindex is " << g.vidtoindex[g.v[i].vid] << endl;
                //g.find_overlap_g(g.v[ind].vid);
                exit(-1); 
            } 
        } */
        vecupdate.clear();
        delete vtxi;
        return 0;
    }
    
    return -1;
}

int attempt_change_edge_type(geometry &g, int heid0, gsl_rng *r)
{
    int newtype = -1;
    int newoptype = -1;
    newtype = gsl_rng_uniform_int(r, 4);

    int heindex0 = g.heidtoindex[heid0];
    int oldtype = g.he[heindex0].type;

    int heopindex0 = g.heidtoindex[g.he[heindex0].opid]; // indexd of opposite edge
    int optype = g.he[heopindex0].type;

    if (oldtype == 3 || oldtype == 0)
    {
        if (g.he[heindex0].din == 1 || g.he[heopindex0].din == 1)
        {
            return -1;
        }
        if (g.he[heindex0].nextid != -1)
        {
            if (g.he[g.heidtoindex[g.he[heindex0].nextid]].din == 1)
                return -1;
        }
        if (g.he[heopindex0].nextid != -1)
        {
            if (g.he[g.heidtoindex[g.he[heopindex0].nextid]].din == 1)
                return -1;
        }
    }

    if (oldtype == newtype)
        return -1;
    if (newtype == 0)
    {
        newoptype = 3;
    }
    else if (newtype == 3)
    {
        newoptype = 0;
    }
    else if (newtype == 2)
    {
        newoptype = 1;
    }
    else if (newtype == 1)
    {
        newoptype = 2;
    }

    //cout << "edge heid0  "  << heid0 <<"heindex0 "<< heindex0 << "opid " << g.he[heindex0].opid <<endl;
    double e1 = 0;

    g.update_half_edge(heid0);
    g.update_half_edge(g.he[heindex0].opid);

    int nextid0 = g.he[heindex0].nextid; // indexd of next of opposite edge

    if (nextid0 != -1)
    {
        int nextindex0 = g.heidtoindex[nextid0];
        int nexttype = g.he[nextindex0].type;
        e1 += g.find_dg(oldtype, nexttype, g.he[nextindex0].din);
        e1 += g.bend_energy(nextindex0);
    }
    //cout << "one  " <<endl;
    int previd0 = g.he[heindex0].previd; // indexd of prev of opposite edge
    if (previd0 != -1)
    {
        int previndex0 = g.heidtoindex[previd0];
        int prevtype = g.he[previndex0].type;
        e1 += g.find_dg(prevtype, oldtype, g.he[heindex0].din);
        e1 += g.bend_energy(previndex0) + g.dimer_bend_energy(previndex0);
    }

    //cout << "222  " <<endl;

    int nextopid0 = g.he[heopindex0].nextid; // id of next of opposite edge
    if (nextopid0 != -1)
    {
        int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
        int opnexttype = g.he[nextopindex0].type;
        e1 += g.find_dg(optype, opnexttype, g.he[nextopindex0].din);
        e1 += g.bend_energy(nextopindex0);
    }
    //cout << "333  " <<endl;
    int prevopid0 = g.he[heopindex0].previd;
    if (prevopid0 != -1)
    {
        int prevopindex0 = g.heidtoindex[prevopid0];
        int opprevtype = g.he[prevopindex0].type;
        e1 += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
        e1 += g.bend_energy(prevopindex0) + g.dimer_bend_energy(prevopindex0);
    }

    e1 += g.bend_energy(heindex0);
    e1 += g.dimer_bend_energy(heindex0) + g.dimer_bend_energy(heopindex0);

    //cout << "333  " <<endl;

    /*if (oldtype == 1 || oldtype == 2)
    {
        if (gsl_rng_uniform(r) < 0.5)
        {
            newtype = 0;
            newoptype = 3;
        }
        else
        {
            newtype = 3;
            newoptype = 0;
        }
    }
    else if (oldtype == 3 || oldtype == 0)
    {
        if (gsl_rng_uniform(r) < 0.5)
        {
            newtype = 1;
            newoptype = 2;
        }
        else
        {
            newtype = 2;
            newoptype = 1;
        }
    }*/

    g.he[heindex0].type = newtype;
    g.he[heopindex0].type = newoptype;

    double e2 = g.bend_energy(heindex0);
    e2 += g.dimer_bend_energy(heindex0) + g.dimer_bend_energy(heopindex0);
    if (nextid0 != -1)
    {
        int nextindex0 = g.heidtoindex[nextid0];
        int nexttype = g.he[nextindex0].type;
        e2 += g.find_dg(newtype, nexttype, g.he[nextindex0].din);
        e2 += g.bend_energy(nextindex0);
    }
    //e2+=g.find_dg(newtype,nexttype,g.he[nextindex0].din);
    if (previd0 != -1)
    {
        int previndex0 = g.heidtoindex[previd0];
        int prevtype = g.he[previndex0].type;
        e2 += g.find_dg(prevtype, newtype, g.he[heindex0].din);
        e2 += g.bend_energy(previndex0) + g.dimer_bend_energy(previndex0);
    }
    //e2+=g.find_dg(prevtype,newtype,g.he[heindex0].din);
    if (nextopid0 != -1)
    {
        int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
        int opnexttype = g.he[nextopindex0].type;
        e2 += g.find_dg(newoptype, opnexttype, g.he[nextopindex0].din);
        e2 += g.bend_energy(nextopindex0);
    }

    if (prevopid0 != -1)
    {
        int prevopindex0 = g.heidtoindex[prevopid0];
        int opprevtype = g.he[prevopindex0].type;
        e2 += g.find_dg(opprevtype, newoptype, g.he[heopindex0].din);
        e2 += g.bend_energy(prevopindex0) + g.dimer_bend_energy(prevopindex0);
    }
    //e2+=g.find_dg(newoptype,opnexttype,g.he[nextopindex0].din);
    //e2+=g.find_dg(opprevtype,newoptype,g.he[heopindex0].din);

    double de = e2 - g.mu[newtype] - (e1 - g.mu[oldtype]);
    double crit = exp(-de / g.T);
    if (gsl_rng_uniform(r) < crit)
    {
        //cout <<" change accepted";
        return newtype;
    }
    else
    {
        g.he[heindex0].type = oldtype; // back to old type
        g.he[heopindex0].type = optype;
        g.update_half_edge(heid0);
        g.update_half_edge(g.he[heindex0].opid);
        //cout << " change_type not accepted"<<endl;
        return -1;
    }
}


int attempt_bind_wedge_dimer(geometry &g, int vid0, gsl_rng *r)
{
    
    //cout << "attempt_bind_wedge_dimer vid0 " << vid0 << endl;
    //cout << " g.Nd is " <<g.Nd<<endl;

    
    int heidprev=g.v[g.vidtoindex[vid0]].hesurfinid;
    int heindexprev=g.heidtoindex[heidprev];
    if (g.he[heindexprev].previd!=-1) { //cout<< " other end is bound in open wedge" <<endl; 
        return -1;}

    //cout <<"bind wedge for heidprev " << heidprev << endl;
    int heidnext=g.next_open_wedge(heidprev);
    

    if (heidnext==-1 ) { //cout <<"no-open-wedge" <<endl; 
    return -1;}
    
    if (g.no_bond_surface(heidnext) < 0) { //cout<< " other end is bound in open wedge" <<endl; 
    return -1;}
   
    int heindexnext=g.heidtoindex[heidnext];
    //cout <<"next open wedge is  " << heidnext << endl;
    
    if (g.he[heindexprev].nextid!=-1 || g.he[heindexnext].previd!=-1) { cout << " wrong geometry in bind-wedge " <<endl; exit(-1);}

    g.he[heindexprev].nextid = heidnext;
    g.he[heindexnext].previd = heidprev;
    if (heidprev==-1) { cout << "what! previndex" <<endl; exit(-1);}
    g.get_normal(heidnext);
    g.get_normal(heidprev);
    //if (g.check_bind_wedge(heidnext)<1 ||g.check_bind_wedge(heidprev)<1){ return -1;}
    //cout << "calculate de " <<endl;
    double de = g.dimer_bend_energy(heindexprev);//+ g.bend_energy(g.heidtoindex[heidnext])+g.bend_energy(g.heidtoindex[heidprev]);
    de += g.find_dg(g.he[heindexprev].type, g.he[heindexnext].type, g.he[heindexnext].din);
    double crit = exp((-de) / g.T);
    
    if (gsl_rng_uniform(r) < crit)
    {
        
        //cout << "bound! heidpre to heid next"  << heidprev  << heidnext << endl;
        return(1);
    }
    else
    {
        //cout << "Did not bound! crit is " << crit <<endl;
        g.he[heindexprev].nextid = -1;
        g.he[heindexnext].previd = -1;
        g.update_normals();
    }
    return(-1);
}


int attempt_unbind_wedge_dimer(geometry &g, int vid0, gsl_rng *r)
{
    //cout << "attempt_unbind_wedge_dimer vid0 " << vid0 << endl;
    if ((g.is_bond_vsurface(vid0))<0) { return -1;}

    int previd=g.v[g.vidtoindex[vid0]].hesurfinid;
    if (previd == -1 )
    {
        cout << "wrong previd attempt_unbind_wedge" << endl;
        exit(-1);
    }
    int previndex=g.heidtoindex[previd];
    int nextid = g.he[previndex].nextid;
    if (nextid == -1 )
    {
        cout << "wrong nextid attempt_unbind_wedge" << endl;
        exit(-1);
    }
    int nextindex=g.heidtoindex[nextid];
    //cout << "attempt_unbind_wedge_dimer" <<endl;
    //cout << " g.Nd is " <<g.Nd<<endl;
    
    double de = -g.find_dg(g.he[previndex].type, g.he[nextindex].type, g.he[nextindex].din);
    de -= g.dimer_bend_energy(previndex)+g.bend_energy(nextindex)+g.bend_energy(previndex);
        //cout <<"de is "<<de<<endl;
        //if (de==0) return 0;
    double crit = exp((-de) / g.T);
    //if (g.check_bind_wedge(previd)<1 ||g.check_bind_wedge(nextid)<1){ crit=1;}
    if (gsl_rng_uniform(r) < crit)
    {
        g.he[previndex].nextid = -1;
        g.he[nextindex].previd = -1;
        return 1;
    } 
    return -1;
}
int attempt_add_drug(geometry &g, int heid0, gsl_rng *r)
{

    //cout << "in attempt adding drug " << heid0 <<endl;
    //cout << " g.Nd is " <<g.Nd<<endl;
    int heindex0 = g.heidtoindex[heid0];
    if (g.is_surface(heid0)>0 ) { return(-1);}
    if (g.he[heindex0].din == 1)
        return 0; //already has drug
    int hetype = g.he[heindex0].type;
    //int nexttype = -1;
    int previd0 = g.he[heindex0].previd;
    if (previd0 == -1 ) { return -1;} // 
    int previndex0 = g.heidtoindex[previd0];
    
    //if (hetype == 1 || hetype == 2 || g.he[previndex0].type==1) //nexttype == 1 || nexttype == 2 ) //
    //{
    //    return -1;
    //}
    double e1 = 0; //(g.gdrug-g.mudrug);
   // int previd0 = g.he[heindex0].previd;

   // if (previd0 != -1)
   // {
        
        
        e1 += g.dimer_bend_energy(previndex0);
        e1 += g.find_dg(g.he[previndex0].type, hetype, 0);
    //}

    g.he[heindex0].din = 1;

    double e2 = 0;
    //if (previd0 != -1)
    //{
        //int previndex0 = g.heidtoindex[previd0];
        e2 += g.dimer_bend_energy(previndex0);
        e2 += g.find_dg(g.he[previndex0].type, hetype, 1);
    //}
    //else
    //{
      //  e2 += (g.gdrug-g.mudrug);
    //}
    double crit = exp((-(e2 - e1)) / g.T);
    if (gsl_rng_uniform(r) < crit)
    {
        //cout << "drug added on edge " << heindex0 <<endl;
        //cout << "prev is " << g.he[heindex0].previd <<endl;
        //cout << "prev is " << g.he[heindex0].previd <<endl;
        //cout<<"edgetype is " << hetype <<endl;

        //cout << " de is is " <<e2-e1<<endl;
        //cout << "crit is" << crit<<endl;
        g.update_half_edge(heid0);
        g.Nd++;
        //cout << " 010 g.Nd is " <<g.Nd<<endl;
        return 1;
    }
    else
    {
        g.he[heindex0].din = 0;
        return -1;
    }
}
int attempt_remove_drug(geometry &g, int heid0, gsl_rng *r)
{
    //cout << "in attempt remove drug " << heid0 <<endl;
    //cout << " g.Nd is " <<g.Nd<<endl;
     /* no drug on surface*/
    int heindex0 = g.heidtoindex[heid0];
    int previd0 = g.he[heindex0].previd;
    if (previd0 == -1) {return -1;}
    //if (previd0==-1 && g.he[heindex0].din!=0) { cout <<"WRONG DRUG POSITION N REMOVE" <<endl; exit(-1);}

    if (g.he[heindex0].din == 0)
        return 0; //no drug
    int hetype = g.he[heindex0].type;
    //int nexttype=-1;
    //if (g.he[heindex0].nextid!=-1) { nexttype=g.he[g.heidtoindex[g.he[heindex0].nextid]].type; }
    //if (hetype == 1 || hetype == 2) //|| nexttype == 1 || nexttype == 2 ) //

    //{
    //    cout << " why drug on other than 0 3 g.he[heindex0].type" << g.he[heindex0].type << endl;
        //exit(-1);
    //}
    double e1 = 0; //(g.gdrug-g.mudrug);

    //if (previd0 != -1)
    //{
        int previndex0 = g.heidtoindex[previd0];
        e1 += g.dimer_bend_energy(previndex0);
        e1 += g.find_dg(g.he[previndex0].type, hetype, 1);
        //cout<< "g.dimer_bend_energy(previndex0); " << g.dimer_bend_energy(previndex0)<<endl;
    /*}
    else
    {
        e1 += (g.gdrug[][]-g.mudrug);
    }*/

    //cout << "g.he[heindex0].din "  <<  g.he[heindex0].din <<endl;
    g.he[heindex0].din = 0;
    //cout << "try removal g.he[heindex0].din "  <<  g.he[heindex0].din <<endl;
    double e2 = 0;
    //if (previd0 != -1)
    //{
        //int previndex0 = g.heidtoindex[previd0];

        e2 += g.dimer_bend_energy(previndex0);
        e2 += g.find_dg(g.he[previndex0].type, hetype, 0);
        //cout<< "g.dimer_bend_energy(previndex0); " << g.dimer_bend_energy(previndex0)<<endl;
    //}
    double crit = exp((-(e2 - e1)) / g.T);
    //cout << " de removal drug is " <<e2-e1<<endl;
    //cout << " 011 g.Nd is " <<g.Nd<<endl;
    if (gsl_rng_uniform(r) < crit)
    {
        //cout << "drug removed from edge " << heindex0 <<endl;
        // cout<<"edgetype is " << hetype <<endl;

        //cout << "crit is" << crit<<endl;
        //cout << "drug removed" <<endl;
        g.Nd--;
        //cout << " 012 g.Nd is " <<g.Nd<<endl;
        return 1;
    }
    else
    {
        g.he[heindex0].din = 1;
        //cout << " did not remove drug " <<g.Nd<<endl;
        //cout << " 013 g.Nd is " <<g.Nd<<endl;

        //cout <<"g.he[g.heidtoindex[heid0)].din"<<g.he[g.heidtoindex[heid0)].din<<endl;
        return 0;
    }
}
void make_initial_triangle(geometry &g)
{
    double xyz0[3];
    xyz0[0] = 0;
    xyz0[1] = 0;
    xyz0[2] = 0;

    for (int i = 0; i < 3; i++)
    {
        g.add_vertex(xyz0);
        xyz0[0] = cos(i * PI / 3);
        xyz0[1] = sin(i * PI / 3);
        xyz0[2] = 0;
    }
    //for (int i=0; i<3; i++) {

    //	int j=i+1;
    //	if (i==2) { j=0;}
    //	if (i>=Nvlast || j>=Nvlast) { cout << "ERROR in make pentamer" << endl; exit(-1); }
    g.add_half_edge_type(g.v[0].vid, g.v[1].vid, 2);
    g.add_half_edge_type(g.v[1].vid, g.v[0].vid, 1);
    g.add_half_edge_type(g.v[1].vid, g.v[2].vid, 0);
    g.add_half_edge_type(g.v[2].vid, g.v[1].vid, 3);
    g.add_half_edge_type(g.v[2].vid, g.v[0].vid, 1);
    g.add_half_edge_type(g.v[0].vid, g.v[2].vid, 2);

    g.set_prev_next(g.he[0].id, g.he[4].id, g.he[2].id);
    g.set_prev_next(g.he[2].id, g.he[0].id, g.he[4].id);
    g.set_prev_next(g.he[4].id, g.he[2].id, g.he[0].id);

    for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); it++)
    {
        cout << "in make pentame updating edge" << it->id << endl;
        ;
        g.update_half_edge(it->id);
        cout << it->id << "in make pentamer  TYPE " << g.he[it->id].type << " opid " << it->opid << " OP TYPE " << g.he[g.heidtoindex[it->opid]].type << endl;
    }

    g.update_surface();
}
void make_initial_pentamer(geometry &g)
{
    double xyz0[3];
    xyz0[0] = 0;
    xyz0[1] = 0;
    xyz0[2] = 0;

    for (int i = 0; i < 3; i++)
    {
        g.add_vertex(xyz0);
        xyz0[0] = cos(i * PI / 3);
        xyz0[1] = sin(i * PI / 3);
        xyz0[2] = 0;
    }
    //for (int i=0; i<3; i++) {

    //	int j=i+1;
    //	if (i==2) { j=0;}
    //	if (i>=Nvlast || j>=Nvlast) { cout << "ERROR in make pentamer" << endl; exit(-1); }
    g.add_half_edge_type(g.v[0].vid, g.v[1].vid, 2);
    g.add_half_edge_type(g.v[1].vid, g.v[0].vid, 1);
    g.add_half_edge_type(g.v[1].vid, g.v[2].vid, 0);
    g.add_half_edge_type(g.v[2].vid, g.v[1].vid, 3);
    g.add_half_edge_type(g.v[2].vid, g.v[0].vid, 1);
    g.add_half_edge_type(g.v[0].vid, g.v[2].vid, 2);

    g.set_prev_next(g.he[0].id, g.he[4].id, g.he[2].id);
    g.set_prev_next(g.he[2].id, g.he[0].id, g.he[4].id);
    g.set_prev_next(g.he[4].id, g.he[2].id, g.he[0].id);

    for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); it++)
    {
        cout << "in make pentame updating edge" << it->id << endl;
        ;
        g.update_half_edge(it->id);
        cout << it->id << "in make pentamer  TYPE " << g.he[it->id].type << " opid " << it->opid << " OP TYPE " << g.he[g.heidtoindex[it->opid]].type << endl;
    }

    g.update_surface();
    //exit(-1);

    cout << " HEREH in pentamer 2" << endl;
    double *vco = new double[3];

    for (int x = 0; x < 3; x++)
    {
        g.update_surface();
        //dump_lammps_data_file(g, frame++);

        cout << " after g.update surface" << endl;
        g.new_vertex(g.Nhe - 1, vco);
        cout << vco[0] << " " << vco[1] << " " << vco[2] << endl;
        //exit(-1);
        g.force_add_dimer(g.Nhe - 1, vco, 0, 1);
        g.update_half_edge(g.Nhelast - 1);
        g.update_half_edge(g.Nhelast - 2);
        g.update_half_edge(g.Nhelast - 3);
        g.update_half_edge(g.Nhelast - 4);
    }
    //exit(-1);
    g.update_surface();
    g.force_add_monomer(1, g.Nhe - 1, 0);
    g.update_half_edge(g.Nhelast - 1);
    g.update_half_edge(g.Nhelast - 2);

    delete[] vco;
    //exit(-1);
    //int vind=-1;
    g.update_surface();
    /*for (vector<VTX>::iterator it = g.v.begin() ; it != g.v.end(); ++it) {
        fprintf(stderr,"\n%li %d %8.3f %8.3f %8.3f\n", distance(g.v.begin(),it), it->vid ,it->co[0], it->co[1], it->co[2]);
        for (vector<int>::iterator ithe = it->hein.begin() ; ithe != it->hein.end(); ithe++) {
            cout << "     hein " << *ithe <<endl;

        }
    }*/
    //vind++;
    /* for (vector<int>::iterator ithe = g.v[vind].hein.begin(); ithe != g.v[vind].hein.end(); ++ithe)
        { //update the geometry
            //if *ithe
            int heindex = g.*ithe;
            cout << " updating edge after move" << *ithe << endl;
            g.update_half_edge(*ithe);
            //cout << " updating edge " << *ithe << endl;
            g.update_half_edge(g.he[heindex].opid);
            if (g.he[heindex].previd!=-1) {

                g.update_half_edge(g.he[heindex].previd);
                g.update_half_edge(g.he[g.heidtoindex[g.he[heindex].previd)].opid);
            }
            //cout << " updating edge " << g.he[heindex].opid << endl;
            //if ( g.check_overlap_he(*ithe)<0 ) { overlapflag=-1;}
        } */
}

void make_seed(geometry &g, gsl_rng *r)
{

    int frame = 0;
    //int monomeradded=0;
    ///int dimeradded=0;
    //int monomerremoved=0;
    //int dimerremoved=0;

    make_initial_pentamer(g);
    double ee = g.compute_energy();

    //dump_lammps_data_file(g, frame++);
    for (int rstep = 0; rstep < 1000; rstep++)
    { //relaxing the shell
        move_vertex(g, r);
        g.update_surface();
    }
    ee = g.compute_energy();
    cout <<ee<<endl;
    //dump_lammps_data_file(g, frame++);
    //fprintf(stderr, "Graph initialized.\n");

    for (int rstep = 0; rstep < 5; rstep++)
    {
       
        g.add_dimer(3 + rstep * 4, r, 3, 3);
        g.update_surface();

        for (int rstep = 0; rstep < 1000; rstep++)
        { //relaxing the shell
            move_vertex(g, r);
            g.update_surface();
        }
    }
    //fprintf(stderr, "Dimers added \n");
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    //fprintf(stderr, "Adding more!  \n");
    cout << endl;
    for (int rstep = 0; rstep < 5; rstep++)
    {
        
        g.add_dimer(21 + rstep * 4, r, 1, 2);
        g.update_surface();
        ee = g.compute_energy();

        for (int rstep = 0; rstep < 1000; rstep++)
        { //relaxing the shell
            move_vertex(g, r);
            g.update_surface();
        }
        ee = g.compute_energy();

        //dump_lammps_data_file(g, frame++);
    }

    //fprintf(stderr, "More Dimers added \n");
    //cout << "TESTING BIND UNBIND" << endl;
    cout << endl;

    for (unsigned int step = 0; step < g.surfheid.size() * 2; step++)
    {
        int ind = gsl_rng_uniform_int(r, g.surfheid.size());

        int e = g.surfheid[ind];
        //cout << "e for binding " << e << endl;
        int tt = attempt_bind_wedge_dimer(g, e, r);
        if (tt > 0)
        {
            cout << "Bound " << endl;
            tt = 0;
        }
        //cout << "try binding"<<endl;
        g.update_surface();
        // dump_lammps_data_file(g, frame++);
        ind = gsl_rng_uniform_int(r, g.surfheid.size());
        //cout <<"ind id for unbinding" << ind<<endl;
        e = g.surfheid[ind];
        cout << "e for unbinding " << e << endl;
        tt = attempt_unbind_wedge_dimer(g, e, r);
        if (tt > 0)
        {
            cout << "UNNNNBound " << endl;
            tt = 0;
        }
        g.update_surface();
        //dump_lammps_data_file(g, frame++);
    }

    cout << "AFTER ALL BINDINGS" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << endl;
    //fprintf(stderr, "Adding monomer???  \n");
    for (int rstep = 0; rstep < 5; rstep++)
    {
        g.add_monomer_dimer(23 + rstep * 4); //change it
        g.update_surface();

        for (int rstep = 0; rstep < 1000; rstep++)
        { //relaxing the shell
            move_vertex(g, r);
            g.update_surface();
        }

        //dump_lammps_data_file(g, frame++);
    }
    dump_lammps_data_file(g, frame++);
    cout << "ALL monomers added" << endl;
    cout << "#########  NHE " << g.Nhe << " #######################" << endl;
    cout << endl;
}

void make_seed_T3(geometry &g, gsl_rng *r)
{
}