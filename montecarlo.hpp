/*
 * MonteCarlo.h
 *
 *  Created on: Apr 25, 2019
 *      Author: farri
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "geometry.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

void move_vertex(geometry &g, gsl_rng *r);
int move_one_vertex(geometry &g, int vid0, gsl_rng *r);
int attempt_add_monomer(geometry &g, int heid0, gsl_rng *r);
int attempt_add_dimer(geometry &g, int heid0, gsl_rng *r);
int attempt_add_monomer_dimer(geometry &g, int heid0, gsl_rng *r);
int attempt_remove_monomer(geometry &g, int heid0, gsl_rng *r);
int attempt_remove_dimer(geometry &g, int heid0, gsl_rng *r);
int attempt_remove_monomer_dimer(geometry &g, int heid0, gsl_rng *r);
int old_attempt_vertex_fusion(geometry &g, int heid0, gsl_rng *r);
int attempt_vertex_fusion(geometry &g, gsl_rng *r);

int attempt_wedge_fusion(geometry &g,  gsl_rng *r);
int attempt_wedge_fission(geometry &g, gsl_rng *r);

int attempt_fusion(geometry &g, gsl_rng *r);
int attempt_fission(geometry &g, gsl_rng *r);

int old_attempt_vertex_fission(geometry &g, int heid0, gsl_rng *r);
int attempt_vertex_fission(geometry &g, gsl_rng *r);
int attempt_change_edge_type(geometry &g, int heid0, gsl_rng *r);
int attempt_change_edge_type_tri(geometry &g, int heid0, gsl_rng *r);
//int add_monomer_dimer(geometry &g, int heid0, gsl_rng *r);
//int remove_monomer_dimer(geometry &g, int heid0, gsl_rng *r);
int old_attempt_bind_wedge_dimer(geometry &g, int heid0, gsl_rng *r);
int attempt_bind_wedge_dimer(geometry &g, int heid0, gsl_rng *r);

int old_attempt_unbind_wedge_dimer(geometry &g, int heid0, gsl_rng *r);
int attempt_unbind_wedge_dimer(geometry &g, int vid0, gsl_rng *r);

int attempt_bind_triangle(geometry &g, int heid0, gsl_rng *r);
int attempt_unbind_triangle(geometry &g, int heid0, gsl_rng *r);

int attempt_add_drug(geometry &g, int heid0, gsl_rng *r);
int attempt_remove_drug(geometry &g, int heid0, gsl_rng *r);
void make_seed(geometry &g, gsl_rng *r);
void make_seed_T3(geometry &g, gsl_rng *r);
void get_dimer_etypes(int etypeheid0, int etypenew1, int etypenew2, gsl_rng *r);
int force_add_monomer_with_next(geometry &g, int heid0, int xid,gsl_rng *r);

#endif
