/*
 * MonteCarlo.h
 *
 *  Created on: Apr 25, 2019
 *      Author: farri
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "Geometry.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace std;

void move_vertex(Geometry &g, gsl_rng *r);
int move_one_vertex(Geometry &g, int vid0, gsl_rng *r);
int attempt_add_monomer(Geometry &g, int heid0, gsl_rng *r);
int attempt_add_dimer(Geometry &g, int heid0, gsl_rng *r);
int attempt_add_monomer_dimer(Geometry &g, int heid0, gsl_rng *r);
int attempt_remove_monomer(Geometry &g, int heid0, gsl_rng *r);
int attempt_remove_dimer(Geometry &g, int heid0, gsl_rng *r);
int attempt_remove_monomer_dimer(Geometry &g, int heid0, gsl_rng *r);
int old_attempt_vertex_fusion(Geometry &g, int heid0, gsl_rng *r);
int attempt_vertex_fusion(Geometry &g, gsl_rng *r);

int attempt_wedge_fusion(Geometry &g,  gsl_rng *r);
int attempt_wedge_fission(Geometry &g, gsl_rng *r);

int attempt_fusion(Geometry &g, gsl_rng *r);
int attempt_fission(Geometry &g, gsl_rng *r);

int old_attempt_vertex_fission(Geometry &g, int heid0, gsl_rng *r);
int attempt_vertex_fission(Geometry &g, gsl_rng *r);
int attempt_change_edge_type(Geometry &g, int heid0, gsl_rng *r);
int attempt_change_edge_type_tri(Geometry &g, int heid0, gsl_rng *r);
//int add_monomer_dimer(Geometry &g, int heid0, gsl_rng *r);
//int remove_monomer_dimer(Geometry &g, int heid0, gsl_rng *r);
int old_attempt_bind_wedge_dimer(Geometry &g, int heid0, gsl_rng *r);
int attempt_bind_wedge_dimer(Geometry &g, int heid0, gsl_rng *r);

int old_attempt_unbind_wedge_dimer(Geometry &g, int heid0, gsl_rng *r);
int attempt_unbind_wedge_dimer(Geometry &g, int vid0, gsl_rng *r);

int attempt_bind_triangle(Geometry &g, int heid0, gsl_rng *r);
int attempt_unbind_triangle(Geometry &g, int heid0, gsl_rng *r);

int attempt_add_drug(Geometry &g, int heid0, gsl_rng *r);
int attempt_remove_drug(Geometry &g, int heid0, gsl_rng *r);
void make_seed(Geometry &g, gsl_rng *r);
void make_seed_T3(Geometry &g, gsl_rng *r);
void get_dimer_etypes(int etypeheid0, int etypenew1, int etypenew2, gsl_rng *r);
int force_add_monomer_with_next(Geometry &g, int heid0, int xid,gsl_rng *r);

#endif
