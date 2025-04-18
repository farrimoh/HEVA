/*
 * montecarlo.hpp
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
