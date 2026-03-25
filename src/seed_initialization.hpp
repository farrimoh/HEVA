#ifndef SEED_INITIALIZATION_HPP
#define SEED_INITIALIZATION_HPP

#include <gsl/gsl_rng.h>
#include <string>

class geometry;

void make_initial_triangle(geometry &g);
void make_initial_pentamer(geometry &g);
void make_initial_hexamer(geometry &g);

bool make_seed(geometry &g, gsl_rng *r, const std::string &seed_config);
void make_seed_T3(geometry &g, gsl_rng *r);

#endif
