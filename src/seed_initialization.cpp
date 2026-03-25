#include "seed_initialization.hpp"

#include "geometry.hpp"

bool make_seed(geometry &g, gsl_rng *r, const std::string &seed_config)
{
    (void)r;

    if (seed_config == "triangle")
    {
        make_initial_triangle(g);
        return true;
    }

    if (seed_config == "pentamer")
    {
        make_initial_pentamer(g);
        return true;
    }

    if (seed_config == "hexamer")
    {
        make_initial_hexamer(g);
        return true;
    }

    return false;
}

void make_seed_T3(geometry &g, gsl_rng *r)
{
    (void)g;
    (void)r;
}
