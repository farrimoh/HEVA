#include "seed_initialization.hpp"

#include "geometry.hpp"

namespace
{
bool add_seed_dimers(geometry &g, gsl_rng *r, int count)
{
    double distance_vector[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < count; ++i)
    {
        g.update_boundary();
        if (g.boundary.empty())
        {
            return false;
        }

        const int heid = g.boundary.back();
        if (g.add_dimer(heid, r, 0, 1, distance_vector) < 0)
        {
            return false;
        }
    }

    g.update_boundary();
    return true;
}
}

bool make_seed(geometry &g, gsl_rng *r, const std::string &seed_config)
{
    if (seed_config == "triangle")
    {
        make_initial_triangle(g);
        return true;
    }

    if (seed_config == "pentamer")
    {
        make_initial_triangle(g);
        return add_seed_dimers(g, r, 2);
    }

    if (seed_config == "hexamer")
    {
        make_initial_triangle(g);
        return add_seed_dimers(g, r, 3);
    }

    return false;
}

void make_seed_T3(geometry &g, gsl_rng *r)
{
    (void)g;
    (void)r;
}
