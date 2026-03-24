#ifndef HEVA_CHECKS_GEOMETRY_INVARIANTS_HPP
#define HEVA_CHECKS_GEOMETRY_INVARIANTS_HPP

#include <string>

class geometry;

bool check_geometry_invariants(const geometry &g, std::string &error_message);

#endif
