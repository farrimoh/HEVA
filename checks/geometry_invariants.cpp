#include "geometry_invariants.hpp"

#include "../src/geometry.hpp"

#include <set>
#include <sstream>

namespace
{
bool is_valid_vertex_id(const geometry &g, int vid)
{
    return vid >= 0 && g.vidtoindex[vid] >= 0 && g.vidtoindex[vid] < static_cast<int>(g.v.size()) && g.v[g.vidtoindex[vid]].vid == vid;
}

bool is_valid_half_edge_id(const geometry &g, int heid)
{
    return heid >= 0 && g.heidtoindex[heid] >= 0 && g.heidtoindex[heid] < static_cast<int>(g.he.size()) && g.he[g.heidtoindex[heid]].id == heid;
}
}

bool check_geometry_invariants(const geometry &g, std::string &error_message)
{
    if (g.Nv != static_cast<int>(g.v.size()))
    {
        error_message = "Nv does not match the vertex container size.";
        return false;
    }

    if (g.Nhe != static_cast<int>(g.he.size()))
    {
        error_message = "Nhe does not match the half-edge container size.";
        return false;
    }

    std::set<int> vertex_ids;
    for (int i = 0; i < static_cast<int>(g.v.size()); ++i)
    {
        const VTX &vertex = g.v[i];
        if (!vertex_ids.insert(vertex.vid).second)
        {
            error_message = "Duplicate vertex id detected.";
            return false;
        }

        if (g.vidtoindex[vertex.vid] != i)
        {
            error_message = "vidtoindex does not map a vertex id back to the correct index.";
            return false;
        }

        for (std::size_t j = 0; j < vertex.hein.size(); ++j)
        {
            if (!is_valid_half_edge_id(g, vertex.hein[j]))
            {
                error_message = "Vertex incoming half-edge list contains an invalid half-edge id.";
                return false;
            }
        }

        for (std::size_t j = 0; j < vertex.vneigh.size(); ++j)
        {
            if (!is_valid_vertex_id(g, vertex.vneigh[j]))
            {
                error_message = "Vertex neighbor list contains an invalid vertex id.";
                return false;
            }
        }
    }

    std::set<int> half_edge_ids;
    for (int i = 0; i < static_cast<int>(g.he.size()); ++i)
    {
        const HE &edge = g.he[i];
        if (!half_edge_ids.insert(edge.id).second)
        {
            error_message = "Duplicate half-edge id detected.";
            return false;
        }

        if (g.heidtoindex[edge.id] != i)
        {
            error_message = "heidtoindex does not map a half-edge id back to the correct index.";
            return false;
        }

        if (!is_valid_vertex_id(g, edge.vin) || !is_valid_vertex_id(g, edge.vout))
        {
            error_message = "Half-edge endpoints reference invalid vertex ids.";
            return false;
        }

        if (!is_valid_half_edge_id(g, edge.nextid) || !is_valid_half_edge_id(g, edge.previd))
        {
            error_message = "Half-edge next/predecessor references are invalid.";
            return false;
        }

        const HE &next_edge = g.he[g.heidtoindex[edge.nextid]];
        const HE &prev_edge = g.he[g.heidtoindex[edge.previd]];
        if (next_edge.previd != edge.id)
        {
            error_message = "Half-edge next/predecessor linkage is not reciprocal.";
            return false;
        }
        if (prev_edge.nextid != edge.id)
        {
            error_message = "Half-edge predecessor/next linkage is not reciprocal.";
            return false;
        }

        if (edge.opid != -1)
        {
            if (!is_valid_half_edge_id(g, edge.opid))
            {
                error_message = "Half-edge opposite reference is invalid.";
                return false;
            }

            const HE &opposite_edge = g.he[g.heidtoindex[edge.opid]];
            if (opposite_edge.opid != edge.id)
            {
                error_message = "Half-edge opposite linkage is not reciprocal.";
                return false;
            }
        }

        if ((edge.nextid_boundary != -1 && !is_valid_half_edge_id(g, edge.nextid_boundary)) ||
            (edge.previd_boundary != -1 && !is_valid_half_edge_id(g, edge.previd_boundary)))
        {
            error_message = "Boundary half-edge linkage contains an invalid reference.";
            return false;
        }
    }

    std::set<int> boundary_ids;
    for (std::size_t i = 0; i < g.boundary.size(); ++i)
    {
        const int heid = g.boundary[i];
        if (!is_valid_half_edge_id(g, heid))
        {
            error_message = "Boundary list contains an invalid half-edge id.";
            return false;
        }

        if (!boundary_ids.insert(heid).second)
        {
            error_message = "Boundary list contains a duplicate half-edge id.";
            return false;
        }
    }

    int total_neighbors = 0;
    for (std::size_t i = 0; i < g.v.size(); ++i)
    {
        total_neighbors += static_cast<int>(g.v[i].vneigh.size());
    }
    if (total_neighbors % 2 != 0)
    {
        error_message = "The total vertex neighbor count is odd.";
        return false;
    }

    error_message.clear();
    return true;
}
