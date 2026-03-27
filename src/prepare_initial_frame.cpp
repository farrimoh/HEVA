#include "geometry.hpp"

#include <cstdlib>
#include <iostream>

namespace
{
constexpr unsigned long kDefaultIndexCapacity = 1000000UL;
}

int main(int argc, char **argv)
{
    if (argc < 3 || argc > 4)
    {
        std::cerr << "usage: ./prepare_initial_frame INPUT_PATH OUTPUT_DIR [INDEX_CAPACITY]" << std::endl;
        return 1;
    }

    unsigned long index_capacity = kDefaultIndexCapacity;
    if (argc == 4)
    {
        char *end = nullptr;
        const unsigned long parsed = std::strtoul(argv[3], &end, 10);
        if (end == argv[3] || (end != nullptr && *end != '\0') || parsed == 0UL)
        {
            std::cerr << "Invalid INDEX_CAPACITY: " << argv[3] << std::endl;
            return 1;
        }
        index_capacity = parsed;
    }

    geometry g;
    g.initialize(4, index_capacity);
    set_output_directory(argv[2]);
    read_initial_frame_compat_data(g, argv[1]);
    dump_restart_lammps_data_file(g, 0);
    std::cout << output_file_path("restart_lammps.dat") << std::endl;
    return 0;
}
