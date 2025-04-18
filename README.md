# GCMC HEVA
**G**rand **C**anonical **M**onte **C**arlo Simulations using **H**alf **E**dge data-structure for **V**irus **A**ssembly

![Hepatitis-B Capsids](HBV.jpg)

## Overview
This repository contains the implementation of a Grand Canonical Monte Carlo (GCMC) simulation framework for studying virus capsid assembly, specifically focusing on the Hepatitis B Virus (HBV). The codebase utilizes a half-edge data structure for efficient geometric computations and assembly simulations.

## Model Description
The simulation framework implements a multiscale model of HBV capsid assembly that:
- Represents capsid proteins and their interactions at a coarse-grained level
- Incorporates thermodynamic and kinetic aspects of assembly
- Accounts for the dimorphic nature of HBV capsids (T=3 and T=4)
- Uses GCMC to sample different assembly states and configurations

### Key Features
- **Half-edge Data Structure**: Efficient representation of dimer subunits using a doubly connected edge list (DCEL)
- **Two Edge Types**: 
  - AB edge (AB and BA half-edges) representing AB dimer conformation
  - CD edge (CD and DC half-edges) representing CD/CC dimer conformation
- **Conformation-Dependent Interactions**: Different binding affinities for various dimer-dimer interactions
- **Elastic Energy Model**: Incorporates stretching, bending, and binding angle moduli

### Energy Model
The total energy of a capsid configuration is given by:
```
H_capsid = H_elastic + H_steric + H_bind + H_conf
```
where:
- `H_elastic`: Elastic energy from edge length, dihedral angle, and binding angle fluctuations
- `H_steric`: Excluded volume interactions
- `H_bind`: Binding energy between dimer pairs
- `H_conf`: Conformational free energy landscape

## Key Findings
Based on the research published in [Multiscale Modeling of Hepatitis B Virus Capsid Assembly and Its Dimorphism](https://pubs.acs.org/doi/10.1021/acsnano.2c02119):

1. **Assembly Pathways**: The model successfully reproduces experimental observations of HBV assembly pathways and products, including T=3 and T=4 icosahedral structures.

2. **Polymorphism Control**: Capsid polymorphism is promoted by the low HBV capsid bending modulus, with key factors being:
   - Conformational energy landscape
   - Protein-protein binding affinities
   - Coupling between protein conformational state and interactions

3. **Ionic Strength Effects**: Higher proportions of T=3 capsids assemble at higher ionic strengths, suggesting that salt concentration shifts the conformational equilibrium toward CD dimer conformations.

## Implementation Details
### Dependencies
- GSL (GNU Scientific Library)
- C++ compiler with C++11 support
- Make build system

### Building the Project

#### Windows (using MSYS2)
1. Install MSYS2 from https://www.msys2.org/
2. Open MSYS2 terminal and install dependencies:
   ```bash
   pacman -Syu
   pacman -S mingw-w64-x86_64-gsl mingw-w64-x86_64-gcc make
   ```

3. Build the project:
   ```bash
   cd src
   make clean
   make assembly
   ```

#### Linux
1. Install dependencies using your package manager:
   ```bash
   # Ubuntu/Debian
   sudo apt-get update
   sudo apt-get install libgsl-dev g++ make

   # Fedora
   sudo dnf install gsl-devel gcc-c++ make

   # Arch Linux
   sudo pacman -S gsl gcc make
   ```

2. Build the project:
   ```bash
   cd src
   make clean
   make assembly
   ```

#### macOS
1. Install dependencies using Homebrew:
   ```bash
   brew update
   brew install gsl gcc make
   ```

2. Build the project:
   ```bash
   cd src
   make clean
   make assembly
   ```

### Code Structure
- `src/geometry.cpp/hpp`: Geometric computations and half-edge data structure
- `src/montecarlo.cpp/hpp`: Monte Carlo simulation implementation
- `src/tri_tri_intersect.cpp/hpp`: Triangle intersection detection
- `src/run_assembly.cpp`: Main simulation driver

### Monte Carlo Moves
The simulation includes several MC moves:
1. Vertex displacement
2. Single edge addition/removal
3. Paired edge addition/removal
4. Simple boundary binding/unbinding
5. Wedge fusion/fission
6. Conformational switch


## Contributing
Contributions are welcome! Please feel free to submit a Pull Request.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation
If you use this code in your research, please cite:
Farzaneh Mohajerani, Botond Tyukodi, Christopher J. Schlicksup, Jodi A. Hadden-Perilla, Adam Zlotnick, and Michael F. Hagan. "Multiscale Modeling of Hepatitis B Virus Capsid Assembly and Its Dimorphism." ACS Nano 16, no. 9 (2022): 13845–13859. https://doi.org/10.1021/acsnano.2c02119.
