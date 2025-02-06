# Chromatin 3D Structure and Gene Expression Simulation

This project reconstructs the 3D structure of chromatin from HiC contact maps using a Bead-Spring Polymer model. It also predicts gene expression levels by analyzing enhancer-promoter (E-P) interactions through a kinetic model.

## Features
- Simulates chromatin 3D conformations using a Bead-Spring Polymer model.
- Computes gene expression based on enhancer-promoter interactions.
- Compatible with LAMMPS for simulations and MATLAB for kinetic modeling.

## Getting Started

### Prerequisites
- [LAMMPS](https://lammps.org) install lammps for md simulation.
- [MATLAB R2022a](https://www.mathworks.com/products/matlab.html) for the kinetic model simulation.



### Workflow

### 1. Coarse-Grain Bead-Spring Polymer Model

The Bead-Spring Polymer model simulates the chromatin structure for the 6MB Sox9-Kcnj2 loci of the mESC wildtype cell line. This is achieved using the LAMMPS simulation package.

### Key Files:

in.run: Input file for LAMMPS simulations.

poly.data: LAMMPS data file.

### Running the Bead-Spring Polymer Model
1. Navigate to the project folder.
2. Execute the simulation using LAMMPS through command line:
   lmp_serial -in in.run > log &

This command initiates the simulation and redirects the output to the log file for monitoring.

### Analyses files:

distance.f90 to calculate 3D distance and cont_prob.f90 to generate the contact map.

### 2. Kinetic Model

The kinetic model calculates gene expression by calculating the cluster size of promoter by modeling the rates of binding and unbinding interactions between enhancers and promoters, as derived from the 3D chromatin structure.

### Key Files:

Kinetic_model_sox9_wildtype: MATLAB script implementing the kinetic model.
rates_binding.csv and rates_unbinding.csv: Enhancer binding and unbinding rates with promoter

### Running the Kinetic Model:

This script is developed in MATLAB R2022a. It can be executed either through the MATLAB GUI or via the command line prompt.

Command Line Execution:

Open the terminal.

Navigate to the directory containing the script.

Run the script using the following command:

matlab -nodisplay -nosplash -r "run('Kinetic_model_sox9_wildtype'); exit;"

### Notes

Ensure LAMMPS and MATLAB are installed and properly configured on your system before running the simulations and kinetic model.

For any issues or inquiries, refer to the documentation or contact the project administrator.
