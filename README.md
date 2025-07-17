# Dissolution in Motion
This is a code repository accompanied with the paper "Dissolution in Motion: A Case Study on Using Computational Physics to Enhance High School Science Learning". This repository includes the code used to run a small-scale simulation of the dissolution of sodium chloride in water, as well as other accompanying tests.

## Dependencies
OpenMP is required to run this program, since it uses parallel-processing to speed up runtime.

## Usage
To edit the initialisation parameters of the program, edit the parameters under the comments labelled "Simulation parameters" and "Setup parameters". The parameters recommended to be modified are `dt`, `rc`, `iterations`, `interval`, `solvtype` (where 0 represents argon and 1 represents water), `tempinit`, `nsaltlen`, `bsizes`, `soluspace` and `solvspaces`. The output data file can also be modified by changing the value of the variable `datafile`.

To compile the program, run the following line in a Terminal with G++:
```bash
g++ -O3 -fopenmp main.cpp -o main.exe
```

The program can then be run by running the executable:

```bash
main.exe
```

The program will write data to the output file from the `datafile` variable in the following format:

```
(Solvent name) (Initial temperature) (Length of salt crystal in atoms) (Spacing between solute particles) (Spacing between solvent particles)
(Box size) (Cutoff radius) (Timestep) (Iterations) (Sampling interval)
(Particle i Iteration j Position, space-separated) (Particle i Iteration j Velocity, space-separated) (Particle i Type), ... | (Energy) | (Temperature)
...
end
```

The data file can then be read by an external program as appropriate.

# License 

This code is licensed under the MIT license. Users are welcome to adapt the code with proper attribution.
