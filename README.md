# Dissolution in Motion
This is a code repository accompanied with the paper "Dissolution in Motion: A Case Study on Using Computational Physics to Enhance High School Science Learning". This repository includes the code used to run a small-scale simulation of the dissolution of sodium chloride in water, as well as other accompanying tests.

## Dependencies
OpenMP is required to run this program, since it uses parallel-processing to speed up runtime.

## Usage
To edit the initialisation parameters of the program, edit the parameters under the comments labelled "Simulation parameters" and "Setup parameters". The parameters recommended to be modified are `dt`, `rc`, `iterations`, `interval`, `solvtype` (where 0 represents argon and 1 represents water), `tempinit`, `nsaltlen`, `bsizes`, `soluspace` and `solvspaces`.

To compile the program, run the following line in a Terminal with G++:
```bash
g++ -O3 -fopenmp main.cpp -o main.exe
```

The program can then be run by running the executable:

```bash
main.exe
```
