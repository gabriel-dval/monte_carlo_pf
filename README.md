# Predicting protein folding using a replica exchange Monte Carlo method

Program for predicting protein folding using a Monte Carlo algorithm.


## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Prerequisite packages

To run the following script, the following modules are necessary. An environment
.yml file (remc_env.yml) has also been provided. This program was written using 
Python v3.12.5

- [numpy](https://numpy.org/), v2.1
- [matplotlib](https://matplotlib.org/stable/users/explain/quick_start.html), v3.9.2
- [tqdm](https://pypi.org/project/tqdm/), v4.66.5

### Installing

If you wish to use the provided environment, you can install it with conda 
using this command.

macOS command line

    conda env create -f remc_env.yml

You can equally create a new environment with the three packages if this is 
simpler.

    conda create -n myenv python=3.12.5 numpy=2.1 matplotlib=3.9.2 tqdm=4.66.5


## Running the program

The program is run from the command line like any python script. Below is an 
example with a sequence of length 20 called 'test'. Although most parameters have
default options, please make sure to input an optimal energy even if it is not 
known. A cutoff time for each run can equally be passed as argument.

Example usage

    python remc.py HPHPPHHPHPPHPHPH test -5

IMPORTANT - make sure the components.py file is in the same directory as remc.py
so all necessary modules can be imported. 

The program outputs a results_log.txt file, containing results for each iteration 
of the simulation, and a FIGURES directory with the final conformation for each run.

### Launch an example proteins

Launch a simulation of protein S1-1

    python remc.py HPHPPHHPHPPHPHHPPHPH S1 -9


## Authors

  - **GabrielDuval** - *Provided README Template* -
    [gabriel-dval](https://github.com/gabriel-dval)


## References

(1) Thachuk, C., Shmygelska, A. & Hoos, H.H. A replica exchange Monte Carlo algorithm for protein folding in the HP model. BMC Bioinformatics 8, 342 (2007). https://doi.org/10.1186/1471-2105-8-342

