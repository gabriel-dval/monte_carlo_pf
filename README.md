# Predicting protein folding using a replica exchange Monte Carlo method

Program for predicting protein folding using a Monte Carlo algorithm.


## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes. See deployment
for notes on deploying the project on a live system.

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


## Running the tests

The program is run from the command line like any python script. Below is an 
example with a sequence of length 20 called 'test'. Although most parameters have
default options, please make sure to input an optimal energy even if it is not 
known. A cutoff time for each run can equally be passed as argument.

Example usage

    python remc.py HPHPPHHPHPPHPHHPPHPH test -9


### Sample Tests

Explain what these tests test and why

    Give an example

### Style test

Checks if the best practices and the right coding style has been used.

    Give an example

## Deployment

Add additional notes to deploy this on a live system

## Built With

  - [Contributor Covenant](https://www.contributor-covenant.org/) - Used
    for the Code of Conduct
  - [Creative Commons](https://creativecommons.org/) - Used to choose
    the license

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code
of conduct, and the process for submitting pull requests to us.

## Versioning

We use [Semantic Versioning](http://semver.org/) for versioning. For the versions
available, see the [tags on this
repository](https://github.com/PurpleBooth/a-good-readme-template/tags).

## Authors

  - **Billie Thompson** - *Provided README Template* -
    [PurpleBooth](https://github.com/PurpleBooth)

See also the list of
[contributors](https://github.com/PurpleBooth/a-good-readme-template/contributors)
who participated in this project.

## License

This project is licensed under the [CC0 1.0 Universal](LICENSE.md)
Creative Commons License - see the [LICENSE.md](LICENSE.md) file for
details

## Acknowledgments

  - Hat tip to anyone whose code is used
  - Inspiration
  - etc

