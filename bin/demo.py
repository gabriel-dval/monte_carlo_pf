'''Test script for structuring and organising ideas

The aim of the project is to succesfully perform a Monte Carlo Replica Exchange.
The original authors did not use any boundary conditions, but it may be interesting
to try and use these once the basic program has been built. 

Method of the authors --> created a lattice at least 4 times larger than the protein
length, so when a boundary was hit, the protein was translated back into the centre.

Steps 

For a single replica, simple Monte Carlo algorithm :

1) Input : Amino acid sequence

- Turn input sequence of amino acids into Hydrophobic/Hydrophilic (very easy)
- Embed this chain into a 2Dlattice (initial position needs to be thought out)

Output : grid/matrix/graph representation of the sequence 


2) Input : Conformation, number of search iterations, search neighbourhood


'''