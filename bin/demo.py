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


# Import of necessary modules

import numpy as np
import math as m


# Experimenting with first class architectures

class HPAminoAcid:
    '''Base class recording the type of amino acid : Hydrophobic (H) or Polar (P)
    
    Attributes 
    ---
    name : corresponds to the identity of the HP amino acid (H2, P9 etc)
    type : type of the amino acid
    coords : np array, containing coordinates of this AA - the initial
            values of this are always empty

    Methods
    ---
    get_coords : access the coordinates of the amino acid
    set_coords : set the coordinates of the amino acid
    '''
    def __init__(self, name):
        self.name = name
        self.type = name[0]
        self.coords = np.array([0, 0])

    @property
    def coords(self):
        return self._coordinates
    
    @coords.setter
    def coords(self, new_coords):
        if new_coords.shape == (2,):
            self._coordinates = new_coords
        else:
            raise ValueError


class Protein:
    '''Represents a protein - this is constructed from multiple HPAminoAcid objects

    Attributes
    ---
    name : name of the protein
    sequence : HP sequence, also constructed as the protein gets longer
    neighbour_map : dictionnary of the contacts between proteins

    Methods
    ---
    add_aa : add a HPAminoAcid to the protein
    calc_length : calculate the length of the protein
    get_sequence : print the composition of the protein
    build_neighbour_dict : fill in the neighbour map - can only be done 
                           if the neighbour map is empty
    '''
    def __init__(self, name):
        self.name = name
        self.sequence = []
        self.neighbour_map = {}
    

    def add_aa(self, aa):
        if isinstance(aa, HPAminoAcid):
            self.sequence.append(aa)


    def calc_length(self):
        return len(self.sequence)
    

    def get_sequence(self):
        return "".join(aa.type for aa in self.sequence)
    

    def build_neighbour_dict(self):
        '''Remember this method will only run if the neighbour map
        if empty - it will tell you to run it again otherwise'''
        if len(self.neighbour_map) == 0:
            for i, aa in enumerate(self.sequence[:-1]):
                self.neighbour_map[aa.name] = {'neighbour' : self.sequence[i+1].name,
                                               'coords' : self.sequence[i+1].coords}
        else:
            print('A neighbour map of this protein has already been built.\nRun delete_neighbour_map and try again')

    
    def delete_neighbour_map(self):
        '''Delete the current neighbour map'''
        self.neighbour_map.clear()
        print(f'Neighbour map deleted successfully for {self.name}')
        

class Lattice2D:
    '''Class representing the lattice on which the protein will be 
    embedded.
    
    Attributes
    ---
    length : how long is the lattice
    width : how wide is the lattice


    Methods
    ---
    border_control : checks that the protein is not too close to the border
    '''
    def __init__(self, length, width):
        self.length = length
        self.width = width
        self.position_manager = np.zeros((length, width))

    
    def border_control(self):
        '''Checks the values present in the lattice - if there is a 1 present
        near the border, the whole structure needs to be re-translated'''
        print('In development')


class Conformation:
    '''This class combines the protein and the lattice into one - it sets their original
    position and also implements all the methods to update the position of the conformation

    Attributes
    ---
    label : this is the id of the conformation - will be used when RE is implemented
    position_matrix
    '''
    def __init__(self, label, protein, lattice):
        self.label = label
        if isinstance(protein, Protein):
            self.protein = protein
        else:
            raise ValueError('Protein is not of class Protein')
        if isinstance(lattice, Lattice2D):
            self.lattice = lattice
        else:
            raise ValueError('Lattice is not of class Lattice2D')
        self.position_manager = np.zeros((lattice.length, lattice.width))

    
    def initialise_horizontal_conformation(self):
        '''Function to initialise the conformation'''
        ref_x = (self.lattice.length / 2) - (self.protein.calc_length() / 2)
        ref_y = self.lattice.width / 2
        for aa in zip(self.protein, zip.protein.neighbour_map):
            aa.coords = np.array([ref_x, ref_y])
            self.position_manager[ref_x, ref_y] = 1
            ref_x += 1



    

if __name__ == "__main__":
    print('test')
    aa1 = HPAminoAcid('H1')
    aa2 = HPAminoAcid('H2')
    aa3 = HPAminoAcid('P3')
    aa4 = HPAminoAcid('H4')

    prot1 = Protein('test')
    prot1.add_aa(aa1)
    prot1.add_aa(aa2)
    prot1.add_aa(aa3)
    prot1.add_aa(aa4)

    print(prot1.neighbour_map)
    length = prot1.build_neighbour_dict()
    print(prot1.neighbour_map)



    

    