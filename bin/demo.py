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
import copy
import random


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
            raise ValueError('Not the coordinate format')


    def is_above(self, aa2):
        '''Check if amino acid is above another'''
        delta = self.coords - aa2.coords
        if delta == np.array([1, 0]):
            return True
        else: return False

    def is_below(self, aa2):
        '''Check if amino acid is below another'''
        delta = self.coords - aa2.coords
        if delta == np.array([-1, 0]):
            return True
        else: return False

    
    def is_left_of(self, aa2):
        '''Check if amino acid is left of another'''
        delta = self.coords - aa2.coords
        if delta == np.array([0, -1]):
            return True
        else: return False

    
    def is_right_of(self, aa2):
        '''Check is amino acid is right of another'''
        delta = self.coords - aa2.coords
        if delta == np.array([0, 1]):
            return True
        else: return False


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

    
    def show_neighbour_map(self):
        '''Display the current neighbour map'''
        print(self.neighbour_map)


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
    def __init__(self, label, protein: Protein, lattice: Lattice2D):
        '''Automatically initialise the first conformation'''
        self.label = label
        self.protein = protein
        self.lattice = lattice
        self.position_manager = np.zeros((lattice.length, lattice.width), dtype = int)
        self.initialise_horizontal_conformation()

    
    def initialise_horizontal_conformation(self):
        '''Function to initialise the conformation'''
        ref_x = int((self.lattice.length / 2) - (self.protein.calc_length() / 2))
        ref_y = int(self.lattice.width / 2)
        i = 1
        for aa in self.protein.sequence:
            aa.coords = np.array([ref_y, ref_x])
            self.position_manager[ref_y, ref_x] = i
            ref_x += 1
            i += 1
        print(self.position_manager)


    def get_protein_sequence(self):
        return self.protein.sequence
    

    def valid_position(self, pos):
        '''Check if a position is within the lattice and not occupied'''
        y, x = pos
        if 0 <= x < self.lattice.length and 0 <= y < self.lattice.width:
            return self.position_manager[y, x] == 0
        return False
    

    def get_neighbours(self, pos: int):
        '''Returns neighboring positions around a given residue position
        This is based around python index not actual residue position.
        '''
        try:
            aa = self.get_protein_sequence()[pos]
        except:
            raise IndexError('Given position out of range')

        if isinstance(aa, HPAminoAcid):
            neighbors = [
                aa.coords + np.array([0, 1]),  # right
                aa.coords + np.array([0, -1]), # left
                aa.coords + np.array([1, 0]),  # down
                aa.coords + np.array([-1, 0])  # up
            ]
            return [n for n in neighbors if self.valid_position(n)]
        else:
            raise ValueError('Argument is not of class HPAminoAcid')

    
    def  calculate_energy(self):
        '''Function to calculate the energy of the current conformation
        
        Iterate through atoms of the protein and check its coordinates
        '''
        energy = 0
        for i, aa1 in enumerate(self.protein.sequence[:-1]):
            for aa2 in self.protein.sequence[i + 1:]:
                delta = np.sort(np.absolute(aa1.coords - aa2.coords))
                if aa1.type == 'P' or aa2.type == 'P':
                    #print('skip P')
                    continue
                elif int(aa2.name[1:]) - int(aa1.name[1:]) < 2:
                    #print('skip neighbours')
                    continue
                elif np.array_equal(delta, np.array([0, 1])) :
                    energy -= 1
                else:
                    continue
                    #print('too far')
        
        return energy


class MonteCarlo:
    '''Run a single series of MonteCarlo permutations on a conformation. 

    Attributes
    ---
    conformation : initial conformation to optimise
    steps : number of optimisation steps
    move_types : which set of moves to use (pull, vhsc or all)

    Methods
    ---
    Implementation of each move is in this class
    '''
    def __init__(self, conformation: Conformation, steps: int, temperature: float):
        self.conformation = conformation
        self.steps = steps
        self.temperature = temperature
    

    def choose_rand_aa(self):
        prot_length = self.conformation.protein.calc_length()
        k = np.random.randint(0, prot_length)
        return k

    
    def try_end_move(self, conf, k):
        '''With the given conformation, try and implement an end
        move at position k.

        Arguments
        ---
        conformation : the input conformation
        k : the amino acid position

        Returns
        ---
        boolean - whether the move was succesful or not
        '''
        seq = conf.get_protein_sequence()    # Get sequence of the protein

        if k == 0 or k == len(seq) - 1:
            if k == 0:
                free_spots = conf.get_neighbours(k + 1)
            else: 
                free_spots = conf.get_neighbours(len(seq) - 2)

            if len(free_spots) != 0:
                new_y, new_x = random.choice(free_spots)  # Pull coords of new spot
                old_y, old_x = seq[k].coords              # Pull coords of old spot

                # Swap values of both spots
                conf.position_manager[new_y, new_x] = conf.position_manager[old_y, old_x]
                conf.position_manager[old_y, old_x] = 0 

                # Update the aa object of that conformation  
                seq[k].coords = np.array([new_y, new_x])

                # Report results
                print('End move successful')
                return True
            else:
                print('End move skipped')
                return False
        
        else:
            print('End move skipped')
            return False
        

    def try_corner_move(self, conf, k):
        '''Test and perform a corner move if possible
        
        Arguments
        ---
        conformation : the input conformation
        k : the amino acid position

        Returns
        ---
        boolean - whether the move was succesful or not
        '''
        seq = conf.get_protein_sequence()    # Get sequence of the protein

        if k != 0 and k != len(seq) - 1:

            # Calculate vectors with neighbours
            past_res = seq[k].coords - seq[k-1].coords
            future_res = seq[k].coords - seq[k+1].coords
            
            # If both vectors are orthogonal
            if np.dot(past_res, future_res) == 0:

                # Check spot is available
                mov = past_res + future_res
                old_y, old_x = seq[k].coords
                if conf.position_manager[old_y - mov[0], old_x - mov[1]] == 0:

                    #Update positions
                    new_y, new_x = old_y - mov[0], old_x - mov[1]
                    conf.position_manager[new_y, new_x] = conf.position_manager[old_y, old_x]
                    conf.position_manager[old_y, old_x] = 0 

                    #Update aa objects
                    seq[k].coords = np.array([new_y, new_x])

                    #Report
                    print('Corner move succesful') 
                    return True

                else:
                    print('Corner move skipped')
                    return False
            else:
                print('Corner move skipped')
                return False
        else:
            print('Corner move skipped')
            return False


    def try_crankshaft_move(self, conf, k):
        '''Try and implement a crankshaft move
        
        Arguments
        ---
        conformation : the input conformation
        k : the amino acid position

        Returns
        ---
        boolean - whether the move was succesful or not
        '''
        seq = conf.get_protein_sequence()    # Get sequence of the protein

        if k > 1 and k < len(seq) - 2:
            
            # Calculate useful vectors
            before_last_res = seq[k-1].coords - seq[k-2].coords
            last_res = seq[k].coords - seq[k-1].coords
            next_res = seq[k].coords - seq[k+1].coords
            next_next_res = seq[k+1].coords - seq[k+2].coords

            # Case 1
            if np.dot(last_res, next_res) == 0 and np.dot(next_res, next_next_res) == 0:
                
                # Check spots are free
                old_y, old_x = seq[k].coords
                next_old_y, next_old_x = seq[k+1].coords
                if conf.position_manager[seq[k-1].coords[0] - last_res[0], 
                                        seq[k-1].coords[1] - last_res[1]] == 0\
                and conf.position_manager[seq[k+2].coords[0] - last_res[0], 
                                        seq[k+2].coords[1] - last_res[1]] == 0:
                    
                    #Update positions
                    new_y, new_x = seq[k-1].coords[0] - last_res[0], seq[k-1].coords[1] - last_res[1]
                    next_new_y, next_new_x = seq[k+2].coords[0] - last_res[0], seq[k+2].coords[1] - last_res[1]
                    conf.position_manager[new_y, new_x] = conf.position_manager[old_y, old_x]
                    conf.position_manager[next_new_y, next_new_x] = conf.position_manager[next_old_y, next_old_x]
                    
                    conf.position_manager[old_y, old_x] = 0
                    conf.position_manager[next_old_y, next_old_x] = 0

                    #Update aa objects
                    seq[k].coords = np.array([new_y, new_x])
                    seq[k+1].coords = np.array([next_new_y, next_new_x])

                    #Report
                    print('Crankshaft move succesful') 
                    return True
                
                else:
                    print('Crankshaft move skipped')
                    return False
  
            
            # Case 2
            if np.dot(before_last_res, last_res) == 0 and np.dot(last_res, next_res) == 0:
                
                # Check spots are free
                old_y, old_x = seq[k].coords
                last_old_y, last_old_x = seq[k-1].coords
                if conf.position_manager[seq[k-2].coords[0] - before_last_res[0], 
                                         seq[k-2].coords[1] - before_last_res[1]] == 0\
                and conf.position_manager[seq[k+1].coords[0] - before_last_res[0], 
                                         seq[k+1].coords[1] - before_last_res[1]] == 0:
                    
                    #Update positions
                    new_y, new_x = seq[k-2].coords[0] - before_last_res[0], seq[k-2].coords[1] - before_last_res[1]
                    last_new_y, last_new_x = seq[k+1].coords[0] - before_last_res[0], seq[k+1].coords[1] - before_last_res[1]
                    conf.position_manager[new_y, new_x] = conf.position_manager[old_y, old_x]
                    conf.position_manager[last_new_y, last_new_x] = conf.position_manager[last_old_y, last_old_x]
                    
                    conf.position_manager[old_y, old_x] = 0
                    conf.position_manager[last_old_y, last_old_x] = 0

                    #Update aa objects
                    seq[k].coords = np.array([new_y, new_x])
                    seq[k-1].coords = np.array([last_new_y, last_new_x])

                    #Report
                    print('Crankshaft move succesful') 
                    return True
                
                else:
                    print('Crankshaft move skipped')
                    return False
            else:
                print('Crankshaft move skipped')
                return False
        else:
            print('Crankshaft move skipped')
            return False




    def choose_move(self, conf, k, move_neighbourhood):
        '''Function which chooses a random move to perform
        
        Arguments
        ---
        conformation : the input conformation
        k : the amino acid position
        '''

        if move_neighbourhood == 'ALL':
            options = ['end', 'corner', 'crankshaft']
        elif move_neighbourhood == 'VSHD':
            options = ['end', 'corner', 'crankshaft']
        elif move_neighbourhood == 'PULL':
            options = ['corner', ]
        else:
            raise ValueError('Not a valid option')
        
        choice = np.random.choice(options, p = [0.4, 0.4, 0.2])

        if choice == 'end':
            print('End move chosen')
            self.try_end_move(conf, k)
        if choice == 'corner':
            print('Corner move chosen')
            self.try_corner_move(conf, k)
        if choice == 'crankshaft':
            print('Crankshaft move chosen')
            self.try_crankshaft_move(conf, k)

        return choice

    
    def run_sim(self):
        '''General function implementing the MC sim
        
        Arguments
        ---
        none

        Returns
        ---
        nothing - modifies the conformation instance in place.
        '''
        current_conformation = copy.deepcopy(self.conformation)
        sequence = self.conformation.get_protein_sequence()

        for s in range(self.steps):
            test_conformation = copy.deepcopy(current_conformation)

            # Choose an amino acid at random
            k = self.choose_rand_aa() 
            print(k)

            # Choose a move at random (for now only one move available)
            move = self.choose_move(test_conformation, k, 'ALL')

            # Calculate energy of current and test
            test = test_conformation.calculate_energy()
            current = current_conformation.calculate_energy()
            delta_e = test - current

            if delta_e <= 0:
                current_conformation = copy.deepcopy(test_conformation)
            else:
                q = np.random.uniform(0, 1)
                if q > np.exp((-delta_e/self.temperature)):
                    current_conformation = copy.deepcopy(test_conformation)
        
        print(current_conformation.position_manager) 
        return current_conformation

        




        





    

if __name__ == "__main__":

    aa1 = HPAminoAcid('H1')
    aa2 = HPAminoAcid('H2')
    aa3 = HPAminoAcid('P3')
    aa4 = HPAminoAcid('H4')
    aa5 = HPAminoAcid('H5')
    aa6 = HPAminoAcid('P6')
    aa7 = HPAminoAcid('P7')
    aa8 = HPAminoAcid('H8')
    aa9 = HPAminoAcid('H8')
    aa10 = HPAminoAcid('H8')
    aa11 = HPAminoAcid('H8')

    prot1 = Protein('test')
    prot1.add_aa(aa1)
    prot1.add_aa(aa2)
    prot1.add_aa(aa3)
    prot1.add_aa(aa4)
    prot1.add_aa(aa5)
    prot1.add_aa(aa6)
    prot1.add_aa(aa7)
    prot1.add_aa(aa8)
    prot1.add_aa(aa9)
    prot1.add_aa(aa10)
    prot1.add_aa(aa11)

    prot2 = Protein('test2')
    prot2.add_aa(aa4)
    prot2.add_aa(aa3)
    prot2.add_aa(aa2)
    prot2.build_neighbour_dict()

    l1 = Lattice2D(24, 24)

    conf1 = Conformation('C1', prot1, l1)
    conf1.calculate_energy()
    
    # Testing simple Monte Carlo
    mc1 = MonteCarlo(conf1, 40, 60)
    mc1.run_sim()



    

    