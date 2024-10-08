'''Class components of the REMC model
------------------------------------------------------

This file contains : 

- Necessary modules for running all classes correctly
- All classes used for the REMC and their associated methods
- All helper functions to run the programs
'''

# Import of necessary modules

import numpy as np
import math as m
import copy
import random
import time
import multiprocessing as mp
from tqdm import tqdm
import matplotlib.pyplot as plt


# REMC Classes

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
    is_above : is this amino acid above another ?
    is_below : is this amino acid below another ?
    is_right_of : is this amino acid right of another ?
    is_left_of : is this amino acid left of another ?
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
    show_neighbour_map : display neighbours of each amino acid
    delete_neighbour_map : delete the neighbour map
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
    embedded. A little redundant but more user-friendly
    
    Attributes
    ---
    length : how long is the lattice
    width : how wide is the lattice


    Methods
    ---
    none
    '''
    def __init__(self, length, width):
        self.length = length
        self.width = width



class Conformation:
    '''This class combines the protein and the lattice into one - it sets their original
    position and also implements all the methods to update the position of the conformation

    Attributes
    ---
    label : this is the id of the conformation - will be used when RE is implemented
    protein : Protein instance
    lattice : Lattice2D instanec
    position_manager : numpy array keeping track of amino acid positions

    Methods
    ---
    initialise_x_conformation : methods to initialise first position of the conformation
    get_protein_sequence : return the amino acids of the protein
    valid_position : check if a position is free
    check_valid_conformation : is the adopted conformation possible ?
    get_neighbours : return free adjacent neighbours
    get_diagonal_neighbbours : return free diagonal neighbours
    calculate_energy : what is the current energy of the conformation ?
    check_border : is the protein about to run off the lattice ?
    view_conformation : strip all the zero rows and columns to view the conformation only
    '''
    def __init__(self, label: str, protein: Protein, lattice: Lattice2D):
        '''Automatically initialise the first conformation'''
        self.label = label
        self.protein = protein
        self.lattice = lattice
        self.position_manager = np.zeros((lattice.length, lattice.width), dtype = int)
        self.initialise_horizontal_conformation()

    
    def initialise_horizontal_conformation(self):
        '''Function to initialise the conformation horizontally'''
        ref_x = int((self.lattice.length / 2) - (self.protein.calc_length() / 2))
        ref_y = int(self.lattice.width / 2)
        i = 1
        for aa in self.protein.sequence:
            aa.coords = np.array([ref_y, ref_x])
            self.position_manager[ref_y, ref_x] = i
            ref_x += 1
            i += 1

    
    def initialise_vertical_conformation(self):
        '''Function to initialise the conformation vertically'''
        ref_x = int(self.lattice.length / 2)
        ref_y = int((self.lattice.width / 2) - (self.protein.calc_length() / 2))
        i = 1
        for aa in self.protein.sequence:
            aa.coords = np.array([ref_y, ref_x])
            self.position_manager[ref_y, ref_x] = i
            ref_y += 1
            i += 1
            
    
    def initialise_diagonal_conformation(self):
        '''Function to initialise the conformation diagonally'''
        ref_x = int((self.lattice.length / 2) - (self.protein.calc_length() / 2))
        ref_y = int((self.lattice.width / 2) - (self.protein.calc_length() / 2))
        i = 1
        for aa in self.protein.sequence:
            aa.coords = np.array([ref_y, ref_x])
            self.position_manager[ref_y, ref_x] = i
            if i % 2 == 0:
                ref_y += 1
            else:
                ref_x += 1
            i += 1

    
    def choose_conformation(self):
        '''Choose initial conformation'''
        c = random.choice([1,2,3])
        if c == 1:
            self.initialise_diagonal_conformation()
        if c == 2:
            self.initialise_horizontal_conformation()
        if c == 3:
            self.initialise_vertical_conformation()


    def get_protein_sequence(self):
        return self.protein.sequence
    

    def valid_position(self, pos):
        '''Check if a position is within the lattice and not occupied'''
        y, x = pos
        if 0 <= x < self.lattice.length and 0 <= y < self.lattice.width:
            return self.position_manager[y, x] == 0
        return False
    

    def check_valid_conformation(self):
        '''Function to make sure that the protein chain is in a valid conformation
        on the lattice. 
        '''
        seq = self.get_protein_sequence()

        # Check that every aa is a distance of one away from each other
        for i in range(len(seq) - 1):
            aa1_x, aa1_y = seq[i].coords
            aa2_x, aa2_y = seq[i + 1].coords

            # Calculate distance
            dist = np.sqrt((aa2_x - aa1_x)**2 + (aa2_y - aa1_y)**2)

            if dist != 1:
                return False
            
        return True


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


    def get_diagonal_neighbours(self, pos: int):
        '''Returns the diagonal neighbours of an amino acid'''
        try:
            aa = self.get_protein_sequence()[pos]
        except:
            raise IndexError('Given position out of range')
        
        if isinstance(aa, HPAminoAcid):
            neighbors = [
                aa.coords + np.array([1, 1]),  # right up
                aa.coords + np.array([-1, 1]), # right down
                aa.coords + np.array([1, -1]),  # left up 
                aa.coords + np.array([-1, -1])  # left down
            ]
            return [n for n in neighbors if self.valid_position(n)]
        else:
            raise ValueError('Argument is not of class HPAminoAcid')

    
    def calculate_energy(self):
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
                elif np.array_equal(delta, np.array([0, 1])):
                    energy -= 1
                else:
                    continue
                    #print('too far')
        
        return energy


    def check_border(self):
        '''Function to verify structure is not too close to the edge of the
        lattice - if this is the case, the whole structure is translated back
        into the centre of the lattice.

        WARNING : if the lattice size is too small, values may be overwritten on
        the position manager - make sure lattice is least twice (but ideally three)
        times the size of the input protein.
        '''
        # Check borders of the lattice to see if any aa is there
        outlier = False

        # Top
        if np.sum(self.position_manager[0,:]) != 0:
            outlier = True
        # Bottom
        if np.sum(self.position_manager[-1,:]) != 0:
            outlier = True
        # Right
        if np.sum(self.position_manager[:,0]) != 0:
            outlier = True
        # Left
        if np.sum(self.position_manager[:,-1]) != 0:
            outlier = True

        # Translate structure if true
        if outlier:
            
            # Calculate vector between central aa and lattice centre
            seq = self.get_protein_sequence()
            centre_aa = seq[int(len(seq)/2)].coords
            centre_lattice = np.array([int(self.lattice.width/2), int(self.lattice.length/2)])

            dif_vec = centre_lattice - centre_aa

            # Move all coordinates by the opposite of this vector
            for aa in seq:
                aa_y, aa_x = aa.coords
                self.position_manager[aa_y + dif_vec[0], aa_x + dif_vec[1]] = self.position_manager[aa_y, aa_x]
                self.position_manager[aa_y, aa_x] = 0

                aa.coords = aa.coords + dif_vec
            
            print('Structure translated to centre')
            time.sleep(5)


    def view_conformation(self):
        '''Function to view the position manager of the conformation.

        Returns a window of the conformation
        '''
        window = copy.deepcopy(self.position_manager)

        # Top
        while np.sum(window[0,:]) == 0:
            window = window[1:,:]
        # Bottom
        while np.sum(window[-1,:]) == 0:
            window = window[:-1,:]
        # Right
        while np.sum(window[:,-1]) == 0:
            window = window[:,:-1]
        # Left
        while np.sum(window[:,1]) == 0:
            window = window[:, 1:]

        return window



class MonteCarlo:
    '''Run a single series of MonteCarlo permutations on a conformation. 

    Attributes
    ---
    conformation : conformation to optimise
    steps : number of optimisation steps
    temperature : temperature of simulation

    Methods
    ---
    choose_rand_aa : select a random position in the sequence
    try_x_move : try a particular MC move - skip if impossible
    choose_move, choose_available_move : different functions to select a random move
    run_sim : launch a single MC sim on the input conformation
    '''
    def __init__(self, conformation: Conformation, steps: int, temperature: float):
        self.conformation = conformation
        self.steps = steps
        self.temperature = temperature
    

    def choose_rand_aa(self):
        prot_length = self.conformation.protein.calc_length()
        
        # Create rand number generator
        rng = np.random.default_rng()
        k = rng.integers(low = 0, high = prot_length)
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
                #print('End move successful')
                return True
            else:
                #print('End move skipped')
                return False
        
        else:
            #print('End move skipped')
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
                    #print('Corner move succesful') 
                    return True

                else:
                    #print('Corner move skipped')
                    return False
            else:
                #print('Corner move skipped')
                return False
        else:
            #print('Corner move skipped')
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
                    #print('Crankshaft move succesful') 
                    return True
                
                else:
                    #print('Crankshaft move skipped')
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
                    #print('Crankshaft move succesful') 
                    return True
                
                else:
                    #print('Crankshaft move skipped')
                    return False
            else:
                #print('Crankshaft move skipped')
                return False
        else:
            #print('Crankshaft move skipped')
            return False


    def try_pull_move(self, conf, k):
        '''Perform a pull move is possible - this can have many different
        effects.
        
        Arguments
        ---
        conformation : the input conformation
        k : the amino acid position

        Returns
        ---
        boolean - whether the move was succesful or not
        '''
        seq = conf.get_protein_sequence()    # Get sequence of the protein

        # Find L

        if k < len(seq) - 1:
            
            current = conf.get_diagonal_neighbours(k)
            next = conf.get_neighbours(k + 1)

            # Convert each numpy array into tuples to make them hashable
            set1 = {tuple(arr) for arr in current}
            set2 = {tuple(arr) for arr in next}

            if len(set1.intersection(set2)) > 0:
                
                # Find the common tuples between the two sets
                common_tuples = set1.intersection(set2)

                # Convert the common tuples back to numpy arrays
                choices = [np.array(tup) for tup in common_tuples]

                l = random.choice(choices)

                # Find C
                vec_with_k = seq[k+1].coords - seq[k].coords
                vec_with_l = seq[k+1].coords - l

                vec_sum = vec_with_k + vec_with_l
                c = np.array([seq[k+1].coords[0] - vec_sum[0], seq[k+1].coords[1] - vec_sum[1]])

                # Case 1 : C is on k - 1
                if np.array_equal(c, seq[k-1].coords):
                    self.try_corner_move(conf, k)
                    return True
                
                # Case 2 : C is on an empty square 
                elif conf.position_manager[c[0], c[1]] == 0:
                    # Make sure not at beginning of chain
                    if k > 0:

                        # Create copy of old positions
                        old_positions = copy.deepcopy(seq)

                        # Save previous coordinates of k and k - 1
                        tmp1 = seq[k].coords
                        tmp2 = seq[k-1].coords

                        # Change coordinates of k and k-1 on position manager
                        conf.position_manager[l[0], l[1]] = conf.position_manager[tmp1[0], tmp1[1]] 
                        conf.position_manager[tmp1[0], tmp1[1]] = 0

                        conf.position_manager[c[0], c[1]] = conf.position_manager[tmp2[0], tmp2[1]] 
                        conf.position_manager[tmp2[0], tmp2[1]] = 0
                        
                        # Update aa objects
                        seq[k].coords = l
                        seq[k-1].coords = c

                        # Implement pull move
                        ref = k
                        conf_status = conf.check_valid_conformation()

                        # While index not out of range and the conformation is invalid
                        while ref - 2 >= 0 and conf_status == False:
                                
                            # Get old ref coords
                            k_old = old_positions[ref].coords

                            # Move k - 2 to k
                            k2_y, k2_x = seq[ref-2].coords
                            conf.position_manager[k_old[0], k_old[1]] = conf.position_manager[k2_y,k2_x] 
                            conf.position_manager[k2_y, k2_x] = 0

                            # Update aa object
                            seq[ref-2].coords = k_old

                            # Check conformation is valid
                            conf_status = conf.check_valid_conformation()

                            # Substract from ref
                            ref -= 1
                        
                        #print('Pull move succesful')
                        return True
                    else:
                        #print('Pull move skipped')
                        return False
                else:
                    #print('Pull move skipped')
                    return False
            else:
                #print('Pull move skipped')
                return False
        else:
            #print('Pull move skipped')
            return False


    def choose_move(self, conf, k, move_neighbourhood):
        '''Function which chooses a random move to perform
        
        Arguments
        ---
        conformation : the input conformation
        k : the amino acid position
        '''

        if move_neighbourhood == 'ALL':
            options = ['end', 'corner', 'crankshaft', 'pull']
        elif move_neighbourhood == 'VSHD':
            options = ['end', 'corner', 'crankshaft', 'corner']
        elif move_neighbourhood == 'PULL':
            options = ['pull', 'pull', 'pull', 'pull']
        else:
            raise ValueError('Not a valid option')
        
        choice = np.random.choice(options, p = [0.2, 0.2, 0.2, 0.4])

        if choice == 'end':
            self.try_end_move(conf, k)
        if choice == 'corner':
            self.try_corner_move(conf, k)
        if choice == 'crankshaft':
            self.try_crankshaft_move(conf, k)
        if choice == 'pull':
            self.try_pull_move(conf, k)

        return choice


    def choose_available_move(self, conf, k, move_neighbourhood):
        '''This time if a move is not possible then another move is tried.
        The order of moves presented is different every iteration
        '''

        if move_neighbourhood == 'ALL':
            options = ['end', 'corner', 'crankshaft', 'pull']
        
        choices = np.random.choice(options, size = 4, replace = False,
                                  p = [0.2, 0.2, 0.2, 0.4])

        for choice in choices:

            if choice == 'end':
                end = self.try_end_move(conf, k)
                if end:
                    break
            if choice == 'corner':
                corner = self.try_corner_move(conf, k)
                if corner:
                    break
            if choice == 'crankshaft':
                crankshaft = self.try_crankshaft_move(conf, k)
                if crankshaft:
                    break
            if choice == 'pull':
                pull = self.try_pull_move(conf, k)
                if pull:
                    break

    
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

        for s in tqdm(range(self.steps), 
                      desc = f'Conformation {self.conformation.label}'):
            
            test_conformation = copy.deepcopy(current_conformation)

            # Choose an amino acid at random
            k = self.choose_rand_aa() 

            # Choose a move at random
            self.choose_available_move(test_conformation, k, 'ALL')

            # Calculate energy of current and test
            test = test_conformation.calculate_energy()
            current = current_conformation.calculate_energy()
            delta_e = test - current

            if delta_e <= 0:
                current_conformation = copy.deepcopy(test_conformation)
            else:
                rng = np.random.default_rng()
                q = rng.random()
                if q > np.exp((-delta_e/self.temperature)):
                    current_conformation = copy.deepcopy(test_conformation)
            
            # Check the current conformation has no outliers
            current_conformation.check_border()
        
        return current_conformation



class REMC:
    '''Implementation of Replica Exchange Monte Carlo

    Attributes
    ---
    chi : list of conformations
    temps : list of temperatures
    steps : number of steps to use for the individual MC sims

    Methods
    ---
    run_remc_sim : run a REMC on the input data
    '''
    def __init__(self, chi: list, temps: list, steps: int):
        '''Initialise class

        chi : list of replicates (aka list of conformation objects)
        temps : list of temperatures
        steps: number of iterations for single MC
        coupling_map : array with conformations and associated temperatures and energies
        '''
        self.chi = chi
        self.temps = temps
        self.steps = steps
        try:
            self.coupling_map = [[c, t, 0] for c, t in zip(self.chi, self.temps)]
        except:
            raise ValueError('Incorrect input format')

    
    def run_remc_sim(self, E_star: int, max_runtime: float):
        '''Main method to run the REMC

        Args
        ---
        E_star : expected optimal energy
        max_runtime : how long can a run take ?
        '''
        # Check input energy
        if E_star > 0:
            raise ValueError('You cannot input a positive energy !!')

        # Initialise energy and offset for pair comparisons
        saved_conformation = None
        model_E = 0
        offset = 0
        iteration = 1

        # Count swaps and non-swaps
        swaps = 0
        non_swaps = 0

        # Runtime condition
        start = time.time()
        runtime = 0
        cutoff = max_runtime

        with open('results_log.txt', 'a') as filin:
                filin.write(f'\nNew Run\n')
        
        print(f'\nNew Run\n')

        # Start loop
        while model_E > E_star and runtime < cutoff:

            print(f'\nIteration : {iteration}\n')
            
            # Run single MC on all conformation/temperature pairs
            for pair in self.coupling_map:

                conf = pair[0]
                temp = pair[1] 

                # Run monte carlo and save new conformations
                sim = MonteCarlo(conf, self.steps, temp)
                res = sim.run_sim()

                # If the energy of conformation is lower than model E, save it
                e = res.calculate_energy()
                pair[2] = e

                if e < model_E:
                    model_E = e
                    saved_conformation = copy.deepcopy(res)

                # Update the coupling map with a copy of the new conformation
                pair[0] = copy.deepcopy(res)
            
            # Log results 
            view = saved_conformation.view_conformation()
            with open('results_log.txt', 'a') as filin:
                filin.write(f'\nIteration : {iteration}\n')
                filin.write(f'Current best model energy : {model_E}\n')
                filin.write('Matrix representation of conformation: \n')
                filin.write(f'{view}\n')
            
            # Now replica exchange
            i = offset
            while i + 1 < len(self.chi):
                
                j = i + 1

                # Recover corresponding conf energies and temperatures
                ei = self.coupling_map[i][2]
                tempi = self.coupling_map[i][1]
                ej = self.coupling_map[j][2]
                tempj = self.coupling_map[j][1]

                # Determine swap conditions
                delta = ((1/tempj) - (1/tempi)) * (ei - ej)
                if delta < 0:
                    self.coupling_map[i][1] = tempj
                    self.coupling_map[j][1] = tempi
                    swaps += 1
                else:
                    rng = np.random.default_rng()
                    q = rng.random()
                    if q <= np.exp(-delta):
                        self.coupling_map[i][1] = tempj
                        self.coupling_map[j][1] = tempi
                        swaps += 1
                    else:
                        non_swaps += 1
                
                # Increment i and runtime
                i += 2
                mark = time.time()
                runtime = mark - start
            
            # Record swap to non swap ratio
            with open('results_log.txt', 'a') as filin:
                filin.write(f'Current number of swaps: {swaps}\n')
                filin.write(f'Current number of non swaps: {non_swaps}\n')
                filin.write(f'Ratio of swaps: {swaps / (swaps + non_swaps)}\n')

            # Increment offset and iteration
            offset = 1 - offset
            iteration += 1
        
        if runtime > cutoff:
            print('\nREMC stopped : Maximum computation time reached')
        else:
            print('REMC stopped : Conformation with input energy reached')
        
        print(f'Runtime : {mark - start}')
        with open('results_log.txt', 'a') as filin:
                filin.write(f'Run finished\n')
                filin.write(f'Runtime: {mark - start}\n')

        return saved_conformation

                        

# Utility functions

def sequence_to_hp(sequence):
    '''Converts a string of amino acids into a numbered hp sequence
    to be passed to the main programme. 

    Arguments
    ---
    sequence : str

    Returns
    ---
    hp_sequence : str
    '''
    ref_dico = {
    'A': 'H',  # Alanine
    'R': 'P',  # Arginine
    'N': 'P',  # Asparagine
    'D': 'P',  # Aspartic acid
    'C': 'H',  # Cysteine
    'Q': 'P',  # Glutamine
    'E': 'P',  # Glutamic acid
    'G': 'P',  # Glycine
    'H': 'P',  # Histidine
    'I': 'H',  # Isoleucine
    'L': 'H',  # Leucine
    'K': 'P',  # Lysine
    'M': 'H',  # Methionine
    'F': 'H',  # Phenylalanine
    'P': 'H',  # Proline
    'S': 'P',  # Serine
    'T': 'P',  # Threonine
    'W': 'H',  # Tryptophan
    'Y': 'P',  # Tyrosine
    'V': 'H'   # Valine
    }

    hp_sequence = []
    for n, aa in enumerate(sequence.upper()):
        x = ref_dico[aa]
        hp_sequence.append(f'{x}{n+1}')

    return hp_sequence


def create_protein_conformations(protein_name, hp_sequence, min_t, max_t, nb_conformations):
    '''Function to create protein object from an input
    HP sequence

    Args
    ---
    protein_name : str
    hp_sequence : list of strings
    min_t : minimum temperature (int)
    max_t : maximum temperature (int)
    nb_conformations : int

    Returns
    ---
    confs : list
    temperatures : list
    ''' 
    confs = []
    for i in range(nb_conformations):
    
        prot = Protein(protein_name)
        
        # Add amino acids
        for aa in hp_sequence:
            prot.add_aa(HPAminoAcid(aa))

        # Create lattice
        lattice = Lattice2D(4 * len(hp_sequence), 4 * len(hp_sequence))

        # Create conformation
        confs.append(Conformation(f'C{i+1}', prot, lattice))


    # Create temperatures - this is the uniform linear function
    def t(k):
        return min_t + ((k - 1) * (max_t - min_t) / (nb_conformations - 1))
    
    
    temps = [t(k + 1) for k in range(nb_conformations)]

    return confs, temps


def plot_final_conformation(conformation: Conformation, figure_path, name):
    '''Function to plot the calculated conformation of the input 
    protein using matplotlib. 
    
    Args
    ---
    conformation : Conformation object
    figure_path : Location to save the figure
    name : Save name of the figure

    Returns
    ---
    Nothing - saves matplotlib plot
    '''
    # Extract matrix form of conformation and recover sequence + energy
    protein_array = conformation.view_conformation()
    seq = conformation.get_protein_sequence()
    energy = conformation.calculate_energy()

    # Find the coordinates and sequence numbers of the residues
    coords = {}
    for i in range(protein_array.shape[0]):
        for j in range(protein_array.shape[1]):
            if protein_array[i, j] != 0:
                coords[protein_array[i, j]] = (i, j)

    # Sort the coordinates by the residue sequence number
    sorted_coords = [coords[key] for key in sorted(coords.keys())]

    # Separate the coordinates into x and y lists for plotting
    x_coords = [coord[1] for coord in sorted_coords]
    y_coords = [-coord[0] for coord in sorted_coords]  # Negate to keep the plot intuitive (higher on top)

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Plot the connections (lines between residues)
    plt.plot(x_coords, y_coords, '-o', color='black')

    # Plot each residue as a circle with a label
    if len(seq) < 30:
        fs = 10
    elif len(seq) < 60:
        fs = 8
    elif len(seq) < 100:
        fs = 6
    else:
        fs = 5

    for i, (x, y, aa) in enumerate(zip(x_coords, y_coords, seq), start=1):
        if aa.type == 'H':
            plt.text(x, y, f"{aa.name}", fontsize=fs, ha='center', va='center', color='white', 
                     bbox=dict(facecolor='blue', edgecolor='black', boxstyle='circle,pad=0.3'))
        else:
            plt.text(x, y, f"{aa.name}", fontsize=fs, ha='center', va='center', color='white', 
                     bbox=dict(facecolor='darkgreen', edgecolor='black', boxstyle='circle,pad=0.3'))

    # Set axis limits and labels
    plt.xlim(-1, protein_array.shape[1])
    plt.ylim(-protein_array.shape[0], 1)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(f"{conformation.protein.name} 2D Representation - E = {energy}")
    plt.xticks([])
    plt.yticks([])
    plt.grid(True)

    # Save the plot
    plt.savefig(f'{figure_path}/{name}_{conformation.protein.name}.png')       






