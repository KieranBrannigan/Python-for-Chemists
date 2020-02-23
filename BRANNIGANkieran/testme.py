# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:06:21 2019

@author: Kiera
"""
import random
import time
import argparse
import os

import numpy as np

from cluster import total_lennard_jones_potential, shortest_distance_of_cluster

def get_cluster_from_file(file):
    """
    takes file as input (position vectors of atoms) and return NumpyArray from that file. 
        input file format is a csv file in form:
          x  y  z
    atom1 
    atom2
    atom3
    etc...
    
    """
    with open(file, 'r') as File:
        # List comprehension
        cluster = [line.split() for line in File]
    # Convert List cluster into a numpy array with each value as a float type.
    cluster = np.array(cluster).astype(np.float) 
    
    return cluster


def task4_1(file,):
    """
    Reads an input file of position vectors of atoms. 
    Returns the total energy of the cluster of atoms from the file.
    Appends the shortest distance between atoms of this cluster to an output file 'shortest_distance.txt'
    
    input file format is a csv file in form:
          x  y  z
    atom1 
    atom2
    atom3
    etc...
    
    """
    
    # gets the NumpyArray for the atoms cluster and assigns to cluster
    cluster = get_cluster_from_file(file) 
    
    # gets the file name if the file that is supplied as an argument
    with open(file,'r') as File: 
        file_name = File.name
    
    # Gets the shortest distance between any pair of atoms.
    # Assigns the distance to shortest_distance and the atom numbers to closest_atoms
    (shortest_distance,closest_atoms)  = shortest_distance_of_cluster(cluster)  
    
    return total_lennard_jones_potential(cluster), shortest_distance, closest_atoms
    
def task4_2(cluster_geometry, cluster_energy, verbose=False):
    """
    Takes as inputs a NumpyArray cluster_geometry (a series of position vectors of atoms) and Float cluster_energy.
    It slightly adjusts one cartesian coordinate of one atom, and recalculates the energy of the new geometry.
    If the new energy is lower than the original it returns the new coordinates and the new energy.
    Else it returns the old coordinates and energy.
    """
    
    # generates a random number for an element (atom) in the NumpyArray (cluster)
    random_atom_number = random.randint(0,len(cluster_geometry)-1)
    
    # Generates a random integer (0,1 or 2) which will refer to a random coordinate axis (x,y or z)
    random_atom_coordinate = random.randint(0,2)
    
    # assigns a random vector from the cluster to var random_coordinate
    random_coordinate = cluster_geometry[random_atom_number][random_atom_coordinate]
    
    # A simple logic block to determine whether the change will be positive or negative
    if random.randint(0,1) == 0:
        change_direction = +1
    else:
        change_direction = -1
    
    # random.random() generates a random float between 0 and 1.  *5 to give a random float between 0 and 0.5
    change_magnitude = random.random() # between 0 and 0.5
    
    # change is now the random float between 0 and 0.5, multiplied by either -1 or +1.
    # this gives a random float between -0.5 and 0.5
    change = change_magnitude * change_direction
    
    # a dictionary to assign the numbers (0,1,2) to their corresponding axes (x,y,z)
    axes = {0:'x',1:'y',2:'z'}
    
    # creates the altered coordinate of the random atom.
    random_coordinate_changed = random_coordinate + change
    
    # logs the change
    if verbose:
        print(f"changed atom number {random_atom_number} coordinate {axes[random_atom_coordinate]} by {change}")
    
    # copies the original cluster geometry to a new copy to be altered. 
    new_cluster_geometry = np.copy(cluster_geometry)
    
    # alters the random coordinate in the new geometry (the copy of the old one)
    new_cluster_geometry[random_atom_number][random_atom_coordinate] = random_coordinate_changed
    
    # this is just to make things more understandble in the code, it doesn't actually have functional meaning
    old_total_energy = cluster_energy    
    
    # calculates the new total energy from the geometry
    new_total_energy = total_lennard_jones_potential(new_cluster_geometry)
    
    # Logic blocks to determine what geometry and energies to return
    # if new energy < old energy return new geometry and energy, 
    # else if old energy < new energy return old geometry and energy,
    # else if they are equal return the new geometry and energy. 
    if new_total_energy < old_total_energy:
        #print("new energy < old energy")
        return new_cluster_geometry, new_total_energy
    elif old_total_energy < new_total_energy:
        #print("old energy < new energy")
        return cluster_geometry, old_total_energy
    elif new_total_energy == old_total_energy:
        #print("The total energies were not changed")
        return new_cluster_geometry, new_total_energy
    
    
def task4_3(cluster, energy, number_of_iterations=10000, verbose=False):
    """
    Inputs: NumpyArray cluster is a matrix for a cluster of atoms
    energy is the initial energy value that we would like to optimise to try and find a minimum energy.
    number_of_iterations defaults to 10000, but can be specified as a keyword argument to set the number of iterations
    that is required before the function assumes a minimum has been reached. 
    This function runs task4_2 on cluster until the lowest energy conformation is determined
    if the new_energy calculated is the old_energy number_of_iterations times then we assume this is a minimum energy conformation.    
    """
    i=0

    # Prepares the arguments to be inputted into the while loop. I figured it would be safer this way rather than inputting
    # the arguments directly into the loop, in case there needed to be modifications in the future.
    old_geometry = cluster
    old_energy = energy
    
    while True:
        #print("looping...")

        # Calculate and assign the new geometry and energy using task4_2 with the "old geometry" and "old energy"
        (new_geometry, new_energy) = task4_2(old_geometry, old_energy)

        # if old_energy is the same as the new energy, that means a new lower energy wasn't found.
        # in this case we iterate our iterator 'i', so that if this occurs a set number of times we can
        # end the function and assume a minimum has been found.
        if old_energy == new_energy:
            i +=1

        # If the returned energy is less than the supplied energy then a new minimum has been found, and
        # we want to restart the iterator on this new minimum, hence i=0.
        # If the verbose argument is passed as True then it will print the new energy, old energy and the
        # difference in the new and old minimums, each time a new minimum is found.
        elif old_energy > new_energy:
            i=0
            if verbose:
                print(f"\nold energy: {old_energy} \nnew energy: {new_energy} \ndifference: {new_energy-old_energy}\n")

        # When the loop has failed to find a new minimum after a set number of iterations, the loop ends and the function
        # assumes we have found a minimum. The function then returns this new geometry and energy of the assumed minimum.
        # The accuracy of the minimum is therefore dependant on the size of the number_of_iterations argument
        # ie. the higher the number_of_iterations argument the more accurate the minimum will be, however this tends to
        # make the function take a much longer time to complete.  
        if i == number_of_iterations:
            return new_geometry, new_energy

        # This assigns the working minimum to old geometry and old energy to be re-entered into the loop for the purpose
        # of comparing the "old" values with the "new" ones.
        old_geometry = new_geometry
        old_energy = new_energy
          
    

# The code within this block will only be run when this is the main file being run. If this file is imported into
# another file, then the code in this block will not be run.   
if __name__ == "__main__":
    
    # assigning our parser object so we can accept arguments from the command line
    # the arguments themselves should be fairly self-explanatory in the help sections.
    
    # I chose not to have a positional argument for the file being read as we should
    # only be reading one file anyway, but I implemented it as a optional argument
    # so that the script could be ran programmatically with many files to automate
    # the generation of many minimums of clusters.
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i","--iteration_maximum", 
                        help="""
                        This is an integer representing the maximum number of 
                        iterations that the simulation will run before assuming 
                        a minimum has been found.
                        A higher number will generally lead to a better minimum energy.""",
                        type=int,
                        default=1000,
                        )
    parser.add_argument("-n","--number_of_simulation_runs",
                        help="""
                        This is an integer which states how many times the
                        simulation will run to obtain assumed minimums.
                        The global minimum will be taken from all the minimum values
                        obtained, therefore a higher number here will generally lead
                        to a better minimum.
                        """,
                        type=int,
                        default=20,
                        )
    parser.add_argument("-v","--verbose",
                        help="""
                        If this optional parameter is supplied then the verbosity
                        will be set to True.  This means the program will print 
                        the progress of the simulation to the console. If it is
                        not given then it will default to false, and it will only
                        print the results of the simulation.
                        """,
                        action="store_true",
                        )
    parser.add_argument("-f", "--files",
                        help="""
                        This is an optional str type argument which will default to 'input.txt'.
                        It is the files to be read which will have their minimum energy
                        calculated.
                        """,
                        type=str,
                        nargs='+',
                        default=['input.txt'])
    args = parser.parse_args()
    
    (number, iterations, verbose,files) = (args.number_of_simulation_runs, args.iteration_maximum,args.verbose, args.files)
    
    # Supplied input file taken from http://www-wales.ch.cam.ac.uk/~jon/structures/LJ/tables.150.html
    
    for file in files:
        
        if not os.path.isfile(file):
            print(f"""
            can't open '{file}': 
            [Errno 2] No such file or directory: '{file}'
            This file will not be run.
            """)
            
        else:
            
            total_energy, shortest_distance, closest_atoms = task4_1(file,)
            
            
            print(
                  f"""
        Task 4.1:
                      
        the total lennard jones potential is: {total_energy}
        
        the shortest distance between two atoms was {shortest_distance}
        
        the closest atoms were {closest_atoms}
                  
                  """)
                
            
            cluster = get_cluster_from_file(file)
            energy = total_lennard_jones_potential(cluster)
            
            (calculated_geometry,calculated_potential) = task4_2(cluster,energy, verbose=verbose)
            
            
            
            print(
                  f"""
        Task 4.2:
                      
        the new lennard jones potential is : {calculated_potential}
        
        the new geometry is : \n\n{calculated_geometry}
                  
                  """)
            
            
            #print(f"old cluster: \n{cluster} \n\nold energy: \n {energy}\n")
            #(new_geometry, new_energy) = task4_2(cluster, energy)
            #print(f"\nnew geometry: \n{new_geometry}\n\nnew energy: \n{new_energy}")
            
        
                
            print("Task 4.3: \n")
            
            # gets the time when the simulation started
            start=time.time()
            
            # I didn't like it this way, because it was calling the function twice in the list
            # comprehensions, hence it wasn't actually just calling the script for the amount 
            # of times specified in number.
        #    calculated_minimums = dict(
        #            zip(
        #                [task4_3(cluster,energy,number_of_iterations=iterations, verbose=verbose)[1] for i in range(number)],
        #                [task4_3(cluster,energy,number_of_iterations=iterations, verbose=verbose)[0] for i in range(number)]
        #                )
        #            )
            
            # A dictionary to store the minimum energies and associated conformations
            # calculated in the simulation. 
            calculated_minimums = {}
            
            # Runs the simulation and calculates an assumed minimum a number of times  
            # The number of times it runs is defined by the int number.
            for i in range(number):
                (new_geometry, new_energy) = task4_3(cluster,energy,
                number_of_iterations=iterations,
                verbose=verbose,)
                
                # The calculated minimum energy and geometry are added to the calculated minimums dictionary
                calculated_minimums[new_energy] = new_geometry
            
            # returns the minimum value obtained from all the assumed minimums
            # calculated in the step above.
            minimum_of_minimums = min(calculated_minimums.keys())
            
            # Gets the associated geometry from the dictionary.
            minimum_geometry = calculated_minimums[minimum_of_minimums]
            
            # Prints some useful output to the console to show what the outcome of
            # the simulations was.
            print(f"number of minimums calculated: {len(calculated_minimums)}\n")
            print(f"calculated minimums :\n{list(calculated_minimums.keys())}\n")
            print(f"minimum energy found: {minimum_of_minimums}\n")
            print(f"conformation associated with minimum: \n{calculated_minimums[minimum_of_minimums]}\n")
            
            # gets the time when the simulation finished.  Taking start from this gives the time
            # taken for the simulation to run.
            finish=time.time()
            print(f"simulation took {finish-start} seconds\n")
            
            # writes the output from the simulations to an output file
            # the output is quite verbose so that the simulation values and
            # parameters might be evaluated at a later stage. 
            with open('minimums.txt','a') as OutputFile:
                OutputFile.write(f"""
minimum for file '{file}' : {minimum_of_minimums}
conformation for this minimum:

{minimum_geometry}

{number} minimums calculated and {iterations} used as maximum iteration.
Time taken was {finish-start} seconds.

----------------------------------------

""")
            # Console directs the user to where the output will be stored for further reference.  
            print("minimum added to output file in current working directory named 'minimums.txt'")
