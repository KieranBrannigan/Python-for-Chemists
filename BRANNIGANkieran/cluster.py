# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 11:31:56 2019

@author: Kieran
"""

def shortest_distance_of_cluster(cluster):
    """
    given an input NumpyArray cluster, it computes the shortest distance between any pair of atoms.
    It returns the shortest_distance and the closest atoms. 
    shortest_distance is a float.
    closest atoms is a formatted string : f"atom n and atom q" where n and q are integers.  
    """
    distances = {} # Will have {atom n and atom q: distance value}
    
    for i in range(len(cluster)-1): # Loop through first till 2nd from last vectors in cluster
        for j in range(i+1, len(cluster)): #  Loop through from i till last vectors in cluster
            distances[f"atom {i+1} and atom {j+1}"] = dist3D(cluster[i],cluster[j])
    
    # Assign the shortest value in distances dictionary to var shortest_distance
    shortest_distance = min(distances.values())
    
    # Loops through the keys and values (as atoms and distance) and if distance is the shortest distance
    # it assigns the atoms (key) to closest_atoms.  This way we know which atoms are the closest.
    for atoms, distance in distances.items():
        if distance == shortest_distance:
            closest_atoms = atoms
    
    return shortest_distance, closest_atoms

def total_lennard_jones_potential(cluster):
    """
    Takes a NumpyArray cluster as input for a cluster of atoms, calculates the total lennard_jones_potential of the cluster
    """
    
    distances = {} # Will have {atom n and atom q: distance value}
    
    for i in range(len(cluster)-1): # Loop through first till 2nd from last vectors in cluster
        for j in range(i+1, len(cluster)): #  Loop through from i till last vectors in cluster
            distances[f"atom {i+1} and atom {j+1}"] = dist3D(cluster[i],cluster[j])
    
    # Creates a List which contains all the lennard jones potentials for each distance of pair of atoms in the cluster
    all_lennard_jones_potentials = [lennard_jones_potential(dist) for dist in distances.values()]
    
    # Computes the total of all the lennard jones potentials in the cluster of atoms
    total_potential = sum(all_lennard_jones_potentials)
    
    return total_potential


def lennard_jones_potential(x):
    """
    calculates the lennard-jones-potential of a given value x
    """
    return 4 * ((x ** -12) - (x ** -6))


def dist3D(start, final):
    """
    calculates the distance between two given vectors.  
    Vectors must be given as two 3-dimensional lists, with the aim of representing
    2 3d position vectors. 
    """
    (x1, y1, z1) = (start[0], start[1], start[2])        
    (x2, y2, z2) = (final[0], final[1], final[2])
    (xdiff,ydiff,zdiff) = (x2-x1,y2-y1,z2-z1)
    
    # Calculate the hypotenuse between the x difference, y difference and the z difference
    # This should give the distance between the two points
    
    distance = (xdiff **2 + ydiff ** 2 + zdiff ** 2) ** (1/2)
    return distance    
    
def centre_of_mass(mass, numpy_array):
    """
    
    takes the mass of the atoms as first argument, then takes a number of vectors as
    arguments and computes the centre of mass of all the atoms.  
    
    """
    
    # assigns all the x positions, y positions and z position values from all the vectors to three individual lists.
    x_values = []
    y_values = []
    z_values = []
    for vector in numpy_array:
        x_values.append(vector[0])
        y_values.append(vector[1])
        z_values.append(vector[2])
        
    # Calculate the centre of mass for the x, y and z direction and return these as a vector. 
    centre_x = mass * sum(x_values)/(mass * len(numpy_array))
    centre_y = mass * sum(y_values)/(mass * len(numpy_array))
    centre_z = mass * sum(z_values)/(mass * len(numpy_array))
    
    return [centre_x, centre_y, centre_z]
    
