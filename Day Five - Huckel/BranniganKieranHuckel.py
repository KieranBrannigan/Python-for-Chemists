import numpy as np


def read_matrix_from_file(file):
    """
    Given input, str file , which corresponds to a matrix, this function will read the file
    and assign it to NumpyArray matrix and return it
    """

    # open the file in read mode
    with open(file, 'r') as File:

        # calls the realines method of File object and appends the line split by
        # whitespace to list, this list is then converted to a numpy array as a int
        # type matrix
        matrix = np.array([line.split() for line in File.readlines()]).astype(np.int)

    return matrix
    

def size_of_matrix(matrix):
        
    matrix_height = len(matrix)
    matrix_width = len(matrix.T)

    return f"{matrix_height} x {matrix_width}"


# This function isn't used in the script because its easier to just use the linalg.eig
# but this serves as a way of putting comments on what all the code does.
def all_eig(matrix):
    """
    Takes the matrix as input.
    Computes the eigenvalues for all the orbitals in the matrix. These are the energy values of each molecular orbital.
    E_k = a + (b * (eigenvalue_k))
    Therefore, the lowest energy orbitals are the ones with the highest eigenvalues. 
    prints the eigenvalues and returns the eigenvalues.
    """
    # Use numpy linalg module to calculate the eigenvalues and eigenvectors of the matrix
    (eigenvalues, eigenvectors) = np.linalg.eig(matrix)

    # Get the indexes for the eigenvalues so we can select the right eigenvectors for the eigenvalues
    sort_indexes = eigenvalues.argsort()

    # Use the sort indexes to grab the eigenvectors assiociated with each eigenvalue
    eigenvectors = eigenvectors[:,sort_indexes]

    return eigenvalues, eigenvectors

def occupied_eig(matrix):
    """
    Takes the matrix as input.
    Computes the eigenvalues for half the orbitals in the matrix. These are the energy values of each molecular orbital.
    E_k = a + (b * (eigenvalue_k))
    Therefore, the lowest energy orbitals are the ones with the highest eigenvalues. 
    prints the eigenvalues and returns the eigenvalues.
    """

    # Use numpy linalg module to calculate the eigenvalues and eigenvectors of the matrix
    (eigenvalues, eigenvectors) = np.linalg.eig(matrix)

    # Couldn't get this to work
    # Sorts the eigenvalues in increasing, then indexes the last half
    #halfway = len(eigenvalues)/2
    #eigenvalues = np.sort(eigenvalues)[halfway:]

    sort_indexes = eigenvalues.argsort()

    # Get the indexes for the eigenvalues so we can select the right eigenvectors for the corresponsponding eigenvalues
    filled_indexes = []
    for idx, val in enumerate(eigenvalues):
        if val >= 0:
            filled_indexes.append(idx)

    # Filter out the negative eigenvalues to just leave the positive eigenvalues
    filled_eigenvalues = eigenvalues[filled_indexes]
    
    #filled_eigenvalues = np.array([i for i in np.sort(eigenvalues) if i >= 0])

    # access the filled eigenvectors by using our indexes 
    filled_eigenvectors = eigenvectors[:,filled_indexes]

    return filled_eigenvalues, filled_eigenvectors
    

def compute_bond_order(n, eigenvectors):
    """
    Takes as inputs: number of electrons int n ; a list of eigenvectors [[eigenvector_i, eigenvector_j],[eigenvector_i,eigenvector_j],...] list list_of_vectors ;

    We calculate the Huckel bond order using equation:

        P_ij = 2 SIGMA_doublyOccupiedOrbitalsk(eigenvector_i * eigenvector_j)

    so:

        using the eigenvectors for the filled orbitals,
        compute the sum of the products of the eigenvectors * number-of-electrons
    
    """

    # Dict containing keys of the bond order values and values of the atoms
    # those values apply to

    # Transpose the eigenvectors so that it matches up with what was in the Huckel theory notes
    eigenvectors = eigenvectors.T

    # This will be a dict that is upadates with keys 'atom {i} and atom {j}' where i and j are the atom numbers
    # and values of bond orders.
    bond_orders = {} 

    # loop through every pair of atoms
    for i in range(len(eigenvectors[0])):
        for j in range(i+1,len(eigenvectors[0])):
            
            # loops through the rows of the matrix and computes the product of atom i and j
            # then appends this product to list products. 
            products = [(eigenvectors[q][i] * eigenvectors[q][j]) for q in range(len(eigenvectors))]

            # These print statements might be useful for debugging purposes
            #print(eigenvectors)
            #print(products)

            # bond order is n (number of electrons) * the sum of all the products
            # here it is rounded to 5 decimal points
            bond_order = round((n * sum(products)), 5)

            # the atoms and bond order value are added to the bond_orders dict
            bond_orders[f'atom {i+1} and atom {j+1}'] = bond_order
                
    
    return bond_orders

def shortest_and_longest_bonds(bond_orders):
    """
    takes an input dict bond_orders , which has keys of 'bond {i} and bond {j}' and values of the bond order between those atoms
    loops through those bond orders and finds which atoms have the smallest and biggest bond orders.
    The atoms with the biggest bond orders are returned as the longest bonds, and the atoms with the lowest bond orders
    are returned as the longest bonds.    
    """

    # loops through the values in our bond orders dict and compares value with minimums
    # if minimums list is in inital state, set the current value to minimum[0]. Same applies to maximums.
    # The first element of each list is used as a working minimum/maximum and each value is compared to that.#
    # If the value is lower than the working minimum, it replaces the minimums list and becomes the new minimum
    # and the loop continues.
    # It replaces the whole list in case there were two non-minimums that were evaluated before a newer minimum
    # was evaluated.
    # If the value is not lower than the minimum but it is equal to the current minimum then it is appended to
    # the minimums list.
    # the same logic is applied to the maximums.
    for value in bond_orders.values():

        # if minimum variable hasn't been initialised yet, set its value to value
        # same applies for maximum variable
        if 'minimum' not in locals():
            minimum = value

        elif value > 0.3 and value < minimum:
            minimum = value

        if 'maximum' not in locals():
            maximum = value

        elif value > maximum:
            maximum = value

    print(minimum,maximum)

    # loops through the dict and compares the value to the minimum and maximum values.
    # if the value is equal to the minimum the key ('atom {i} and atom {j}) is appended to
    # longest_bonds list.
    # if the value is equal to the maximum the key is appended to the shortest_bonds list.
    longest_bonds = []
    shortest_bonds = []
    for key,val in bond_orders.items():
        if val == minimum:
            longest_bonds.append(key)
        if val == maximum:
            shortest_bonds.append(key)


    return shortest_bonds, longest_bonds

output = ""
if __name__ == "__main__":
    
    # Run our code

    # Task 4.1
    # Define our input files
    input_files = ['input2.txt']

    # All the outputs will be stored here to be written to an
    # output file at the end of the script.
    total_output = ""

    # run the code for each input file
    for input_file in input_files:
    
        # Read the file
        matrix = read_matrix_from_file(input_file)

        # Task 4.2
        # Get the matrix size
        matrix_size = size_of_matrix(matrix)

        # Task 4.3
        (eigenvalues, eigenvectors) = occupied_eig(matrix)

        bond_orders = compute_bond_order(2, eigenvectors)

        shortest_bonds, longest_bonds = shortest_and_longest_bonds(bond_orders)

        (filled_eigenvalues,filled_eigenvectors) = (eigenvalues,eigenvectors)
        (values,vectors) = np.linalg.eig(matrix)
        newline = ',\n'
        
        output = (f"""
input file : {input_file}
matrix : 
{matrix}
matrix size : {matrix_size}
all eigenvalues and eigenvectors : 
{values}
{vectors}
eigenvalues and eigenvectors for filled orbitals :
{filled_eigenvalues}
{filled_eigenvectors}
bond_orders :
{str(bond_orders).replace(', ',newline)}
shortest bonds :
{shortest_bonds}
longest bonds :
{longest_bonds}
""")
        
        # print the output as the code for each file as it is ran
        print(output)
        total_output+=output
    # Write all the outputs to an output file.
    # It is opened in write mode so the output file will get overwritten each time the script is ran.
    with open('BranniganKieranHuckel-output.txt', 'w') as File:
            File.write(total_output)
