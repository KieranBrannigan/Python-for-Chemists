# Day 6 Python for chemists

import os
import sys
import csv

import numpy as np
import matplotlib.pyplot as plt

def read_file(path):
    """
    Arguments is the path to the file
    reads a text file input of x and y coordinates and returns the values as a 
    2 dimensional array
    """
    
    # open the file in read mode and create a numpy array of a list
    # which is each line of the file, split by whitespace
    # this is then converted to a numpy float type matrix.
    with open(path,'r') as ReadFile:
        matrix = np.array(
                [line.split() for line in ReadFile.readlines()]
                ).astype(np.float)
        
    return matrix

def show_graph_from_matrix(matrix):
    """
    Arguments is a 2 dimensional matrix of coordinates x and y
    [[x1,y1]
    [x2,y2]
    ...]
    The 2d matrix will be displayed using matplotlib.pyplot
    """
    
    plt.plot(matrix.T[0],matrix.T[1],)
    

def read_data_files(path):
    """
    Given path to the data files directory,
    returns the list of file paths.
    """
    
    # returns the list of all file names in the given directory
    file_names = os.listdir(path)
    
    # convert the file_names to full (not-relative) paths
    file_paths = [os.path.join(path,name) for name in file_names]
    
    return file_paths


def find_maxima(matrix):
    """
    given input of a 2d matrix of x,y coordinates,
    this will find all y-maxima in the matrix and return
    these maxima as a list.

    2 dimensional matrix of coordinates x and y
    [[x1,y1]
    [x2,y2]
    ...]

    """
    
    maxima = []

    # loop through the indexes of the y coordinates of the matrix matrix.T[1] 
    # from the 2nd index until the 2nd from last index
    for i in range(1,len(matrix.T[1])-1):
        # check if the y value is greater than the value before and after
        # if this is the case it is a maximum, so append the x,y coordinates
        # to the 
        if matrix.T[1][i-1] < matrix.T[1][i] > matrix.T[1][i+1]:
            maxima.append(matrix[i])

    # convert maxima into a numpy array
    maxima = np.array(maxima)          

    return maxima
    


def definite_integral(matrix):
    """
    Given matrix(**) of x,y coordinates compute the definite integral between 
    x=0 and x=1
    
    (**) 2 dimensional matrix of coordinates x and y
    [[x1,y1]
    [x2,y2]
    ...]
    """
    
    # Calculates the integral of the supplied matrix, using the supplied x values
    # the dx for the x values is consistently 0.01 so accuracy shouldn't be 
    # a big issue.
    integral = np.trapz(matrix.T[1], matrix.T[0])
    
    return integral

def compare_files(og_maxima,new_maxima, compare_file, until=100, divisor=1000):
    """
    given input of the maxima of a graph, compare it to the maxima from data100.txt

    maxima will be a series of x,y coordinates corresponding to the x,y values of a maximum from a file.

    First see if there is a maxima with the same x value as data100.txt, if there is not expand the x value ranges
    until a maximum is found. Find out what this dx is for the new file.  
    Note do it for all the peaks of data100.txt at once, so that if it finds a peak for the 2nd peak of data100.txt,
    it doesn't also assign this to the first peak as well.
    
    kewyword arguments until and divisor:
        for the dx loop the loop will increase dx from 0 until until/divisor in steps of 1/divisor
        eg for default values until=100 and divisor=1000,
        it will increase dx from 0 until 100/1000 (=0.1) in steps of 1/1000 (=0.001)

        changing these arguments will lead to more or less peak matching, which could 
        affect the results of the calculation significantly.

    """
    

    if compare_file == 'data100.txt':
        return None

    
    # Whenever there is a match we will iterate this, so that we can compare
    #this at the end?
    number_of_matches = 0

    # Initiate two lists to contain all the dx and dy values for each peak that
    # is matched by the code.
    dx_values = []
    dy_values = []
    
    # Loop through the original maxima list (supplied as an argument) 
    # and also loop through the maxima from the file being compared.
    for og_idx,og_val in enumerate(og_maxima.T[0]):
        for idx,val in enumerate(new_maxima.T[0]):

            
            #this will loop dx from 0 to (until)/divisor in steps of 1/divisor
            for x in range(until+1):
                dx = x/divisor
                # For the current value of dx see if there is a matching
                # peak between the data100.txt file and the file being compared.
                # There is a match if the val from the compare_file is within the range
                # of the original peak x value +/- the dx value.
                if og_val - dx <= val <= og_val + dx:
                    #if there is a match print some logging information to the console.
                    print(f"Peak Match : index {og_idx} from data100.txt and {idx} from {compare_file}")
                    print(f"values are {og_val} and {val} respectively")

                    # iterate the number of peak matches between the two files being compared.
                    number_of_matches+=1
                    # append the current dx value to our running list which will keep track
                    # of the dx values for all the matched peaks
                    dx_values.append(dx)
                    
                    # Get the absolute value of the difference in y values (dy)
                    dy = abs(og_maxima.T[1][og_idx] - new_maxima.T[1][idx])
                    dy_values.append(dy)
                    #breaks us out of the "for x in range" loop
                    break 
            
            # If the for loop (for x in range ...) isn't terminated by a break statement 
            # I.E. we didn't get a match
            else:
                "move onto next peak in new_maxima"
                continue

            # If the for loop does get terminated by the break statement 
            # I.E. we get a match
            """compare next peak in og_maxima, IE break the new_maxima loop and move onto
            next in the original maxima list"""
            break
        

    # Calculate the absolute value of the difference in number of peaks 
    # between the two data files
    different_no_peaks = abs(len(new_maxima) - len(og_maxima))

    return [dx_values, dy_values, number_of_matches, different_no_peaks]



# If this file is imported by another file, then this code won't be run.
if __name__ == "__main__":
    
    # Specify the path to the directory containing the data files
    data = 'C:/Users/Ale/Data/'

    # Get the paths for all the files contained within the data directory.
    files = read_data_files(data)
    
    # Initialise list which will contain files that have one maximum with
    # x value in the range 0 to 0.05
    maxima_in_range = []

    # Calculate the matrix for data100.txt
    data100_matrix = read_file(files[0])

    data100_maxima = find_maxima(data100_matrix)

    comparisons = []

    for path in files:
            
        file_name = os.path.basename(path)
        
        matrix = read_file(path)
        
        #show_graph_from_matrix(matrix)
        #plt.savefig(os.path.join('plots',f"{file_name.split('.')[0]}.png"))
        #plt.clf()

        # Log which file is being processed
        print(file_name)

        # Calculate the peaks of the data.
        maxima = find_maxima(matrix)
    
        print(f"maxima for {file_name} : \n{maxima}")
        
        # Compute the integral for the file
        integral = definite_integral(matrix)
        
        print(integral)
        
        normalisation = 1/integral
        
        # Convert the matrix to a numpy array of itself where all the y values are
        # multiplied by the normalisation constant
        # basically matrix = np.array(xvalues,yvalues * norm const)
        # where the x values are the first element of the transpose of the matrix
        # and the y values are the second element.
        # The numpy array is then transposed at the end to restore it to its original
        # orientation.
        matrix = np.array([matrix.T[0],matrix.T[1]*normalisation]).T
        
        # Integral of the normalised matrix (should be equal to 1)
        normal_integral = definite_integral(matrix)
        print(f"normalised integral : {normal_integral}")
        #show_graph_from_matrix(matrix)
        #plt.show()
        #plt.clf()

        # creates a list of maximums where the x value is between 0 and 0.05
        # if this list length is exactly 1 then it appends the file name to the
        # maxima_in_range list.
        if 2 > len(
            [maximum for maximum in maxima if 0 < maximum[0] < 0.05]
            ) > 0 : 
            maxima_in_range.append(file_name)

        # Calls the compare files function on the current working file and data100.txt
        # returns the dx values, dy values, number of peak matches, difference in number 
        # of peaks
        comparison_data = compare_files(data100_maxima, maxima, file_name)

        # comparison data returns None if it is comparing data100.txt and data100.txt
        if comparison_data is not None:
            # Calculate the comparison value for the current file
            comparison_value = sum(comparison_data[0])*1000 + sum(comparison_data[1])/100 + abs(6-comparison_data[2])*1000 + comparison_data[3]/10

            # Add the file name, comparison data and comparison values to the comparisons list
            comparisons.append([file_name]+comparison_data + [comparison_value])

            print(comparison_value)
        else:
            print(f"COMPARISON_DATA IS NONE : {file_name}")
        
    print("data that contains one maximum with x value in range 0-0.05 : ")
    
    for i in maxima_in_range:print(i)
    
    # Sort the files by their comparison values (the comparison value is stored as the last element
    # of each row)
    comparisons = sorted(comparisons, key=lambda row: row[-1])


    output=f"""
We used an algorithm to match peaks from all the files in the dataset to the peaks from data100.txt
This is due to peaks having the most chemical significance for absorption spectra.  All data was 
normalised so that integral = 1 to reduce any differences that could arise because of differences
in preparation of sample or machine used.
If a peak was matched then the dx and dy values of this peak was noted. Using this data an comparison
score was calculated.  Most emphasis was placed on how close the number of matched peaks was to 6 (there
were 6 peaks present in data100.txt), after that the sum of the dx and dy values were summed and also the
difference in the number of peaks was acounted for.  
After sorting the files by these metrics we found the file that most closely matched data100.txt
was {comparisons[0][0]}
"""


    print(output)




    # Code for outputting to a file all the compared files with their comparison values
    """

    with open('output_values.csv', 'w',) as WriteFile:
        writer = csv.writer(WriteFile,dialect='excel')
        writer.writerow(
            ['file_name'
            ,'dx values'
            ,'dy values'
            ,'number of peak matches'
            ,'difference in number of peaks'
            , 'comparison value']
            )
        for row in comparisons:
            writer.writerow(row)
            
    """

    """

    # Code to check the maximum x values for every file
    x_maxima = []
    for file in files:
        matrix = read_file(file)
        x_values = matrix.T[0]
        x_maxima.append(max(x_values))
    
    # Check all maxima are 0.99
    print(x_maxima)

    """
"""

# For reference
maxima values for data100.txt
[[0.1        1.00067571]
 [0.24       2.00482677]
 [0.35       3.30382273]
 [0.4        2.89867253]
 [0.57       4.00386628]
 [0.8        6.15272174]]
    
DONE
see what graphs are similar in terms of having peaks in similar 
positions to these maxima.
Is this specific enough? if not see ideas below for further refinement

 
- do above but for the minimum values (this doesn't make sense chemically because
the troughs are not as chemically significant as the peaks, if at all significant.)
- try to find peaks that are in the ranges and see what the dx and dy values 
are for these similar peaks, then try to quantify this to see which has the
lowest dx values and lowest dy values.

- run the code so that we end up with a number for each file
- we can then have a loop that reruns the compare_matrix function
with smaller numbers for the dx range until we refine it down?
- *** that sounds more consuming than just starting low and continuing until
we have a single file. 


""" 
