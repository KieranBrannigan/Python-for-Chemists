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
    Given path to the data files,
    returns the list of file paths.
    """
    
    file_names = os.listdir(path)
    
    # convert the file_names to full (not-relative) paths
    file_paths = [os.path.join(path,name) for name in file_names]
    
    return file_paths


def find_maxima(matrix):
    """
    given input of a 2d matrix of x,y coordinates,
    this will find all y-maxima in the matrix and return
    these maxima as a list.
    """
    
    # loop through the indexes of the y coordinates of the matrix matrix.T[1]
    maxima = []
    for i in range(1,len(matrix.T[1])-1):
        
        if matrix.T[1][i-1] < matrix.T[1][i] > matrix.T[1][i+1]:
            maxima.append(matrix[i])

    maxima = np.array(maxima)            
            
    return maxima
    


def definite_integral(matrix):
    """
    Given matrix of x,y coordinates compute the definite integral between 
    x=0 and x=1
    """
    
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

    """
    

    if compare_file == 'data100.txt':
        return None

    
    # Whenever there is a match we will iterate this, so that we can compare
    #this at the end?
    number_of_matches = 0

    dx_values = []
    dy_values = []
    
    
    
    for og_idx,og_val in enumerate(og_maxima.T[0]):
        for idx,val in enumerate(new_maxima.T[0]):
            #print(f"val = {val} og_val = {og_val}")
            #this will loop dx from 0 to (until-1)/divisor in steps of 1/divisor
            until = until + 1
            for x in range(until):
                dx = x/divisor
                if og_val - dx <= val <= og_val + dx:
                    print(f"we have a match index {og_idx} from og and {idx} from comparefile")
                    print(f"values are {og_val} and {val}")
                    number_of_matches+=1
                    dx_values.append(dx)
                    
                    # Get the absolute value of the difference in y values
                    dy = abs(og_maxima.T[1][og_idx] - new_maxima.T[1][idx])
                    dy_values.append(dy)
                    #breaks us out of the for x in range loop
                    break 
            
            # If the for loop (for x in range ...) doesn't run IE we got a match
            else:
                "move onto next peak in og_maxima"
                continue
            break
            #If the for loop does run IE we don't get a match
            "compare next peak in new_maxima, IE continue"
        
        #If the for loop (for idx,val ...) didn't run completely IE we got a match
        

    # TODO how about instead see what the number of peak matches is?
    # a peak could have the same number of peaks but 
    different_no_peaks = abs(len(og_maxima) - len(new_maxima))

    #what if they only have 2 peaks that match

    return [dx_values, dy_values, number_of_matches, different_no_peaks]




if __name__ == "__main__":
    
    
    data = 'Data'
    
    files = read_data_files(data)
    
    maxima_in_range = []
    
    file100_matrix = read_file(files[0])

    file100_maxima = find_maxima(file100_matrix)

    comparisons = []

    comparison_values = {} # [file_name:value]

    for path in files:
            
        file_name = os.path.basename(path)
        
        matrix = read_file(path)
        
        show_graph_from_matrix(matrix)
        
        plt.axis([0,1,0,10])
        plt.savefig(os.path.join('normalised-plots',f"{file_name.split('.')[0]}.png"))
        plt.clf()
        
        print(file_name)
        
        maxima = find_maxima(matrix)
    
        #print(maxima)
        
        integral = definite_integral(matrix)
        
        #print(integral)
        
        normalisation = 1/integral
        
        matrix = np.array([matrix.T[0],matrix.T[1]*normalisation]).T
        

        # Integral of the normalised matrix (should be equal to 1)
        normal_integral = definite_integral(matrix)
        print(f"normalised integral : {normal_integral}")
        #show_graph_from_matrix(matrix)
        #plt.show()
        #plt.clf()
        for maximum in maxima:
            if 0 < maximum[0] < 0.05:
                maxima_in_range.append(file_name)

        comparison_data = compare_files(file100_maxima, maxima, file_name)
        if comparison_data is not None:
            print(f"adding {file_name} to comparison_values dictionary")
            comparison_value = sum(comparison_data[0])*1000 + sum(comparison_data[1])/100 + abs(6-comparison_data[2])*1000 + comparison_data[3]/10
            comparisons.append([file_name]+ comparison_data + [comparison_value])
            print(comparison_value)
            comparison_values[file_name] = {comparison_value}
        else:
            print(f"COMPARISON_DATA IS NONE : {file_name}")
        


    #print("data that contains maxima with x value in range 0-0.05 : ")
    
    #for i in maxima_in_range:print(i)
    

    #TODO sort the output by their comparison_values
    
    comparisons = sorted(comparisons, key=lambda row: row[-1])

    with open('output_values_v2.csv', 'w',) as WriteFile:
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
            writer.writerow([sum(row[1])*1000, sum(row[2])/100, abs(6-row[3])*1000, row[4]/10])


    
    """
    x_maxima = []
    for file in files:
        matrix = read_file(file)
        x_values = matrix.T[0]
        x_maxima.append(max(x_values))
    
    # all maxima are 0.99
    print(x_maxima)

    """
"""


maxima values for data100.txt
[[0.1        1.00067571]
 [0.24       2.00482677]
 [0.35       3.30382273]
 [0.4        2.89867253]
 [0.57       4.00386628]
 [0.8        6.15272174]]
    
TODO
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
