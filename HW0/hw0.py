# Note: this is the Python syntax for importing modules
# We sometimes use package as a synonym for module, and formally
# a package is a module with a path on the file system.
# Modules may have submodules, which are delimited with .
# i.e. module.submodule
import argparse, os, sys

# Note: 'as' provides an alias for a module so that you can
# access its functions with syntax like np.mean() instead of numpy.mean()
import numpy as np

np.random.seed(776)

# Note: this is the syntax for Python functions. You do not need to
# specify a return type in advance. Multiple values can be returned
# using tuples such as (value1, value2, value3).
def read_file(input_filename):
	# Note: informally, 'with' statements are used for file reading and
	# writing, among other things. They guarantee that the file is properly
	# closed regardless of whether the code block runs normally or throws
	# an exception. The line below opens the file in read mode and creates
	# a file object.
	with open(input_filename, "r") as in_file:
		# Read the file contents.
		# For normal NumPy data, you can load by np.load.
		# But for practice, you need to read them from plain text here.
		matrix = []
		for line in in_file:
			# https://stackoverflow.com/a/27222113
			# map() handles \n fine
			matrix.append( list( map( float, line.split(" ") ) ) )
	return np.array(matrix)
	
def print_dict(dict, out):
	for i in dict.items():
		print(i[0] + ":" + str(i[1]), end=" ", file=out)
	print("\n", file=out)
	
def print_matrix(matrix, out):
	for y in matrix:
		print(*y, file=out)
	print(file=out)

def hw0_test(input_file, output_file):
	# Note: out is the output file you can write to with the print function.
	with open(output_file, "w") as out:
		#####
		#Part 1: strings
		#####
		print("Part 1:", file=out)
		# Create variables first_name and last_name with your name.
		first_name = "Noah"
		last_name = "Cohen Kalafut"

		# Create full_name by concatenating the variables with a tab
		# separating them and print it.
		full_name = first_name + "\t" + last_name
		print(full_name, end="\n\n", file=out)

		# Transform the full_name string to all uppercase letters and
		# print it.
		print(full_name.upper(), end="\n\n", file=out)
		

		#####
		# Part 2: lists
		#####
		print("\n\nPart 2:", file=out)
		# Initialize a list x with four zeros.
		x = [0] * 4
		
		# Prepend one 1 to the head and append one 1 to the tail. Print x.
		x = [1] + x + [1]
		print(x, end="\n\n", file=out)

		# Set y to be the first four elements of x.
		y = x[:4]
		
		# Change the second to last element of y to 2. Print y.
		y[-2] = 2
		print(y, end="\n\n", file=out)

		# Calculate the product of the elements in y.
		# Pass (skip over) the element if it is 0.  Print the result.
		prod = 1
		for e in y:
			if e != 0:
				prod *= e
		print(prod, end="\n\n", file=out)
		
		print(file=out)
		# Note: Python strings can be indexed in the same manner as lists.
		course_str = "Advanced Bioinformatics"
		
		# Print the index of 'B'.
		B = course_str.find("B")
		print(B, end="\n\n", file=out)

		# Print the number of occurrences of "i" and assign it to variable k.
		k=course_str.count("i")
		print(k, end="\n\n", file=out)

		# Print the substring of length k starting with 'B'.
		print(course_str[B:B+k], end="\n\n", file=out)
		

		#####
		# Part 3: dictionaries and sets
		#####
		print("\n\nPart 3:", file=out)
		# Note: This is set syntax.
		keys = {"a", "b", "c", "d"}
		# Create a dictionary called hash_map.
		hash_map = {}
		
		# Map the chararacters a-d to 1-4. Save the mapping in hash_map.
		for l,n in zip("abcd",range(1,5)):
			hash_map[l] = n

		# Check if "e" exists in hash_map. If not, map it to 5.
		if "e" not in hash_map:
			hash_map["e"] = 5
		
		# Print all key-value pairs in format <key:value> like "a:1".
		print_dict(hash_map, out)
		# Is this outside of point 5 on the requirements and hints?

		# Map "e" to 6.
		hash_map["e"] = 6
		
		# Print all key-value pairs in format <key:value> like "a:1".
		print_dict(hash_map, out)
		

		#####
		# Part 4: NumPy arrays
		# 
		# Note: you may write a function print_matrix(matrix, output_file) 
		# to print matrices in the desired format.
		#####
		print("\n\nPart 4:", file=out)
		u = [1, 2, -3]
		v = [3, -2, 1]
		
		# Convert u and v into NumPy arrays.
		u = np.array(u)
		v = np.array(v)

		# Calculate a, the inner product of u and v. Print a.
		# Hint: a is a scalar.
		a = u.dot(v)
		print(a, end="\n\n", file=out)

		# Calculate B, the outer product of u and v. Print B.
		# Hint: B is a 3 by 3 matrix.
		B = np.outer(u, v)
		print_matrix(B, out)

		# Calculate C = a * B. Print C.
		C = a * B
		print_matrix(C, out)

		# Create R, a 3 by 3 random matrix. Print R.
		R = np.random.rand(3,3)
		print_matrix(R, out)

		# Calculate the matrix product RC. Print the result.
		print_matrix(R.dot(C), out)

		# Calculate the matrix product CR. Print the result.
		print_matrix(C.dot(R), out)

		# Calculate the elementwise product of R and C. Print the result.
		print_matrix(R * C, out)
		

		#####
		# Part 5: NumPy sorting and binning
		#####
		bins = 5
		
		print("\n\nPart 5:", file=out)
		# Complete the read_file function and call it to
		# read the matrix D from input_file into a NumPy array.
		D = read_file(input_file)

		# Sort the first column of D in ascending order.
		# Assuming that this means the column by itself
		D_sort = np.sort(D[:,0])
		
		# Print the first three values of the sorted list.
		print(*D_sort[:3], end="\n\n", file=out)

		print("\nEqual Width:", file=out)
		# Assign the sorted values into five bins with equal width.
		# The left edge of the first bin is the min value.
		# The right edge of the last bin is the max value.
		# For each bin, print
		#  1) its left and right edges, and
		#  2) the values that fall inside the bin.
		_, bounds = np.histogram(D_sort, bins)
		
		iter_r = 0
		for i in range(1, len(bounds)):
			l,r = bounds[i-1],bounds[i]
			iter_l = iter_r
		
			# Values
			# https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
			# Bins are half-open on the right-hand side for histogram()
			# https://numpy.org/doc/stable/reference/generated/numpy.digitize.html
			# Bins are *default* half-open open on the right for digitize() as well
			print("Values: ", end="", file=out)
			if i == len(bounds)-1:
				print(D_sort[iter_l:], file=out)
			else:
				while D_sort[iter_r] < r:
					iter_r += 1
				print(D_sort[iter_l:iter_r], file=out)
			
			# Edges
			print("Left Edge: " + str(l), file=out)
			print("Right Edge: " + str(r), file=out)
			print(file=out)
		

		print("\nEqual Size:", file=out)
		# Assign the sorted values into five bins such that each bin 
		# contains the same number of elements.
		# Print the values that fall inside each bin.
		increment = len(D_sort) / bins
		
		# An arbitrary bin b can be accessed with [round(b*increment):round((b+1)*increment)]
		# Tries to work for incompatible bin numbers too using round()
		for b in range(bins):
			print("Values: ", end="", file=out)
			print(D_sort[round(b*increment):round((b+1)*increment)], file=out)
		print(file=out)
		
		print(file=out)
		# Sort the last row of D in descending order.
		D_sort_desc = -np.sort(-D[-1])
		
		# Print the last three values of the sorted list.
		print(*D_sort_desc[-3:], file=out)

def main(args):
	input_file = args.inputfile
	output_file = args.outputfile

	hw0_test(input_file, output_file)


# Note: this syntax checks if the Python file is run as the main program
# and will not execute if the module is imported into a different module.
if __name__ == "__main__":
	# Note: this example shows named command line arguments. See the argparse
	# documentation for positional arguments and other examples.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	# Note: you can use ' or " for strings.
	parser.add_argument('--inputfile',
						help='input data file path.',
						type=str,
						default='data.txt')
	parser.add_argument('--outputfile',
						help='output file path.',
						type=str,
						default='out.txt')

	args = parser.parse_args()
	# Note: this simply calls the main function above, which we could have
	# given any name.
	main(args)
