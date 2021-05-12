import argparse, os, sys
import numpy as np
# You can choose to write classes in other python files
# and import them here.

# The output process is deprecated on newer versions of numpy, but the
# output remains valid.
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

# Read the given file
def read_exp(filename):
	data = []
	with open(filename, 'r') as f:
		# Get rid of title line
		f.readline()
		
		# Read data
		for line in f:
			line_split = line.split('\t')
			if len(line_split) < 1:
				break
			elif line_split[-1][-1] == '\n':
				line_split[-1] = line_split[-1][:-1]
			data.append(line_split)
			
	# Return data
	return np.array(data).astype(float)

def mi_out(mi, filename):
	with open(filename,'w') as f:
		for e in mi:
			print(str(e[0]).replace(' ',''), end="\t", file=f)
			print("%.3f" % round(e[1],3), file=f)
	

# Discretize the data
def bin(data, bin_str="uniform", bin_num=5):
	counts = np.zeros(data.shape)
	
	# Bin for each gene
	for g in range(data.shape[1]):
		gene_data = data[:,g]
		
		# Bins are left-side inclusive
		# Execute chosen binning method
		if bin_str.upper() == "UNIFORM":
			min = np.amin(gene_data)
			max = np.amax(gene_data)
			interval = (max-min)/bin_num
			
			bins = min + (1+np.arange(bin_num)) * interval
		elif bin_str.upper() == "DENSITY":
			percentiles = (1+np.arange(bin_num)) * 100 / bin_num
			bins = np.percentile(gene_data,percentiles)
		else:
			raise Exception("Binning method not found")
		
		# Calculate for each entry
		temp_counts = np.zeros(gene_data.shape)
		for i in range(bin_num):
			if i == bin_num-1:
				temp_counts += gene_data <= bins[i]+.1
				break
			temp_counts += gene_data < bins[i]
		temp_counts = bin_num - temp_counts
		
		# Append to discrete matrix
		counts[:,g] = temp_counts
		
	return counts.astype(int)
	

# Generate our count matrix
def c_matrix(disc,gene_i,gene_j,pseudo=.1):
	unique = np.sort(np.unique(disc))
	
	count = np.zeros((unique.shape[0],unique.shape[0]))
	for i in range(disc.shape[0]):
		count[disc[i,gene_i],disc[i,gene_j]] += 1
		
	return pseudo+count


# Calculate mutual information for one count matrix
def MI_Sub(count):
	dim = count.shape[0]
	
	sum = 0
	for i in range(dim):
		for j in range(dim):
			# Calculate probabilities
			x = np.sum(count,axis=1)
			Px = x[i] / np.sum(x)
			y = np.sum(count,axis=0)
			Py = y[j] / np.sum(y)
			Pxy = count[i,j] / np.sum(count)
			
			sum += Pxy * np.log2( Pxy / ( Px * Py ) )
	
	return sum
	
	
# Calculate mutual information for all combinations
def MI(disc):
	mi_data = np.zeros((0,2))
	for gene_i in range(disc.shape[1]):
		for gene_j in range(gene_i+1,disc.shape[1]):
			count = c_matrix(disc,gene_i,gene_j)
			mutual = MI_Sub(count)
			mi_data = np.vstack((mi_data,[(gene_i+1,gene_j+1),mutual]))
	
	return mi_data[np.flip(mi_data[:,1].argsort(), axis=0),:]
	
	
def main(args):
	# Parse input arguments
	data_file_path = args.dataset
	bin_num = args.bin_num
	bin_str = args.bin_str
	output_file_path = args.out
	# asdf verify?
	
	data = read_exp(data_file_path)[:,1:]
	disc = bin(data, bin_str=bin_str, bin_num=bin_num)
	
	mutual = MI(disc)
	
	mi_out(mutual, output_file_path)

# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module
if __name__ == "__main__":
	# Note: this example shows named command line arguments.  See the argparse
	# documentation for positional arguments and other examples.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('dataset',
						help='input gene expression data file path',
						type=str)
	parser.add_argument('--bin_num',
						help='number of bins',
						type=int,
						default=5)
	parser.add_argument('--bin_str',
						help='binning strategy',
						type=str,
						choices={'uniform', 'density'},
						default='uniform')
	parser.add_argument('--out',
						help='MI output file path',
						type=str,
						default='uniform.txt')

	args = parser.parse_args()
	
	main(args)
