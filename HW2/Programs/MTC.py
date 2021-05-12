import argparse, os, sys
import numpy as np
from scipy import stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# random seed set for reproducibility
# DO NOT CHANGE.
np.random.seed(2020)

def t_statistics(g1, g2):
	n1, n2 = g1.shape[1], g2.shape[1]
	mu1 = np.mean(g1, axis=1, keepdims=True)
	mu2 = np.mean(g2, axis=1, keepdims=True)
	var1 = np.sum((g1 - mu1)**2, axis=1, keepdims=True) / (n1 - 1)
	var2 = np.sum((g2 - mu2)**2, axis=1, keepdims=True) / (n2 - 1)
	t = (mu2 - mu1) / np.sqrt(var1/n1 + var2/n2)
	return t

# This is the main function provided for you.
# Define additional functions to implement mutual information.
def main(args):
	# Parse input arguments.
	# It is generally good practice to validate the input arguments, e.g.,
	# verify that the input and output filenames are provided.
	input_file_path = args.counts
	procedure = args.procedure
	alpha = args.alpha
	lamb = args.lamb
	fig_path = args.fig

	brca1, brca2 = 7, 8 # two groups: BRCA1 tumors vs. BRCA2 tumors
	C = brca1 + brca2   # number of samples
	B = 100			 # number of permutations

	# Construct input file path.
	dirname = os.path.dirname(__file__)
	fname = os.path.join(dirname, input_file_path)
	
	# Loads log2-counts from input file.
	data = np.loadtxt(fname, 
		delimiter=',', skiprows=1, usecols=range(1, C+1))
	N = data.shape[0]   # number of genes

	# Performs permutation test to obtain p-values.
	# (Permutation test is a nonparametric, distribution-free test.)
	t = t_statistics(data[:, :brca1], data[:, brca1:])
	abs_t = np.abs(t)
	counts = np.zeros((N,))
	for i in range(B):
		per_data = data[:, np.random.permutation(C)]
		per_t = t_statistics(per_data[:, :brca1], per_data[:, brca1:])
		abs_per_t = np.reshape(np.abs(per_t), (N,))
		counts += np.sum(abs_per_t >= abs_t, axis=1)
	pvals = counts / (N * B)

	# Sorts p-values in ascending order.
	pvals = pvals[np.argsort(pvals)]

	# TODO: implement the Bonferroni correction procedure.
	def bonferroni(pvals, alpha=0.05):
		# Calculate beta
		beta = 1-(1-alpha)**(1/len(pvals))
		
		# Return indices
		return np.argwhere(pvals <= beta)

	# TODO: implement the Benjamini-Hochberg correction procedure.
	def benjamini_hochberg(pvals, alpha=0.05):
		# Calculate max k
		# Note that, even if pval > (k/len(pvals))*alpha, we do not 'break'
		k_max = -1
		for k, pval in enumerate(pvals):
			if pval <= (k/len(pvals))*alpha:
				k_max = k
				
		return range(k_max+1)

	# TODO: implement the Storey-Tibshirani correction procedure.
	def storey_tibshirani(pvals, lamb, alpha=0.05):
		# Calculate q values
		# Note, we can use k+1 in the denominator.  Since q-value decreases with
		# the index and we only track the maximal p-value threshold, the calculation
		# #{P_i <= t} is still done and is the minimal q-value for any p-value
		# grouping.
		k_max = -1
		for k,pval in enumerate(pvals):
			qval = (lamb*len(pvals)*pval)/(k+1)
			if qval <= alpha:
				k_max = k
		
		return range(k_max+1)

	# Call the user-specified correction procedure to find significant genes.
	if procedure == 'bf':
		pos_idx = bonferroni(pvals, alpha)
	elif procedure == 'bh':
		pos_idx = benjamini_hochberg(pvals, alpha)
	else:
		##### COMMENT OUT THIS CODE AFTER LAMBDA HAS BEEN ESTIMATED #####
		# Plots p-value histogram for estimation lambda
		"""
		fig, pl = plt.subplots(1, 1)
		nbins = 20
		(weights, bins) = np.histogram(pvals, bins=nbins)
		pl.hist(bins[:-1], bins=bins, weights=weights/(len(pvals)/nbins), color='white', ec='black')
		plt.xlabel("p-values")
		plt.ylabel("Sample uniform bin counts")
		plt.savefig(fig_path, dpi=512)
		exit()
		"""
		##### COMMENT OUT THIS CODE AFTER LAMBDA HAS BEEN ESTIMATED #####
		
		pos_idx = storey_tibshirani(pvals, lamb, alpha)

	def plot():
		# Colors the markers
		facecolors = np.array(['white'] * N)
		facecolors[pos_idx] = 'red'

		caps = [sum(pvals<0.002), sum(pvals<0.01)] # number of points plotted
		gaps = [cap//5 for cap in caps]   # gap size between x-ticks
		margins = [1e-4, 1e-3]

		fig, axes = plt.subplots(2, 1)
		for i in range(2):
			cap, gap, margin = caps[i], gaps[i], margins[i]
			x = list(range(1, cap+1))   # x-axis: indices
			y = pvals[:cap]			 # y-axis: p-values
			
			# Generates a scatter plot
			axes[i].scatter(x, y, s=20, 
				facecolors=facecolors[:cap], edgecolors='gray')
			axes[i].set_xticks(list(range(0, cap+1, gap)))
			axes[i].set_ylim([-margin, y.max()+margin])
			axes[i].set_ylabel('p-value')

		plt.xlabel('gene index')
		plt.savefig(fig_path, dpi=512)
		plt.show()

	plot()

# Note: this syntax checks if the Python file is being run as the main program
# and will not execute if the module is imported into a different module.
if __name__ == "__main__":
	# Note: this example shows named command line arguments.  See the argparse
	# documentation for positional arguments and other examples.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('counts',
						help='log2 gene counts from two biological conditions',
						type=str)
	parser.add_argument('--procedure',
						help='multiple testing correction procedure',
						type=str,
						choices={'bf', 'bh', 'st'},
						default='bf')
	parser.add_argument('--alpha',
						help='significance threshold',
						type=float,
						default=0.05)
	parser.add_argument('--lamb',
						help='estimated lambda used by Storey-Tibshirani',
						type=float,
						default=1.0)
	parser.add_argument('--fig',
						help='path to save the plot',
						type=str,
						default='fig.png')

	args = parser.parse_args()
	# Note: this simply calls the main function above, which we could have
	# given any name.
	main(args)
