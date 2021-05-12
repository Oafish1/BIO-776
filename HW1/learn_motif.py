# Interesting Notes:
# Without the exhaustive subsequencing approach, example2
#	matches nearly exactly -- but shifted 1 to the right.
#	Perhaps the genome sequence in example 2 is actually
#	12+ width?

import argparse, os, sys
import numpy as np

np.random.seed(42)

def readSubsequences(f_name):
	with open(f_name, 'r') as f:
		sequences_raw = f.read().split('\n')
	length = len(sequences_raw[0])
	
	sequences = []
	for i in range(len(sequences_raw)):
		if len(sequences_raw[i]) == length:
			sequences.append(list(sequences_raw[i]))
	return np.array(sequences)
	

# Average expected value for given subsequences, p, and Z
# Notes:
# AEV is calculated in EM to save computation time
# subsequences is assumed to be a 2D numpy array
def AEV(sequences, Z, p, dic, W):
	possible = len(sequences[0])-W+1
	pxt_new = 0
	for i in range(len(sequences)):
		pxt_new_row = 0
		for j in range(possible):
			prod = 1
			# Prior Background
			for k in range(j):
				prod *= p[dic[sequences[i][k]],0]
			# Sequence
			for k in range(W):
				prod *= p[dic[sequences[i][j+k]],k+1]
			# Posterior Background
			for k in range(len(sequences[0])-(j+W)):
				prod *= p[dic[sequences[i][j+W+k]],0]
			
			pxt_new_row += prod * Z[i,j] # Calculate Previous P(X|theta)
		# We take the average instead of multiplying because the computer doesn't like numbers that small
		pxt_new += pxt_new_row 
	pxt_new /= len(sequences)
	
	return pxt_new
	

# Run the EM algorithm
def EM(sequences, W, pseudo=1):
	# Initialize Z and n
	types = np.unique(sequences); possible = len(sequences[0])-W+1
	Z = np.ones([len(sequences),possible]) / possible
	n = np.zeros([len(types),W+1]) # +1 for background
	
	# Embed
	dic = {}
	for i,c in enumerate(types):
		dic[c] = i
	
	# Initialize p
	# Note: Exhaustive subsequencing is done on only a
	# 	subset of the sequences.  Otherwise, the time is
	#	unmanageable for larger W.
	print("Choosing starting p...")
	
	# p_max defaults to all equal values
	p_max=np.ones([len(types),W+1]) / len(types);
	pi=.7;l_max=0;tries={}
	if len(sequences) > 10:
		subset = np.random.choice(range(len(sequences)),10)
		subset = sequences[subset]
	else:
		subset = sequences
	
	for s_num,s in enumerate(subset):
		print(str(1+s_num) + "/" + str(len(subset)), end='\r')
		for st in range(possible):
			if str(s[st:st+W]) not in tries:
				p=np.zeros([len(types),W+1])
				for i,e in enumerate(s[st:st+W]):
					p[dic[e],i+1] = 1
				p = pi*p + (1-pi)*(1-p)/(len(types)-1)
				p[:,0] = np.ones([len(types)]) / len(types)
				l = AEV(sequences,Z,p,dic,W)
				if l > l_max:
					p_max=p;l_max=l
				tries[str(s[st:st+W])] = 1
		print(end='')
	p=p_max
	
	# Main Loop
	print("Running EM algorithm...")
	
	t=0; pxt=999; change=999; goal=1e-40
	while change > goal:
		# E step
		pxt_new = 0
		for i in range(len(sequences)):
			pxt_new_row = 0
			for j in range(possible):
				prod = 1
				# Prior Background
				for k in range(j):
					prod *= p[dic[sequences[i][k]],0]
				# Sequence
				for k in range(W):
					prod *= p[dic[sequences[i][j+k]],k+1]
				# Posterior Background
				for k in range(len(sequences[0])-(j+W)):
					prod *= p[dic[sequences[i][j+W+k]],0]
				
				pxt_new_row += prod * Z[i,j] # Calculate Previous P(X|theta)
				Z[i,j] = prod
			# We take the average instead of multiplying because the computer doesn't like numbers that small
			pxt_new += pxt_new_row 
		pxt_new /= len(sequences)
		
		Z = Z.T / np.sum(Z,1); Z=Z.T
		
		# M step
		for i,c in enumerate(types):
			matches = sequences == c
			for k in range(W,-1,-1):
				if k==0:
					n[i,k] = sum(sum(matches)) - sum(n[i,1:])
				else:
					if k != 1:
						temp = np.roll(matches,-k+1,1)
						temp[:,-k+1:] = False
					temp = temp[:,:possible]
					n[i,k] = sum(sum(temp * Z))
		p = (n+pseudo) / (np.sum(n+pseudo,0))
		
		# Calculate Probability Change
		# Note: This actually does an extra iteration because pxt_new
		# is calculated on p_(t-1).  This is done to ease the
		# computation requirements but is easily changeable with the
		# following line:
		# pxt_new = AEV(Z,p,dic)
		change = abs(pxt_new - pxt)
		pxt = pxt_new
		t += 1
	
	print("Done!")
	return dic,p,Z
	
	
def main(args):
	# Parse input arguments
	# It is generally good practice to validate the input arguments
	seq_file_path = args.sequences_filename
	W = args.width
	model_file_path = args.model
	position_file_path = args.positions
	subseq_file_path = args.subseqs
	
	# Read the File
	sequences = readSubsequences(seq_file_path)
	
	# Run EM
	dic,p,Z = EM(sequences, W)
	
	# Save to Files
	# model
	# https://stackoverflow.com/questions/483666/reverse-invert-a-dictionary-mapping
	rev_dic = {v:k for k,v in dic.items()}
	with open(model_file_path, 'w') as f:
		#print(dic.,file=f)
		for i,row in enumerate(p):
			print(rev_dic[i],end='\t',file=f)
			for e in row:
				print(round(e,3),end='\t',file=f)
			print(file=f)
	
	# positions
	positions = np.argmax(Z,1)
	with open(position_file_path, 'w') as f:
		#print(dic.,file=f)
		for pos in positions:
			print(pos,file=f)
	
	# subsequences
	with open(subseq_file_path, 'w') as f:
		#print(dic.,file=f)
		for i,pos in enumerate(positions):
			print(''.join(sequences[i][pos:pos+W]),file=f)
	

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('sequences_filename',
						help='sequences file path.',
						type=str)
	parser.add_argument('--width',
						help='width of the motif.',
						type=int,
						default=6)
	parser.add_argument('--model',
						help='model output file path.',
						type=str,
						default='model.txt')
	parser.add_argument('--positions',
						help='position output file path.',
						type=str,
						default='positions.txt')
	parser.add_argument('--subseqs',
						help='subsequence output file path.',
						type=str,
						default='subseqs.txt')

	args = parser.parse_args()
	
	main(args)
