import gzip
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def get_arguments():
	parser = argparse.ArgumentParser(description="Arguments the user wishes to pass, including file and k-mer length")
	parser.add_argument( "-f", "--filename", help="The name of the file you wish to input",\
						required=True, type=str)
	return parser.parse_args()

args = get_arguments()
file = args.filename
file2 = file[38:-9]
index_indentify = file[38:]

def convert_phred(letter):
	"""Converts a single character into a phred score"""
	phred = ord(letter) - 33
	return(phred)

def populate_array(file):
	'''Opens file which loops through every FASTQ record and coverts
	Phred score to corresponding number and then sums the Phred score for each base pair.
	Postion 0 will be the average Phred score of the first nucleotide in the read and so on until poisiton 101
	'''
	LN = 0
	with gzip.open(file, "rt") as fh:
		if index_indentify == "1294_S1_L008_R2_001.fastq.gz" or index_indentify == "1294_S1_L008_R3_001.fastq.gz":
			mean_scores = np.zeros(8, dtype = float) #index array
			for line in fh:
				LN += 1
				if LN%4==0: # Quality score line
					for char in range(len(line)-1): #takes char and converts to phred scroe along each position of line
						mean_scores[char] += convert_phred(line[char])

		else:
			mean_scores = np.zeros(101, dtype = float) # seqeucne array
			for line in fh:
				LN += 1
				if LN%4==0: # Quality score line
					for char in range(len(line)-1): #takes char and converts to phred scroe along each position of line
						mean_scores[char] += convert_phred(line[char])
		return mean_scores, LN


mean_scores, LN = populate_array(file)

for score in range(len(mean_scores)):
	mean_scores[score] /= (LN/4)

plt.plot(mean_scores)
plt.xlabel("Position")
plt.ylabel("Mean Quality Score")
plt.savefig(file2 + ".png")
