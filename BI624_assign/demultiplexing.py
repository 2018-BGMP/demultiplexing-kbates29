import argparse
import numpy as np
import gzip


def get_arguments():
	parser = argparse.ArgumentParser(description= "Arguments the user wishes to pass")
	parser.add_argument( "-f", "--file",
	help="The name of the files you wish to demultiplex as well as the file containing all indexes used. Files must be ordered 'Read 1' 'Index 1' 'Index 2' 'Read 2' 'All_Indexes'",
	required=True, type=str, nargs='+')
	return parser.parse_args()

args = get_arguments()
files = args.file
read_1 = files[0]
index_1 = files[1]
index_2 = files[2]
read_2 = files[3]
all_index = files[4]

def reverse_complement(seq):
	'''
	Takes sequence and returns the reverse complement
	'''
	seq = seq[::-1]
	intab = "ATCG"
	outtab = "TAGC"
	transtab = str.maketrans(intab, outtab)
	seq = seq.translate(transtab)
	return seq

def convert_phred(letter):
	"""Converts a single character into a phred score"""
	phred = ord(letter) - 33
	return(phred)

index_hopping = 0 # keeps track amount of index hopping
index_dict = {} # dictionary which tracks indexes and the number of reads which are properly sorted
low_qs = 0 # counts number of indexes with low quality scores
count = 0 # counts number of lines in file
index_list = []

#looks through file to find indexes and adds to index_list
with open(all_index, "r") as indexes:
	for line in indexes:
		line = line.strip("\n").split("\t")
		index_list.append(line[4])

index_list = index_list[1:] # removes header from list
files_dict = {}
# loops through index_list to create and open good forward and reverse files for each index
for i in index_list:
	files_dict[i] = open(i + "_Forward_Read.fastq", "a")
	files_dict[reverse_complement(i)] = open(i + "_Reverse_Read.fastq", "a")


with gzip.open(read_1, "rt") as r1:
	with gzip.open(read_2, "rt") as r2:
		with gzip.open(index_1, "rt") as i1:
			with gzip.open(index_2, "rt") as i2:
				with open("Bad_Forward_Read.fastq", "a") as bad_f:
					with open("Bad_Reverse_Read.fastq", "a") as bad_r:
						with open("Bad_QC_Forward.fastq", "a") as bad_qc_f:
							with open("Bad_QC_Reverse.fastq", "a") as bad_qc_r:
								while True:
									header1 = r1.readline()
									seq1 = r1.readline()
									plus1 = r1.readline()
									qs1 = r1.readline()
									header2 = i1.readline()
									seq2 = i1.readline()
									plus2 = i1.readline()
									qs2 = i1.readline()
									header3 = i2.readline()
									seq3 = i2.readline()
									plus3 = i2.readline()
									qs3 = i2.readline()
									header4 = r2.readline()
									seq4 = r2.readline()
									plus4 = r2.readline()
									qs4 = r2.readline()
									if len(header1)==0: # used to signify that end of  fastq file has been reached and will break out of loop
										break
									header1 = header1.strip("\n") # stips new line and appends index to end of header line
									header1 = header1 + ":" + seq2
									header4 = header4.strip("\n")
									header4 = header4 + ":" + seq3
									if "N" in seq2 or "N" in seq3: # removes reads that have N in index. Considered index hopping
										bad_f.write(header1)
										bad_f.write(seq1)
										bad_f.write(plus1)
										bad_f.write(qs1)
										bad_r.write(header4)
										bad_r.write(seq4)
										bad_r.write(plus3)
										bad_r.write(qs4)
										index_hopping += 1
									else:
										#finds mean quality score for index read and then filters if below 30
										# write to bad quality control file
										phred_r1 = np.zeros(len(qs2)-1, dtype = float)
										phred_r2 = np.zeros(len(qs3)-1, dtype = float)
										for char in range(len(qs2)-1):
											phred_r1[char] = convert_phred(qs2[char])
										for char in range(len(qs2)-1):
											phred_r2[char] = convert_phred(qs3[char])
										phred_r1_mean = np.mean(phred_r1)
										phred_r2_mean = np.mean(phred_r2)
										if phred_r1_mean < 30 or phred_r2_mean < 30:
											bad_qc_f.write(header1)
											bad_qc_f.write(seq1)
											bad_qc_f.write(plus1)
											bad_qc_f.write(qs1)
											bad_qc_r.write(header4)
											bad_qc_r.write(seq4)
											bad_qc_r.write(plus4)
											bad_qc_r.write(qs4)
											low_qs += 1
										else:
											# compares indexes to one another
											seq3_reverse = reverse_complement(seq3[:-1])
											if seq3_reverse == seq2[:-1] and seq2[:-1] in index_list: # ensures that each sequence is in index list
												files_dict[seq2[:-1]].write(header1)
												files_dict[seq2[:-1]].write(seq1)
												files_dict[seq2[:-1]].write(plus1)
												files_dict[seq2[:-1]].write(qs1)
												files_dict[seq3[:-1]].write(header4)
												files_dict[seq3[:-1]].write(seq4)
												files_dict[seq3[:-1]].write(plus4)
												files_dict[seq3[:-1]].write(qs4)
												if seq2[:-1] in index_dict: # keeps track of the number of reads which enter a given index file
													index_dict[seq2[:-1]] += 1
												else:
													index_dict[seq2[:-1]] = 1
											else:
												# writes to bad file if indexes are not in index_list
												# ensures that matching does not occur due to sequening error
												# considered index hopping
												bad_f.write(header1)
												bad_f.write(seq1)
												bad_f.write(plus1)
												bad_f.write(qs1)
												bad_r.write(header4)
												bad_r.write(seq4)
												bad_r.write(plus3)
												bad_r.write(qs4)
												index_hopping +=1
									count += 4
for i in files_dict: # closes all files
	files_dict[i].close()

#prints quality metrics including number of times index hopping occured and % reads that were removed due bad QC and % of reads in each file
print("Overall index hopping: " + str(index_hopping))
print("Percentage of reads removed due to low mean quality score: " + str((low_qs/(count/4))*100))
for key in index_dict:
	print("Percentage of reads from " + key + " : " + str((index_dict[key]/(count/4))*100))
