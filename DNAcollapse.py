#code to collapse and quantify exact match peptides in DNA space

from _collections import defaultdict

# define collapsing function

def collapse(strings, threshold):
	sums = defaultdict(int)
	for string in strings:
		sums[string] += 1
	out = {key:sum for key, sum in sums.items() if sums[key] >= threshold}
	return out

print 'Building'

# initialize lists

lib =[]
seq_unsort = []
count =[]

# open DNA fasta file

with open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK pepdna.fasta','r') as myfile:
	lines = myfile.readlines()
	for i,line in enumerate(lines):
		# isolate peptide lines
		if i%2 == 1:
			lib.append(line)

	# collapse all peptides

	print 'Collapsing'
	sequences = collapse(lib,1)

# write collapsed values to lists

for sequence, sum in sequences.items():
	var = sequence.replace('\n','')
	seq_unsort.append(var)
	count.append(sum)

sort_inds = sorted(range(len(count)), key=lambda k: count[k], reverse = True)

seq = []

for i in range(len(seq_unsort)):
	seq.append(seq_unsort[sort_inds[i]])

count.sort(reverse = True)

print 'Sorted'
print len(seq)

# define the Hamming distance

def hamdist(str1, str2):
	diffs = 0
	for ch1, ch2 in zip(str1, str2):
		if ch1 != ch2:
			diffs += 1
	return diffs

seq_d = []
count_d = []

# Identify sequences that are hamming distance 1 away from a more prevalent sequence, 2 away from a sequence 100 times
# more prevalent than itself or 3 away from 10k greater

for i in range(len(seq)):
	if count[i] > 1:
		if i%100 == 0:
			print i

		for ii in range((i+1),len(seq)):
			ratio = float(count[i])/count[ii]
			if ratio > 1:
				hd = hamdist(seq[i],seq[ii])
				if hd < 4:
					if hd == 1:
						if seq[ii] not in seq_d:
							seq_d.append(seq[ii])
							count_d.append(count[ii])
					elif hd == 2 and ratio >= 100:
						if seq[ii] not in seq_d:
							seq_d.append(seq[ii])
							count_d.append(count[ii])
					elif hd == 3 and ratio >= 10000:
						if seq[ii] not in seq_d:
							seq_d.append(seq[ii])
							count_d.append(count[ii])

print len(seq_d)
print count_d

seq_out = []
count_out = []

for i in range(len(seq)):
	if seq[i] not in seq_d:
		seq_out.append(seq[i])
		count_out.append(count[i])

# output to CSV and fasta files:

print len(seq_out)

fastafile = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect pepdna grouped.fasta','w')

for i,line in enumerate(seq_out):
	fastafile.write(">Clone {} count:{}\n".format((i+1),count_out[i]))
	fastafile.write(seq_out[i] + '\n')

expfastafile = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect pepdna.fasta','w')

for i in range(len(seq_out)):
	for ii in range(count_out[i]):
		expfastafile.write(">Clone {} copy {}\n".format((i+1),(ii+1)))
		expfastafile.write(seq_out[i] + '\n')

file = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect pepdna grouped.csv','w')

for i,line in enumerate(seq_out):
	file.write("{}\t{}\t{}\n".format((i+1),seq_out[i],count_out[i]))
