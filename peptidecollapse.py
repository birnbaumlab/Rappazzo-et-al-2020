from _collections import defaultdict

def collapse(strings, threshold) :
	sums = defaultdict(int)
	for string in strings:
		sums[string] += 1
	out = {key:sum for key, sum in sums.items() if sums[key] >= threshold}
	return out

print 'Building'

lib =[]
seq =[]
count =[]

with open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect peptides filtered.fasta','r') as myfile:
	lines = myfile.readlines()
	for i,line in enumerate(lines):
		if i%2 == 1:
			lib.append(line)
	print 'Collapsing'
	sequences = collapse(lib,1)

for sequence, sum in sequences.items():
	var = sequence.replace('\n','')
	seq.append(var)
	count.append(sum)

fastafile = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect peptides filtered grouped.fasta','w')

for i,line in enumerate(seq):
	fastafile.write(">Clone {} count:{}\n".format((i+1),count[i]))
	fastafile.write(seq[i] + '\n')

print len(seq)

file = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect peptides filtered grouped.csv','w')

for i,line in enumerate(seq):
	file.write("{}\t{}\t{}\n".format((i+1),seq[i],count[i]))