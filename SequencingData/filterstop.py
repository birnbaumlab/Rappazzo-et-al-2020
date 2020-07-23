out = []
clones = []

with open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect peptides.fasta','r') as myfile:
	lines = myfile.readlines()
	for i,line in enumerate(lines):
		if i%2 == 1:
			if '*' not in line:
				out.append(lines[i-1])
				out.append(lines[i])

file = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK PCRcorrect peptides filtered.fasta','w')

for i in out:
	file.write("%s" % i)

print len(out)/2
