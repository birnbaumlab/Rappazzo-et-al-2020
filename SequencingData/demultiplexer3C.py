#Code to demultiplex DR401 NNK1Reg Library with strict filters for length, intact 3C sequence, peptide flanks, phred score, and NNK codons

# Establish strict parameter values
ID_3C = 'CTGGAAGTTCTGTTCCAGGGGCCC'
flank1 = 'GCTGCC'
flank2 = 'TGGGAAGAAGGT'
correct_len = 221

# Establish accepted Phred (>19) scores and barcodes
BClist = ['ATCACG','CAGATC','ACTTGA','GATCAG','TAGCTT','GGCTAC','CGATGT','TTAGGC','TGACCA','ACAGTG','GCCAAT']
Qlist = ['5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K']

# Initialize peptide groupings
peptides0 = []
peptides1 = []
peptides2 = []
peptides3 = []
peptides4 = []
peptides5 = []
peptides6 = []
peptides7 = []
peptides8 = []
peptides9 = []
peptides10 = []

count_pep = 0
count_bc0 = 0
count_bc1 = 0
count_bc2 = 0
count_bc3 = 0
count_bc4 = 0
count_bc5 = 0
count_bc6 = 0
count_bc7 = 0
count_bc8 = 0
count_bc9 = 0
count_bc10 = 0

# Demultiplexing

with open('171201 DR401 NNK1Reg overlaps.fastq','r') as myfile:
	print 'Im open!'
	lines = myfile.readlines()
	print 'Ive read'
	for i,line in enumerate(lines):
		# Single out line containing DNA sequence
		if i%4 == 1:
			# Filter for peptides with correct length and intact 3C sequence
			if len(line) == correct_len and ID_3C in line:

				# Establish peptide, barcode, and score for read
				pep = line[41:86]
				BC = line[8:14]
				score_line = lines[i+2]
				scores = score_line[41:86]

				# reset count variables
				count_nnk = 0
				count_phred = 0

				# count number of codons that follow NNK
				for ii in range(9):
					if pep[3*ii+8] == 'G' or pep[3*ii+8] == 'T':
						count_nnk += 1

				# count number of bases with phred >19
				for ii in range(len(scores)):
					if scores[ii] in Qlist:
						count_phred +=1

				# only collect peptides with correct flanks and with correct NNK and phred counts

				if pep[0:6] == flank1 and pep[33:45] == flank2 and count_nnk == 9 and count_phred == 45 and BC in BClist:
					# group peptides 
					count_pep += 1
					if BC == 'ATCACG':
						count_bc0 += 1
						peptides0.append('>' +lines[i-1])
						peptides0.append(pep + '\n')
					elif BC == 'CAGATC':
						count_bc1 += 1
						peptides1.append('>' +lines[i-1])
						peptides1.append(pep + '\n')
					elif BC == 'ACTTGA':
						count_bc2 += 1
						peptides2.append('>' +lines[i-1])
						peptides2.append(pep + '\n')
					elif BC == 'GATCAG':
						count_bc3 += 1
						peptides3.append('>' +lines[i-1])
						peptides3.append(pep + '\n')
					elif BC == 'TAGCTT':
						count_bc4 += 1
						peptides4.append('>' +lines[i-1])
						peptides4.append(pep + '\n')
					elif BC == 'GGCTAC':
						count_bc5 += 1
						peptides5.append('>' +lines[i-1])
						peptides5.append(pep + '\n')
					elif BC == 'CGATGT':
						count_bc6 += 1
						peptides6.append('>' +lines[i-1])
						peptides6.append(pep + '\n')
					elif BC == 'TTAGGC':
						count_bc7 += 1
						peptides7.append('>' +lines[i-1])
						peptides7.append(pep + '\n')
					elif BC == 'TGACCA':
						count_bc8 += 1
						peptides8.append('>' +lines[i-1])
						peptides8.append(pep + '\n')
					elif BC == 'ACAGTG':
						count_bc9 += 1
						peptides9.append('>' +lines[i-1])
						peptides9.append(pep + '\n')
					elif BC == 'GCCAAT':
						count_bc10 += 1
						peptides10.append('>' +lines[i-1])
						peptides10.append(pep + '\n')

# write peptides out

pepfile0 = open('DR401 NNK1Reg Library R0 No UMI phredandNNK pepdna.fasta','w')

for i in peptides0:
	pepfile0.write("%s" % i)

pepfile1 = open('DR401 NNK1Reg Library R1 3C Only No UMI phredandNNK pepdna.fasta','w')

for i in peptides1:
	pepfile1.write("%s" % i)

pepfile2 = open('DR401 NNK1Reg Library R2 3C Only No UMI phredandNNK pepdna.fasta','w')

for i in peptides2:
	pepfile2.write("%s" % i)

pepfile3 = open('DR401 NNK1Reg Library R3 3C Only No UMI phredandNNK pepdna.fasta','w')

for i in peptides3:
	pepfile3.write("%s" % i)

pepfile4 = open('DR401 NNK1Reg Library R4 3C Only No UMI phredandNNKpepdna.fasta','w')

for i in peptides4:
	pepfile4.write("%s" % i)

pepfile5 = open('DR401 NNK1Reg Library R5 3C Only No UMI phredandNNK pepdna.fasta','w')

for i in peptides5:
	pepfile5.write("%s" % i)

pepfile6 = open('DR401 NNK1Reg Library R1 3CDM No UMI phredandNNK pepdna.fasta','w')

for i in peptides6:
	pepfile6.write("%s" % i)

pepfile7 = open('DR401 NNK1Reg Library R2 3CDM No UMI phredandNNK pepdna.fasta','w')

for i in peptides7:
	pepfile7.write("%s" % i)

pepfile8 = open('DR401 NNK1Reg Library R3 3CDM No UMI phredandNNK pepdna.fasta','w')

for i in peptides8:
	pepfile8.write("%s" % i)

pepfile9 = open('DR401 NNK1Reg Library R4 3CDM No UMI phredandNNK pepdna.fasta','w')

for i in peptides9:
	pepfile9.write("%s" % i)

pepfile10 = open('DR401 NNK1Reg Library R5 3CDM No UMI phredandNNK pepdna.fasta','w')

for i in peptides10:
	pepfile10.write("%s" % i)

print BC
print pep
print scores
print count_nnk
print count_phred
print count_pep
print count_bc0
print count_bc1
print count_bc2
print count_bc3
print count_bc4
print count_bc5
print count_bc6
print count_bc7
print count_bc8
print count_bc9
print count_bc10
