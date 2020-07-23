## Generate expression-matched decoy peptides
### Generate a set of decoy peptides whose lengths match the distribution of positive example peptides and come from the source proteins of the positive example peptides.

# Initialize:
import pandas as pd
from Bio import SeqIO
import random
import numpy as np

#Import MS data:
msdata = pd.read_csv('401msdata.csv')

#Import source protein sequences
## Dictionary of the uniprot IDs and their corresponding protein sequences
prot_dict = SeqIO.to_dict(SeqIO.parse('Source Proteins 401.fasta', "fasta"))

## Keep only central 6 unit Uniprot ID in this dictionary
prot_dict_short = dict((key[3:9], value) for (key, value) in prot_dict.items())


# Make a list of possible decoy peptides that can serve as peptide sources:
## -Tile over a protein, with subsequent peptides overlapping by 8 amino acids
## -Select each peptide with a length randomly from length_dist_MS

def generate_decoys(length_dist_MS, seq):
    #'length_dist_MS': lists lengths of peptides in MS data (for size-matching)
    #'seq': sequence of protein working with right now
    
    #protein_length: length of protein in question
    protein_length=len(seq)

    ## randomly select lengths from "length_dist_MS" distribution
    rand_len = np.random.choice(length_dist_MS,protein_length) #(generate extra random_len values)

    ## Tile over a single protein, starting at the first position, where the overlap between subsequent peptides is 8 amino acids
    ### List starting positions of them here.
    ### If last peptide dangles over end, then frame adjust whole thing randomly within protein

    index = 0 #increase by one each time; corresponds to current index position in 'rand_len'
    pointer = 0 #place within protein
    tile_starting_positions = [0] #start at position 0
    tile_ending_positions = []
    end=0
    while pointer < protein_length and end < protein_length:
        end = pointer+rand_len[index]
        tile_ending_positions.append(end)
        pointer = end - 8 #max overlap is 8 amino acids
        tile_starting_positions.append(pointer) 
        index +=1
        
    
    #Remove last positions which are over the edge of the protein:
    if len(tile_ending_positions) > 1: #check that at least 1 peptide found
        tile_starting_positions = tile_starting_positions[0:-2] #remove last pointer position (before 'while' stopped) and 2nd to last (over edge)
        tile_ending_positions = tile_ending_positions[0:-1]

        
        #Is the last peptide starting after end of peptide?
        ## If so, randomly shift window in protein (amount: "rand_start")
        if tile_ending_positions[-1]>=protein_length:
            window = protein_length - tile_ending_positions[-1] #extra number of amino acids that can shift peptides

            rand_start = np.random.choice(range(0,window+1)) #random place to start the peptides

            tile_starting_positions += rand_start
            tile_ending_positions += rand_start

    else:
        tile_starting_positions = []
        tile_ending_positions = []

    #generate decoy peptides
    n_decoys = len(tile_starting_positions)
    decoys = ['']*n_decoys
    for i in range(0,n_decoys):
        decoys[i] = str(seq[tile_starting_positions[i]:tile_ending_positions[i]])

    return decoys


#Call 'generate_decoys', generating possible decoys from the source protein
possible_decoys = list()
for i in range(0,len(msdata['uniprot'])):
    #seq:            sequence of the source protein
    #ms_seq:         sequence of the ms peptide
    #length_dist_MS: length distribution of MS peptides
    
    seq = prot_dict_short[msdata['uniprot'][i]].seq
    
    length_dist_MS = msdata['sequence'].str.len()
    
    new_decoys = generate_decoys(length_dist_MS, seq)
    possible_decoys = possible_decoys + new_decoys

possible_decoys_df = pd.DataFrame({'decoys':possible_decoys})
possible_decoys_df['len']=possible_decoys_df['decoys'].str.len()


#Select peptides from the list of possible decoys 
## -Start with the MS-identified minimum epitope cores 
## -Select sequence such that decoys being added don't have overlapping 9mer sequences with peptides already on the list

final_set = msdata[['sequence','uniprot']]

# Check if 'singlet' contains any 9mers which overlap with the peptides in df (a list)
def check_9mers(singlet, df):
    #generate all 9mers from singlet:
    all_9mers = []
    n_9mers = len(singlet) - 8
    if n_9mers > 0:
        for x in range(0,n_9mers):
            all_9mers.append(singlet[x:x+9])
            
    #check all 9mers against df:
    overlap = 0
    for i in all_9mers:
        if df.str.contains(i).any():
            overlap +=1
    
    if overlap == 0:
        found_flag = 1 #this singlet works! (no overlap in 9mers)
    else: 
        found_flag = 0
    
    return found_flag


#randomly select one decoy to length match, making sure it doesn't contain overlapping 9mers:
fold_decoys = [1,18] #won't reset - so is 1, 1+18
for x in fold_decoys:
    for y in range(0,x):
        for i in length_dist_MS:
            found_flag = 0
            while found_flag == 0:
                rand_one = np.random.choice(possible_decoys_df[possible_decoys_df['len']==i]['decoys'],1)[0]
                found_flag = check_9mers(rand_one, final_set['sequence'])
            new= pd.DataFrame({'sequence':[rand_one], 'uniprot':'decoy'})
            final_set = final_set.append(new)
    fold_name = np.sum(fold_decoys[0:x])
    final_set.to_csv('401ms+'+str(fold_name)+'xdecoys.csv')
        



