# HLA-DR401
This repository is a collection of scripts utilized in "Repertoire-scale determination of class II MHC peptide binding via yeast display improves antigen prediction".

## Processing Sequencing Data
The in-house scripts comprising our pipeline are:
- Demultiplexer3C.py - Identify, filter, and demultiplex peptide-containing reads. 
	- Ensure 3C fidelity and length
	- Accept peptide DNA with only constant flanking regions, all scores Phred >= 20
	- Accept peptide DNA with NNK format
- DNAcollapse.py - Group and enumerate peptide DNA sequences, account for PCR and sequencing errors
	- Remove sequences too close to more prevalent sequences as likely daughter sequences
- Filterstop.py - Remove peptides with stop codons
- Peptidecollapse.py - Group and enumerate peptides
	- Groups peptides with same amino acid sequence but different peptide DNA encoding

## NNAlign-Trained Models
We trained the NNAlign-2.0 algorithm (Morten & Andreatta, 2017) on yeast display (YD) and mass spectrometry (MS) data. The resulting model files can be uploaded to NNAlign to perform additional predictions. These model files are included in this repository:
- 401_9mer_YD-trained_model.txt
- 401_13mer_YD-trained_model.txt
- 401_MS-trained_model.txt
- 402_9mer_YD-trained_model.txt
- 402_MS-trained_model.txt

## Generating Test Sets
We generated test sets of peptides which included expression-matched decoy peptides. An example script for generating these tests sets is included in this repository:
- Expression_Matched_Decoy_Peptides_Generator.py

## Code Authors
Brooke D. Huisman & C. Garrett Rappazzo, Birnbaum Lab
