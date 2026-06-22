# CSE190 Lehigh in Ireland Research

This document is meant for navigation and notes throughout this research project.

## Phase 1

In this phase of the project, we read research papers and navigated nih.gov through the nucleotide database
in order to get the percentage of As, Ts, Gs, Cs and the separate percentages of AT and GC content in 
S. aureus and P. aeringuosa genomic sequeces. 

## Phase 2

In this phase, we are looking at antibiotic resistance genes to see if they contain certain genomes,
and if they do, what is the GCAT sequence, the first and last 3 bases, and the percentages of AT and GC ratio of 
that specific sequence. 

Specifically, we are looking at:
- hpt
- apt
- adk
- tdk
- pstc

## Phase 3

We have now moved on to analysing and utilizing data that the lab has collected. We were given a list of genes and asked to run
our code on it. This will hopefully be a check to see if the code written previously is working correctly (aka if the bases change
from ATG an TAA). The code written in this phase will be tested on *one* S.a. / P.a. organism to begin with, looking at multiple genes.
It is important to note that Zoe was tasked with looking at S. aureus and Reid at P. aeringuosa.

## Which Files are Which
- **genome_id_fetcher.py** -> gets the ids of both SA and PA from the nucleotide database and outputs PA_nucleotide_genomes.txt and SA_nucleotide_genomes.txt
- **read_PA_seq.py** -> using the nucleotide database and PA_nucleotide_genomes.txt, outputs PA_local_ids_summary.csv which contains: Accession,Description,Total_Bases,A_pct,T_pct,G_pct,C_pct,GC_pct,AT_pct
- **read_SA_seq.py** -> using nucleotide and SA_nucleotide_genomes.txt, outputs SA_local_ids_summary.csv which contains: Accession,Description,Total_Bases,A_pct,T_pct,G_pct,C_pct,GC_pct,AT_pct
- **SA_specific_genes.py** -> takes in antibiotic resistance genes MANUALLY and 

