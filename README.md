# IBMScore
 IAP Binding Motif (IBM) scoring algorithm for identifying IBMs present in protein and peptide sequences.

## IBMScore Description
IBMScore using a position weighted scoring matrix derived from data in Kurakin et. al. J. Mol Recognit (2007), 
Eckelman et. al. Cell. Death Differ. (2008), and Chen et. al. biorxiv (2023). A full description of the methods
used to develop this scoring algorithm can be found in XXXXXX.

## Installation Guide
__Programming Language:__ Python
This code was specifically written in Python 3.7.5

__Required Python Packages:__
None

## Running IBMScore
IBMScore can be run in the command-line, where inputs are specified using an arguement parser. An example run is:
```
run IBMScore.py -in IgA.fasta
```
Each of the parser flags are described below and can also be displayed using -h:
```
-in   --Input_FASTA_File  Required  Name of the text file containing the FASTA sequence of the protein of interest. Carot should not be in same line as sequence, UniProt format preferred.
-out  --Output_File       Optional  Name of the text file that will contain the results. Use when you want the output saved to a file instead of printing to the terminal.
-cut  --Cutoff_Value      Optional  Cutoff value above which a peptide is considered to be an IBM. The default is 0.5.
-pwm  --PWSM              Optional  Supply a specific position weighted score matrix for score calculation.
-all  --Return_All        Optional  Use of this flag will return all scored peptides instead of just peptides above the cutoff.
```