# DFLpred - High-throughput prediction of disordered flexible linker regions in protein sequences  
DFLpred is designed for the prediction of disordered flexible linker (DFL). The method generates numeric score for each residue in the input protein sequence that quantifies putative propensity to form a DFL. Larger value of propensity denotes higher likelihood to form DFL. It also provides putative binary annotations (a given residue is predicted either as a DFL or not a DFL) based on false positive rate = 0.05 using threshold on the propensity score = 0.18). Residues with propensity > 0.18 are assumed to form DFLs and otherwise they are assumed not to form DFLs.  

## System requirement
Java Runtime Environment 1.5 or later  
Any ANSI C compiler e.g. GUN C compiler to compile "myiupred.c" (if "myiupred" is not runnable on your  computer)

## Usage (command line)
java -jar DFLpred.jar [input_fasta_file] [output_file]  
Accepts up to 5000 FASTA sequences in a single fasta file

## Example
java -jar DFLpred.jar examples.fasta results.txt

## Format of results
Line 1: >protein ID  
Line 2: protein sequence using 1-letter amino acid encoding. Lower/upper case indicates a residue is/is not predicted as a flexible disordered linker (DFL) by the score threshold 0.18  
Line 3: comma-separated propensity score of each residue to be a DFL

## DFLpred webserver
http://biomine-ws.ece.ualberta.ca/DFLpred

## Citation
Meng, F. and Kurgan, L. DFLpred: High-throughput prediction of disordered flexible linker regions in protein sequences. Bioinformatics 2016;32(12):i341-i350

## Acknowledgments
We acknowledge with thanks the following software used as a part of DFLpred:  
IUPred (http://iupred.enzim.hu/) - Prediction of Intrinsically Unstructured Proteins  
AAindex (http://www.genome.jp/aaindex/) - Amino acid indices, substitution matrices and pair-wise contact potentials 

## Disclaimer and license agreement
http://biomine-ws.ece.ualberta.ca/Disclaimer.html

## Biomine home page
http://biomine.ece.ualberta.ca/