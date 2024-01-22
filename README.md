# Guide to the REP_finder.m (REP-like) and Gene_Location.m (general) Programs
---
## REP_finder.m

![Program pipeline and workflow of parsing and analyzing REP sequences](https://i.postimg.cc/xTMV9YM5/Program-Flowchart.png)

The REP_finder.m is a program to search for REP-like repeats (palindromic inverted repeats with a defined 6-base pair palindrome region based on REP consensus sequences) requires the following additional MATLAB-based programs:
1.	GetChromoAccNumfromAssemSearch.m: a program that obtains an accession number list that will be used in the REP_finder.m. The list is obtained from the NCBI website based on search terms input in the “AssAcc” command line. Requires the following programs to run:
    a.	getAccNumfromSearch.m
    b.	getAssemblyInfo2.m
    c.	getAssNamefromAccNum.m
    d.	getPhylum.m
2.  downloadGenomeKEH.m: a program that downloads the genome sequence and feature table based on the input accession number.
3.	palindrome_flip.m: a function that produces an inverted complimentary sequence based on an input sequence. The function will first invert the order of the sequence and then replace the inverted sequence with the inverted complimentary sequence. This function is used in the regular expression to produce the palindromic sections. 
4.	nestedSortStruct.m: a function provided by Jake Hughey in the MATLAB file exchange that sorts data structures (used to organize data obtained in this program). The data structure rows are organized in ascending or descending order of a parameter of the structure, which will be input into the function. Requires: 
    a.	sortStruct.m
5.	nucleotide2dec.m: a function that converts a gene sequence into a decimal number code. The conversion is done by first converting each base in the sequence into a base-4 digit by the following rules: G = 0; A = 1; T = 2; C = 3. The resultant base-4 number will then be converted into the base-10 number code used to represent the original sequence. This function is used to define the types of palindromic sequences present in matches. It allows easier display of the types of palindrome sequences in the one-line table.
The REP_finder.m first obtains the accession number (genome) list for the program run. For each genome, the gene sequence and feature table are then downloaded. This program allows for continuation of previous runs, in case of any interruptions. The program will load any preexisting saved data from a previous run to continue the run. To keep track, a REP_Counter file is built, which contains the number of genomes processed in a previous run. If this number is below the number of genomes in the accession number list, the following number can be input to continue the program run. (For example, if the counter file displays ‘42’, then when continuing the run, the number ‘43’ should be input alongside the program name. The program can also run without inputting a genome number; the program will just start from the beginning of the genome list.) The program will then proceed to find matches based on the regular expression:
['(G[CT]C[CGT]GA).{' num2str(Min_Middle_Section) ',' num2str(Max_Middle_Section) '}(??@palindrome_flip($1))'].
Two parameters can be input for this expression, with both affecting the number base pairs that can be allowed in the region of the matched sequence between the palindrome regions (termed the ‘middle section’.) Min_Middle_Section defines the minimum number of base pairs in the middle section; Max_Middle_Section defines the maximum.

After finding the matches, the following parameters are included in an excel file (.csv and .xlsx due to potential bugs) for each genome:
1.	Accession Number
2.	Genome Name
3.	StartIndex: base pair number of genome where matched sequence starts
4.	EndIndex: base pair number of genome where matched sequence ends
5.	Total_Seq: matched sequence
6.	Total_Length: base pair length of matched sequence (last row displays average base pair length of matched sequences in the genome)
7.	Palindrome_Seq: sequence of palindrome section
8.	Palindrome_Length: base pair length of palindrome section
9.	Reg_Expression: regular expression used (as denoted above)
10.	Genome_Length: length of genome
11.	Total_Hits: number of matched sequences
12.	Mean_Distance: average distance (in bp) between each matched sequence in genome
(If there are no matched sequences found, the table will display “No REP-like matches found”.)

The Excel file will be saved as: ‘[Accession Number]REP_matches’ in the folder Genome Repeats Data/REP-like.
After obtaining the data for each genome on the list, a one-line table consisting of all data from all genomes is developed. The one-line table will have the following information (with each line corresponding to one genome):
1.	Accession Number
2.	Genome Name
3.	Genome_Length
4.	REP_Hits: number of matched sequences (The last two rows displays number of genomes with matched sequences (upper row) and the number of genomes without matched sequences (lower row))
5.	Expected_Random_REP_Hits: the number of matched sequences that would appear by random chance, determined by genome length / 4^6. 6 is the length of the palindrome section.
6.	REP_Mean_Length: average length of matched sequences in the genome (The last row displays the average of the average length of matched sequences for all genomes with matched sequences.)
7.	REP_General_Mean_Distance: average distance between each matched sequence in genome (The last row displays the average of the average distances between matched sequences for all genomes with matched sequences.)
8.	CodeFreq: collection of columns that display first the decimal code for the palindrome section followed by the number of times it has appeared in matched sequences of the genome. The decimal code is defined in the nucleotide2dec.m function. 
(This file will be saved as ‘REP_One_Line_All_Genomes’ under the same REP-like folder. A documentation file containing similar and specific information on certain parameters will be made, labeled as ‘REP_One_Line_All_Genomes_Doc’. Note that sometimes, information can be missing in the xlsx file, so a csv file containing the same format is also constructed to contain all information. The documentation for xlsx is found under the Documentation sheet, while the documentation for csv will be in a separate file.)

The general repeat program follows a similar workflow except that the regular expression utilized is '(.{' num2str(Palindrome_Number) '}).{' num2str(Min_Middle_Section) ',' num2str(Max_Middle_Section) '}(??@palindrome_flip($1))'. The palindrome sequence is no longer a pre-specified sequence, and the length of the sequence is dictated by the parameter 'Palindrome number'. 

The nestedSortStruct.m and sortStruct.m were obtained from https://www.mathworks.com/matlabcentral/fileexchange/28573-nested-sort-of-structure-arrays?focused=5166120&tab=function.
