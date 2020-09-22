# Gene Prediction by Similarity Approach (GPSA)
 This project, Gene Prediction by Similarity Approach (GPSA), aims to provide a graphical tool to identify the genomic region of an unknown DNA sequence based on well-known genes by similarity approach. The program allows users to select an unknown DNA sequence, a template or mRNA (well-known gene), and specify scoring scheme (for a match, mismatch or gap) and the threshold of the alignment score. Local Alignment algorithm will be used to align the two sequences. A putative exon is an alignment according to the given threshold found in the Local Alignment process. Exon Chaining algorithm will be used to find the maximum chain of non-overlapping putative genomic regions.

![Image of GPSA](https://i.imgur.com/4AFztxY.jpg)

1. Files paths: browse to select Unknown DNA Sequence and DNA template accordingly. It filters to select *.txt or *.fasta (ignore line that starts with character '>') files.
1. Scoring scheme: 
   1. Match: spin box of integer ranges from 1 to 99
   1. Mismatch and Gap: spin box of integer range from -99 to -1
   1. Threshold: textfield that can only accept integer
1. View Alignment Scores: button to run the Local Alignment algorithm on given sequences
1. Find Exons: button to run the Exon Chaining algorithm on the result of Local Alignment
1. Table: to display the result of Local Alignment and Exon Chaining
   1. Start and End: format as index_of_template-index_of_unknown_sequence
   1. Score: score of the alignment from the scoring scheme
   1. View Alignment: double click to display the alignment in the plain textfield below
1. Alignment: to display the selected alignment
