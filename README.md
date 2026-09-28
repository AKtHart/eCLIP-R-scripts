**R scripts to analyze eCLIP datasets**

_Scripts used for datasets originating from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86040_

"8-mer enrichment calculation from bed file with extension.R"
Analyzes the enrichment of each possible 8-mer in the significantly enriched peaks in the bed-file. Peaks are extended 20 nts upstream of the crosslink site, while 10 nts of the downstream peak body are retained. 8-mer enrichment is tested against three variations of dinucleotide shuffled sequence using Fisher's Exact Test.

"Search overlapping 8mers"
Analyzes whether a specific 8mer overlaps with the next according to its genomic location to assemble longer sequences bound by the target protein. Uses the .csv file produced by the previous script and the bed file of the eCLIP analysis as input. Only enriched sequences (FDR < 0.05 & Odds_ratio > 1.5) are used. 

_Scripts used for datasets originating from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1031_65

Same as above scripts but with _MCF7 extension. The only difference is the different handling of statistics as this dataset reports its peak significant in -log10. 

_Scripts used to compare analyzed datasets_

"TRRR bar graph comparison between cell types.R"
Uses the bed files from different experiments and calculates the enrichment of all possible TRRR sequences and plots them in a bar graph. 
