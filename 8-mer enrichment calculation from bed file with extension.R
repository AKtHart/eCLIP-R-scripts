library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19) # Load human genome hg19
library(universalmotif)

# ==========================================
# 1. Setup Interactive File Selection
# ==========================================
message("Please select your .bed file in the pop-up window...")
bed_file_path <- file.choose()

working_dir <- dirname(bed_file_path)
input_filename <- basename(bed_file_path) 
output_filename <- "eCLIP_8mers_ranked_ext20_10_MN3.csv"
output_file_path <- file.path(working_dir, output_filename)

# Define the number of background shuffles to average over
num_iterations <- 3

message("Using file: ", bed_file_path)
message(sprintf("Running analysis using %d averaged background shuffles...", num_iterations))

# ==========================================
# 2. Load and Process eCLIP Sequences with Crosslink Extension
# ==========================================
message("Loading and processing BED file peaks...")
bed_df <- read.table(bed_file_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)

# Filter peaks based on the p-value column (Column 5 <= 0.05)
bed_df <- bed_df[bed_df$V5 <= 0.05, ]

peaks_gr <- GRanges(
  seqnames = bed_df$V1,
  ranges = IRanges(start = bed_df$V2 + 1, end = bed_df$V3),
  strand = bed_df$V6
)

genome <- BSgenome.Hsapiens.UCSC.hg19
peaks_gr <- peaks_gr[seqnames(peaks_gr) %in% seqnames(genome)]

# Strand-aware extension focusing on the 5' crosslink site (RT stop):
# 'upstream = 20' extends 20 bp upstream of the 5' crosslink site on RNA.
# 'downstream = 10' retains 10 bp downstream into the peak body.
peaks_gr <- promoters(peaks_gr, upstream = 20, downstream = 10)

# Clamp ranges so they don't exceed chromosome boundaries
peaks_gr <- trim(peaks_gr)

peak_seqs <- getSeq(genome, peaks_gr)

# ==========================================
# 3. Calculate Real Pool Frequencies
# ==========================================
message("Counting all 8-mers in real peak pool...")
real_counts <- colSums(oligonucleotideFrequency(peak_seqs, width = 8))
total_real_positions <- sum(width(peak_seqs) - 8 + 1)

# ==========================================
# 4. Multi-Iteration Shuffled Control Pool
# ==========================================
# Initialize a vector of zeros to accumulate control counts across iterations
accumulated_shuffled_counts <- numeric(length(real_counts))
names(accumulated_shuffled_counts) <- names(real_counts)

for (j in 1:num_iterations) {
  message(sprintf(" -> Generating and counting background shuffle iteration %d/%d...", j, num_iterations))
  
  # Shuffle the sequence pool using the Eulerian path dinucleotide method
  shuffled_seqs <- shuffle_sequences(peak_seqs, k = 2, method = "euler")
  
  # Accumulate the counts
  iter_counts <- colSums(oligonucleotideFrequency(shuffled_seqs, width = 8))
  accumulated_shuffled_counts <- accumulated_shuffled_counts + iter_counts
}

# Calculate the mean (average) expected count for each 8-mer
mean_shuffled_counts <- accumulated_shuffled_counts / num_iterations
total_shuffled_positions <- sum(width(peak_seqs) - 8 + 1) # Opportunities remain identical per run

# Build an integrated analysis dataframe
results_df <- data.frame(
  kmer = names(real_counts),
  Real_Count = as.numeric(real_counts),
  Avg_Shuffled_Count = as.numeric(mean_shuffled_counts),
  stringsAsFactors = FALSE
)

# ==========================================
# 5. Mass Statistical Testing (Fisher's Exact Test)
# ==========================================
message("Performing Fisher's Exact Tests on all 8-mers using averaged backgrounds...")

p_values <- numeric(nrow(results_df))
odds_ratios <- numeric(nrow(results_df))

# Loop row-by-row to evaluate every 8-mer
for (i in 1:nrow(results_df)) {
  a <- results_df$Real_Count[i]
  
  # Round the averaged background count to the nearest whole number 
  # because Fisher's Exact Test strictly requires integers (whole counts)
  b <- round(results_df$Avg_Shuffled_Count[i])
  
  c <- total_real_positions - a
  d <- total_shuffled_positions - b
  
  # Construct the 2x2 table
  contingency_matrix <- matrix(c(a, b, c, d), nrow = 2)
  
  # Compute one-tailed Fisher's exact test for Enrichment
  ft <- fisher.test(contingency_matrix, alternative = "greater")
  
  p_values[i] <- ft$p.value
  odds_ratios[i] <- as.numeric(ft$estimate)
}

# Add stats to data frame
results_df$p_value <- p_values
results_df$Odds_Ratio <- odds_ratios

# Adjust for multiple testing errors (Benjamini-Hochberg False Discovery Rate)
results_df$FDR_adjusted_p <- p.adjust(results_df$p_value, method = "BH")

# ==========================================
# 6. Sorting and Saving Output
# ==========================================
message("Sorting and saving the data table...")

# Sort exclusively by p-value (most strongly enriched 8-mers rise to the top)
results_df <- results_df[order(results_df$p_value), ]

# Save as a standard CSV format
write.csv(results_df, file = output_file_path, row.names = FALSE)

message("Analysis complete! The averaged, stable results table has been saved.")