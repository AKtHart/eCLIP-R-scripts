library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(ggplot2)

# ==========================================
# 1. Load Files & Set Up Environment
# ==========================================
message("Select your peak .bed file...")
bed_path <- file.choose()

message("Select your ranked k-mer analysis .csv output...")
kmer_path <- file.choose()

genome <- BSgenome.Hsapiens.UCSC.hg19

# Parse eCLIP peaks (apply the -log10(p) >= 1.301 filter and 20 nt 5' extension)
bed_df <- read.table(bed_path, sep = "", header = FALSE, stringsAsFactors = FALSE, comment.char = "#", fill = TRUE)
bed_df$V5 <- as.numeric(bed_df$V5)
bed_df <- bed_df[!is.na(bed_df$V5) & bed_df$V5 <= 0.05, ]

peaks_gr <- GRanges(
  seqnames = bed_df$V1,
  ranges = IRanges(start = bed_df$V2 + 1, end = bed_df$V3),
  strand = bed_df$V6
)
peaks_gr <- peaks_gr[seqnames(peaks_gr) %in% seqnames(genome)]
peaks_gr <- promoters(peaks_gr, upstream = 20, downstream = 10)
peaks_gr <- trim(peaks_gr)

# Extract peak sequences
peak_seqs <- getSeq(genome, peaks_gr)

# Load ranked k-mers and select top hits
kmer_df <- read.csv(kmer_path)
top_kmers <- kmer_df %>%
  filter(FDR_adjusted_p < 0.05 & Odds_Ratio > 1.5) %>%
  pull(kmer)

message(sprintf("Loaded %d peaks and %d significantly enriched k-mers.", length(peaks_gr), length(top_kmers)))

# ==========================================
# 2. Map Top k-mers to Peak Coordinates
# ==========================================
message("Locating k-mer positions within peak sequences...")

hit_ranges_list <- list()

for (kmer in top_kmers) {
  # Ensure k-mer uses DNA alphabet (T instead of U) to match BSgenome
  kmer_dna <- gsub("U", "T", kmer)
  
  # vmatchPattern searches across all sequences in peak_seqs (XStringSet) at once
  matches_list <- vmatchPattern(kmer_dna, peak_seqs)
  
  # Find indices of peaks that contain at least one match
  matching_peaks <- which(elementNROWS(matches_list) > 0)
  
  for (i in matching_peaks) {
    peak_i <- peaks_gr[i]
    local_ir <- matches_list[[i]] # IRanges of matches within peak sequence i
    
    # Map local relative offset to absolute genomic coordinates
    if (as.character(strand(peak_i)) == "+") {
      g_start <- start(peak_i) + start(local_ir) - 1
      g_end <- start(peak_i) + end(local_ir) - 1
    } else {
      # Handle minus strand coordinate orientation
      g_start <- end(peak_i) - end(local_ir) + 1
      g_end <- end(peak_i) - start(local_ir) + 1
    }
    
    gr_hit <- GRanges(
      seqnames = seqnames(peak_i),
      ranges = IRanges(start = g_start, end = g_end),
      strand = strand(peak_i),
      kmer = kmer,
      peak_id = i
    )
    hit_ranges_list[[length(hit_ranges_list) + 1]] <- gr_hit
  }
}

all_hits_gr <- unlist(GRangesList(hit_ranges_list))

# ==========================================
# 3. Collapse Overlapping Hits into Contigs
# ==========================================
message("Merging overlapping k-mer matches...")

# reduce() merges overlapping or directly adjacent ranges on the same strand
assembled_contigs_gr <- reduce(all_hits_gr, with.revmap = TRUE)

# Extract the reconstructed longer genomic sequence for each assembled region
assembled_seqs <- getSeq(genome, assembled_contigs_gr)

# Build a summary data frame of assembled binding regions
contigs_df <- data.frame(
  chrom = seqnames(assembled_contigs_gr),
  start = start(assembled_contigs_gr),
  end = end(assembled_contigs_gr),
  strand = strand(assembled_contigs_gr),
  width = width(assembled_contigs_gr),
  sequence = as.character(assembled_seqs),
  num_overlapping_kmers = sapply(mcols(assembled_contigs_gr)$revmap, length),
  stringsAsFactors = FALSE
)

# ==========================================
# 4. Summary & Diagnostic Plotting
# ==========================================
cat("\n=== Assembly Summary ===\n")
cat("Total individual k-mer instances mapped:", length(all_hits_gr), "\n")
cat("Total unique collapsed contigs formed:", nrow(contigs_df), "\n")
cat("Average length of assembled sites:", round(mean(contigs_df$width), 1), "bp\n")
cat("Max length of assembled site:", max(contigs_df$width), "bp\n\n")

# Distribution of assembled sequence lengths
p_len <- ggplot(contigs_df, aes(x = width)) +
  geom_histogram(binwidth = 1, fill = "#2b5c8f", color = "black", alpha = 0.8) +
  theme_classic(base_size = 12) +
  labs(
    title = "Length Distribution of Assembled K-mer Clusters",
    x = "Assembled Sequence Length (bp)",
    y = "Frequency"
  )

print(p_len)

# Save assembled sequence regions to CSV
output_csv <- file.path(dirname(kmer_path), "assembled_kmer_contigs.csv")
write.csv(contigs_df, output_csv, row.names = FALSE)
message("Saved assembled sequence regions to: ", output_csv)