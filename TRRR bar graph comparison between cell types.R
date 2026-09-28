library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(ggplot2)
library(tools)

# 0. Check and load required packages
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)

genome <- BSgenome.Hsapiens.UCSC.hg19

# ==========================================
# 1. Load Target Sequences
# ==========================================
message("Select your sequence text file (one sequence per line)...")
seq_file_path <- file.choose()

target_seqs_raw <- readLines(seq_file_path)
target_seqs_raw <- trimws(target_seqs_raw)
target_seqs_raw <- target_seqs_raw[target_seqs_raw != ""] # Remove blank lines
target_seqs_dna <- gsub("U", "T", toupper(target_seqs_raw))

# ==========================================
# 2. Function to Process a Single BED File
# ==========================================
process_bed_file <- function(bed_path, target_seqs_dna, target_seqs_raw, genome) {
  bed_df <- read.table(
    bed_path, 
    sep = "", 
    header = FALSE, 
    stringsAsFactors = FALSE, 
    comment.char = "#", 
    fill = TRUE
  )
  
  bed_df$V5 <- as.numeric(bed_df$V5)
  
  # Detect score scaling and filter
  if (max(bed_df$V5, na.rm = TRUE) <= 1.0) {
    bed_df <- bed_df[!is.na(bed_df$V5) & bed_df$V5 <= 0.05, ]
  } else {
    bed_df <- bed_df[!is.na(bed_df$V5) & bed_df$V5 >= 1.301, ]
  }
  
  total_sig_peaks <- nrow(bed_df)
  if (total_sig_peaks == 0) stop(paste("No peaks survived in:", bed_path))
  
  # Build GRanges with 20 nt 5' extension
  peaks_gr <- GRanges(
    seqnames = bed_df$V1,
    ranges = IRanges(start = bed_df$V2 + 1, end = bed_df$V3),
    strand = bed_df$V6
  )
  peaks_gr <- peaks_gr[seqnames(peaks_gr) %in% seqnames(genome)]
  peaks_gr <- promoters(peaks_gr, upstream = 20, downstream = 10)
  peaks_gr <- trim(peaks_gr)
  
  peak_seqs <- getSeq(genome, peaks_gr)
  
  # Count occurrences per target sequence
  res_list <- list()
  for (i in seq_along(target_seqs_dna)) {
    seq_dna <- target_seqs_dna[i]
    seq_orig <- target_seqs_raw[i]
    
    matches_list <- vmatchPattern(seq_dna, peak_seqs)
    peaks_with_match <- sum(elementNROWS(matches_list) > 0)
    pct_peaks <- (peaks_with_match / total_sig_peaks) * 100
    
    res_list[[i]] <- data.frame(
      Sequence = seq_orig,
      Percent_Peaks = pct_peaks,
      stringsAsFactors = FALSE
    )
  }
  return(bind_rows(res_list))
}

# ==========================================
# 3. GUI File Selection & Cell Type Naming
# ==========================================
# Check for svDialogs package to ensure reliable GUI prompt
if (!requireNamespace("svDialogs", quietly = TRUE)) install.packages("svDialogs")
library(svDialogs)

# Cell Type A Selection
message("\n--- CELL TYPE A ---")
message("Select BED File for Cell Type A - Replicate 1...")
file_A1 <- file.choose()
message("Select BED File for Cell Type A - Replicate 2...")
file_A2 <- file.choose()

# Cell Type B Selection
message("\n--- CELL TYPE B ---")
message("Select BED File for Cell Type B - Replicate 1...")
file_B1 <- file.choose()
message("Select BED File for Cell Type B - Replicate 2...")
file_B2 <- file.choose()

# GUI Popups for Cell Type Names
cell_type_A_name <- dlgInput("Enter name for Cell Type A (e.g., K562):", default = "Cell Type A")$res
if (length(cell_type_A_name) == 0 || cell_type_A_name == "") cell_type_A_name <- "Cell Type A"

cell_type_B_name <- dlgInput("Enter name for Cell Type B (e.g., HepG2):", default = "Cell Type B")$res
if (length(cell_type_B_name) == 0 || cell_type_B_name == "") cell_type_B_name <- "Cell Type B"

# ==========================================
# 4. Process Files & Aggregate Replicates
# ==========================================
message("\nProcessing BED files...")

df_A1 <- process_bed_file(file_A1, target_seqs_dna, target_seqs_raw, genome) %>% mutate(Cell_Type = cell_type_A_name, Replicate = "Rep1")
df_A2 <- process_bed_file(file_A2, target_seqs_dna, target_seqs_raw, genome) %>% mutate(Cell_Type = cell_type_A_name, Replicate = "Rep2")
df_B1 <- process_bed_file(file_B1, target_seqs_dna, target_seqs_raw, genome) %>% mutate(Cell_Type = cell_type_B_name, Replicate = "Rep1")
df_B2 <- process_bed_file(file_B2, target_seqs_dna, target_seqs_raw, genome) %>% mutate(Cell_Type = cell_type_B_name, Replicate = "Rep2")

# Combine all counts
all_reps_df <- bind_rows(df_A1, df_A2, df_B1, df_B2)

# Calculate Mean and SD across replicates
summary_df <- all_reps_df %>%
  group_by(Cell_Type, Sequence) %>%
  summarise(
    Mean_Percent = mean(Percent_Peaks),
    SD_Percent = sd(Percent_Peaks),
    .groups = "drop"
  )

# Preserve exact input sequence order & set group order
summary_df <- summary_df %>%
  mutate(
    Sequence = factor(Sequence, levels = target_seqs_raw),
    Cell_Type = factor(Cell_Type, levels = c(cell_type_A_name, cell_type_B_name))
  )

# ==========================================
# 5. Save Summary CSV & Plot SVG
# ==========================================
output_dir <- dirname(seq_file_path)
write.csv(summary_df, file.path(output_dir, "cell_type_comparison_summary.csv"), row.names = FALSE)

output_svg <- file.path(output_dir, "cell_type_comparison_bargraph.svg")

# Generate Grouped Bar Graph with Error Bars
p_comp <- ggplot(summary_df, aes(x = Sequence, y = Mean_Percent, fill = Cell_Type)) +
  # Side-by-side grouped bars
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.3) +
  
  # Add Replicate SD Error Bars
  geom_errorbar(
    aes(ymin = pmax(0, Mean_Percent - SD_Percent), ymax = Mean_Percent + SD_Percent),
    position = position_dodge(width = 0.75),
    width = 0.2,
    linewidth = 0.3
  ) +
  
  # Custom Colors
  scale_fill_manual(values = setNames(c("#2b5c8f", "#d95f02"), c(cell_type_A_name, cell_type_B_name))) +
  
  scale_y_continuous(
    limits = c(0, max(summary_df$Mean_Percent + summary_df$SD_Percent, na.rm = TRUE) * 1.15),
    expand = c(0, 0)
  ) +
  
  # Set base font size to 8 pt with Arial family
  theme_classic(base_size = 8, base_family = "Arial") +
  labs(
    title = "Sequence Occurrence Across Cell Types",
    x = "Sequence Motif",
    y = "Peaks Containing Sequence (%)",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5, family = "Arial"),
    axis.title = element_text(size = 8, family = "Arial"),
    axis.text.x = element_text(size = 8, family = "Arial", angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 8, family = "Arial"),
    axis.line = element_line(linewidth = 0.4),
    legend.text = element_text(size = 7, family = "Arial"),
    legend.position = "top",
    legend.key.size = unit(0.3, "cm")
  )

# Save SVG Output (8 cm Total Width, editable text)
ggsave(
  filename = output_svg, 
  plot = p_comp, 
  width = 8, 
  height = 7.5, 
  units = "cm", 
  device = svglite::svglite
)

cat("\nSummary table and grouped bar graph saved to:\n", output_dir, "\n")