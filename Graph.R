#!/usr/bin/env

suppressPackageStartupMessages({
  library(ComplexUpset)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggvenn)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Graph.R <pangene_matrix.tsv>")
}
input_file <- args[1]
output_prefix <- sub("\\.tsv$", "", input_file)
output_png <- paste0(output_prefix, "_upset_plot.png")
output_venn <- paste0(output_prefix, "_venn_plot.png")
output_barplot <- paste0(output_prefix, "_validated_gene_counts_barplot.png")
summary_file <- sub("pangene_matrix.tsv$", "pangene_summary.tsv", input_file)

# ----------------------
# Load pangene matrix
# ----------------------
pangene_df <- read_tsv(input_file)
binary_cols <- names(pangene_df)[-1]
pangene_df <- pangene_df %>%
  mutate(across(all_of(binary_cols), ~ .x == 1))

# ----------------------
# UpSet Plot
# ----------------------
set_sizes <- pangene_df %>%
  summarise(across(all_of(binary_cols), sum)) %>%
  pivot_longer(everything(), names_to = "Genome", values_to = "Count")

p <- upset(
  pangene_df,
  intersect = binary_cols,
  name = "Pangene Sets",
  width_ratio = 0.2,
  min_size = 1,
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(size = 3),
      bar_number_threshold = 1
    )
  ),
  set_sizes = (
    upset_set_size(
      aes(x=Count),
      text = list(size = 3)
    ) +
      scale_x_continuous(expand = c(0, 0)) +
      labs(x = "Gene Families per Genome")
  ),
  themes = upset_default_themes(text_scale = 1.2)
)

ggsave(output_png, plot = p, width = 10, height = 6, dpi = 300)
cat(paste0("UpSet plot saved to: ", output_png, "\n"))

# ----------------------
# Venn Diagram (if ≥2 genomes)
# ----------------------
if (length(binary_cols) >= 2) {
  venn_list <- lapply(binary_cols, function(col) {
    pangene_df %>% filter(!!sym(col)) %>% pull(GeneFamilyID)
  })
  names(venn_list) <- binary_cols
  
  venn_plot <- ggvenn(venn_list, show_percentage = FALSE, stroke_size = 0.3, set_name_size = 4)
  ggsave(output_venn, plot = venn_plot, width = 6, height = 5, dpi = 300)
  cat(paste0("Venn diagram saved to: ", output_venn, "\n"))
}

# ----------------------
# Barplot: Final Validated Genes Per Genome
# ----------------------
if (file.exists(summary_file)) {
  summary_lines <- read_lines(summary_file)
  idx <- which(summary_lines == "Final Validated Gene Counts per Genome")
  if (length(idx) > 0 && length(summary_lines) > idx + 1) {
    validated_df <- read_tsv(
      paste(summary_lines[(idx+2):length(summary_lines)], collapse="\n"),
      col_names = TRUE
    )
    barplot <- ggplot(validated_df, aes(x=Genome, y=`Validated Genes`, fill=Genome)) +
      geom_bar(stat="identity", width=0.7) +
      geom_text(aes(label=`Validated Genes`), vjust=-0.3) +
      theme_minimal(base_size = 13) +
      labs(title = "Final Validated Genes per Genome", y = "Gene Count", x = "Genome") +
      theme(legend.position="none")
    
    ggsave(output_barplot, plot=barplot, width=7, height=5, dpi=300)
    cat(paste0("Validated gene barplot saved to: ", output_barplot, "\n"))
  }
}
