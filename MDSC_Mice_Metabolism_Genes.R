# Load necessary library
library(readxl)

# Get a list of all .xlsx files in the current working directory
file_list <- list.files(pattern = "*.xlsx")

# Loop over the file list, read each file into a data frame and name it based on the file name
for (file in file_list) {
  # Create a valid variable name from the file name (removing the .xlsx extension and any invalid characters)
  df_name <- gsub("[^[:alnum:]_]", "", tools::file_path_sans_ext(file))
  
  # Read the file into a data frame
  assign(df_name, read_excel(file))
}

# Load required libraries
library(biomaRt)
library(dplyr)
library(openxlsx)

# 1. Convert "Gene symbol" column in RM1MetabolismGeneList from human to mouse annotation
human_to_mouse <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

# Get mouse homologs for the human gene symbols
human_genes <- RM1MetabolismGeneList$GeneSymbol
mouse_genes <- getLDS(attributes = c("hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values = human_genes,
                      mart = human_to_mouse,
                      attributesL = c("mgi_symbol"),
                      martL = mouse_mart)

# Add Mouse_gene_symbol to RM1MetabolismGeneList
RM1MetabolismGeneList <- RM1MetabolismGeneList %>%
  left_join(mouse_genes, by = c("Gene Symbol" = "HGNC.symbol"))

# 2. Rename columns named "...1" as "Gene symbol" in the three DE data frames
dfs <- list(
  DE_gene_per_fine_celltype_agegroup_young_vs_old,
  DE_gene_per_fine_celltype_sepsisgroup_sepsis_vs_naive,
  DE_gene_per_fine_celltype_sexgroup_male_vs_female
)

dfs_renamed <- lapply(dfs, function(df) {
  colnames(df)[colnames(df) == "...1"] <- "Gene symbol"
  df
})

# Reassign the renamed dataframes
DE_gene_per_fine_celltype_agegroup_young_vs_old <- dfs_renamed[[1]]
DE_gene_per_fine_celltype_sepsisgroup_sepsis_vs_naive <- dfs_renamed[[2]]
DE_gene_per_fine_celltype_sexgroup_male_vs_female <- dfs_renamed[[3]]

###################
#GSEA##############
###################

# Load the necessary libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation
library(msigdbr)       # Gene sets
library(dplyr)

# Prepare ranked gene list for GSEA (without filtering)
prepare_gsea_input <- function(df) {
  # Arrange the genes by log2FoldChange (descending order)
  ranked_genes <- df %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::select(GeneSymbol, log2FoldChange)
  
  # Create a named vector for GSEA (gene names as names, log2FoldChange as values)
  ranked_vector <- setNames(ranked_genes$log2FoldChange, ranked_genes$GeneSymbol)
  return(ranked_vector)
}

# Prepare GSEA inputs for each DE dataframe
gsea_inputs <- lapply(dfs_renamed, prepare_gsea_input)

# Perform GSEA on each DE dataframe
gsea_results <- lapply(gsea_inputs, function(ranked_vector) {
  GSEA(geneList = ranked_vector,
       TERM2GENE = cbind("Metabolism_Genes",na.omit(RM1MetabolismGeneList$MGI.symbol)), # Define custom gene set
       pvalueCutoff = 0.05,
       minGSSize = 10,  # Minimum size of gene sets
       maxGSSize = 500)
})

# Create a README dataframe with descriptions for each column
readme_df <- data.frame(
  "Column Name" = c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment"),
  "Description" = c(
    "ID of the gene set",
    "Description of the gene set or pathway",
    "Number of genes in the gene set",
    "Enrichment score, a running sum statistic calculated by GSEA",
    "Normalized Enrichment Score (adjusted for differences in gene set size)",
    "P-value for the enrichment of the gene set",
    "Adjusted p-value (corrected for multiple testing)",
    "False Discovery Rate (FDR) q-value",
    "Rank at which the maximum enrichment score is obtained",
    "Summary of core (leading edge) genes driving the enrichment",
    "List of core (leading edge) genes that drive the enrichment"
  )
)


# Create a workbook for the results
wb <- createWorkbook()
# Add the README sheet as the first sheet in the workbook
addWorksheet(wb, sheetName = "README")

# Write the README dataframe to the first sheet
writeData(wb, sheet = "README", readme_df)

# Loop through GSEA results and save each to a separate Excel sheet
for (i in seq_along(gsea_results)) {
  # Extract the GSEA result as a dataframe
  gsea_df <- as.data.frame(gsea_results[[i]]@result)
  
  # Extract the core enrichment (leading-edge) genes for each term
  core_enrichment_genes <- lapply(gsea_df$core_enrichment, function(x) {
    strsplit(x, "/")[[1]]  # Split the core.enrichment string by '/'
  })
  
  # Get the corresponding DE dataframe to match leading edge genes with their log2FoldChange and p-values
  de_df <- dfs_renamed[[i]]
  
  # Prepare the worksheet name, removing the first part before the underscore
  sheet_name <- paste0("GSEA_", sub("^[^_]*_", "", df_names[i]))
  addWorksheet(wb, sheetName = sheet_name)
  
  # Write the GSEA result to the worksheet
  writeData(wb, sheet = sheet_name, gsea_df)
  
  # Write the leading-edge genes, log2FoldChange, and p-values under the GSEA result for each term
  for (j in seq_along(core_enrichment_genes)) {
    # Create a data frame for leading-edge genes and their log2FoldChange, p-value from DE dataframe
    leading_edge_df <- data.frame(
      "Gene" = core_enrichment_genes[[j]],
      "log2FoldChange" = de_df$log2FoldChange[match(core_enrichment_genes[[j]], de_df$GeneSymbol)],
      "padj-value" = de_df$padj[match(core_enrichment_genes[[j]], de_df$GeneSymbol)]
    )
    
    # Write a label for the leading edge genes section
    writeData(wb, sheet = sheet_name, startRow = nrow(gsea_df) + 5 * j, 
              x = data.frame("Leading edge genes for GSEA term"))
    
    # Write the leading edge genes along with log2FoldChange and p-value
    writeData(wb, sheet = sheet_name, startRow = nrow(gsea_df) + 5 * j + 1, 
              x = leading_edge_df)
  }
}

# Save the workbook
saveWorkbook(wb, "GSEA_results_with_leading_edge_genes.xlsx", overwrite = TRUE)

#############################################
#############################################
# Perform GO enrichment for each DE dataframe
#############################################
#############################################
# Calculate the overlap between ranked gene list and RM1MetabolismGeneList$MGI.symbol
calculate_overlap <- function(ranked_vector, gene_set) {
  # Get the gene names from the ranked vector
  ranked_genes <- names(ranked_vector)
  
  # Identify the overlap
  overlap <- intersect(ranked_genes, gene_set)
  
  return(overlap)
}

# Apply the overlap calculation to each GSEA input
overlaps <- lapply(gsea_inputs, calculate_overlap, na.omit(RM1MetabolismGeneList$MGI.symbol))

# Perform the hypergeometric test
perform_hypergeometric_test <- function(overlap, ranked_vector, gene_set, total_genes) {
  # Number of genes in the ranked list
  M <- length(ranked_vector)  # Total genes in the DE list
  # Number of genes in the gene set
  N <- length(gene_set)  # Genes in RM1MetabolismGeneList
  # Number of genes in both (overlap)
  k <- length(overlap)
  # Total number of genes in the genome/universe
  T <- total_genes  # e.g., all genes in the genome (use appropriate value)
  
  # Perform the hypergeometric test
  p_value <- phyper(k - 1, N, T - N, M, lower.tail = FALSE)
  return(p_value)
}

# Total number of genes (you can adjust this based on the known size of your background)
total_genes_in_genome <- 20000

p_values <- lapply(overlaps, function(overlap, ranked_vector) {
  p_value <- perform_hypergeometric_test(overlap, ranked_vector, na.omit(RM1MetabolismGeneList$MGI.symbol), total_genes_in_genome)
  format(p_value, scientific = TRUE)  # Format the p-value here
}, ranked_vector = gsea_inputs)

# Create a README dataframe with descriptions for each column
readme_hyper_df <- data.frame(
  "Column Name" = c("DE_Dataframe", "Overlap_with_RM1", "Total_DE_Genes", "Total_RM1_Genes", "P_Value"),
  "Description" = c(
    "The name of the DE dataframe used in the overlap analysis",
    "Number of overlapping genes between the DE dataframe and the RM1MetabolismGeneList",
    "Total number of genes in the DE dataframe",
    "Total number of genes in the RM1MetabolismGeneList",
    "The p-value from the hypergeometric test for the significance of the overlap"
  )
)

# Create a workbook for the results
wb <- createWorkbook()
# Add the README sheet for hypergeometric overlap results
addWorksheet(wb, sheetName = "README_Overlap")

# Write the README dataframe to the first sheet for overlap results
writeData(wb, sheet = "README_Overlap", readme_hyper_df)

# Loop through each overlap and save results
for (i in seq_along(overlaps)) {
  # Retrieve the overlap and corresponding p-value
  overlap <- overlaps[[i]]
  p_value <- perform_hypergeometric_test(overlap, gsea_inputs[[i]], na.omit(RM1MetabolismGeneList$MGI.symbol), total_genes_in_genome)
  formatted_p_value <- format(p_value, scientific = TRUE)  # Format the p-value
  
  # Create a summary dataframe to hold overlap and p-value information
  summary_df <- data.frame(
    "DE_Dataframe" = df_names[i],
    "Overlap_with_RM1" = length(overlap),
    "Total_DE_Genes" = length(gsea_inputs[[i]]),
    "Total_RM1_Genes" = length(na.omit(RM1MetabolismGeneList$MGI.symbol)),
    "P_Value" = formatted_p_value
  )
  
  # Add a new worksheet with an informative name
  sheet_name <- df_names[i]
  addWorksheet(wb, sheetName = sheet_name)
  
  # Write the summary information about the overlap and p-value
  writeData(wb, sheet = sheet_name, summary_df)
  
  # Get the corresponding DE dataframe to match overlapping genes with their log2FoldChange and p-values
  de_df <- dfs_renamed[[i]]
  
  # Create a data frame for overlapping genes and their log2FoldChange, p-value from DE dataframe
  overlap_df <- data.frame(
    "Overlapping_Genes" = overlap,
    "log2FoldChange" = de_df$log2FoldChange[match(overlap, de_df$GeneSymbol)],
    "padj-value" = de_df$padj[match(overlap, de_df$GeneSymbol)]
  )
  
  # Write the list of overlapping genes with log2FoldChange and p-value below the summary
  writeData(wb, sheet = sheet_name, startRow = 6, x = overlap_df)
}

# Save the workbook
saveWorkbook(wb, "Overlap_P_Values_with_log2FC_and_p.xlsx", overwrite = TRUE)

###########
#GSEA plots
library(enrichplot)
library(ggplot2)

# Define a function to generate GSEA enrichment plots
create_gsea_plot <- function(gsea_result, gene_set_id, output_filename) {
  # Check if the gene set is significant
  if (gsea_result@result$p.adjust[gene_set_id] < 0.05) {
    # Create the enrichment plot
    p <- gseaplot(gsea_result, geneSetID = gene_set_id)
    
    # Save the plot to a file
    ggsave(output_filename, plot = p, width = 8, height = 6)
    
    # Optionally, display the plot in R
    print(p)
  }
}

# Loop through each GSEA result and generate plots for significant gene sets
for (i in seq_along(gsea_results)) {
  gsea_result <- gsea_results[[i]]
  significant_gene_sets <- which(gsea_result@result$p.adjust < 0.05)
  
  # Loop through significant gene sets and create plots
  for (gene_set_id in significant_gene_sets) {
    # Prepare the output filename (use the dataframe name and gene set ID)
    output_filename <- paste0("GSEA_", df_names[i], "_GeneSet_", gene_set_id, ".png")
    
    # Create the GSEA plot for the significant gene set
    create_gsea_plot(gsea_result, gene_set_id, output_filename)
  }
}
# GSEA plot for significant gene set
p1 <- gseaplot(gsea_results[[1]], geneSetID = 1,title = "Enrichment of Metabolism Gene Set (933 genes) in Young vs. Old log2FC Ranked DE Genes")
p2 <- gseaplot(gsea_results[[2]], geneSetID = 1,title = "Enrichment of Metabolism Gene Set (933 genes) in Sepsis vs. Naive log2FC Ranked DE Genes")  # Replace 1 with the actual significant gene set ID

print(p)  # Display the plot


