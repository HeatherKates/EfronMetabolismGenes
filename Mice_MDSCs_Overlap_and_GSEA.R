# Load necessary libraries
library(readxl)
library(dplyr)
library(openxlsx)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse gene annotation
library(msigdbr)       # Gene sets
library(enrichplot)
library(ggplot2)

# Load metabolism genes and convert "Gene symbol" column from human to mouse annotation
RM1MetabolismGeneList <- read_excel("Data/RM1-MetabolismGeneList.xlsx")
RM1MetabolismGeneList <- RM1MetabolismGeneList  %>% dplyr::rename("GeneSymbol"="Gene Symbol")
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
  left_join(mouse_genes, by = c("GeneSymbol" = "HGNC.symbol"))

# Get a list of all .xlsx files in the current working directory
file_list <- list.files(pattern = "Data/*.xlsx")
file_IDs <- c("Age", "Sepsis", "Sex", "None")  # Define file IDs

############
# Functions#
############

# Calculate the overlap between ranked gene list and RM1MetabolismGeneList$MGI.symbol
calculate_overlap <- function(ranked_vector, gene_set) {
  overlap <- intersect(names(ranked_vector), gene_set)
  return(overlap)
}

# Perform the hypergeometric test
perform_hypergeometric_test <- function(overlap, ranked_vector, gene_set, total_genes) {
  M <- length(ranked_vector)  # Total genes in the DE list
  N <- length(gene_set)       # Genes in RM1MetabolismGeneList
  k <- length(overlap)        # Number of genes in both (overlap)
  T <- total_genes            # Total genes in the genome/universe
  
  p_value <- phyper(k - 1, N, T - N, M, lower.tail = FALSE)
  return(p_value)
}

# Prepare ranked gene list for GSEA
prepare_gsea_input <- function(df) {
  ranked_genes <- df %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::select(GeneSymbol, log2FoldChange)
  
  ranked_vector <- setNames(ranked_genes$log2FoldChange, ranked_genes$GeneSymbol)
  return(ranked_vector)
}

# Abbreviate the sheet names using file_IDs and input sheet name
abbreviate_sheet_name <- function(test_type, file_index, input_sheet_name, file_IDs) {
  file_short <- file_IDs[file_index]
  sheet_short <- paste0(test_type, "_", file_short, "_", input_sheet_name)
  sheet_short <- substr(sheet_short, 1, 31)  # Ensure within Excel's 31-character limit
  return(sheet_short)
}

###############
# README Setup #
###############

# Combine the GSEA and Hypergeometric README descriptions
readme_combined <- rbind(
  data.frame(
    "Result_Type" = "GSEA",
    "Column_Name" = c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment"),
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
  ),
  data.frame(
    "Result_Type" = "Hypergeo",
    "Column_Name" = c("DE_Dataframe", "Overlap_with_RM1", "Total_DE_Genes", "Total_RM1_Genes", "P_Value"),
    "Description" = c(
      "The name of the DE dataframe used in the overlap analysis",
      "Number of overlapping genes between the DE dataframe and the RM1MetabolismGeneList",
      "Total number of genes in the DE dataframe",
      "Total number of genes in the RM1MetabolismGeneList",
      "The p-value from the hypergeometric test for the significance of the overlap"
    )
  )
)

#############
# Main Loop #
#############

# Loop over each Excel file
for (i in seq_along(file_list[1:3])) {
  file <- file_list[i]
  sheet_names <- excel_sheets(file)  # Get sheet names
  
  # Create a workbook for the current file
  wb <- createWorkbook()
  
  # Add the combined README sheet
  addWorksheet(wb, sheetName = "README")
  writeData(wb, sheet = "README", readme_combined)
  
  for (sheet in sheet_names) {
    df <- read_excel(file, sheet = sheet)
    colnames(df)[colnames(df) == "...1"] <- "GeneSymbol"  # Rename columns if needed
    
    # Prepare GSEA input
    ranked_vector <- prepare_gsea_input(df)
    
    # Perform GSEA with error handling
    gsea_result <- tryCatch({
      GSEA(
        geneList = ranked_vector,
        TERM2GENE = cbind("Metabolism_Genes", na.omit(RM1MetabolismGeneList$MGI.symbol)),
        pvalueCutoff = 0.05,
        minGSSize = 10,
        maxGSSize = 500
      )
    }, error = function(e) {
      message(paste("GSEA failed for:", sheet, "in file:", file))
      return(NULL)
    })
    
    # If GSEA failed or no terms enriched, skip the sheet
    if (is.null(gsea_result) || nrow(gsea_result@result) == 0) {
      next
    }
    
    # Abbreviate sheet names
    gsea_sheet_name <- abbreviate_sheet_name("GSEA", i, sheet, file_IDs)
    hypergeo_sheet_name <- abbreviate_sheet_name("Hypergeo", i, sheet, file_IDs)
    
    # Write GSEA results to workbook
    gsea_df <- as.data.frame(gsea_result@result)
    addWorksheet(wb, sheetName = gsea_sheet_name)
    writeData(wb, sheet = gsea_sheet_name, gsea_df)
    
    # Extract and write leading edge genes
    core_enrichment_genes <- lapply(gsea_df$core_enrichment, function(x) strsplit(x, "/")[[1]])
    leading_edge_df <- data.frame(
      "Gene" = core_enrichment_genes[[1]],
      "log2FoldChange" = df$log2FoldChange[match(core_enrichment_genes[[1]], df$GeneSymbol)],
      "padj-value" = df$padj[match(core_enrichment_genes[[1]], df$GeneSymbol)]
    )
    writeData(wb, sheet = gsea_sheet_name, leading_edge_df, startRow = nrow(gsea_df) + 5)
    
    # Perform hypergeometric test and write results
    overlap <- calculate_overlap(ranked_vector, na.omit(RM1MetabolismGeneList$MGI.symbol))
    p_value <- perform_hypergeometric_test(overlap, ranked_vector, na.omit(RM1MetabolismGeneList$MGI.symbol), 20000)
    
    summary_df <- data.frame(
      "DE_Dataframe" = paste(tools::file_path_sans_ext(file), sheet, sep = "_"),
      "Overlap_with_RM1" = length(overlap),
      "Total_DE_Genes" = length(ranked_vector),
      "Total_RM1_Genes" = length(na.omit(RM1MetabolismGeneList$MGI.symbol)),
      "P_Value" = format(p_value, scientific = TRUE)
    )
    
    addWorksheet(wb, sheetName = hypergeo_sheet_name)
    writeData(wb, sheet = hypergeo_sheet_name, summary_df)
    
    # Write overlap genes
    overlap_df <- data.frame(
      "Overlapping_Genes" = overlap,
      "log2FoldChange" = df$log2FoldChange[match(overlap, df$GeneSymbol)],
      "padj-value" = df$padj[match(overlap, df$GeneSymbol)]
    )
    writeData(wb, sheet = hypergeo_sheet_name, overlap_df, startRow = 6)
  }
  
  # Save the workbook once all sheets are processed
  saveWorkbook(wb, paste0("Results/Results_", file_IDs[i], ".xlsx"), overwrite = TRUE)
}

# Load necessary libraries
library(readxl)

# Initialize an empty data frame to store the summary
summary_df <- data.frame(
  DE = character(),
  CellType = character(),
  NumberofOverlap = numeric(),
  Significance_of_overlap = numeric(),
  NES = numeric(),
  Significance_of_GSEA = numeric(),
  stringsAsFactors = FALSE
)

# Get the list of all result files
result_files <- list.files(pattern = "Results/Results_.*\\.xlsx")

# Loop through each workbook
for (file in result_files) {
  # Get all sheet names from the workbook
  sheet_names <- excel_sheets(file)
  
  # Filter out the README sheet
  analysis_sheets <- sheet_names[!sheet_names %in% "README"]
  
  # Loop through each sheet in the workbook
  for (sheet in analysis_sheets) {
    # Extract DE and CellType from the sheet name
    DE <- strsplit(sheet, "_")[[1]][2]  # Extract DE from the sheet name
    CellType <- strsplit(sheet, "_")[[1]][3]  # Extract Cell Type from the sheet name
    
    # Initialize placeholders
    NumberofOverlap <- NA
    Significance_of_overlap <- NA
    NES <- NA
    Significance_of_GSEA <- NA
    
    # Check if the sheet contains Hypergeo results
    if (grepl("Hypergeo", sheet)) {
      # Read the Hypergeometric results
      hypergeo_data <- read_excel(file, sheet = sheet)
      
      # Extract the number of overlaps and significance (p-value) for the overlap
      NumberofOverlap <- hypergeo_data$Overlap_with_RM1[1]
      Significance_of_overlap <- hypergeo_data$P_Value[1]
      
      # Check if this DE-CellType combination already exists in summary_df
      existing_row <- which(summary_df$DE == DE & summary_df$CellType == CellType)
      
      if (length(existing_row) > 0) {
        # If it exists, update the overlap information
        summary_df$NumberofOverlap[existing_row] <- NumberofOverlap
        summary_df$Significance_of_overlap[existing_row] <- Significance_of_overlap
      } else {
        # If it doesn't exist, create a new row with only the Hypergeo data
        summary_df <- rbind(
          summary_df,
          data.frame(
            DE = DE,
            CellType = CellType,
            NumberofOverlap = NumberofOverlap,
            Significance_of_overlap = Significance_of_overlap,
            NES = NA,
            Significance_of_GSEA = NA,
            stringsAsFactors = FALSE
          )
        )
      }
      
    } else if (grepl("GSEA", sheet)) {
      # Read the GSEA results
      gsea_data <- read_excel(file, sheet = sheet)
      
      # Extract the NES and significance (adjusted p-value) for the GSEA results
      if (!is.null(gsea_data$NES)) {
        NES <- gsea_data$NES[1]  # NES (first entry)
      }
      if (!is.null(gsea_data$p.adjust)) {
        Significance_of_GSEA <- gsea_data$p.adjust[1]  # Adjusted p-value for GSEA (first entry)
      }
      
      # Check if this DE-CellType combination already exists in summary_df
      existing_row <- which(summary_df$DE == DE & summary_df$CellType == CellType)
      
      if (length(existing_row) > 0) {
        # If it exists, update the GSEA information
        summary_df$NES[existing_row] <- NES
        summary_df$Significance_of_GSEA[existing_row] <- Significance_of_GSEA
      } else {
        # If it doesn't exist, create a new row with only the GSEA data
        summary_df <- rbind(
          summary_df,
          data.frame(
            DE = DE,
            CellType = CellType,
            NumberofOverlap = NA,
            Significance_of_overlap = NA,
            NES = NES,
            Significance_of_GSEA = Significance_of_GSEA,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}

# View the summary dataframe
print(summary_df)

# Optionally, write the summary dataframe to an Excel file
write.xlsx(summary_df, "Results/Summary_Results.xlsx", overwrite = TRUE)
