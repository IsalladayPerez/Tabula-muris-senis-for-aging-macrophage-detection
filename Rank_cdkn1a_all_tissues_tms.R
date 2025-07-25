#### Project MSen enrichment across all tissues and age in tabula muris senis data set 09062023 ###
#### A new bioconductor package dropped September 2023 Authored by Soneson C. et al 2023 on bioconductor ###
### This new package provides acccess to RNA-seq and scRNAseq from the Tabula muris senis data set so it bypassess doanloading all data loacally to computer###

### Install packages needed for the analysis 
library(ggplot2);
library(TabulaMurisSenisData);
library(SingleCellExperiment);
library(pheatmap);
library(RColorBrewer);
library(patchwork);
library(ggplotify);

### import data with Tabulamurissenis package from bioconductor####
tms_droplet <- TabulaMurisSenisDroplet(tissues = "Liver", processedCounts = TRUE);
tms_droplet <- TabulaMurisSenisDroplet(tissues = "All", processedCounts = TRUE)$All;

# tissue colors
tissue_cols <- c(Skin                = "#e6550d", 
                 Pancreas            = "#3182bd", 
                 Limb_Muscle         = "#d6616b",
                 Heart               = "#e7ba52", 
                 Spleen              = "#fd8d3c", 
                 Diaphragm           = "#8c6d31", 
                 Trachea             = "#636363", 
                 Tongue              = "#a1d99b", 
                 Thymus              = "#31a354", 
                 `Brain_Non-Myeloid` = "#cedb9c", 
                 Brain_Myeloid       = "#b5cf6b", 
                 Bladder             = "#637939", 
                 Large_Intestine     = "#843c39", 
                 BAT                 = "#9c9ede", 
                 GAT                 = "#bd9e39", 
                 MAT                 = "#a55194", 
                 SCAT                = "#6baed6", 
                 Lung                = "#7b4173", 
                 Liver               = "#e7969c", 
                 Marrow              = "#de9ed6", 
                 Kidney              = "#e7cb94", 
                 Aorta               = "#393b79", 
                 Mammary_Gland       = "#ce6dbd"
                 );

# get dataset with all tissues
se <- tms_droplet$Liver
# prepare data set for ggplot
ds <- as.data.frame(reducedDim(se, "UMAP"));
ds <- cbind(ds, tissue = colData(se)$tissue);
ds <- cbind(ds, Cell_type = colData(se)$cell);
ds <- cbind(ds, Age = colData(se)$age);
ds <- cbind(ds, Sex = colData(se)$sex);

# plot umap for all cells in all tissues

ggplot(ds, 
      aes(x = UMAP1, y = UMAP2, color = Cell_type)) + 
      geom_point(size = 0.05) + 
      theme_test() + 
      labs(title = "Tabula Muris Senis",
           color = "Cell Type")

ggplot



## isolate all genes that match list with in genes_of_interest
genes_of_interest <- c("Cdkn1a",
                       "Cdkn2a",
                       "Trp53",
                       "Cdkn2b")

## query all counts data 
log_genes_all_cells <- as.data.frame(assays(se)$counts);
