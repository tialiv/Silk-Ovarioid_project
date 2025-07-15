setwd("~/Desktop/PhD/Ovarioid/Silk-Ovarioid_project/codes/")

source("Deconvolution_functions.R")

#get gene description
mart <- useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl")
output <- getBM(attributes=c("ensembl_gene_id","external_gene_name","entrezgene_id","description"), mart = mart)

# Load bulk data
bulk <- counts(dataset[["dds_Cortex"]], normalized = F) %>% as.data.frame()
bulk <- counts(dataset[["dds_Medulla"]], normalized = F) %>% as.data.frame()
bulk <- edgeR::cpm(bulk, normalized.lib.sizes = TRUE) %>% as.data.frame()
bulk$ENSEMBL <- rownames(bulk)
temp1 <- output[which(output$ensembl_gene_id %in% bulk$ENSEMBL),]
temp1 <- temp1[,1:2] 
colnames(temp1)[1] <- "ENSEMBL"

# Keep the genes with the gene names
bulk <- bulk[which(bulk$ENSEMBL %in% temp1$ENSEMBL),]
bulk_df <- base::merge(bulk, temp1, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = T, all.y = F)
bulk_df <- bulk_df[!duplicated(bulk_df),]
bulk_df <- bulk_df[-which(duplicated(bulk_df$external_gene_name)),]
bulk_df <- bulk_df[,-1]

# Relocate gene column, gene name must be in the first column of the bulk data matrix
bulk_df <- bulk_df %>% relocate(!!"external_gene_name")
colnames(bulk_df)[1] <- "gene_id"

# Load scRNAseq data
#sc <- readRDS("~/Desktop/PhD/Ovarioid/data/Unsorted_Manuscript figures.rds")
sc_in <- readRDS("~/Desktop/PhD/Ovarioid/data/Integrated_Manuscript figures.rds")

# Check single cell labels
#Idents(sc) <- sc@active.ident
#sc <- AddMetaData(sc, sc@active.ident, col.name = "clusters")
#sc_v5 <- UpdateSeuratObject(sc)

#DefaultAssay(sc_v5) <- "RNA"
#sc_sub <- subset(x = sc_v5, idents = c("2_immune", "3_gran", "4_endo", "5_pv", "6_stroma"))

# For integrated data
Idents(sc_in)
sc_in_v5 <- UpdateSeuratObject(sc_in)
sc_in_v5 <- AddMetaData(sc_in_v5, Idents(sc_in), col.name = "clusters")
DimPlot(sc_in_v5, group.by = "clusters", repel = T, label = T)
sc_in_v5 <- RenameIdents(sc_in_v5, "0" = "Stroma_theca", "1" = "Stroma_theca", "2" = "Stroma_theca", "9" = "Stroma_theca",
                         "13" = "Perivascular", "6" = "Perivascular", "3" = "Endothelial",
                         "12" = "Endothelial", "7" = "T_NK", "11" = "Granulosa",
                         "10" = "Monocytes", "4" = "Granulosa", "5" = "Granulosa", "8" = "Granulosa")
sc_in_v5 <- AddMetaData(sc_in_v5, sc_in_v5@active.ident, col.name = "integrated_clusters")
DimPlot(sc_in_v5, group.by = "integrated_clusters", repel = T, label = T)

DefaultAssay(sc_in_v5) <- "RNA"
#sc_sub <- subset(x = sc_in, idents = c("2_immune", "3_gran", "4_endo", "5_pv", "6_stroma"))
sc_sub <- sc_in_v5

# Get scRNAseq matrix
sc.mtx <- sc_sub[["RNA"]]$counts 
sc.meta <- sc_sub@meta.data %>% as.data.frame()

# Randomly sample 1000 cells in each group from each cluster. 
# If cluster has cells < 1000, get the cell id
cellid <- ""
gro <- unique(sc.meta$integrated_clusters)

for (i in gro) {
  temp1 <- sc.meta %>% filter(integrated_clusters == i)
  if (nrow(temp1) > 1000) {
    inTest <- sample(1:nrow(temp1), 1000, replace = F)
    id <- rownames(temp1)[inTest]
  } else if (nrow(temp1) < 1000) {
    id <- rownames(temp1)
  }
  cellid <- c(cellid, id)
}

# Subset the matrix, rownames of the sc.mtx should be the gene name
sc.mtx <- sc.mtx[, which(colnames(sc.mtx) %in% cellid)] %>% as.matrix()
sc.meta$cellid <- rownames(sc.meta)
sc.meta <- sc.meta[which(sc.meta$cellid %in% cellid),]

# Check the order and cell id are the same in labels and matrix
all(sc.meta$cellid == colnames(sc.mtx))
labels <- sc.meta$integrated_clusters

# Deconvolution
Signature <- buildSignatureMatrixMAST(sc.mtx, labels, "results")

allCounts_DWLS <- NULL
#allCounts_OLS <- NULL
#allCounts_SVR <- NULL
for (j in 1:(dim(bulk_df)[2]-1)) {
  S <- Signature
  Bulk <- bulk_df[,j+1]
  names(Bulk) <- bulk_df[,1]
  Genes <- intersect(rownames(S), names(Bulk))
  B <- Bulk[Genes]
  S <- S[Genes,]
  #solOLS <- solveOLS(S,B)
  solDWLS <- solveDampenedWLS(S,B)
  #solSVR <- solveSVR(S,B)
  
  allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
  #allCounts_OLS <- cbind(allCounts_OLS, solOLS)
  #allCounts_SVR <- cbind(allCounts_SVR, solSVR)
}

# Assign the column name to the final result matrix
colnames(allCounts_DWLS) <- colnames(bulk_df)[-1]
df <- as.data.frame(allCounts_DWLS)
df$celltype <- rownames(df)
df <- df %>% tidyr::pivot_longer(!celltype, names_to = "sample", values_to = "percent")  

info <- read.csv("/Users/tili/Desktop/PhD/Ovarioid/data/processed/info_Cortex.csv") %>% dplyr::select(name, culture, group)
info <- read.csv("/Users/tili/Desktop/PhD/Ovarioid/data/processed/info_Medulla.csv") %>% dplyr::select(name, culture, group)

anno <- rbind(info, info, info, info, info, info)

colnames(df)[2] <- "name"
df1 <- cbind(df, anno)
colnames(df1)[4] <- "sample"
all(df1$name == df1$sample)

ggplot(df1, aes(x = name, y = percent, fill = celltype)) + geom_col(position = "stack") + 
  facet_grid(~culture, scales = "free") +
  theme_bw() + #facet_zoom(ylim = c(0.98, 1)) + 
  xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stroma_theca" = "grey", "Granulosa" = "#00B6EB", "T_NK" = "#A58AFF",
                               "Monocytes" = "#f4a582",
                               "Endothelial" = "#53B400", "Perivascular" = "#C49A00"))

df1 <- df1 %>% dplyr::filter(culture != "Tissue")
df1 %>% dplyr::filter(culture == "3D") %>% 
  ggplot(aes(x = name, y = percent, fill = celltype)) + geom_col(position = "stack") + #facet_grid(~culture) +
  theme_bw() + facet_zoom(ylim = c(0.97, 1)) + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stroma_theca" = "grey", "Granulosa" = "#00B6EB", "T_NK" = "#A58AFF",
                               "Monocytes" = "#f4a582",
                               "Endothelial" = "#53B400", "Perivascular" = "#C49A00"))

df1 %>% dplyr::filter(culture == "2D") %>% 
  ggplot(aes(x = name, y = percent, fill = celltype)) + geom_col(position = "stack") + #facet_grid(~culture) +
  theme_bw() + facet_zoom(ylim = c(0.98, 1)) + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Stroma_theca" = "grey", "Granulosa" = "#00B6EB", "T_NK" = "#A58AFF",
                               "Monocytes" = "#f4a582",
                               "Endothelial" = "#53B400", "Perivascular" = "#C49A00"))
