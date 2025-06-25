# Figure 1: Single-cell Analysis of Integrated Ventricles -----------------

# 1. Setup ----------------------------------------------------------------
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Set working directory and random seed
set.seed(1)


# 2. Load Data ------------------------------------------------------------
integrated_ventricles <- readRDS("MAID_Integrated_PV_Region_MS_vs_CTRL.RDS")


# 3. Figure 1A: UMAP by CellType ------------------------------------------
DimPlot(
  object    = integrated_ventricles,
  reduction = "umap",
  label     = TRUE,
  pt.size   = 0.5,
  group.by  = "CellType"
)


# 4. Figure 1B: Define & Plot CellGroup ----------------------------------

## 4.1. Define mapping from CellType → broader CellGroup
cell_group_mapping <- c(
  ASTROCYTE               = "Astrocytes",
  NEURONS_1               = "Neurons",
  NEURONS_2               = "Neurons",
  GABA_NEURONS            = "Neurons",
  OLIGOS_1                = "Oligodendrocytes & OPCs",
  OLIGOS_2                = "Oligodendrocytes & OPCs",
  OLIGOS_3                = "Oligodendrocytes & OPCs",
  OPC                     = "Oligodendrocytes & OPCs",
  MACROPHAGE              = "Immune Cells",
  B_CELL                  = "Immune Cells",
  CYTOTOXIC_T             = "Immune Cells",
  IMMUNE_1                = "Immune Cells",
  IMMUNE_2                = "Immune Cells",
  IMMUNE_3                = "Immune Cells",
  IMMUNE_4                = "Immune Cells",
  IMMUNE_5                = "Immune Cells",
  IMMUNE_6                = "Immune Cells",
  EPENDYMAL               = "Ependymal Cells",
  ENDOTHELIAL             = "Vascular Cells",
  PERICYTE_SMOOTH_MUSCLE  = "Vascular Cells"
)

## 4.2. Assign CellGroup
integrated_ventricles$CellType <- as.character(integrated_ventricles$CellType)
integrated_ventricles$CellGroup <- unname(cell_group_mapping[integrated_ventricles$CellType])

# Verify assignment
table(integrated_ventricles$CellGroup)


## 4.3. UMAP with custom colors
custom_colors <- c(
  "Astrocytes"               = "#E41A1C",
  "Neurons"                  = "#377EB8",
  "Oligodendrocytes & OPCs"  = "#4DAF4A",
  "Immune Cells"             = "#984EA3",
  "Ependymal Cells"          = "#FF7F00",
  "Vascular Cells"           = "#FFFF33"
)

umap_plot <- DimPlot(
  object    = integrated_ventricles,
  reduction = "umap",
  label     = TRUE,
  pt.size   = 0.5,
  group.by  = "CellGroup"
) +
  scale_color_manual(values = custom_colors)

# Save result
ggsave(
  filename = "UMAP_CellGroup_CustomColors.png",
  plot     = umap_plot,
  width    = 8, height = 6, dpi = 300
)


# 5. Figure 1C: DotPlot of Marker Genes -----------------------------------

features <- c(
  "AQP4",   "SLC1A2",  # Astrocytes
  "FOXJ1",  "PIFO",    # Ependymal
  "C1QA",   "CD68",    # Immune
  "SYT1",   "SYNPR",   # Neurons
  "CNP",    "MOG",     # Oligodendrocytes
  "PECAM1", "MCAM"     # Vascular
)

DefaultAssay(integrated_ventricles) <- "RNA"

dotplot <- DotPlot(
  object   = integrated_ventricles,
  features = features,
  group.by = "CellGroup"
) +
  scale_color_gradient(low = "yellow", high = "red") +
  scale_size(range = c(2, 8)) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title   = element_blank()
  ) +
  labs(color = "Avg. Expr.", size = "Pct. Expr.")

print(dotplot)


# 6. Figure 1D: CellGroup Proportions per Patient -------------------------

## 6.1. Summarize counts & proportions
cell_counts <- integrated_ventricles@meta.data %>%
  group_by(Patient, CellGroup) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Patient) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

## 6.2. Rename patients
sample_mapping <- c(
  P251 = "HSP1",  P252 = "ALS1", P257 = "ALS2", P261 = "ALS3",
  P253 = "MS1",   P259 = "MS2",  P276 = "MS3",  P280 = "MS4"
)
cell_counts$Patient <- recode(cell_counts$Patient, !!!sample_mapping)

## 6.3. Stacked bar plot
bar_plot <- ggplot(cell_counts, aes(x = Patient, y = Proportion, fill = CellGroup)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 20) +
  labs(
    title = "Proportion of Each Cell Group Per Patient",
    x     = "Patient",
    y     = "Proportion",
    fill  = "Cell Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 22),
    axis.text.y = element_text(size = 22),
    axis.title  = element_text(size = 26, face = "bold"),
    legend.text = element_text(size = 22),
    legend.title= element_text(size = 26, face = "bold"),
    plot.title  = element_text(size = 30, face = "bold", hjust = 0.5)
  )

ggsave("CellGroup_Proportion_Per_Patient.tiff",
       plot   = bar_plot,
       width  = 10, height = 8, dpi = 300)


# 7. Reclustering Ependymal Cells -----------------------------------------

# 7.1. Create new Seurat object from raw counts
raw_counts <- GetAssayData(integrated_ventricles, slot = "counts")
new_obj    <- CreateSeuratObject(counts = raw_counts)

# 7.2. Transfer metadata
meta_cols <- c("seurat_clusters", "CellType", "Patient", "Group")
for (col in meta_cols) {
  new_obj[[col]] <- integrated_ventricles@meta.data[col, drop = FALSE][colnames(new_obj), ]
}

# 7.3. Subset to ependymal cells
ependymal <- subset(new_obj, subset = CellType == "EPENDYMAL")

# 7.4. Standard preprocessing & clustering
ependymal <- ependymal %>%
  NormalizeData(norm.method = "LogNormalize", scale.factor = 1e4) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.)) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.75) %>%
  RunUMAP(dims = 1:30)

# 7.5. UMAP plots
DimPlot(ependymal, reduction = "umap", label = TRUE, pt.size = 1, group.by = "Group")
ggsave("ependymal_umap.png", width = 12, height = 8, dpi = 300)


# 8. QC Violin Plots ------------------------------------------------------

## 8.1. Compute percent mito
mito_genes <- grep("^MT-", rownames(ependymal), value = TRUE)
ependymal$percent.mito <- 
  Matrix::colSums(GetAssayData(ependymal, slot = "counts")[mito_genes, ]) /
  Matrix::colSums(GetAssayData(ependymal, slot = "counts")) * 100

meta_data <- ependymal@meta.data

## 8.2. Plotting function
plot_violin <- function(df, var, label, file) {
  p <- ggplot(df, aes_string(x = "Sample", y = var, fill = "Sample")) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white", alpha = 0.7) +
    theme_minimal() +
    labs(title = paste(label, "Per Patient"), x = "Patient", y = label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3")
  
  ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
}

# Generate QC plots
plot_violin(meta_data, "percent.mito",    "Mitochondrial Percentage", "violin_mito_percentage.png")
plot_violin(meta_data, "nCount_RNA",      "UMI Count",                "violin_UMI_count.png")
plot_violin(meta_data, "nFeature_RNA",    "nFeature RNA Count",       "violin_nFeature_RNA.png")
