# scRNA-seq Embryo Reference Atlas + Query Mapping
# --------------------------------------------------
# 1. Builds an integrated human pre-implantation embryo reference atlas
#    from THREE published datasets:
#      - Yan et al. (2013)        PMID: 23934149
#      - Petropoulos et al. (2016) PMID: 27062923
#      - Mazid et al. (2022)       PMID: 35314832
#    using Seurat's anchor-based integration
#    (Stuart, Butler et al., Cell 2019).
# 2. Projects our scRNA-seq data onto the reference using MapQuery.
# 3. Produces three UMAP plots (atlas only, atlas + query, atlas grey + query).
#
# Note on parameters: nfeatures used in BOTH FindVariableFeatures (per dataset)
# AND anchor.features (in FindIntegrationAnchors) is 3000, matching the
# manuscript text.

library(Seurat)
library(dplyr)
library(ggplot2)

set.seed(42)

# ============================================================
# 1. Load datasets
# ============================================================

# ---- Yan et al. 2013 ----
yan_counts <- read.csv("path/to/Yan_2013_PMID_23934149.counts.gz",
                       row.names = 1, check.names = FALSE)
yan_meta <- read.table("path/to/Yan_2013_PMID_23934149.meta.tsv",
                       sep = "\t", header = TRUE)
yan_meta_r <- yan_meta
rownames(yan_meta_r) <- yan_meta_r$cell
yan_meta_r$cell <- NULL

yan_seu <- CreateSeuratObject(counts    = as.matrix(yan_counts),
                              meta.data = yan_meta_r,
                              project   = "Yan2013")
yan_seu$dataset <- "Yan_2013"
yan_seu$stage   <- yan_seu$devTime

# ---- Petropoulos et al. 2016 ----
pet_counts <- read.csv("path/to/Petropoulos_2016_PMID_27062923.counts.gz",
                       row.names = 1, check.names = FALSE)
pet_meta <- read.table("path/to/Petropoulos_2016_PMID_27062923.meta.tsv",
                       sep = "\t", header = TRUE)
pet_meta_r <- pet_meta
rownames(pet_meta_r) <- pet_meta_r$cell
pet_meta_r$cell <- NULL

pet_seu <- CreateSeuratObject(counts    = as.matrix(pet_counts),
                              meta.data = pet_meta_r,
                              project   = "Petropoulos2016")
pet_seu$dataset <- "Petropoulos_2016"
pet_seu$stage   <- pet_seu$raw_annotation

# ---- Mazid et al. 2022 ----
mazid_counts <- read.csv("path/to/Mazid_2022_PMID_35314832.counts.gz",
                         row.names = 1, check.names = FALSE)
mazid_meta <- read.table("path/to/Mazid_2022_PMID_35314832.meta.tsv",
                         sep = "\t", header = TRUE)
mazid_meta_r <- mazid_meta
rownames(mazid_meta_r) <- mazid_meta_r$cell
mazid_meta_r$cell <- NULL

mazid_seu <- CreateSeuratObject(counts    = as.matrix(mazid_counts),
                                meta.data = mazid_meta_r,
                                project   = "Mazid2022")
mazid_seu$dataset <- "Mazid_2022"
mazid_seu$stage   <- mazid_seu$raw_annotation

# ---- Our query data ----
obj <- readRDS("path/to/prbjn.rds")
obj$dataset <- "Our_data"
obj$stage   <- obj$sample

cat("Yan:",         dim(yan_seu),   "\n")
cat("Petropoulos:", dim(pet_seu),   "\n")
cat("Mazid:",       dim(mazid_seu), "\n")
cat("Our data:",    dim(obj),       "\n")

# ============================================================
# 2. Common genes + normalization (for the reference)
# ============================================================
common_genes_ref <- Reduce(intersect, list(
  rownames(yan_seu),
  rownames(pet_seu),
  rownames(mazid_seu)
))
cat("Reference common genes:", length(common_genes_ref), "\n")

# nfeatures = 3000 (matches the manuscript text and anchor.features below)
yan_r   <- yan_seu[common_genes_ref, ]   %>%
            NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
pet_r   <- pet_seu[common_genes_ref, ]   %>%
            NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)
mazid_r <- mazid_seu[common_genes_ref, ] %>%
            NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)

# ============================================================
# 3. Build reference atlas (Seurat anchor-based integration)
# ============================================================
cat("Computing integration anchors...\n")
anchors_ref <- FindIntegrationAnchors(
  object.list     = list(yan_r, pet_r, mazid_r),
  dims            = 1:20,
  anchor.features = 3000   # consistent with the per-dataset nfeatures
)

atlas3 <- IntegrateData(anchorset = anchors_ref, dims = 1:30)

DefaultAssay(atlas3) <- "integrated"
atlas3 <- ScaleData(atlas3)
atlas3 <- RunPCA(atlas3, npcs = 30)

# ---- Petropoulos devTime correction (with sanity check) ----
pet_devtime <- setNames(pet_meta$devTime, pet_meta$cell)
pet_cells   <- rownames(atlas3@meta.data[atlas3@meta.data$dataset == "Petropoulos_2016", ])
pet_cells_clean <- gsub("_[0-9]+$", "", pet_cells)

# Sanity check: every cleaned cell ID should map back to the metadata
unmatched_pet <- pet_cells_clean[!pet_cells_clean %in% names(pet_devtime)]
if (length(unmatched_pet) > 0) {
  warning("Some Petropoulos cell IDs did not match metadata after stripping _<int>$. ",
          "First few: ", paste(head(unmatched_pet), collapse = ", "))
}

atlas3@meta.data$stage_corrected <- atlas3@meta.data$stage
atlas3@meta.data[pet_cells, "stage_corrected"] <- pet_devtime[pet_cells_clean]

# ---- Yan raw_annotation correction (with sanity check) ----
yan_meta_map    <- setNames(yan_meta$raw_annotation, yan_meta$cell)
yan_cells       <- rownames(atlas3@meta.data[atlas3@meta.data$dataset == "Yan_2013", ])
yan_cells_clean <- gsub("_[0-9]+$", "", yan_cells)

unmatched_yan <- yan_cells_clean[!yan_cells_clean %in% names(yan_meta_map)]
if (length(unmatched_yan) > 0) {
  warning("Some Yan cell IDs did not match metadata after stripping _<int>$. ",
          "First few: ", paste(head(unmatched_yan), collapse = ", "))
}

atlas3@meta.data[yan_cells, "stage_corrected"] <- yan_meta_map[yan_cells_clean]

# ---- UMAP + clustering ----
atlas3 <- RunUMAP(atlas3,
                  dims          = 1:30,
                  min.dist      = 1.1,
                  spread        = 3,
                  n.neighbors   = 70,
                  return.model  = TRUE,
                  seed.use      = 42)

atlas3 <- FindNeighbors(atlas3, dims = 1:30)
atlas3 <- FindClusters(atlas3, resolution = 0.5)
cat("Atlas built:", dim(atlas3), "\n")

# ============================================================
# 4. Map our query data onto the reference (MapQuery)
# ============================================================
obj_query <- obj[intersect(common_genes_ref, rownames(obj)), ] %>%
  NormalizeData() %>% FindVariableFeatures(nfeatures = 3000)

cat("Computing transfer anchors...\n")
anchors_query <- FindTransferAnchors(
  reference           = atlas3,
  query               = obj_query,
  dims                = 1:30,
  reference.reduction = "pca"
)

obj_mapped <- MapQuery(
  anchorset           = anchors_query,
  query               = obj_query,
  reference           = atlas3,
  refdata             = list(stage = "stage_corrected"),
  reference.reduction = "pca",
  reduction.model     = "umap"
)
cat("Query mapped\n")

cat("Predicted stage distribution:\n")
print(table(obj_mapped$predicted.stage, obj_mapped$sample))

# ============================================================
# 5. Plotting
# ============================================================

# ---- Atlas UMAP data ----
umap3 <- as.data.frame(Embeddings(atlas3, "umap"))
umap3$dataset         <- atlas3@meta.data$dataset
umap3$stage_corrected <- atlas3@meta.data$stage_corrected

umap3$stage_label <- case_when(
  umap3$stage_corrected == "8 cell"   ~ "8C",
  umap3$stage_corrected == "2-4 cell" ~ "2-4C",
  umap3$stage_corrected == "Morula"   ~ "Morula",
  umap3$stage_corrected == "E3"       ~ "E3",
  umap3$stage_corrected == "E4"       ~ "E4",
  umap3$stage_corrected == "E5"       ~ "E5",
  umap3$stage_corrected == "E6"       ~ "E6",
  umap3$stage_corrected == "E7"       ~ "E7",
  umap3$stage_corrected == "8CLC"     ~ "Sorted 8CLC",
  umap3$stage_corrected == "4CL"      ~ "4CL-D12",
  umap3$stage_corrected == "e4CL"     ~ "e4CL-D5",
  umap3$stage_corrected == "Primed"   ~ "Primed",
  umap3$stage_corrected == "hES_P10"  ~ "hESC P10",
  umap3$stage_corrected == "hES_P0"   ~ "hESC P0",
  umap3$stage_corrected %in% c("Epiblast", "Hypoblast", "TE") ~ "Late blastocyst",
  umap3$stage_corrected %in% c("Oocyte", "Zygote") ~ NA_character_,
  TRUE ~ umap3$stage_corrected
)

umap3_plot <- umap3[!is.na(umap3$stage_label), ]
umap3_plot$umap_2_flip <- -umap3_plot$umap_2

# ---- Query UMAP data ----
query_umap <- as.data.frame(Embeddings(obj_mapped, "ref.umap"))
colnames(query_umap) <- c("umap_1", "umap_2")
query_umap$stage_label <- obj_mapped$sample
query_umap$umap_2_flip <- -query_umap$umap_2

# ---- Color palette ----
stage_colors <- c(
  "Sorted 8CLC"     = "#F5C400",
  "e4CL-D5"         = "#7FCDCD",
  "4CL-D12"         = "#4DA6D4",
  "Primed"          = "#2166AC",
  "E3"              = "#D73027",
  "E4"              = "#E86530",
  "E5"              = "#F4A442",
  "E6"              = "#FDCC5C",
  "E7"              = "#9DC96A",
  "8C"              = "#4DAF4A",
  "2-4C"            = "#B7E4C7",
  "Morula"          = "#1B7837",
  "Late blastocyst" = "#084594",
  "hESC P10"        = "#003218",
  "hESC P0"         = "#1B4332",
  "DMSO"            = "#2D6A4F",
  "PRBJN"           = "#D62828"
)

# ---- Label positions ----
label_pos_atlas <- umap3_plot %>%
  group_by(stage_label) %>%
  summarise(umap_1 = median(umap_1), umap_2_flip = median(umap_2_flip))

label_pos_query <- query_umap %>%
  group_by(stage_label) %>%
  summarise(umap_1 = median(umap_1), umap_2_flip = median(umap_2_flip))

label_pos_all <- rbind(label_pos_atlas, label_pos_query)

# ============================================================
# Plot 1: Reference atlas only
# ============================================================
p_atlas <- ggplot(umap3_plot,
                  aes(x = umap_1, y = umap_2_flip,
                      color = stage_label, shape = dataset)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_text(data = label_pos_atlas,
            aes(x = umap_1, y = umap_2_flip, label = stage_label),
            color = "black", size = 3.5, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = stage_colors, name = "Samples") +
  scale_shape_manual(
    values = c("Yan_2013" = 18, "Petropoulos_2016" = 15, "Mazid_2022" = 17),
    name   = "Dataset"
  ) +
  labs(title = "Reference Atlas", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12) +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text     = element_text(size = 9)) +
  guides(color = guide_legend(override.aes = list(size = 3), order = 1),
         shape = guide_legend(override.aes = list(size = 3), order = 2))

# ============================================================
# Plot 2: Atlas + query (colored)
# ============================================================
p_full <- ggplot() +
  geom_point(data = umap3_plot,
             aes(x = umap_1, y = umap_2_flip,
                 color = stage_label, shape = dataset),
             size = 1.2, alpha = 0.7) +
  geom_point(data = query_umap,
             aes(x = umap_1, y = umap_2_flip, color = stage_label),
             shape = 16, size = 1.5, alpha = 0.9) +
  geom_text(data = label_pos_all,
            aes(x = umap_1, y = umap_2_flip, label = stage_label),
            color = "black", size = 3.5, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = stage_colors, name = "Samples") +
  scale_shape_manual(
    values = c("Yan_2013" = 18, "Petropoulos_2016" = 15, "Mazid_2022" = 17),
    name   = "Dataset"
  ) +
  labs(title = "Reference Atlas + Our Data (MapQuery)",
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12) +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text     = element_text(size = 9)) +
  guides(color = guide_legend(override.aes = list(size = 3), order = 1),
         shape = guide_legend(override.aes = list(size = 3), order = 2))

# ============================================================
# Plot 3: Atlas grey + query (colored)
# ============================================================
p_ourdata <- ggplot() +
  geom_point(data = umap3_plot,
             aes(x = umap_1, y = umap_2_flip),
             color = "grey80", size = 0.8, alpha = 0.5) +
  geom_point(data = query_umap,
             aes(x = umap_1, y = umap_2_flip, color = stage_label),
             shape = 16, size = 1.8, alpha = 0.9) +
  geom_text(data = label_pos_query,
            aes(x = umap_1, y = umap_2_flip, label = stage_label),
            color = "black", size = 4, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_manual(values = stage_colors, name = "Samples") +
  labs(title = "Our Data on Reference Atlas",
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 12) +
  theme(legend.key.size = unit(0.5, "cm"),
        legend.text     = element_text(size = 9))

# ============================================================
# Save plots and atlas object
# ============================================================
ggsave("atlas_reference_only.png",   p_atlas,   width = 11, height = 7, dpi = 300)
ggsave("atlas_mapquery_full.png",    p_full,    width = 11, height = 7, dpi = 300)
ggsave("atlas_mapquery_ourdata.png", p_ourdata, width = 11, height = 7, dpi = 300)
cat("3 plots saved\n")

saveRDS(atlas3, "atlas3_reference.rds")
cat("Atlas saved: atlas3_reference.rds\n")
