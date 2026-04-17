set.seed(1234)
polychrome_26 <- DiscretePalette_scCustomize(num_colors = 26, palette = "polychrome", shuffle_pal = T)
donor_palette <- c(
  "#E41A1C", # strong red
  "#377EB8", # blue
  "#4DAF4A", # green
  "#984EA3", # purple
  "#FF7F00", # orange
  "#FFFF33", # yellow (deep, not pastel)
  "#276419", # Dark Olive Green
  "#F781BF", # pink
  "#999933", # olive
  "#66C2A5", # teal
  "#FC8D62", # salmon
  "#8DA0CB", # periwinkle
  "#E78AC3", # magenta
  "#A6D854", # lime
  "#FFD92F", # golden yellow
  "#7E1137",  # Dark Wine Red
  "#1B9E77", # dark teal
  "#D95F02", # burnt orange
  "#7570B3", # indigo
  "#66A61E", # olive green
  "#E6AB02", # mustard
  "#B15928", # reddish brown
  "#3288BD", # cyan-blue
  "#D53E4F", # crimson
  "#5E3C99"  # deep violet
)



pal26 <- c(
  # Blues
  "#060580",  # navy
  "#0021CF",  # royal blue
  "#3D63B7",  # blue
  "#0085FF",  # azure
  
  # Cyan / Teal
  "#00C6CF",  # cyan
  "#07848D",  # teal
  "#84D6C9",  # aqua
  "#4BB580",  # green-teal
  
  # Greens
  "#32953F",  # green
  "#1A4C20",  # dark green
  #"#8FF547",  # neon green
  "#7A8F00",  # olive
  
  # Yellows
  "#E6C800",  # yellow
  "#FFEB61",  # light yellow
  
  # Oranges
  "#A63A00",  # burnt orange
  "#F56D00",  # vivid orange
  "#FF9138",  # orange
  
  # Reds
  "#CF142B",  # red
  "#8F0032",  # wine red
  "#650B39",  # maroon
  
  # Pinks / Magentas
  "#E7298A",  # pink
  "#FF94A5",  # light pink
  "#B80093",  # magenta
  "#F000FF",  # neon magenta
  
  # Purples
  "#9A6DCF",  # lavender
  "#6D009E",   # purple
  "gray9" #dark gray for a special last category
  )

#PalettePlot(new_pal)

dim_palette <- c("#D55E00", # Dark Orange
                 "#0072B2", # Dark Blue
                 "#009E73", # Dark Bluish Green
                 "#276419", # Dark Olive Green
                 "#CC79A7", # Dark Reddish Purple
                 "#E69F00", # Dark Gold
                 "#56B4E9", # Dark Sky Blue
                 "#004D47", # Dark Green
                 "#882E72", # Dark Magenta
                 "#1965B0", # Dark Cerulean
                 "#DC050C", # Dark Crimson
                 "#7E1137"  # Dark Wine Red
)


pt_ctl <- c( "#007dbc", "#ffab40")
mf <-  c("#900C3F", "#237873")

colors_ident <- c("#000000","#dbdbdb", "#f71c05","#ff63ce","#49fa44", "#3792f2")
polychrome_pal <- sample(DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome"))


###TO GET FEATURES WITHOUT .##
rename_features <- function(data_path){
  if(file.exists(file.path(data_path, "features_original.tsv"))){
    return()
  }
  
  file.copy(file.path(data_path, "features.tsv"), file.path(data_path, "features_original.tsv"))
  df <- read_delim(file.path(data_path, "features.tsv"), delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) %>%
    dplyr::rename(ensembl_id = X1) %>%
    mutate(ensembl_id = str_remove(ensembl_id, "\\.\\d+(?=_PAR_Y|$)")) %>%
    relocate(ensembl_id, .before = 1)
  
  write_tsv(df, col_names = FALSE, file.path(data_path, "features.tsv"))
}

custom_read_solo <- function(filemetadata, .village_id) {
  file_paths <- filemetadata %>%
    dplyr::filter(village_id == .village_id) # if you don't have a village short id, you can use any unique identifier
  ReadSTARsolo(file_paths$path_to_star_solo[[1]], feature.column = 1, unique.features = TRUE)
}

NPCvsOther <- read_csv("/u/project/mfwells/SHARED/SCRIPTS/scrnaseq/NPCvsOther_genelist.csv") %>% arrange(celltype_or_region)
NPCshortlist <- read_csv("/u/project/mfwells/SHARED/SCRIPTS/scrnaseq/NPCmarkers_shortgenelist.csv")
nehme <- read_csv("/u/project/mfwells/SHARED/SCRIPTS/scrnaseq/nehme_markers_ensembl.csv") %>%
  dplyr::rename(gene_name = gene_symbol, 
                celltype_or_region = class,
                gene_id = ensembl_id)

ensg_to_symbol <- read_tsv("/u/project/mfwells/SHARED/SCRIPTS/scrnaseq/ensembl_to_symbol.tsv")
master_genelist <- read_csv("/u/project/mfwells/SHARED/SCRIPTS/scrnaseq/master_genelist.csv")


#----------------------------------------------------------------------------------------------------------------------------
#FUNCTIONS FOR COUNTING AND PLOTTING QC METRICS
#----------------------------------------------------------------------------------------------------------------------------
#Function: cells_per_condition_barchart
#Args: gex, seurat object
#      label, string with save info
#      var, string variable to separate by
#      ylimit, int for y limit if needed, default is NA so auto limit
cells_per_condition_barchart <- function(gex, label, var, colorscheme = polychrome_pal, ylimit = NULL, save = TRUE){
  meta_data <- Fetch_Meta(object = gex)
  prefilter_ncells <-
    meta_data %>%
    ggplot(aes(x = orig.ident, fill = !!sym(var))) +
    geom_bar() +
    geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = colorscheme) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("Number of Cells")
  
  if (!is.null(ylim)) {
    prefilter_ncells <- prefilter_ncells + expand_limits(y = ylimit)
  }
  
  print(prefilter_ncells)
  
  if(save){
  savename <- paste0("QC_Plots_NumCellsPerCondition_", label, ".png")
  ggsave(here(output_dir, savename), width = 4, height = 4)}
}

#Function: qc_violins
#Args: gex, seurat object
#      label, string with save info
qc_violins <- function(gex, label, save = TRUE){
  p1 <- QC_Plots_Genes(seurat_object = gex)
  p2 <- QC_Plots_UMIs(seurat_object = gex, y_axis_log = TRUE)
  p3 <- QC_Plots_Mito(seurat_object = gex)
  p4 <- QC_Plots_Complexity(seurat_object = gex, high_cutoff = 1)
  composite <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2, nrow = 2)
  print(composite)
  
  if(save) {
  savename <- paste0("QC_Plots_Genes_UMIs_Mito_Complex_", label, ".png")
  ggsave(here(output_dir, savename), width = 8.5, height = 11)}
}

#Function: qc_scatters
#Args: gex, seurat object
#      label, string with save info
qc_scatters <- function(gex, label, colors = colors_ident, save = TRUE){
  p1 <- FeatureScatter(object = gex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", shuffle = TRUE) +
    labs(title = "Genes vs UMIs per Cell", subtitle = "", x = "UMIs per Cell", y = "Genes per Cell") +
    scale_color_manual(values = colors)
  p2 <- QC_Plot_UMIvsGene(seurat_object = gex, meta_gradient_name = "percent_mito")  +
    labs(title = "Genes vs UMIs per Cell", subtitle = "", x = "UMIs per Cell", y = "Genes per Cell")
  p3 <- FeatureScatter(object = gex, feature1 = "nCount_RNA", feature2 = "percent_mito", group.by = "orig.ident", shuffle = TRUE) +
    labs(title = "UMIs vs % Mito. Expr. per Cell", x = "UMIs per Cell", y = "Mitochondrial Read Percentage") +
    scale_color_manual(values = colors)
  p4 <- FeatureScatter(object = gex, feature1 = "log10GenesPerUMI", feature2 = "percent_mito", group.by = "orig.ident", shuffle = TRUE) +
    labs(title = "Complexity vs % Mito. Expr. per Cell", x = "log10GenesPerUMI", y = "Mitochondrial Read Percentage") +
    scale_color_manual(values = colors)
  composite <- patchwork::wrap_plots(p1, p2, p3, p4,  ncol = 2, nrow = 2)
  print(composite)
  
  if(save) {
  savename <- paste0("QC_Plots_ComprehensiveFeatureScatter_", label, ".png")
  ggsave(here(output_dir, savename), width = 9, height = 11)}
}

#----------------------------------------------------------------------------------------------------------------------------
#FUNCTIONS FOR DECIDING AND APPLYING FILTERING THRESHOLDS
#----------------------------------------------------------------------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

#Note: all of the following functions use vectors called "thresholds" with the EXACT following format:
#c(high_feature, low_umi, high_umi, low_ribo, high_ribo, high_mito, complexity)

#Function: threshold_scatter_density
#Args: gex, seurat object
#      label, string with save info
#      thresholds, vector with the format low_feature, high_feature, low_umi, high_umi, low_ribo, high_ribo, high_mito, complexity
threshold_scatter_density <- function(gex, label, thresholds, save = TRUE){
  low_feature <- thresholds[["low_feature"]]  
  high_feature <- thresholds[["high_feature"]] 
  low_umi <- thresholds[["low_umi"]] 
  high_umi <- thresholds[["high_umi"]] 
  low_ribo <- thresholds[["low_ribo"]] 
  high_ribo <- thresholds[["high_ribo"]] 
  high_mito <- thresholds[["high_mito"]]
  complexity <- thresholds[["complexity"]]
  
  meta_data <- gex@meta.data
  #Examine the potential filtering thresholds and where they fall
  p1 <- QC_Plot_UMIvsGene(seurat_object = gex, meta_gradient_name = "percent_mito", 
                          low_cutoff_gene = low_feature, high_cutoff_gene = high_feature,
                          low_cutoff_UMI = low_umi, high_cutoff_UMI = high_umi)  +
    labs(title = "Genes vs UMIs per Cell", subtitle = "", x = "UMIs per Cell", y = "Genes per Cell")
  
  p2 <- QC_Plot_UMIvsGene(seurat_object = gex, meta_gradient_name = "percent_mito", 
                          low_cutoff_gene = low_feature, high_cutoff_gene = high_feature,
                          low_cutoff_UMI = low_umi, high_cutoff_UMI = high_umi)  +
    scale_x_log10()+
    scale_y_log10()+
    labs(title = "Genes vs UMIs per Cell", subtitle = "", x = "log10(UMIs per Cell)", y = "log10(Genes per Cell)")
  
  p3 <- meta_data %>%
    ggplot(aes(x = nFeature_RNA, fill = orig.ident)) +
    geom_density(alpha = 0.2) +
    labs(title = "Features (Genes) detected", x = "Features detected per cell", y = "Density") +
    geom_vline(xintercept = low_feature, linetype = 2) + 
    geom_vline(xintercept = high_feature, linetype = 2) + 
    theme_minimal()
  
  p4 <- meta_data %>%
    ggplot(aes(x = nCount_RNA, fill = orig.ident)) +
    geom_density(alpha = 0.2) +
    labs(title = "UMI Counts per cell", x = "nCount_RNA", y = "Density") +
    scale_x_log10() +
    geom_vline(xintercept = low_umi, linetype = 2) + 
    geom_vline(xintercept = high_umi, linetype = 2) + #choose a line that represents your possible nCount UMI threshold
    theme_minimal() 
  
  p5 <- meta_data %>%
    ggplot(aes(x = percent_mito, fill = orig.ident)) +
    geom_density(alpha = 0.2) +
    labs(title = "Percent Mitochondrial Expression") +
    scale_x_log10()+
    geom_vline(xintercept = high_mito, linetype = 2) + #choose a line that represents your possible mitochondrial threshold
    theme_minimal() 
  
  p6 <- meta_data %>%
    ggplot(aes(x = percent_ribo, fill = orig.ident)) +
    geom_density(alpha = 0.2) +
    labs(title = "Percent Ribosomal Expression") +
    #scale_x_log10()+
    geom_vline(xintercept = low_ribo, linetype = 2) + #choose a line that represents your possible lower ribosomal threshold
    geom_vline(xintercept = high_ribo, linetype = 2) + #choose a line that represents your possible upper ribosomal threshold
    theme_minimal() 
  
  composite <- patchwork::wrap_plots(p1, p2, p3, p4,p5,p6, ncol = 2, nrow = 3)
  print(composite)
  
  if(save) {
  savename <- paste0("QC_Plots_FeatureScatter&DensityPlots_ExamineThresholds_", label, ".png")
  ggsave(here(output_dir, savename), width = 9, height = 11)}
}

#Function: threshold_qc_violin
#Args: gex, seurat object
#      label, string with save info
#      thresholds, vector with the format low_feature, high_feature, low_umi, high_umi, low_ribo, high_ribo, high_mito, complexity
threshold_qc_violin <- function(gex, label, thresholds, save = TRUE){
  low_feature <- thresholds[["low_feature"]]  
  high_feature <- thresholds[["high_feature"]] 
  low_umi <- thresholds[["low_umi"]] 
  high_umi <- thresholds[["high_umi"]] 
  low_ribo <- thresholds[["low_ribo"]] 
  high_ribo <- thresholds[["high_ribo"]] 
  high_mito <- thresholds[["high_mito"]]
  complexity <- thresholds[["complexity"]]
  
  p1 <- QC_Plots_Genes(seurat_object = gex , low_cutoff = low_feature, high_cutoff = high_feature) 
  p2 <- QC_Plots_UMIs(seurat_object = gex , low_cutoff = low_umi, high_cutoff = high_umi, y_axis_log = TRUE)
  p3 <- QC_Plots_Mito(seurat_object = gex, high_cutoff = high_mito)
  p4 <- QC_Plots_Complexity(seurat_object = gex , high_cutoff = complexity)
  composite <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2, nrow = 2)
  print(composite)
  
  if(save) {
  savename <- paste0("QC_Plots_Genes_UMIs_Mito_Complex_ExamineThresholds_", label, ".png")
  ggsave(here(output_dir, savename), width = 8.5, height = 8.5)}
}

#Function: count_remaining_cells
#Args: gex, seurat object
#      thresholds, vector with the format low_feature, high_feature, low_umi, high_umi, low_ribo, high_ribo, high_mito, complexity
count_remaining_cells <- function(gex, thresholds){
  low_feature <- thresholds[["low_feature"]]  
  high_feature <- thresholds[["high_feature"]] 
  low_umi <- thresholds[["low_umi"]] 
  high_umi <- thresholds[["high_umi"]] 
  low_ribo <- thresholds[["low_ribo"]] 
  high_ribo <- thresholds[["high_ribo"]] 
  high_mito <- thresholds[["high_mito"]]
  complexity <- thresholds[["complexity"]]
  
  print("PRE-THRESHOLDING:")
  print(table(Idents(gex)))
  
  print("AFTER THRESHOLDING:")
  print(table(Idents(subset(x = gex, 
                            subset = (nFeature_RNA >= low_feature) &
                              (nFeature_RNA <= high_feature) &
                              (nCount_RNA >= low_umi) &
                              (nCount_RNA <= high_umi) &
                              (percent_ribo > low_ribo) &
                              (percent_ribo < high_ribo) &
                              (percent_mito < high_mito) &
                              (log10GenesPerUMI >= complexity))))
  )
}

#Function: apply_thresholds
#Args: gex, seurat object
#      thresholds, vector with the format low_feature, high_feature, low_umi, high_umi, low_ribo, high_ribo, high_mito, complexity
apply_thresholds <- function(gex, thresholds){
  low_feature <- thresholds[["low_feature"]]  
  high_feature <- thresholds[["high_feature"]] 
  low_umi <- thresholds[["low_umi"]] 
  high_umi <- thresholds[["high_umi"]] 
  low_ribo <- thresholds[["low_ribo"]] 
  high_ribo <- thresholds[["high_ribo"]] 
  high_mito <- thresholds[["high_mito"]]
  complexity <- thresholds[["complexity"]]
  return(subset(x = gex, 
                subset = (nFeature_RNA >= low_feature) &
                  (nFeature_RNA <= high_feature) &
                  (nCount_RNA >= low_umi) &
                  (nCount_RNA <= high_umi) &
                  (percent_ribo > low_ribo) &
                  (percent_ribo < high_ribo) &
                  (percent_mito < high_mito) &
                  (log10GenesPerUMI >= complexity)
  )
  )
}

ensg_to_symbol <- read_tsv("/u/project/mfwells/SHARED/SCRIPTS/scrnaseq/ensembl_to_symbol.tsv")

ens2sym <- setNames(ensg_to_symbol$gene_name, ensg_to_symbol$gene_id)


# 2) Core relabeler for ONE ggplot
label_with_symbols <- function(p, mapping = ens2sym, strip_version = TRUE) {
  map_one <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x0 <- if (strip_version) sub("\\.\\d+$", "", x) else x
    y  <- mapping[x0]
    y[is.na(y)] <- x[is.na(y)]
    as.character(y)
  }
  map_vec <- function(v) vapply(v, map_one, character(1))
  
  p$labels$title <- map_one(p$labels$title)
  p$labels$x     <- map_one(p$labels$x)
  p$labels$y     <- map_one(p$labels$y)
  
  p <- p +
    scale_x_discrete(labels = function(x) map_vec(x)) +
    scale_y_discrete(labels = function(x) map_vec(x)) +
    scale_color_discrete(labels = function(x) map_vec(x)) +
    scale_fill_discrete(labels = function(x) map_vec(x)) +
    scale_shape_discrete(labels = function(x) map_vec(x)) +
    scale_linetype_discrete(labels = function(x) map_vec(x))
  
  if (inherits(p$facet, "Facet")) {
    p$facet$params$labeller <- labeller(.default = function(strings) map_vec(strings))
  }
  p
}

# 3) Universal wrapper
symbols_wrap <- function(x, ..., mapping = ens2sym, ncol = 3, strip_version = TRUE) {
  relabel_one <- function(p) label_with_symbols(p, mapping = mapping, strip_version = strip_version)
  
  # Case A: x is a plotting function (e.g., VlnPlot)
  if (is.function(x)) {
    f <- x
    args <- list(...)
    # If the function supports `combine`, default to FALSE so we can relabel each panel
    fmls <- tryCatch(names(formals(f)), error = function(e) character())
    if ("combine" %in% fmls && is.null(args$combine)) args$combine <- FALSE
    plots <- do.call(f, args)
    
    # Case B: x is an already-built ggplot or a list of ggplots
  } else {
    plots <- x
  }
  
  # Single ggplot
  if (inherits(plots, "ggplot")) return(relabel_one(plots))
  
  # List of ggplots (multi-feature outputs)
  if (is.list(plots) && length(plots) > 0 && all(vapply(plots, inherits, logical(1), "ggplot"))) {
    return(wrap_plots(lapply(plots, relabel_one), ncol = ncol))
  }
  
  stop("Input must be a ggplot, a list of ggplots, or a plotting function returning one of those.")
}

VlnPlot_symbols_clean <- function(
    seurat_obj,
    features,
    ...,
    mapping = ens2sym,
    ncol = 3,
    strip_version = TRUE,
    y_label = NULL,          # NULL = keep Seurat default, "" = remove, string = custom
    x_label = NULL,          # NULL = keep Seurat default, "" = remove, string = custom
    drop_legend = TRUE,
    rotate_x_ticks = NULL    # NULL = keep default; numeric angle (e.g., 45 or 90) to rotate
){
  plist <- VlnPlot(seurat_obj, features = features, combine = FALSE, ...)
  
  plist <- lapply(plist, function(p){
    p <- label_with_symbols(p, mapping = mapping, strip_version = strip_version)
    if (drop_legend) p <- p + theme(legend.position = "none")
    
    # y-axis
    if (!is.null(y_label)) {
      p <- p + if (identical(y_label, "")) ylab(NULL) else ylab(y_label)
    }
    # x-axis
    if (!is.null(x_label)) {
      p <- p + if (identical(x_label, "")) xlab(NULL) else xlab(x_label)
    }
    # rotate x tick labels (only if requested)
    if (!is.null(rotate_x_ticks)) {
      # simple sensible defaults for angled labels
      hjust <- if (rotate_x_ticks == 0) 0.5 else 1
      vjust <- if (rotate_x_ticks == 0) 0.5 else 1
      p <- p + theme(axis.text.x = element_text(angle = rotate_x_ticks, hjust = hjust, vjust = vjust))
    }
    
    p
  })
  
  wrap_plots(plist, ncol = ncol)
}


VlnPlot_symbols_shared_axes <- function(
    seurat_obj,
    features,
    ...,
    mapping = ens2sym,
    ncol = 3,
    strip_version = TRUE,
    shared_y_label = NULL,    # e.g., "Expression"; NULL = no shared Y label
    shared_x_label = NULL,    # e.g., "Clusters";   NULL = no shared X label
    drop_legend   = TRUE,
    rotate_x_ticks = NULL     # e.g., 45 or 90; NULL = default
){
  # Build individual panels
  plist <- VlnPlot(seurat_obj, features = features, combine = FALSE, ...)
  
  # Per-panel tweaks: relabel titles, drop legends, drop axis titles, rotate ticks
  plist <- lapply(plist, function(p){
    p <- label_with_symbols(p, mapping = mapping, strip_version = strip_version)
    if (drop_legend) p <- p + theme(legend.position = "none")
    p <- p + xlab(NULL) + ylab(NULL)
    if (!is.null(rotate_x_ticks)) {
      hjust <- if (rotate_x_ticks == 0) 0.5 else 1
      vjust <- if (rotate_x_ticks == 0) 0.5 else 1
      p <- p + theme(axis.text.x = element_text(angle = rotate_x_ticks, hjust = hjust, vjust = vjust))
    }
    p
  })
  
  # Arrange panels
  grid <- wrap_plots(plist, ncol = ncol)
  
  # Add shared axis labels (if provided)
  if (is.null(shared_x_label) && is.null(shared_y_label)) {
    return(grid)
  }
  
  out <- ggdraw(grid)
  if (!is.null(shared_y_label) && nzchar(shared_y_label)) {
    out <- out + draw_label(shared_y_label, x = 0.02, y = 0.5, angle = 90, vjust = 0.5, hjust = 0.5)
  }
  if (!is.null(shared_x_label) && nzchar(shared_x_label)) {
    out <- out + draw_label(shared_x_label, x = 0.5, y = 0.02, vjust = 0.5, hjust = 0.5)
  }
  out
}


ens2sym <- setNames(ensg_to_symbol$gene_name, ensg_to_symbol$gene_id)


# 2) Core relabeler for ONE ggplot
label_with_symbols <- function(p, mapping = ens2sym, strip_version = TRUE) {
  map_one <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x0 <- if (strip_version) sub("\\.\\d+$", "", x) else x
    y  <- mapping[x0]
    y[is.na(y)] <- x[is.na(y)]
    as.character(y)
  }
  map_vec <- function(v) vapply(v, map_one, character(1))
  
  p$labels$title <- map_one(p$labels$title)
  p$labels$x     <- map_one(p$labels$x)
  p$labels$y     <- map_one(p$labels$y)
  
  p <- p +
    scale_x_discrete(labels = function(x) map_vec(x)) +
    #scale_y_discrete(labels = function(x) map_vec(x)) +
    scale_color_discrete(labels = function(x) map_vec(x)) +
    scale_fill_discrete(labels = function(x) map_vec(x)) +
    scale_shape_discrete(labels = function(x) map_vec(x)) +
    scale_linetype_discrete(labels = function(x) map_vec(x))
  
  if (inherits(p$facet, "Facet")) {
    p$facet$params$labeller <- labeller(.default = function(strings) map_vec(strings))
  }
  p
}

# 3) Universal wrapper
symbols_wrap <- function(x, ..., mapping = ens2sym, ncol = 3, strip_version = TRUE) {
  relabel_one <- function(p) label_with_symbols(p, mapping = mapping, strip_version = strip_version)
  
  # Case A: x is a plotting function (e.g., VlnPlot)
  if (is.function(x)) {
    f <- x
    args <- list(...)
    # If the function supports `combine`, default to FALSE so we can relabel each panel
    fmls <- tryCatch(names(formals(f)), error = function(e) character())
    if ("combine" %in% fmls && is.null(args$combine)) args$combine <- FALSE
    plots <- do.call(f, args)
    
    # Case B: x is an already-built ggplot or a list of ggplots
  } else {
    plots <- x
  }
  
  # Single ggplot
  if (inherits(plots, "ggplot")) return(relabel_one(plots))
  
  # List of ggplots (multi-feature outputs)
  if (is.list(plots) && length(plots) > 0 && all(vapply(plots, inherits, logical(1), "ggplot"))) {
    return(wrap_plots(lapply(plots, relabel_one), ncol = ncol))
  }
  
  stop("Input must be a ggplot, a list of ggplots, or a plotting function returning one of those.")
}

VlnPlot_symbols_clean <- function(
    seurat_obj,
    features,
    ...,
    mapping = ens2sym,
    ncol = 1,
    strip_version = TRUE,
    y_label = NULL,          # NULL = keep Seurat default, "" = remove, string = custom
    x_label = NULL,          # NULL = keep Seurat default, "" = remove, string = custom
    drop_legend = TRUE,
    rotate_x_ticks = NULL    # NULL = keep default; numeric angle (e.g., 45 or 90) to rotate
){
  plist <- VlnPlot(seurat_obj, features = features, combine = FALSE, ...)
  
  plist <- lapply(plist, function(p){
    p <- label_with_symbols(p, mapping = mapping, strip_version = strip_version)
    if (drop_legend) p <- p + theme(legend.position = "none")
    
    # y-axis
    if (!is.null(y_label)) {
      p <- p + if (identical(y_label, "")) ylab(NULL) else ylab(y_label)
    }
    # x-axis
    if (!is.null(x_label)) {
      p <- p + if (identical(x_label, "")) xlab(NULL) else xlab(x_label)
    }
    # rotate x tick labels (only if requested)
    if (!is.null(rotate_x_ticks)) {
      # simple sensible defaults for angled labels
      hjust <- if (rotate_x_ticks == 0) 0.5 else 1
      vjust <- if (rotate_x_ticks == 0) 0.5 else 1
      p <- p + theme(axis.text.x = element_text(angle = rotate_x_ticks, hjust = hjust, vjust = vjust))
    }
    
    p
  })
  
  wrap_plots(plist, ncol = ncol)
}

DotPlot_symbols_clean <- function(
    seurat_obj,
    features,
    ...,
    mapping = ens2sym,   # either a named char vector or a 2-col data.frame
    strip_version = TRUE,
    y_label = NULL,
    x_label = NULL,
    drop_legend = FALSE,
    rotate_x_ticks = 45
){
  # --- normalize `mapping` to a named character vector: names = ENSG, values = SYMBOL
  to_named_vec <- function(m) {
    if (is.null(m)) return(setNames(character(), character()))
    if (is.character(m) && !is.null(names(m))) return(m)
    if (is.data.frame(m)) {
      id_col  <- intersect(names(m), c("gene_id","ensembl","ensg","id"))[1]
      sym_col <- intersect(names(m), c("gene_name","symbol","hgnc_symbol","name"))[1]
      stopifnot(length(id_col) == 1L, length(sym_col) == 1L)
      ids <- as.character(m[[id_col]])
      if (strip_version) ids <- sub("\\.\\d+$","", ids)
      return(setNames(as.character(m[[sym_col]]), ids))
    }
    stop("`mapping` must be a named character vector or a 2-column data.frame.")
  }
  map_vec <- function(v, mvec) {
    v <- as.character(v)
    key <- if (strip_version) sub("\\.\\d+$","", v) else v
    out <- unname(mvec[key])
    out[is.na(out) | out == ""] <- v[is.na(out) | out == ""]
    out
  }
  mvec <- to_named_vec(mapping)
  
  p <- DotPlot(seurat_obj, features = features, ...)
  
  # --- relabel factor levels directly (most robust)
  lev <- levels(p$data$features.plot)
  new_lev <- map_vec(lev, mvec)
  p$data$features.plot <- factor(p$data$features.plot, levels = lev, labels = new_lev)
  
  # (Optional) axis labels
  if (!is.null(y_label)) p <- p + if (identical(y_label, "")) ylab(NULL) else ylab(y_label)
  if (!is.null(x_label)) p <- p + if (identical(x_label, "")) xlab(NULL) else xlab(x_label)
  
  # (Optional) rotate x ticks
  if (!is.null(rotate_x_ticks)) {
    hjust <- if (rotate_x_ticks == 0) 0.5 else 1
    vjust <- if (rotate_x_ticks == 0) 0.5 else 1
    p <- p + theme(axis.text.x = element_text(angle = rotate_x_ticks, hjust = hjust, vjust = vjust))
  }
  
  if (drop_legend) p <- p + theme(legend.position = "none")
  p
}

FeaturePlot_symbols_clean <- function(
    seurat_obj,
    features,
    ...,
    mapping       = ens2sym,
    ncol          = 1,
    strip_version = TRUE,
    y_label       = NULL,   # NULL keep, "" remove, or string
    x_label       = NULL,   # NULL keep, "" remove, or string
    drop_legend   = TRUE,
    rotate_x_ticks = NULL
){
  # map helper (same behavior as your vln version)
  map_one <- function(x) {
    if (is.null(x)) return(NULL)
    x <- as.character(x)
    x0 <- if (strip_version) sub("\\.\\d+$", "", x) else x
    y  <- unname(mapping[x0])
    ifelse(is.na(y), x, y)
  }
  map_vec <- function(v) vapply(v, map_one, character(1), USE.NAMES = FALSE)
  
  # build list so we can retitle
  plist <- Seurat::FeaturePlot(seurat_obj, features = features, combine = FALSE, ...)
  
  # precompute display names
  disp <- map_vec(features)
  
  plist <- Map(function(p, ttl){
    # set panel title to symbol
    p <- p + ggplot2::labs(title = ttl)
    
    if (isTRUE(drop_legend)) p <- p + ggplot2::theme(legend.position = "none")
    
    # y-axis
    if (!is.null(y_label)) p <- p + if (identical(y_label, "")) ggplot2::ylab(NULL) else ggplot2::ylab(y_label)
    # x-axis
    if (!is.null(x_label)) p <- p + if (identical(x_label, "")) ggplot2::xlab(NULL) else ggplot2::xlab(x_label)
    
    # optional x tick rotation
    if (!is.null(rotate_x_ticks)) {
      hjust <- if (rotate_x_ticks == 0) 0.5 else 1
      vjust <- if (rotate_x_ticks == 0) 0.5 else 1
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = rotate_x_ticks, hjust = hjust, vjust = vjust))
    }
    
    p
  }, plist, disp)
  
  patchwork::wrap_plots(plist, ncol = ncol)
}
