

#################################
##         Calculation          
#################################

#' Seurat_FindMarkers
#'
#' @param obj seurat
#' @param groupby cell_type2
#' @param outdir out dir
#'
#' @export
#'
Seurat_FindMarkers <- function(
    obj=NULL, 
    groupby='cell_type2',
    ident1='', 
    min_pct=0.3,
    logfc=1,
    recorrect_umi=TRUE)
{
    markers <- Seurat::FindMarkers(obj, ident.1 = ident1, group.by = groupby, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc, recorrect_umi=recorrect_umi)
    
    return(markers)
}




#' ConverGeneToMouse
#'
#' @param gene gene name
#'
#' @return gene
#' @export
#'
ConverGeneToMouse <- function(gene){
    mouse_gene <- gprofiler2::gorth(gene, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    return(mouse_gene)
}



#' SeuratAddScore
#'
#' @param obj seurat object 
#' @param name score name  
#' @param features gene set  
#'
#' @return seurat obj
#' @export
#'
SeuratAddScore <- function(obj=NULL, name='', features=''){
    all_gene <- rownames(obj[['SCT']]@data)
    gene_intersect <- intersect(features, all_gene)

    obj <- AddModuleScore(
        object = obj,
        assay = 'SCT',
        slot = "data",
        features = list(gene_intersect),
        name = name
    )

    return(obj)
}



#' Seurat_EstimateCycling
#'
#' @param obj seurat object 
#' @param ref human/mouse   
#'
#' @return data frame
#' @export
#'
Seurat_EstimateCycling <- function(obj=NULL, ref='human')
{
    # cell cycle
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes

    if (ref == 'mouse'){
        s.genes <- ConverGeneToMouse(cc.genes.updated.2019$s.genes)
        g2m.genes <- ConverGeneToMouse(cc.genes.updated.2019$g2m.genes)
    }
    
    # CellCycleScoring
    obj <- CellCycleScoring(
        object = obj,
        s.features = s.genes,
        g2m.features = g2m.genes
    )

    # CellCycleScoring2 for all cell cycle genes
    cc.genes <- c(s.genes, g2m.genes)
    obj <- SeuratAddScore(obj=obj, name='cc_score', features=cc.genes)

    df <- obj@meta.data[,c('cc_score1', 'S.Score', 'G2M.Score', 'Phase', 'seurat_clusters', 'cell_type2')]
    df$group <- 'Non-cycling'
    
    # define cycling
    df$group[df$cc_score1 > 0.2 & (df$S.Score > 0.2 | df$G2M.Score > 0.2) ] <- 'Cycling'

    return(df)
}



#################################
##           Plot          
#################################

#' Score_UMAPPlot
#'
#' @param obj seurat object 
#' @param umap umap 
#' @param name score name 
#'
#' @return plot
#' @export
#'
Score_UMAPPlot <- function(obj=NULL, umap='umap', name='')
{
    library(ggplot2)
    #library(RColorBrewer)

    obj$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    obj$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]

    df <- obj@meta.data[,c('UMAP_1', 'UMAP_2', paste0(name,'1'))]
    colnames(df) <- c('UMAP_1', 'UMAP_2', 'signal')

    df <- df[order(df$signal, decreasing=FALSE), ]

    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=signal)) + 
            geom_point(size=0.2) +
            theme_classic(base_line_size=0.2) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            theme(legend.key.width=unit(5,"mm")) +
            theme(legend.title = element_text(size=9)) +
            labs(x = "UMAP 1", y="UMAP 2")

    
    #p <- p + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
    color <- c("#2166AC", "#F5F9CD", "#B2182B")
    p <- p + scale_colour_gradientn(colours = color)


    return(p)
}




#' Score_UMAPPlot
#'
#' @param data data frame
#' @param group group 
#' @param color new color
#'
#' @return plot
#' @export
#'
ScatterPlot2D <- function(
    data=NULL, 
    x='differentiation_score',
    y='progenitor_score',
    group='group', 
    x_lab='Differentiation Score',
    y_lab='Progenitor Score',
    color=NULL)
{
    # lineage_score, progenitor_score, ...
    data <- data[order(data$group, decreasing = TRUE),]

    #colors <- c("Cycling"="#E64130", "Non-cycling"="#B3B3B3")

    x_min <- min(data[[x]])
    x_max <- max(data[[x]])
    if (abs(x_min) > abs(x_max)){
        x_max <- abs(x_min)
    }

    y_min <- min(data[[y]])
    y_max <- max(data[[y]])

    p <- ggplot(data, aes_string(x=x, y=y, color=group)) + 
            geom_point(shape=16, size=0.4, alpha=0.8) +  
            theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(text=element_text(size=8), plot.title=element_text(size=9, hjust=0.5)) +
            labs(title='', x=x_lab, y=y_lab) +
            xlim(x_min, x_max) +
            theme(legend.position="top", legend.title = element_blank(),
                legend.text = element_text(), legend.spacing.x = unit(0, 'mm')) +
            guides(color = guide_legend(override.aes = list(size = 2)))

    if (length(color)>0){
        p <- p + scale_color_manual(values=color)
    }

    p <- p + geom_segment(aes(x=x_min, y=y_min/3.5, xend=x_min, yend=y_max/5), arrow=arrow(length=unit(1, 'mm'), type='closed', ends = "both"), color='black', size = 0.2) +
        annotate('text', x=x_min, y=y_max/1.3, label='Undifferentiated', size=2.5, angle='90') +
        annotate('text', x=x_min, y=y_min/1.4, label='Differentiated', size=2.5, angle='90')

    return(p)
}





#################################
##          Main pipe          
#################################

#' RunTriangleAnalysis
#'
#' @param seu Seurat object
#' @param ref human/mouse  
#' @param root root cell type 
#' @param armL left arm cell type
#' @param armR  right arm cell type
#' @param min_pct_shared pct cutoff for shared
#' @param logfc_shared logfc cutoff for shared 
#' @param min_pct_armL pct cutoff for armL
#' @param logfc_armL logfc cutoff for armL
#' @param min_pct_armR pct cutoff for armR
#' @param logfc_armR logfc cutoff for armR 
#' @param top_n_genes cutoff top genes 
#' @param expr_ratio cutoff ratio of expressed cells 
#' @param out_dir out dir.
#'
#' @return data frame
#' @export
#'
RunTriangleAnalysis <- function(
    seu=NULL, ref = 'human', out_dir = 'out_triangelplot', 
    root='CycProg', armL='Astrocyte', armR='Neuron',
    min_pct_shared = 0.8, logfc_shared = 0.05,
    min_pct_armL = 0.5, logfc_armL = 1,
    min_pct_armR = 0.5, logfc_armR = 0.5,
    top_n_genes = 50, expr_ratio = 0.8
) {
  
  # Load required libraries
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  
  # Create output directory
  dir.create(out_dir, showWarnings = FALSE)
  
  # Subset to target cell types
  seu <- subset(x = seu, subset = cell_type2 %in% c(root, armL, armR))
  
  # Estimate cycling cells
  d_cc <- Seurat_EstimateCycling(obj = seu, ref = ref)
  
  # Get expression matrix
  sct_mtx <- as.data.frame(as.matrix(seu[['SCT']]@data))
  d_data <- as.matrix(seu[['SCT']]@data)
  all_genes <- rownames(d_data)
  
  # Filter mitochondrial or ribosomal genes
  if (ref == 'human') {
    filtered_genes <- all_genes[!grepl("^RPS|^RPL|^MT-", all_genes)]
  } else {
    filtered_genes <- all_genes[!grepl("Rps|Rpl|mt-", all_genes)]
  }
  
  # Find expressed genes (in >80% cells)
  exp_genes <- rownames(d_data[rowSums(d_data > 0) > ncol(d_data) * expr_ratio,])
  exp_genes <- intersect(exp_genes, filtered_genes)
  cat("\nNumber of expressed genes:", length(exp_genes), "\n")
  
  #--- 1.root cell type and shared markers
  # Find shared markers for root cells
  markers_root <- Seurat_FindMarkers(obj = seu, ident1 = root, 
                                    min_pct = min_pct_shared, logfc = logfc_shared, 
                                    recorrect_umi = FALSE)
  markers_root <- filter(markers_root, avg_log2FC < 0.2)
  markers_root <- markers_root[order(markers_root$avg_log2FC),]
  
  if (nrow(markers_root) > 30) {
    shared_genes <- intersect(rownames(markers_root), exp_genes)
  } else {
    shared_genes <- rownames(markers_root)
  }
  #cat("Number of shared genes:", length(shared_genes), "\n")
  
  # Get top shared genes
  if (length(shared_genes) > top_n_genes) { 
        shared_genes <- shared_genes[1:top_n_genes] 
  }
  
  # Calculate average expression for shared genes
  df_shared_mean <- as.data.frame(sapply(sct_mtx[shared_genes, ], FUN = mean))
  colnames(df_shared_mean) <- 'avgExp_shared'
  

  #--- 2. non-root cell type and specific markers
  # Find markers for armL-Astro cells
  markers_armL <- Seurat_FindMarkers(obj = seu, ident1 = armL, 
                                     min_pct = min_pct_armL, logfc = logfc_armL, 
                                     recorrect_umi = FALSE)
  markers_armL <- markers_armL[order(-markers_armL$avg_log2FC),]
  gene_armL <- intersect(rownames(markers_armL), filtered_genes)
  if (length(gene_armL) > top_n_genes) { 
    gene_armL <- gene_armL[1:top_n_genes] 
  }
  
  # Calculate average expression for cell type markers
  df_mean_armL <- as.data.frame(sapply(sct_mtx[gene_armL, ], FUN = mean))
  colnames(df_mean_armL) <- 'avgExp_armL'
  

  # Find markers for armR-Neuron cells
  markers_armR <- Seurat_FindMarkers(obj = seu, ident1 = armR, 
                                     min_pct = min_pct_armR, logfc = logfc_armR, 
                                     recorrect_umi = FALSE)
  markers_armR <- markers_armR[order(-markers_armR$avg_log2FC),]
  gene_armR <- intersect(rownames(markers_armR), filtered_genes)
  if (length(gene_armR) > top_n_genes) { 
    gene_armR <- gene_armR[1:top_n_genes] 
  }
  #cat("Number of Neuron markers:", length(gene_armR), "\n")
  df_mean_armR <- as.data.frame(sapply(sct_mtx[gene_armR, ], FUN = mean))
  colnames(df_mean_armR) <- 'avgExp_armR'
  
  # Save gene lists
  writeLines(shared_genes, paste0(out_dir, '/shared.gene'))
  writeLines(gene_armL, paste0(out_dir, '/', armL, '.gene'))
  writeLines(gene_armR, paste0(out_dir, '/', armR, '.gene'))
  

  # Merge all scores
  df_merged <- df_shared_mean
  df_merged <- cbind(df_merged, df_mean_armL)
  df_merged <- cbind(df_merged, df_mean_armR)
  df_merged <- merge(df_merged, d_cc, by = 'row.names')
  
  # Calculate progenitor and differentiation scores
  df_merged$progenitor_score <- df_merged$avgExp_shared - pmax(df_merged$avgExp_armR, df_merged$avgExp_armL)
  df_merged$differentiation_score <- df_merged$avgExp_armR - df_merged$avgExp_armL
  
  # Save results
  write.table(df_merged, paste0(out_dir, '/progenitor_score.out'), sep = '\t', quote = FALSE, row.names = FALSE)
  
  cat("\nAnalysis completed. Results saved to:", out_dir, "\n")
  
  return(df_merged)
}









