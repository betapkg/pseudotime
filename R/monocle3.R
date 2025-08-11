

#' Monocle3MakeCDS
#'
#' @param obj seurat object
#' @param assay start root
#' @param umap out dir
#'
#' @export
#'
Monocle3MakeCDS <- function(
    obj=NULL, 
    assay='SCT', 
    umap='wnn.umap'
){  
    Seurat::DefaultAssay(obj) <- assay

    obj[["UMAP"]] <- obj[[umap]]

    # make cds
    cds <- SeuratWrappers::as.cell_data_set(obj)
    cds <- monocle3::cluster_cells(cds = cds, reduction_method = "UMAP")
    cds <- monocle3::learn_graph(cds, use_partition = TRUE)
    
    return(cds)
}




#' RunMonocle3Pipe
#'
#' @param seu seurat object
#' @param root_cells start root
#' @param outdir out dir
#'
#' @export
#'
RunMonocle3Pipe <- function(
    seu=NULL,
    root_cells=NULL,
    umap='umap',
    outdir='.'
){
    library(monocle3)
    library(ggplot2)

    dir.create(outdir)

    # 1.make cds
    cds <- Monocle3MakeCDS( obj=seu, assay='SCT', umap=umap)
    saveRDS(cds, paste0(outdir, '/1.monocle3_cds.rds'))

    # 2.order cells
    #meta_filtered <- meta %>% filter(orig.ident=='E12' & cell_type2=='RGC')
    #root_cells <- rownames(meta_filtered)
    cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

    # 3. export pseudotime
    pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
    d_pseudotime <- data.frame(pseudotime)
    write.table(d_pseudotime, paste0(outdir, '/2.pseudotime.out'), sep='\t', quote=F, col.names=NA)

    # 4. plot cells
    p <- plot_cells(
      cds = cds,
      color_cells_by = "pseudotime",
      show_trajectory_graph = TRUE,
      label_leaves=FALSE,
      label_branch_points=FALSE
    )
    ggsave(paste0(outdir, '/2a.pseudotime.umap.pdf'), w=6.5, h=5)
    

    # 5. make matrix for heatmap
    modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
    # In Monocle3, "Morans_I" refers to a statistical measure called Moran's I, which is used to assess the spatial autocorrelation of gene expression along a cell trajectory, essentially indicating how similar the expression levels of a gene are between neighboring cells on the trajectory; a high Moran's I value suggests that nearby cells tend to have similar expression levels for that gene
    genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.5))

    # used pseudotime
    pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
    #Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
    #normalized_counts(cds, norm_method = "log")

    pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
    pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
    rownames(pt.matrix) <- genes

    write.table(pt.matrix, paste0(outdir,'/3.pt_matrix.tsv'), sep='\t', quote=F, row.names=F)

}






