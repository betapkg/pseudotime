

#' RunMonocle2
#'
#' @param seu seurat
#' @param root start root
#' @param outdir out dir
#'
#' @export
#'
RunMonocle2 <- function(
    seu=NULL,
    root=NULL,
    outdir='.'
){
    library(monocle)
    
    data <- as(as.matrix(seu@assays$SCT@data), 'sparseMatrix')
    pd <- new('AnnotatedDataFrame', data = seu@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)

    HSMM <- monocle::newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = uninormal())

    HSMM$cells <- colnames(data)  # same as rownames(HSMM@phenoData)

    # add Size_Factor
    HSMM <- estimateSizeFactors(HSMM) 
    #HSMM <- estimateDispersions(HSMM)


    # Find DEG genes using 'M3DropFeatureSelection' (it's good enough)
    norm <- M3Drop::M3DropConvertData(data, is.counts=FALSE)
    d_diff_exp_m3 <- M3Drop::M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=0.001)
    write.table(d_diff_exp_m3, paste0(outdir,'/out_m3drop.featureSelection.xls'), sep='\t', quote=F, row.names=F)

    # only use top 2000 genes (best practise)
    # effect.size order will get similar results with marker gene based results, but q-val order is different
    diff_ordering_genes <- row.names(d_diff_exp_m3)
    diff_ordering_genes <- row.names(d_diff_exp_m3[1:2000,])


    #Step 1: choosing genes that define progress
    #Step 2: reducing the dimensionality of the data
    #Step 3: ordering the cells in pseudotime

    # Don't use VariableFeatures(seu) as diff_ordering_genes. It will get messy results.

    HSMM <- monocle::setOrderingFilter(HSMM, ordering_genes=diff_ordering_genes)
    # https://rdrr.io/bioc/monocle/man/reduceDimension.html
    # If you don't want any transformation at all, set norm_method to "none" and pseudo_expr to 0. This maybe useful for single-cell qPCR data, or data you've already transformed yourself in some way.

    ## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
    # max_components=2 can get the best result
    HSMM <- monocle::reduceDimension(HSMM, norm_method="none", 
                            reduction_method="DDRTree",
                            max_components=2,  
                            scaling=TRUE,
                            verbose=TRUE,
                            pseudo_expr=0)

    HSMM <- orderCells(HSMM)
    saveRDS(HSMM, file=paste0(outdir, '/monocle2_ordercell.rds'))


    #------- set 'root' cell type, according to the above results 'plot_cell_trajectory'

    ## ordering cells by assigning root nodes
    GM_state <- function(cds=NULL, root_celltype=NULL){
        if (length(unique(cds$State)) > 1){
            T0_counts <- table(cds$State, cds$cell_type2)[, root_celltype]
            return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
        } else {
            return (1)
        }
    }

    # set root cell
    #start <- 'CycProg-Like'
    if (length(root) > 0){
        HSMM <- orderCells(HSMM, root_state=GM_state(HSMM, root))
        saveRDS(HSMM, file=paste0(outdir, '/monocle2_ordercell.root.rds'))
    }
    

    d_cds <- data.frame(cells=HSMM$cells, pseudotime=HSMM$Pseudotime, cell_type2=HSMM$cell_type2)
    write.table(d_cds, paste0(outdir, '/monocle2_ordercell.root.pseudotime.xls'), sep='\t', quote=F, row.names=F)
    
    
    library(igraph)
    library(ggplot2)
    
    # save final plot_cell_trajectory
    p <- monocle::plot_cell_trajectory(HSMM, 
                         color_by = "cell_type2",  # seurat_clusters
                         theta = -15,
                         show_branch_points = TRUE,
                         show_tree = TRUE,
                         cell_link_size = 0.2,
                         cell_size = 0.1)  +
        theme(legend.title=element_blank(), legend.spacing.x = unit(0, 'cm'))+
        guides(color = guide_legend(override.aes = list(size = 2)))

    ggsave(paste0(outdir, '/plot_cell_trajectory.pdf'), width=3, height=3, useDingbats=FALSE)



    # save final plot_cell_trajectory with Pseudotime
    p <- monocle::plot_cell_trajectory(HSMM, color_by = "Pseudotime", cell_link_size = 0.2, cell_size = 0.1) +
        scale_color_viridis_c() + 
        theme(legend.key.height=unit(4, 'mm'))

    ggsave(paste0(outdir, '/plot_cell_trajectory.withPseudotime.pdf'), width=3, height=3, useDingbats=FALSE)

}






