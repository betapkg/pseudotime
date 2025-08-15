
#' RunSlingshotSeurat
#'
#' @param seu seurat
#' @param assay assay
#' @param reduction reducedDim
#' @param start_root root name
#' @param ident ident labels
#'
#' @export
#'
RunSlingshotSeurat <- function(
    seu = NULL,
    assay = 'SCT',
    ident = 'seurat_clusters',
    reduction = 'PCA',
    start_root = NULL,
    allow_breaks = TRUE,
    save_out = FALSE,
    outdir = 'out_slingshot'
){
    dir.create(outdir)

    # seurat to sce
    sce <- as.SingleCellExperiment(seu, assay=assay)
    
    # UMAP or RNA.UMAP
    # set root cluster 
    sce_slingshot <- slingshot(sce, reducedDim=reduction, start.clus=start_root, clusterLabels=ident, allow.breaks = allow_breaks)

    summary_lineages <- slingshot::SlingshotDataSet(sce_slingshot)
    print(summary_lineages@lineages)


    if (save_out){
        print('save ...')
        # save summary lineage
        library(conflicted)
        library(tidyverse)
        s_list <- summary_lineages@lineages
        #d_lineages <- as.data.frame(summary_lineages@lineages)
        d_lineages <- map(s_list, ~ .x %>% `length<-`(max(lengths(s_list)))) %>% bind_rows() %>% as.data.frame()
        write.table(d_lineages, file=paste0(outdir, '/1.slingshot.summary_lineages.xls'), sep='\t', quote=F, col.names=NA)

        # save slingshot object
        saveRDS(sce_slingshot, paste0(outdir, '/1.slingshot.rds'))

        # save pesudotime
        pseudotimeED <- slingshot::slingPseudotime(sce_slingshot, na=FALSE)
        write.table(pseudotimeED, paste0(outdir, '/2.slingshot.pseudotime.xls'), sep='\t', quote=F, col.names=NA)

        # plot
        PlotSlingshot(sce_slingshot=sce_slingshot, reduction=reduction, lineage=1, w=4.5, h=3, outdir=outdir)
    } else{
        return(sce_slingshot)
    }
    
}






#' Plot-Slingshot
#'
#' @param sce_slingshot sce_slingshot
#' @param lineage lineage n
#'
#' @import RColorBrewer 
#' @export
#'
PlotSlingshot <- function(
    sce_slingshot = NULL,
    reduction = NULL,
    lineage = 1,
    w=4.5,
    h=3,
    outdir = 'out_slingshot'
){  
    library(RColorBrewer)
    # plot
    colors <- colorRampPalette(rev(brewer.pal(11, 'Spectral'))[-6])(100)
    # library(grDevices)
    lineage_name <- paste0('slingPseudotime_', lineage)
    plotcol <- colors[cut(sce_slingshot[[lineage_name]], breaks=100)]

    pdf(paste0(outdir, '/', lineage_name, '.pdf'), width=w, height=h, useDingbats=FALSE)
    plot(slingshot::reducedDims(sce_slingshot)[[reduction]], col=plotcol, pch=16, asp=0.5)
    lines(slingshot::SlingshotDataSet(sce_slingshot), lwd=1, col='black')
    dev.off()
}






#' MakePTMatrix
#'
#' @param slingX sling object
#' @param counts count matrix
#' @param n total genes
#' @param lineage lineage model
#'
#' @return pt-matrix
#'
#' @export
#' 
MakePTMatrix <- function(
    slingX = 'sce_slingshot',
    counts = NULL,
    nknots = 10,
    n = 200,
    lineage = 1
){
    pseudotimeED <- slingPseudotime(slingX, na=FALSE)
    cellWeightsED <- slingCurveWeights(slingX)

    # define top genes
    library(tradeSeq)
    X_slingshot <- fitGAM(counts = as.matrix(counts), pseudotime=pseudotimeED, cellWeights=cellWeightsED, nknots=nknots, verbose=T)
    ATres <- associationTest(X_slingshot)
    topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:n]

    # make matrix
    lineage_model <- paste0('slingPseudotime_', lineage)
    pst.ord <- order(slingX[[lineage_model]], na.last=NA)
    pt.matrix <- assays(slingX)$counts[topgenes, pst.ord]
    pt.matrix <- t(apply(pt.matrix, 1, function(x){smooth.spline(x, df=3)$y}))
    pt.matrix <- t(apply(pt.matrix, 1, function(x){(x-mean(x))/sd(x)}))

    return(pt.matrix)
}





#' RunComplexHeatmap
#'
#' @param pt_mtx pt-matrix
#' @param km km
#'
#' @return ht
#'
#' @export
#' 
RunComplexHeatmap <- function(pt_mtx='pt.matrix', km=4, outdir='.') {
    htkm <- ComplexHeatmap::Heatmap(
          pt_mtx,
          name                 = "z-score",
          col                  = circlize::colorRamp2(seq(from=-2,to=2,length=11),rev(RColorBrewer::brewer.pal(11, "Spectral"))),
          show_row_names       = TRUE,
          show_column_names    = FALSE,
          row_names_gp         = grid::gpar(fontsize = 6),
          km                   = km,
          row_title_rot        = 0,
          cluster_rows         = TRUE,
          cluster_row_slices   = FALSE,
          cluster_columns      = FALSE
        )


    HM <- ComplexHeatmap::draw(htkm)  #Show the heatmap

    r.dend <- ComplexHeatmap::row_dend(HM)  #If needed, extract row dendrogram
    rcl.list <- ComplexHeatmap::row_order(HM)  #Extract clusters (output is a list)
      
    lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

    library(magrittr) # needed to load the pipe function '%%'
     
    clu_df <- lapply(names(rcl.list), function(i){
                out <- data.frame(GeneID = rownames(pt_mtx[rcl.list[[i]],]),
                         Cluster = paste0("cluster", i),
                         stringsAsFactors = FALSE)
                return(out)
            }) %>% do.call(rbind, .)

    #export
    write.table(clu_df, file= paste0(outdir, '/3.gene_clusters.kmean_heatmap.txt'), sep="\t", quote=F, row.names=FALSE)

    # save plot
    pdf(paste0(outdir, '/heatmap.pseudotime.pdf'), w=3.5, h=4, useDingbats=FALSE)
    HM
    dev.off()
}





#' RunSlingshotPipe_Seurat
#'
#' @param seu pt-matrix
#' @param d_marker data frame for markers
#'
#' @export
#' 
RunSlingshotPipe_Seurat <- function(
    seu = NULL,
    assay = 'SCT',
    reduction = 'PCA',
    start_root = NULL,
    allow_breaks=TRUE,
    ident = 'seurat_clusters',
    outdir = 'out_slingshot',
    w=4.5,
    h=3,
    d_marker='',
    avg_log2FC=0.585,
    pct_1=0.3,
    diff_pct=0.2,
    n = 200,
    nknots = 10,
    lineage = 1,
    km=4
){
    print('1.pseudotime')
    RunSlingshotSeurat(seu=seu, assay=assay, reduction=reduction, start_root=root, ident=ident, allow_breaks=allow_breaks, outdir=outdir, save_out=TRUE)  
    sce_slingshot <- readRDS(paste0(outdir, '/1.slingshot.rds'))

    print('2-1.diff markers')
    library(dplyr)
    d_marker <- d_marker %>% filter(p_val_adj<0.001) %>% filter(avg_log2FC > avg_log2FC) %>% filter(pct.1 > pct_1) %>% filter(diff_pct > diff_pct)
    diff.genes <- unique(d_marker$gene)


    counts <- as.matrix(seu@assays$SCT@counts[diff.genes,])

    print('2-2.matrix')
    pt.matrix <- MakePTMatrix(slingX=sce_slingshot, counts=counts, n=n, nknots=nknots, lineage=lineage)
    # save pt.matrix
    saveRDS(pt.matrix, paste0(outdir, '/2.matrix.lineage_', lineage, '.rds'))

    print('3.heatmap')
    RunComplexHeatmap(pt.matrix, km, outdir)
}











