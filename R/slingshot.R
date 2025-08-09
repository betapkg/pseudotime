
#' RunSlingshotSeurat
#'
#' @param seu seurat
#' @param assay assay
#' @param reduction reducedDim
#' @param root root name
#' @param clus cluster labels
#'
#' @export
#'
RunSlingshotSeurat <- function(
    seu = NULL,
    assay = 'SCT',
    reduction = 'RNA.UMAP',
    root = '5',
    clus = 'seurat_clusters',
    w=4.5,
    h=3,
    outdir = 'out_slingshot'
){
    dir.create(outdir)

    # seurat to sce
    sce <- as.SingleCellExperiment(seu, assay='SCT')
    
    # UMAP or RNA.UMAP
    # set root cluster 
    sce_slingshot <- slingshot(sce, reducedDim=reduction, start.clus=root, clusterLabels=sce[[clus]])

    # save slingshot object
    saveRDS(sce_slingshot, paste0(outdir, '/slingshot.rds'))

    # call lineage
    SlingshotDataSet(sce_slingshot)

    # save pesudotime
    pseudotimeED <- slingPseudotime(sce_slingshot, na=FALSE)
    write.table(pseudotimeED, paste0(outdir, '/slingshot.pseudotime.out'), sep='\t', quote=F, col.names=NA)


    # plot
    # set color
    library(RColorBrewer)
    colors <- colorRampPalette(rev(brewer.pal(11, 'Spectral'))[-6])(100)
    # library(grDevices)
    plotcol <- colors[cut(sce_slingshot$slingPseudotime_1, breaks=100)]

    pdf(paste0(outdir, '/sligshort.pseudotime_1.pdf'), width=w, height=h, useDingbats=FALSE)
    plot(reducedDims(sce_slingshot)[[reduction]], col=plotcol, pch=16, asp=0.5)
    lines(SlingshotDataSet(sce_slingshot), lwd=1, col='black')
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
    n = 200,
    lineage = 1
){
    pseudotimeED <- slingPseudotime(slingX, na=FALSE)
    cellWeightsED <- slingCurveWeights(slingX)

    library(tradeSeq)
    slingX <- fitGAM(counts = as.matrix(counts), pseudotime=pseudotimeED, cellWeights=cellWeightsED, nknots=5, verbose=T)
    ATres <- associationTest(slingX)
    topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:n]

    lineage_model <- paste0('slingPseudotime_', lineage)
    pst.ord <- order(slingX[[lineage_model]], na.last=NA)
    pt.matrix <- assays(slingX)$counts[topgenes, pst.ord]
    pt.matrix <- t(apply(pt.matrix, 1, function(x){smooth.spline(x,df=3)$y}))
    pt.matrix <- t(apply(pt.matrix, 1, function(x){(x-mean(x))/sd(x)}))

    pt.matrix
}





#' RunComplexHeatmap
#'
#' @param mtx pt-matrix
#' @param km km
#'
#' @return ht
#'
#' @export
#' 
RunComplexHeatmap <- function(mtx='pt.matrix', km=4) {
    htkm <- ComplexHeatmap::Heatmap(
          mtx,
          name                 = "z-score",
          col                  = circlize::colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
          show_row_names       = TRUE,
          show_column_names    = FALSE,
          row_names_gp         = gpar(fontsize = 6),
          km                   = km,
          row_title_rot        = 0,
          cluster_rows         = TRUE,
          cluster_row_slices   = FALSE,
          cluster_columns      = FALSE
        )

    htkm
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
    reduction = 'RNA.UMAP',
    root = '5',
    clus = 'seurat_clusters',
    outdir = 'out_slingshot',
    w=4.5,
    h=3,
    d_marker='',
    n = 200,
    lineage = 1,
    km=4
){
    RunSlingshotSeurat(seu=seu, assay=assay, reduction=reduction, root=root, clus=clus, w=w, h=h, outdir=outdir)

    library(dplyr)
    d_marker <- d_marker %>% filter(p_val_adj<0.001) %>% filter(avg_log2FC >1) %>% filter(pct.1 >0.5) %>% filter(diff_pct >0.3)
    diff.genes <- unique(d_marker$gene)

    counts <- as.matrix(seu@assays$SCT@counts[diff.genes,])

    sce_slingshot <- readRDS(paste0(outdir, '/slingshot.rds'))
    pt.matrix <- MakePTMatrix(slingX = sce_slingshot, counts = counts, n = n, lineage = lineage)

    ht <- RunComplexHeatmap(pt.matrix, km)

}











