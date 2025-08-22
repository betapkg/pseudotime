

#####################################
#       Pseudotime - Heatmap
#####################################

#' MakePTMatrix - Slingshot
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
MakePTMatrixSlingshot <- function(
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
    write.table(clu_df, file= paste0(outdir, '/pseudotime_heatmap.gene_clusters.txt'), sep="\t", quote=F, row.names=FALSE)

    # save plot
    pdf(paste0(outdir, '/pseudotime_heatmap.pdf'), w=3.5, h=4, useDingbats=FALSE)
    HM
    dev.off()
}









