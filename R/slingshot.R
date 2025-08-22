

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
    plot(reducedDims(sce_slingshot)[[reduction]], col=plotcol, pch=16, asp=0.5)
    lines(SlingshotDataSet(sce_slingshot), lwd=1, col='black')
    dev.off()
}






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
    reduction = 'PHATE',
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
    sce_slingshot <- slingshot(data=sce, reducedDim=reduction, start.clus=start_root, clusterLabels=ident, allow.breaks = allow_breaks)

    summary_lineages <- SlingshotDataSet(sce_slingshot)
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
        pseudotimeED <- slingPseudotime(sce_slingshot, na=FALSE)
        write.table(pseudotimeED, paste0(outdir, '/2.slingshot.pseudotime.xls'), sep='\t', quote=F, col.names=NA)

        # plot
        PlotSlingshot(sce_slingshot=sce_slingshot, reduction=reduction, lineage=1, w=4.5, h=3, outdir=outdir)
    } else{
        return(sce_slingshot)
    }
    
}













