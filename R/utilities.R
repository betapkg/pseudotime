


#' RunM3DropFeatureSelection
#'
#' @param seu seurat
#' @param root start root
#' @param outdir out dir
#'
#' @export
#'
RunM3DropFeatureSelection <- function(
    seu=NULL,
    assay='SCT',
    is_counts=FALSE,
    outdir='.'
){
    if (is_counts){
        data <- as(as.matrix(seu@assays[[assay]]@counts), 'sparseMatrix')
    } else {
        data <- as(as.matrix(seu@assays[[assay]]@data), 'sparseMatrix')
    }
    

    # Find DEG genes using 'M3DropFeatureSelection' (it's good enough)
    norm <- M3Drop::M3DropConvertData(data, is.counts=is_counts)
    if (is_counts){
        norm <- M3Drop::M3DropConvertData(log2(norm+1), is.log=TRUE, pseudocount=1)
    }
    d_diff_exp_m3 <- M3Drop::M3DropFeatureSelection(norm, mt_method="fdr", mt_threshold=0.001)
    write.table(d_diff_exp_m3, paste0(outdir,'/out_m3drop.featureSelection.xls'), sep='\t', quote=F, row.names=F)

    # only use top 2000 genes (best practise)
    # effect.size order will get similar results with marker gene based results, but q-val order is different
    diff_ordering_genes <- row.names(d_diff_exp_m3)
    diff_ordering_genes <- row.names(d_diff_exp_m3[1:2000,])

    writeLines(diff_ordering_genes, paste0(outdir,'/out_m3drop.DEG2k.gene'))

}




#' Find Root Cell
#'
#' @param seu seurat
#' @param gene gene
#' @param celltype cell type
#' @param sample sample name
#'
#' @export
#'
FindRootCell <- function(
    seu = NULL,
    gene = 'Pax6',
    celltype = 'RGC',
    sample = NULL
){
    d <- as.data.frame(seu@assays$SCT@data[gene,])
    colnames(d) <- gene

    d$cell_type2 <- seu$cell_type2[match(rownames(d), rownames(seu@meta.data))]
    d$sample <- seu$orig.ident[match(rownames(d), rownames(seu@meta.data))]
    d2 <- d[order(d[[gene]], decreasing=T),]

    if (length(sample) > 0){
        rc <- head(rownames(d2[d2$sample==sample & d2$cell_type2==celltype,]), n=1)
        } else {
            rc <- head(rownames(d2[d2$cell_type2==celltype,]), n=1)
        }

    return(rc)
}









