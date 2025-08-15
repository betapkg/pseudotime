
#' RunPhate
#'
#' @param seurat_object seurat
#' @param use_py python path
#' @param assay seurat assay
#'
#' @return seurat object
#' @export
#'
RunPhate <- function(
    seu=NULL,
    use_py='/Users/hua/opt/miniconda3/envs/phate/bin/python',
    assay='SCT',
    seed=42
){
    library(reticulate)
    use_python(use_py)
    library(phateR)

    DefaultAssay(seurat_object) <- assay

    # Calculate PHATE
    oupPhate <- phate(t(GetAssayData(seurat_object)[VariableFeatures(seurat_object), ]), knn=30, npca=50, seed=seed)
    oupDR = oupPhate$embedding
    oupDR = oupDR / 10^(floor(log10(diff(range(oupDR)))))
    rownames(oupDR) = colnames(seurat_object)
    colnames(oupDR) = c("PHATE_1", "PHATE_2")

    seurat_object[["phate"]] <- CreateDimReducObject(embeddings = oupDR, key = "PHATE_", assay = assay)

    #Idents(seurat_object) <- 'cell_type2'
    #DimPlot(seurat_object, reduction = "phate", pt.size = 0.1, label = TRUE)

    return(seurat_object)
}











