
#' RunPhate
#'
#' @param seu_obj seurat
#' @param use_py python path
#' @param assay seurat assay
#'
#' @return seurat object
#' @export
#'
RunPhate <- function(
    seu_obj=NULL,
    use_py='/Users/hua/opt/miniconda3/envs/phate/bin/python',
    assay='SCT',
    seed=42
){
    library(reticulate)
    use_python(use_py)
    library(phateR)

    DefaultAssay(seu_obj) <- assay


    mat <- as.matrix(GetAssayData(seu_obj, slot = "data")[VariableFeatures(seu_obj), ])

    # Calculate PHATE
    oupPhate <- phate(t(mat), knn=30, npca=50, seed=seed)
    oupDR = oupPhate$embedding
    oupDR = oupDR / 10^(floor(log10(diff(range(oupDR)))))
    rownames(oupDR) = colnames(seu_obj)
    colnames(oupDR) = c("PHATE_1", "PHATE_2")

    seu_obj[["phate"]] <- CreateDimReducObject(embeddings = oupDR, key = "PHATE_", assay = assay)

    #Idents(seu_obj) <- 'cell_type2'
    #DimPlot(seu_obj, reduction = "phate", pt.size = 0.1, label = TRUE)

    return(seu_obj)
}











