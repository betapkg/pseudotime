
#' RunPhate
#'
#' @param seu seurat
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

    # Calculate PHATE
    oupPhate <- phate(t(GetAssayData(seu)[VariableFeatures(seu), ]), knn=30, npca=50, seed=seed)
    oupDR = oupPhate$embedding
    oupDR = oupDR / 10^(floor(log10(diff(range(oupDR)))))
    rownames(oupDR) = colnames(seu)
    colnames(oupDR) = c("PHATE_1", "PHATE_2")

    seu[["phate"]] <- CreateDimReducObject(embeddings = oupDR, key = "PHATE_", assay = assay)

    #Idents(seu) <- 'cell_type2'
    #DimPlot(seu, reduction = "phate", pt.size = 0.1, label = TRUE)

    return(seu)
}











