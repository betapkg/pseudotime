

#' Run URD
#'
#' @param seu seurat
#' @param root start root
#' @param cluster cluster
#'
#' @return urd object
#' @export
#'
RunURD <- function(
    seu=NULL,
    cluster='cell_type2',
    root='RGC',
    seed=42
) {

    # 2.RNA count
    metadata <- seu@meta.data
    # Seurat 5
    #matrix_count <- seu[["RNA"]]$counts  # must use $RNA@counts rather than $SCT@counts
    matrix_count <- seu@assays$SCT@counts
    
    # 3.Create URD object
    axial <- createURD(count.data = matrix_count, meta = metadata, min.cells=3, min.counts=3)


    # 4.var genes
    var.genes <- seu@assays$SCT@var.features
    axial@var.genes <- var.genes

    # 5.Calculate PCA and consider those PCs that with standard deviation 2x expected by noise as significant
    axial <- calcPCA(axial, mp.factor = 2)

    # pca cutoff
    pcSDPlot(axial)

    # 6.Calculate tSNE
    set.seed(seed)
    axial <- calcTsne(object = axial)


    # 7. DM (dim map)
    # n_pcs depend on 'pcSDPlot(axial)' cutoff
    #axial <- calcDM(axial, knn = 100, sigma=16)
    axial <- calcDM(axial, knn=NULL, sigma.use=NULL)

    # diffusion map
    plotDimArray(axial, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map", label=cluster, plot.title="", legend=F)

    # tsne
    # plotDim(axial, "cell_type2", transitions.plot = 10000, plot.title="cell type2 (with transitions)")



    # 6.Calculate pseudotime
    # Here we use all cells from the first stage as the root
    axial@group.ids$cell_type2 <- as.character(axial@meta[rownames(axial@group.ids), cluster])
    root.cells <- cellsInCluster(axial, cluster, root)


    # create axial.floods
    axial.floods <- floodPseudotime(axial, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)
    # saveRDS(axial.floods, 'axial.floods.rds')

    # The we process the simulations into a pseudotime
    axial <- floodPseudotimeProcess(axial, axial.floods, floods.name="pseudotime")
    # saveRDS(axial, 'axial.rds')

    pseudotimePlotStabilityOverall(axial)

    # show pseudotime with tSNE
    plotDim(axial, "pseudotime")
    # show pseudotime with density
    plotDists(axial, "pseudotime", cluster, plot.title="Pseudotime by cell type")

    return(axial)

}



