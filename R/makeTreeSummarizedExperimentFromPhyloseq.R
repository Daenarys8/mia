#' Coerce a \code{phyloseq} object to a \code{TreeSummarizedExperiment}
#'
#' \code{makeTreeSEFromPhyloseq} converts \code{phyloseq}
#' objects into \code{TreeSummarizedExperiment} objects.
#'
#' All data stored in a \code{phyloseq} object is transferred.
#'
#' @param x a \code{phyloseq} object
#'
#' @return An object of class \code{TreeSummarizedExperiment}
#'
#' @importFrom S4Vectors SimpleList DataFrame make_zero_col_DFrame
#' @importFrom SummarizedExperiment colData colData<-
#'
#' @export
#'
#' @name makeTreeSEFromPhyloseq
#' @seealso
#' \code{\link[=makeTreeSEFromBiom]{makeTreeSEFromBiom}}
#' \code{\link[=makeTreeSEFromDADA2]{makeTreeSEFromDADA2}}
#' \code{\link[=importQIIME2]{importQIIME2}}
#' \code{\link[=importMothur]{importMothur}}
#'
#' @examples
#' if (requireNamespace("phyloseq")) {
#'     data(GlobalPatterns, package="phyloseq")
#'     makeTreeSEFromPhyloseq(GlobalPatterns)
#'     data(enterotype, package="phyloseq")
#'     makeTreeSEFromPhyloseq(enterotype)
#'     data(esophagus, package="phyloseq")
#'     makeTreeSEFromPhyloseq(esophagus)
#' }
makeTreeSEFromPhyloseq <- function(x) {
    # input check
    .require_package("phyloseq")
    if(!is(x,"phyloseq")){
        stop("'x' must be a 'phyloseq' object")
    }
    #
    # Get the assay
    counts <- x@otu_table@.Data
    # Check the orientation, and transpose if necessary
    if( !x@otu_table@taxa_are_rows ){
        counts <- t(counts)
    }
    # Create a list of assays
    assays <- SimpleList(counts = counts)
    
    if(!is.null(x@tax_table@.Data)){
        rowData <- DataFrame(data.frame(x@tax_table@.Data))
    } else{
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(assays$counts))
        rownames(rowData) <- rownames(assays$counts)
    }
    if(!is.null(x@sam_data)){
        colData <- DataFrame(data.frame(x@sam_data))
    } else{
        colData <- S4Vectors::make_zero_col_DFrame(ncol(assays$counts))
        rownames(colData) <- colnames(assays$counts)
    }
    if(!is.null(x@phy_tree)){
        rowTree <- x@phy_tree
    } else {
        rowTree <- NULL
    }
    if (!is.null(x@refseq)) {
        referenceSeq <- x@refseq
    } else {
        referenceSeq <- NULL
    }
    TreeSummarizedExperiment(assays = assays,
                            rowData = rowData,
                            colData = colData,
                            rowTree = rowTree,
                            referenceSeq = referenceSeq)
}

####################### makeTreeSummarizedExperimentFromPhyloseq #######################
#' @rdname makeTreeSEFromPhyloseq
#' @export
setGeneric("makeTreeSummarizedExperimentFromPhyloseq", signature = c("x"),
    function(x)
        standardGeneric("makeTreeSummarizedExperimentFromPhyloseq"))

#' @rdname makeTreeSEFromPhyloseq
#' @export
setMethod("makeTreeSummarizedExperimentFromPhyloseq", signature = c(x = "ANY"),
    function(x){
        makeTreeSEFromPhyloseq(x)
    })
