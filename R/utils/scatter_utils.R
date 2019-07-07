#' add label to colnames
#' @param add.colnames add colnames infront of barcode
#' @param genes.tsv read gene names from genes.tsv column. Default is 1
read10xCounts.1 <- function (samples, col.names = TRUE, type = c("auto", "sparse","HDF5"),
                             group = NULL, add.colnames = NULL, genes.tsv = 1) 
{
        nsets <- length(samples)
        full_data <- vector("list", nsets)
        gene_info_list <- vector("list", nsets)
        cell_info_list <- vector("list", nsets)
        type <- match.arg(type)
        for (i in seq_len(nsets)) {
                run <- samples[i]
                type <- .type_chooser(run, type)
                if (type == "sparse") {
                        info <- .read_from_sparse(run)
                }
                else {
                        info <- .read_from_hdf5(run, group = group)
                }
                full_data[[i]] <- info$mat
                gene_info_list[[i]] <- info$gene.info
                cell.names <- info$cell.names
                cell_info_list[[i]] <- DataFrame(Sample = rep(run, length(cell.names)), 
                                                 Barcode = cell.names, row.names = NULL)
        }
        if (nsets > 1 && length(unique(gene_info_list)) != 1L) {
                stop("gene information differs between runs")
        }
        gene_info <- gene_info_list[[1]]
        colnames(gene_info) <- c("ID", "Symbol")
        rownames(gene_info) <- gene_info[,genes.tsv]
        full_data <- do.call(cbind, full_data)
        cell_info <- do.call(rbind, cell_info_list)
        if (col.names && nsets == 1L) {
                if(is.null(add.colnames)) {
                        colnames(full_data) <- cell_info$Barcode
                }
                else {
                        colnames(full_data) <- paste0(add.colnames,
                                                      "_",cell_info$Barcode)
                }
        }
        SingleCellExperiment(list(counts = full_data), rowData = gene_info, 
                             colData = cell_info)
}

#' Extract delimiter information from a string.
#'
#' Parses a string (usually a cell name) and extracts fields based on a delimiter
#'
#' @param string String to parse.
#' @param field Integer(s) indicating which field(s) to extract. Can be a vector multiple numbers.
#' @param delim Delimiter to use, set to underscore by default.
#'
#' @return A new string, that parses out the requested fields, and (if multiple), rejoins them with the same delimiter
#'
#' @export
#'
#' @examples
#' ExtractField(string = 'Hello World', field = 1, delim = '_')
#'
ExtractField <- function(string, field = 1, delim = "_") {
        fields <- as.numeric(x = unlist(x = strsplit(x = as.character(x = field), split = ",")))
        if (length(x = fields) == 1) {
                return(strsplit(x = string, split = delim)[[1]][field])
        }
        return(paste(strsplit(x = string, split = delim)[[1]][fields], collapse = delim))
}

.read_from_sparse <- function (path) 
{
        barcode.loc <- file.path(path, "barcodes.tsv")
        gene.loc <- file.path(path, "genes.tsv")
        features.loc <- file.path(path, "features.tsv.gz")
        matrix.loc <- file.path(path, "matrix.mtx")
        pre_ver_3 <- file.exists(gene.loc)
        if (!pre_ver_3) {
                addgz <- function(s) {
                        return(paste0(s, ".gz"))
                }
                barcode.loc <- addgz(s = barcode.loc)
                matrix.loc <- addgz(s = matrix.loc)
        }
        mat <- readMM(matrix.loc)
        cell.names = readLines(barcode.loc)
        if (all(grepl(pattern = "\\-1$", x = cell.names))) {
                cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                                    FUN = ExtractField, field = 1, delim = "-")))
        }
        gene.info <- read.delim(file = ifelse(test = pre_ver_3, 
                                                  yes = gene.loc, no = features.loc), 
                                header = FALSE, colClasses = "character",stringsAsFactors = FALSE,
                                quote = "", comment.char = "")
        list(mat=mat, cell.names=cell.names, gene.info=gene.info)
}

.read_from_hdf5 <- function (path, group = NULL) 
{
        if (is.null(group)) {
                available <- h5ls(path, recursive = FALSE)
                available <- available[available$otype == "H5I_GROUP", 
                                       ]
                if (nrow(available) > 1L) {
                        to.see <- head(available$name, 3)
                        if (length(to.see) == 3L) {
                                to.see[3] <- "..."
                        }
                        stop("more than one available group (", paste(to.see, 
                                                                      collapse = ", "), ")")
                }
                else if (nrow(available) == 0L) {
                        stop("no available groups")
                }
                group <- available$name
        }
        mat <- TENxMatrix(path, group)
        dimnames(mat) <- NULL
        list(mat = mat, cell.names = as.character(h5read(path, paste0(group, 
                                                                      "/barcodes"))), gene.info = data.frame(as.character(h5read(path, 
                                                                                                                                 paste0(group, "/genes"))), as.character(h5read(path, 
                                                                                                                                                                                paste0(group, "/gene_names"))), stringsAsFactors = FALSE))
}

.type_chooser <- function(path, type) {
        if (type=="auto") {
                type <- if (grepl("\\.h5", path)) "HDF5" else "sparse"
        }
        type
}
