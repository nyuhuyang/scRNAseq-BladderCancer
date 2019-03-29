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

.read_from_sparse <- function (path) 
{
        barcode.loc <- file.path(path, "barcodes.tsv")
        gene.loc <- file.path(path, "genes.tsv")
        matrix.loc <- file.path(path, "matrix.mtx")
        list(mat = as(readMM(matrix.loc), "dgCMatrix"), cell.names = readLines(barcode.loc), 
             gene.info = read.delim(gene.loc, header = FALSE, colClasses = "character", 
                                    stringsAsFactors = FALSE, quote = "", comment.char = ""))
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
