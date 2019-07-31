library(dplyr)
library(magrittr)
library(Seurat)

#' .AddMetaColor: prepare meta.data data frame to store color code
#' @param mat factor at column 1, rownames is cell.names
#' @param colors vector of hex colors
#' @param df_colors data frame to add to meta data directly
#' @example 
# singlerDF = .AddMetaColor(mat = singler$singler[[2]]$SingleR.single$labels,
#                        colors = singler_colors[3:26])
.AddMetaColor <- function(mat, colors){
    if(class(mat) != "data.frame") mat = as.data.frame(mat)
    mat$cell.names = rownames(mat)
    mat$index <- as.numeric(as.factor(mat[,1]))
    if(length(colors)<length(unique(mat$index))) {
        stop(paste("Not enough colors! Provide at least", 
                   length(unique(mat$index)),"different colors"))}
    df_colors = data.frame(colors[1:length(unique(mat$index))],
                           "index" = 1:length(unique(mat$index)))
    mat_colors <- dplyr::full_join(mat, df_colors, by = "index")
    # remove NA color if there is any
    mat_colors <- mat_colors[(mat_colors$cell.names %in% rownames(mat)),]
    rownames(mat_colors) = mat_colors$cell.names
    df_colors <- data.frame(mat_colors$colors,
                            row.names = mat_colors$cell.names)
    colnames(df_colors)[1] = paste(colnames(mat)[1],"colors",sep =".")
    
    return(df_colors)
}


#' AddMetaColor: convert one MetaData label into color scheme and store into MetaData
#' @object seurat object
#' @label colname in metadata
#' @colors vector of hex colors
# MCL <- AddMetaColor(object = MCL, label= "singler1sub", colors = singler_colors)
AddMetaColor<- function(object, label = NULL, colors = NULL){
    
    if(is.null(label)) label <- FindIdentLabel(object)
    if(is.null(colors)) colors <- SingleR:::singler.colors
    mat = data.frame(object@meta.data[,label],
                     row.names = object@cell.names)
    colnames(mat) = get("label")
    newMetaData = .AddMetaColor(mat = mat, colors = colors)
    object <- AddMetaData(object = object, metadata = newMetaData)
    
    return(object)
}


# keep enrich.name as colname
# add option for only.pos = FALSE
.AddModuleScore <- function (object, genes.list = NULL, genes.pool = NULL, n.bin = 25, 
          seed.use = 1, ctrl.size = 100, use.k = FALSE, enrich.name = "Cluster", 
          random.seed = 1, only.pos = F) 
{
    set.seed(seed = random.seed)
    genes.old <- genes.list
    if (use.k) {
        genes.list <- list()
        for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
            genes.list[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                                                   i))
        }
        cluster.length <- length(x = genes.list)
    }
    else {
        if (is.null(x = genes.list)) {
            stop("Missing input gene list")
        }
        genes.list <- lapply(X = genes.list, FUN = function(x) {
            return(intersect(x = x, y = rownames(x = object@data)))
        })
        cluster.length <- length(x = genes.list)
    }
    if (!all(Seurat:::LengthCheck(values = genes.list))) {
        warning(paste("Could not find enough genes in the object from the following gene lists:", 
                      paste(names(x = which(x = !Seurat:::LengthCheck(values = genes.list)))), 
                      "Attempting to match case..."))
        genes.list <- lapply(X = genes.old, FUN = CaseMatch, 
                             match = rownames(x = object@data))
    }
    if (!all(Seurat:::LengthCheck(values = genes.list))) {
        stop(paste("The following gene lists do not have enough genes present in the object:", 
                   paste(names(x = which(x = !Seurat:::LengthCheck(values = genes.list)))), 
                   "exiting..."))
    }
    if (is.null(x = genes.pool)) {
        genes.pool = rownames(x = object@data)
    }
    data.avg <- Matrix::rowMeans(x = object@data[genes.pool, 
                                                 ])
    data.avg <- data.avg[order(data.avg)]
    data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg)/n.bin)))
    names(x = data.cut) <- names(x = data.avg)
    ctrl.use <- vector(mode = "list", length = cluster.length)
    for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        for (j in 1:length(x = genes.use)) {
            ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                                                                                      data.cut[genes.use[j]])], size = ctrl.size, 
                                                               replace = FALSE)))
        }
    }
    ctrl.use <- lapply(X = ctrl.use, FUN = unique)
    ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                          ncol = ncol(x = object@data))
    for (i in 1:length(ctrl.use)) {
        genes.use <- ctrl.use[[i]]
        ctrl.scores[i, ] <- Matrix::colMeans(x = object@data[genes.use, 
                                                             ])
    }
    genes.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                           ncol = ncol(x = object@data))
    for (i in 1:cluster.length) {
        genes.use <- genes.list[[i]]
        data.use <- object@data[genes.use, , drop = FALSE]
        genes.scores[i, ] <- Matrix::colMeans(x = data.use)
    }
    genes.scores.use <- genes.scores - ctrl.scores
    rownames(x = genes.scores.use) <- enrich.name
    genes.scores.use <- as.data.frame(x = t(x = genes.scores.use))
    if(only.pos) genes.scores.use <- genes.scores.use - apply(genes.scores.use,2,min)
    rownames(x = genes.scores.use) <- colnames(x = object@data)
    object <- AddMetaData(object = object, metadata = genes.scores.use, 
                          col.name = colnames(x = genes.scores.use))
    gc(verbose = FALSE)
    return(object)
}


#' find Alias gene names
#' @param df marker data frame with Alias notation
#' @param gene gene name
#' @example Alias(df = df_markers, gene = "CD19")
Alias1 <- function(df, gene = "PDCD1"){
        
        if(!all(grep("Alias|alias",colnames(df)) %%2 ==0)) stop("Check Alias Names!")
        df = as.data.frame(df)
        new_df <- lapply(0:(ncol(df)/2-1), function(i) df[,c(i*2+1, i*2+2)])
        for(i in 1:length(new_df)) colnames(new_df[[i]]) = c("gene","alias")
        new_df <- bind_rows(new_df) %>% 
                .[complete.cases(.),] %>% 
                .[!duplicated(.[,1]),] %>% 
                as.data.frame
        rownames(new_df) = new_df$gene
        if(is.na(new_df[gene,"alias"])) {
                return(NULL)
        } else if(new_df[gene,"alias"] != gene) 
                return(paste0(" (",new_df[gene,"alias"],")"))
}

#' find Alias gene names
#' @param df marker data frame with Alias notation
#' @param gene gene name
#' @example Alias(df = df_markers, gene = "Cd19")
Alias <- function(df, gene = "HLA-DRB1"){
        
        df = as.data.frame(df)
        df_remove_alias = df[,-grep("Alias",colnames(df))]
        species <- CheckSpecies(gene)
        gene = toupper(gene)
        ind <- which(df_remove_alias == gene, arr.ind = TRUE)
        if(nrow(ind) == 0) return(NULL) # no record
        # unique gene name has alias
        ind <- which(df == gene, arr.ind = TRUE)
        if(nrow(ind) == 1) {
                cell_type <- colnames(df)[ind[2]]
                gene.alias <-  df[ind[1],paste0(cell_type,".Alias")]
        } else {
                # duplicate gene name and has no alias
                if(nrow(ind) > 1 & (ind[2,2]-ind[1,2]==1)) return(NULL)
        
                # duplicate gene name and has alias
                if(nrow(ind) > 1 & (ind[2,2]-ind[1,2]!=1)) {
                        cell_type <- colnames(df)[ind[1,2]]
                        gene.alias = df[ind[1,1],paste0(cell_type,".Alias")]
                }
        }
        if(CheckSpecies(gene)=="Mouse") {
            gene.alias = paste0(" (",Hmisc::capitalize(tolower(gene.alias)),")")
        }
        if(CheckSpecies(gene)=="Human") {
            gene.alias = paste0(" (",toupper(tolower(gene.alias)),")")
        }
        return(gene.alias)
        }


#' check sepcies
#' @param subject could be Seurat object, or gene names
#' @example CheckSpecies("CD8A")
#' @example CheckSpecies(object)
CheckSpecies <- function(subject){
        if(class(subject)=="seurat"){
                subject = rownames(subject@raw.data)[1]    
        }
        # Human gene
        if(subject == toupper(subject))
                return("Human")
        # Mouse gene
        if(subject == Hmisc::capitalize(tolower(subject)))
                return("Mouse") 

}

# Combine and print multiple PNG,
#' @param ... PNG path
#' @param ncol grid.arrange argments
#' @param do.print TRUE/FALSE print in device
#' @param do.save TRUE/FALSE save in output/date folder
#' @param bottom_text grid.arrange argments
#' @example CombPngs(GSEA.plots.path.list, ncol = 3)
CombPngs <- function(...,ncol = 2, do.print = FALSE, do.save = TRUE, bottom_text = NULL){
    Img <- lapply(..., function(x) grid::rasterGrob(as.raster(png::readPNG(x)),interpolate = TRUE))
    
    path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    if(do.print) {
        do.call(gridExtra::grid.arrange ,c(Img, ncol = ncol, bottom_text = NULL))
    } else if(do.save) {
        jpeg(paste0(path,deparse(substitute(...)),"_CombPngs.jpeg"), units="in", width=10, height=7,res=600)
        do.call(gridExtra::grid.arrange ,c(Img, ncol = ncol, bottom_text = NULL))
        dev.off()
    }
}


# TRUE Global Environment
#' Produce a contingency table with certain gene expression at different ident
#'
#' @param object seurat object
#' @param subset.name Parameter to subset on. Eg, the name of a gene, PC1, ... 
#' Any argument that can be retreived using FetchData
#' @export Cells_ident nX2 dataframe with orig.ident at column1, Freq at column2
#' @examples
#' CountsbyIdent(SSCs,"Gfra1")
CountsbyIdent <- function(object,subset.name,...){
    "Generate nX2 dataframe with orig.ident at column1, Freq at column2"
    cells.use <- WhichCells(object=object,subset.name=subset.name,...)
    Cells.use <- data.frame("cells"=cells.use,"ident"=sub('_.*', '', cells.use))
    Cells_ident <- as.data.frame(table(Cells.use$ident))
    colnames(Cells_ident)[2] <- subset.name
    return(Cells_ident)
}


#' a supprting function for SingleFeaturePlot.1 and FeatureHeatmap.1
#' Change ggplot color scale to increase contrast gradient
#' #https://github.com/satijalab/seurat/issues/235
#' @param p ggplot object
#' @param alpha.use Define transparency of points
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
ChangeColorScale <- function(p1, alpha.use = 1,legend.title = "log(UMI)",
                             gradient.use = c("yellow", "red"),
                             scaled.expression.threshold = NULL) {
    # Order data by scaled gene expresion level
    # Compute maximum value in gene expression
    if (length(p1$data$scaled.expression)>0){                # FeatureHeatmap.1
        p1$data = p1$data[order(p1$data$scaled.expression),]
        max.scaled.exp <- max(p1$data$scaled.expression)
    } else if (length(p1$data$gene)>0){                   # SingleFeaturePlot.1 
        p1$data = p1$data[order(p1$data$gene),] 
        max.scaled.exp <- max(p1$data$gene) 
    }
    
    # Define lower limit of scaled gene expression level
    if (!is.null(scaled.expression.threshold)) {
        scaled.expression.threshold <- scaled.expression.threshold
    } else if (is.null(scaled.expression.threshold)) {
        if (length(p1$data$scaled.expression)>0){
            scaled.expression.threshold <- min(p1$data$scaled.expression)+0.0001
        } else if (length(p1$data$gene)>0) {
            scaled.expression.threshold <- min(p1$data$gene)+0.0001
        }
    }
    
    # Fill points using the scaled gene expression levels
    p1$layers[[1]]$mapping$fill <- p1$layers[[1]]$mapping$colour
    
    # Define transparency of points
    p1$layers[[1]]$mapping$alpha <- alpha.use
    
    # Change fill and colour gradient values
    p1 = p1 + guides(colour = FALSE)
    p1 = p1 + scale_colour_gradientn(colours = gradient.use, guide = F,
                                   limits = c(scaled.expression.threshold,
                                              max.scaled.exp),
                                   na.value = "grey") +
        scale_fill_gradientn(colours = gradient.use,
                             name = legend.title,
                             limits = c(scaled.expression.threshold,
                                        max.scaled.exp),
                             na.value = "grey") +
        scale_alpha_continuous(range = alpha.use, guide = F)
    
    # Return plot
    return(p1)
}


#' Convert data frame to list
#'
#' This function will convert a data frame to a list, even if they are unequal length
#'
#' @param df
#' @export
#' @examples
#' library(GSVAdata)
#' data(brainTxDbSets)
#' brainTxDbSets_df <- list2df(brainTxDbSets)
#' genelist <- df2list(brainTxDbSets_df)
df2list <- function(df){
    if(is.matrix(df)) df <- as.data.frame(df)
    list <- lapply(df, as.vector) # as.vector! not as.character
    list <- lapply(list, function(x) x[!is.na(x)])
    list <- lapply(list, function(x) x[!(x == "")])
    list <- lapply(list, function(x) gsub("\\s","", x)) #remove space
    names(list) <- names(df)
    return(list)
}


DimElbowPlot.1 <- function (object, reduction.type = "pca", dims.plot = 20,slot = "sdev",
                            xlab = "", ylab = "", title = "") 
{
    data.use <- GetDimReduction(object = object, reduction.type = reduction.type, 
                                slot = slot)
    if (length(data.use) == 0) {
        stop(paste("No standard deviation info stored for", 
                   reduction.type))
    }
    if (length(x = data.use) < dims.plot) {
        warning(paste("The object only has information for", 
                      length(x = data.use), "PCs."))
        dims.plot <- length(x = data.use)
    }
    data.use <- data.use[1:dims.plot]
    dims <- 1:length(x = data.use)
    data.plot <- data.frame(dims, data.use)
    plot <- ggplot(data = data.plot, mapping = aes(x = dims, 
                                                   y = data.use)) + geom_point()
    if (reduction.type == "pca") {
        plot <- plot + labs(y = "Standard Deviation of PC", 
                            x = "PC", title = title)
    }
    else if (reduction.type == "ica") {
        plot <- plot + labs(y = "Standard Deviation of IC", 
                            x = "IC", title = title)
    }
    else {
        plot <- plot + labs(y = ylab, x = xlab, title = title)
    }
    return(plot)
}


#' DoHeatmap.1, automatically group_top by cluster, order by Time_points
#example    top <- group_top_mutate(df = All_markers, major_cells)
#           top %>% head(20) %>% kable() %>% kable_styling()
#           DoHeatmap.1(SSCs,top,Top_n = 15, 
#                   group.order = major_cells,ident.use = "all cell types",
#                   group.label.rot = T,cex.row = 5,remove.key =T)
#
DoHeatmap.1 <- function(object, marker_df, add.genes = NULL, gene.sets = NULL,
                        Top_n = 10, group.order = NULL,ident.use,
                        group.label.rot =T,cex.row = 8,remove.key =T,use.scaled = T,
                        group.cex = 13,title.size = 14,do.print = FALSE,...){
    if (!missing(x = marker_df)) {
        colnames(marker_df)[8] = "cluster"
        top <-  marker_df %>% 
            group_by(cluster) %>% 
            top_n(Top_n, avg_logFC)
        } else { 
            top <- list()
            top$gene = NULL
            }
    if(!is.null(add.genes)) {
        top_gene = c(as.character(top$gene),add.genes)
    } else top_gene <- top$gene
    heatmap <- DoHeatmap(object = object, genes.use = top_gene, #gene.sets = gene.sets,
                         group.order = group.order, use.scaled = use.scaled,
                         slim.col.label = TRUE, remove.key = remove.key,cex.row = cex.row,
                         group.cex = group.cex, rotate.key = T,group.label.rot = group.label.rot,
                         title = paste("Top",Top_n,
                                       "differentially expressed genes in",
                                       ident.use),...)
    heatmap = heatmap+ theme(plot.title = element_text(size=title.size))
    if(do.print){
            path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
            if(!dir.exists(path)) dir.create(path, recursive = T)
            jpeg(paste0(path,ident.use,"_.jpeg"), units="in", width=10, height=7,res=600)
            print(heatmap)
            dev.off()
    } else return(heatmap)
}


#' Modified Seurat::DimPlot function to add ggrepel::geom_text_repel and ggrepel::geom_label_repel
#' A supporting function for TSNEPlot.1
#' @param object Seurat object
#' @param text.repel Adds text directly to the plot.
#' @param label.repel Draws a rectangle underneath the text, making it easier to read.
#' @param force Force of repulsion between overlapping text labels. Defaults to 1.
#' @param ... all other paramethers are the same as Seurat::DimPlot
#' @export p/p3 ggplot object
#' @example DimPlot.1(SSCs, reduction.use = "tsne")
DimPlot.1 <- function (object, reduction.use = "pca", dim.1 = 1, dim.2 = 2,
                       cells.use = NULL, alpha = 1,pt.size = 1, do.return = TRUE, do.bare = FALSE,do.print =FALSE,
                       colors.use = NULL,cols.use = colors.use, group.by = "ident", pt.shape = NULL, do.hover = FALSE,
                       data.hover = "ident", do.identify = FALSE, do.label = FALSE,
                       label.size = 4, no.legend = FALSE, coord.fixed = FALSE,
                       no.axes = FALSE, dark.theme = FALSE, plot.order = NULL,
                       cells.highlight = NULL, cols.highlight = "red", sizes.highlight = 1,
                       plot.title = NULL, vector.friendly = FALSE, png.file = NULL,
                       png.arguments = c(10, 10, 100), na.value = "grey50",
                       text.repel = TRUE, label.repel = FALSE,force=1, ...)
{
    if (vector.friendly) {
        previous_call <- blank_call <- png_call <- match.call()
        blank_call$pt.size <- -1
        blank_call$do.return <- TRUE
        blank_call$vector.friendly <- FALSE
        png_call$no.axes <- TRUE
        png_call$no.legend <- TRUE
        png_call$do.return <- TRUE
        png_call$vector.friendly <- FALSE
        png_call$plot.title <- NULL
        blank_plot <- eval(blank_call, sys.frame(sys.parent()))
        png_plot <- eval(png_call, sys.frame(sys.parent()))
        png.file <- Seurat:::SetIfNull(x = png.file, default = paste0(tempfile(),
                                                                      ".png"))
        ggsave(filename = png.file, plot = png_plot, width = png.arguments[1],
               height = png.arguments[2], dpi = png.arguments[3])
        to_return <- AugmentPlot(plot1 = blank_plot, imgFile = png.file)
        file.remove(png.file)
        if (do.return) {
            return(to_return)
        }
        else {
            print(to_return)
        }
    }
    embeddings.use <- GetDimReduction(object = object, reduction.type = reduction.use,
                                      slot = "cell.embeddings")
    if (length(x = embeddings.use) == 0) {
        stop(paste(reduction.use, "has not been run for this object yet."))
    }
    cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = object@data))
    dim.code <- GetDimReduction(object = object, reduction.type = reduction.use,
                                slot = "key")
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(x = embeddings.use)
    cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
    data.plot <- data.plot[cells.use, dim.codes]
    ident.use <- as.factor(x = object@ident[cells.use])
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(object = object,
                                             vars.all = group.by)[cells.use, 1])
    }
    data.plot$ident <- ident.use
    data.plot$x <- data.plot[, dim.codes[1]]
    data.plot$y <- data.plot[, dim.codes[2]]
    data.plot$pt.size <- pt.size
    if (!is.null(x = cells.highlight)) {
        if (is.character(x = cells.highlight)) {
            cells.highlight <- list(cells.highlight)
        }
        else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
            cells.highlight <- as.list(x = cells.highlight)
        }
        cells.highlight <- lapply(X = cells.highlight, FUN = function(cells) {
            cells.return <- if (is.character(x = cells)) {
                cells[cells %in% rownames(x = data.plot)]
            }
            else {
                cells <- as.numeric(x = cells)
                cells <- cells[cells <= nrow(x = data.plot)]
                rownames(x = data.plot)[cells]
            }
            return(cells.return)
        })
        cells.highlight <- Filter(f = length, x = cells.highlight)
        if (length(x = cells.highlight) > 0) {
            if (!no.legend) {
                no.legend <- is.null(x = names(x = cells.highlight))
            }
            names.highlight <- if (is.null(x = names(x = cells.highlight))) {
                paste0("Group_", 1L:length(x = cells.highlight))
            }
            else {
                names(x = cells.highlight)
            }
            sizes.highlight <- rep_len(x = sizes.highlight,
                                       length.out = length(x = cells.highlight))
            cols.highlight <- rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
            highlight <- rep_len(x = NA_character_, length.out = nrow(x = data.plot))
            if (is.null(x = cols.use)) {
                cols.use <- "black"
            }
            cols.use <- c(cols.use[1], cols.highlight)
            size <- rep_len(x = pt.size, length.out = nrow(x = data.plot))
            for (i in 1:length(x = cells.highlight)) {
                cells.check <- cells.highlight[[i]]
                index.check <- match(x = cells.check, rownames(x = data.plot))
                highlight[index.check] <- names.highlight[i]
                size[index.check] <- sizes.highlight[i]
            }
            plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
            plot.order[is.na(x = plot.order)] <- "Unselected"
            highlight[is.na(x = highlight)] <- "Unselected"
            highlight <- as.factor(x = highlight)
            data.plot$ident <- highlight
            data.plot$pt.size <- size
            if (dark.theme) {
                cols.use[1] <- "white"
            }
        }
    }
    if (!is.null(x = plot.order)) {
        if (any(!plot.order %in% data.plot$ident)) {
            stop("invalid ident in plot.order")
        }
        plot.order <- rev(x = c(plot.order, setdiff(x = unique(x = data.plot$ident),
                                                    y = plot.order)))
        data.plot$ident <- factor(x = data.plot$ident, levels = plot.order)
        data.plot <- data.plot[order(data.plot$ident), ]
    }
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y, colour =ident)) +
        geom_point(mapping = aes(colour = factor(x = ident),
                                 size = pt.size),
                   alpha = alpha)
    if (!is.null(x = pt.shape)) {
        shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use,
                                                                     1]
        if (is.numeric(shape.val)) {
            shape.val <- cut(x = shape.val, breaks = 5)
        }
        data.plot[, "pt.shape"] <- shape.val
        p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) +
            geom_point(mapping = aes(colour = factor(x = ident),
                                     shape = factor(x = pt.shape), size = pt.size),
                       alpha = alpha)
    }
    if (!is.null(x = cols.use)) {
        p <- p + scale_colour_manual(values = cols.use, na.value = na.value)
    }
    if (coord.fixed) {
        p <- p + coord_fixed()
    }
    p <- p + guides(size = FALSE)
    p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) +
        scale_size(range = c(min(data.plot$pt.size), max(data.plot$pt.size)))
    p3 <- p2 + Seurat:::SetXAxisGG() + Seurat:::SetYAxisGG() +
        Seurat:::SetLegendPointsGG(x = 6) +
        Seurat:::SetLegendTextGG(x = 12) + Seurat:::no.legend.title + theme_bw() +
        Seurat:::NoGrid()
    if (dark.theme) {
        p <- p + DarkTheme()
        p3 <- p3 + DarkTheme()
    }
    p3 <- p3 + theme(legend.title = element_blank())
    if (!is.null(plot.title)) {
        p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
    }
    if (no.legend) {
        p3 <- p3 + theme(legend.position = "none")
    }
    if (do.label) {
        centers <- data.plot %>% dplyr::group_by(ident) %>%
            dplyr::summarize(x = median(x), y = median(y))
        p3 = p3 + geom_point(data = centers, aes(x = x, y = y),
                             size = 0, alpha = alpha)
        if (label.repel == TRUE) {
            p3 = p3 + ggrepel::geom_label_repel(data = centers,
                                                aes(label = ident),
                                                size = label.size,
                                                force = force)
        }
        else if (text.repel == TRUE){
            p3 = p3 + ggrepel::geom_text_repel(data = centers,
                                               aes(label = ident),
                                               size = label.size,
                                               force = force,
                                               color = "black")
        }
        if (label.repel == FALSE & text.repel == FALSE) {
            p3 = p3 + geom_text(data = centers,
                                aes(label = ident),
                                size = label.size,
                                color = "black")
        }
        p3 = p3 + guides(colour = FALSE)
        x.range = layer_scales(p)$x$range$range
        add_to_x = sum(abs(x.range)) * 0.03
        p3 = p3 + xlim(x.range[1] - add_to_x, x.range[2] + add_to_x)
    }
    if (no.axes) {
        p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(),
                         axis.text.y = element_blank(), axis.ticks = element_blank(),
                         axis.title.x = element_blank(), axis.title.y = element_blank(),
                         panel.background = element_blank(), panel.border = element_blank(),
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                         plot.background = element_blank())
    }
    if (do.identify || do.hover) {
        if (do.bare) {
            plot.use <- p
        }
        else {
            plot.use <- p3
        }
        if (do.hover) {
            if (is.null(x = data.hover)) {
                features.info <- NULL
            }
            else {
                features.info <- FetchData(object = object,
                                           vars.all = data.hover)
            }
            return(HoverLocator(plot = plot.use, data.plot = data.plot,
                                features.info = features.info, dark.theme = dark.theme))
        }
        else if (do.identify) {
            return(FeatureLocator(plot = plot.use, data.plot = data.plot,
                                  dark.theme = dark.theme, ...))
        }
    }
    if(do.print){
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        s <- unique(object@meta.data[,FindIdentLabel(object)])
        jpeg(paste0(path,paste(s[1:min(length(s),4)],collapse = "_"),
        ".jpeg"), units="in", width=10, height=7,res=600)
        print(p3)
        dev.off()
    }
    if (do.return) {
        if (do.bare) {
            return(p)
        }
        else {
            return(p3)
        }
    }
    if (do.bare) {
        print(p)
    }
}


#' ExtractMetaColor: extract color code from meta.data
#' @object seurat object, meta.data slot must have "color"
#' @label corresponding label
#' @example 
# ExtractMetaColor(object = MCL, color_index = "singler1sub.colors")
ExtractMetaColor <- function(object,color_index = NULL){
    meta.data =object@meta.data
    if(is.null(color_index)) {
        color_index <- paste0(FindIdentLabel(object),".colors")
    }
    if(length(color_index) == 0) {
        return(NULL)
    }
     else
        if(!any(color_index %in% colnames(meta.data))) {
            return(NULL)
            } else {
                label = sub("colors","",color_index)
                label = sub("[.]$","",label)
                meta.data = meta.data[,c(label,color_index)]
                meta.data$index <- as.numeric(as.factor(meta.data[,1]))
                df_colors = meta.data[!duplicated(meta.data$index),]
                df_colors = df_colors[order(df_colors$index),]
            }
    
    return(as.character(df_colors[,color_index]))
}




#' use Findmark results to generate eulerr data frame for venn diagram
#' @param df Seruat::FindMarkers results.
#' @param shape shape of venn diagram
#' @param key Names of venn diagram catergroy, which must be the sub components of df$cluster.
#'  defaul NULL means all components of df$cluster.
#' @param cut_off choose within c("p_val","p_val_adj","avg_logFC")
#' @param cut_off_value corrsponding cut off value.
#' @param eulerr::euler table
#' @param do.legend TRUE/FALSE
#' @param do.return TRUE/FALSE return plot
#' @param do.print print figures
#' @export g plot object
#' @example eulerr(T_cells_markers,shape =  "ellipse",cut_off = "avg_logFC", cut_off_value = 0.01)
eulerr <- function(df, shape =  "circle", key = NULL,cut_off = "avg_logFC",
                   cut_off_value = 0.05, do.lenged = TRUE,do.return = TRUE, do.print = FALSE,...){
        file.name = deparse(substitute(df))
        df$cluster <- as.vector(df$cluster)
        df$gene <- as.vector(df$gene)
        if(!is.null(key)) df <- df[(df$cluster %in% key),]
        df_list <- split(df,df$cluster)
        
        if(cut_off == "avg_logFC"){
                pos.share_genes <- sapply(df_list, function(df) df[(df$avg_logFC > -cut_off_value),"gene"])
        }  
        if(any(cut_off %in% c("p_val","p_val_adj"))){
                pos_genes <- lapply(df_list, function(df) df[(df$avg_logFC > 0),"gene"])
                shared_genes <- lapply(df_list, function(df) df[(abs(df[,cut_off]) > cut_off_value),"gene"])
                pos.share_genes <- mapply(function(x,y) unique(c(x,y)), pos_genes, shared_genes)
        }
        euler_df <- eulerr::euler(pos.share_genes,shape = shape,...)
        
        g <- plot(euler_df, quantities = TRUE, lty = 1:6,
                  legend = do.lenged, main = paste(cut_off," : ",cut_off_value))
        if(do.print) {
                path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
                if(!dir.exists(path)) dir.create(path, recursive = T)
                jpeg(paste0(path,"Venn_",file.name,"_",cut_off,
                            "_",cut_off_value,".jpeg"), units="in", width=10, height=7,res=600)
                print(g)
                dev.off()
        }
        if(do.return) return(g)
        
}


#' Modified Seurat::FeatureHeatmap function to increase contrast gradient
#' 
#' #https://github.com/satijalab/seurat/issues/235
#' @param object Seurat object
#' @param features.plot gene name
#' @param gradient.use Change fill and colour gradient values
#' @param scaled.expression.threshold Define lower limit of scaled gene expression level
#' @export p ggplot object
#' @example FeatureHeatmap.1(SSCs,"Gfra1")
FeatureHeatmap.1 <- function(object, features.plot, mouse.genes = TRUE,
                             group.by ="orig.ident",  alpha.use = 1,legend.title = "log(UMI)",
                             pt.size = 0.5, scaled.expression.threshold = 0.1,
                             gradient.use = c("grey50", "red4"),
                             pch.use = 20) {
    
    if (mouse.genes) features.plot = FilterGenes(object =  object, unique = T,
                                                marker.genes = features.plot)
    else features.plot = FilterGenes(object =  object, unique = T,
                                    marker.genes = features.plot)
    
    x <- FeatureHeatmap(object = object, features.plot = features.plot,
                        group.by = group.by, sep.scale = T, pt.size = pt.size, 
                        cols.use = c("gray98", "red"), pch.use = pch.use, do.return = T)
    p <- ChangeColorScale(x, alpha.use = alpha.use,legend.title=legend.title,
                          scaled.expression.threshold = scaled.expression.threshold,
                          gradient.use = gradient.use)
    return(p)
    
}


.FeaturePlot <- function (object, features.plot, min.cutoff = NA, max.cutoff = NA, 
                         dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1, cols.use = c("yellow", 
                                                                                           "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE, 
                         data.hover = "ident", do.identify = FALSE, reduction.use = "tsne", 
                         use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE, 
                         coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE, 
                         vector.friendly = FALSE, png.file = NULL, png.arguments = c(10, 
                                                                                     10, 100)) 
{
    cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = object@data))
    if (is.null(x = nCol)) {
        nCol <- 2
        if (length(x = features.plot) == 1) {
            nCol <- 1
        }
        if (length(x = features.plot) > 6) {
            nCol <- 3
        }
        if (length(x = features.plot) > 9) {
            nCol <- 4
        }
    }
    num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 
        1
    if (overlay | do.hover) {
        num.row <- 1
        nCol <- 1
    }
    par(mfrow = c(num.row, nCol))
    dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                slot = "key")
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(GetCellEmbeddings(object = object, 
                                                 reduction.type = reduction.use, dims.use = c(dim.1, 
                                                                                              dim.2), cells.use = cells.use))
    x1 <- paste0(dim.code, dim.1)
    x2 <- paste0(dim.code, dim.2)
    data.plot$x <- data.plot[, x1]
    data.plot$y <- data.plot[, x2]
    data.plot$pt.size <- pt.size
    names(x = data.plot) <- c("x", "y")
    data.use <- t(x = FetchData(object = object, vars.all = features.plot, 
                                cells.use = cells.use, use.imputed = use.imputed))
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
                                                            ]), no = cutoff)
    }, cutoff = min.cutoff, feature = features.plot)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
                                                            ]), no = cutoff)
    }, cutoff = max.cutoff, feature = features.plot)
    check_lengths = unique(x = vapply(X = list(features.plot, 
                                               min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check_lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    if (overlay) {
        pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot, 
                                data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                                cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff, 
                                max.cutoff = max.cutoff, coord.fixed = coord.fixed, 
                                no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
    }
    else {
        pList <- mapply(FUN = SingleFeaturePlot.1, feature = features.plot, 
                        min.cutoff = min.cutoff, max.cutoff = max.cutoff, 
                        coord.fixed = coord.fixed,
                        MoreArgs = list(data.use = data.use, 
                                        data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                                        cols.use = cols.use, dim.codes = dim.codes, 
                                        no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme, 
                                        vector.friendly = vector.friendly, png.file = png.file, 
                                        png.arguments = png.arguments), SIMPLIFY = FALSE)
    }
    if (do.hover) {
        if (length(x = pList) != 1) {
            stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
        }
        if (is.null(x = data.hover)) {
            features.info <- NULL
        }
        else {
            features.info <- FetchData(object = object, vars.all = data.hover)
        }
        return(HoverLocator(plot = pList[[1]], data.plot = data.plot, 
                            features.info = features.info, dark.theme = dark.theme, 
                            title = features.plot))
    }
    else if (do.identify) {
        if (length(x = pList) != 1) {
            stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
        }
        return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, 
                              dark.theme = dark.theme))
    }
    else {
        print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
    }
    ResetPar()
    if (do.return) {
        return(pList)
    }
}


#' Combine FindAllMarkers and calculate average UMI
#' Modified Seurat::FindAllMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' #' @export gde.all data frame
#' @example FindAllMarkers.UMI(SSCs)
FindAllMarkers.UMI <- function (object, genes.use = NULL, logfc.threshold = 0.25, 
                                test.use = "MAST", min.pct = 0.1, min.diff.pct = -Inf, 
                                print.bar = TRUE, only.pos = TRUE, max.cells.per.ident = Inf, 
                                return.thresh = 0.01, do.print = FALSE, random.seed = 1, 
                                min.cells.gene = 3, min.cells.group = 3, latent.vars = "nUMI",
                                assay.type = "RNA", get.slot = "data",...)
{
    data.1 <- GetAssayData(object = object, assay.type = assay.type, 
                           slot = "data")
    genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.1))
    if ((test.use == "roc") && (return.thresh == 0.01)) {
        return.thresh = 0.7
    }
    idents.all <- sort(x = unique(x = object@ident))
    genes.de <- list()
    for (i in 1:length(x = idents.all)) {
        genes.de[[i]] <- tryCatch({
            FindMarkers.UMI(object = object, assay.type = assay.type, 
                        ident.1 = idents.all[i], ident.2 = NULL, genes.use = genes.use, 
                        logfc.threshold = logfc.threshold, test.use = test.use, 
                        min.pct = min.pct, min.diff.pct = min.diff.pct, 
                        print.bar = print.bar, min.cells.gene = min.cells.gene, 
                        min.cells.group = min.cells.group, latent.vars = latent.vars, 
                        max.cells.per.ident = max.cells.per.ident, ...)
        }, error = function(cond) {
            return(NULL)
        })
        if (do.print) {
            message(paste("Calculating cluster", idents.all[i]))
        }
    }
    gde.all <- data.frame()
    for (i in 1:length(x = idents.all)) {
        if (is.null(x = unlist(x = genes.de[i]))) {
            next
        }
        gde <- genes.de[[i]]
        if (nrow(x = gde) > 0) {
            if (test.use == "roc") {
                gde <- subset(x = gde, subset = (myAUC > return.thresh | 
                                                     myAUC < (1 - return.thresh)))
            }
            else {
                gde <- gde[order(gde$p_val, -gde$avg_logFC), 
                           ]
                gde <- subset(x = gde, subset = p_val < return.thresh)
            }
            if (nrow(x = gde) > 0) {
                gde$cluster <- idents.all[i]
                gde$gene <- rownames(x = gde)
            }
            if (nrow(x = gde) > 0) {
                gde.all <- rbind(gde.all, gde)
            }
        }
    }
    if ((only.pos) && nrow(gde.all) > 0) {
        return(subset(x = gde.all, subset = avg_logFC > 0))
    }
    rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
    if (nrow(gde.all) == 0) {
        warning("No DE genes identified.")
    }
    return(gde.all)
}


#' FindIdentLabel: Find identical label between ident and metadata
#' @object seurat object
#' @label colname in metadata
FindIdentLabel <- function(object){
    ident.label <- as.vector(object@ident)
    labels <- sapply(object@meta.data,
                     function(x) all(ident.label == x)) %>% .[.] %>% .[!is.na(.)]
    return(names(labels[labels])[1])
}


#' Calculate average UMI and attach to FindMarkers results 
#' Modified Seurat::FindMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
#' @example FindMarkers.UMI(SSCs,ident.1 = "1")
FindMarkers.UMI <- function (object, ident.1, ident.2 = NULL, genes.use = NULL,
                             logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1,
                             min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
                             max.cells.per.ident = Inf, random.seed = 1, latent.vars = "nUMI",
                             min.cells.gene = 3,min.cells.group=3, pseudocount.use = 1, 
                             assay.type = "RNA",...)
{
    genes.de <- FindMarkers(object = object, assay.type = assay.type, 
                            ident.1 = ident.1, ident.2 = ident.2, genes.use = genes.use, 
                            logfc.threshold = logfc.threshold/log2(exp(1)), test.use = test.use, 
                            min.pct = min.pct, min.diff.pct = min.diff.pct, 
                            print.bar = print.bar, only.pos = only.pos, min.cells.gene = min.cells.gene, 
                            min.cells.group = min.cells.group, latent.vars = latent.vars, 
                            max.cells.per.ident = max.cells.per.ident,...)
    genes.de$avg_logFC = log2(exp(1)) * genes.de$avg_logFC
    ave_UMI.1 <- Matrix::rowMeans(expm1(x = object@data[, WhichCells(object = object,
                                                            ident = ident.1)]))
    if(is.null(ident.2)) {
        ave_UMI.2 <- Matrix::rowMeans(expm1(x = object@data[, WhichCells(object = object,
                                                                        ident.remove = ident.1)]))
    } else ave_UMI.2 <- Matrix::rowMeans(expm1(x = object@data[, WhichCells(object = object,
                                                                       ident = ident.2)]))
    
    avg_UMI <-data.frame(ave_UMI.1, ave_UMI.2)
    genes.de <- cbind(genes.de,avg_UMI[match(rownames(genes.de),rownames(avg_UMI)),])

    return(genes.de)
}


#' find marker across by conditions
#' Modified FindMarkers.UMI function, compare the same ident across conditions
#' @param ident.1 dent.1 list
#' @param ident.2 dent.2 list
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
#' @export save.path folder to save
#' @example FindPairMarkers(object, ident.1 = 1:8, ident.2 = c(5:8,1:4))
FindPairMarkers <- function(object, ident.1, ident.2 = NULL, genes.use = NULL,return.thresh = 0.05,
                            logfc.threshold = 0.05, test.use = "MAST", min.pct = 0.1,
                            min.diff.pct = -Inf, print.bar = TRUE, only.pos = FALSE,
                            max.cells.per.ident = Inf, random.seed = 1, latent.vars = "nUMI",
                            min.cells.gene = 3,min.cells.group=3, pseudocount.use = 1, 
                            assay.type = "RNA",save.path = NULL,save.files = TRUE,...){
    #prepare save folder name
    if(class(ident.1) == "numeric" & class(ident.2) == "numeric") {
        ident1="";ident2=""
    } else {
        ident1 <- unique(gsub('\\_.*', '', ident.1))
        ident2 <- unique(gsub('\\_.*', '', ident.2))
    }
    if(is.null(save.path)){
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        save.path <- paste0(path,ident1,"_vs_",ident2,"/")
    }
    gde <- list()
    for(i in 1:length(ident.1)) {
        ident.1vs2 <- paste(ident.1[i], ident.2[i], sep = ".vs.")
        print(ident.1vs2)
        num1 <- unique(gsub('.*\\_', '', ident.1[i]))
        num2 <- unique(gsub('.*\\_', '', ident.2[i]))
        if(num1 == num2) num1 = ""
        gde[[i]] <- FindMarkers.UMI(object = object, ident.1 = ident.1[i],
                                    ident.2 = ident.2[i], assay.type = assay.type, 
                                    genes.use = genes.use, 
                                    logfc.threshold = logfc.threshold, test.use = test.use, 
                                    min.pct = min.pct, min.diff.pct = min.diff.pct, 
                                    print.bar = print.bar, only.pos = only.pos, min.cells.gene = min.cells.gene, 
                                    min.cells.group = min.cells.group, latent.vars = latent.vars, 
                                    max.cells.per.ident = max.cells.per.ident)
        gde[[i]] <- gde[[i]][order(-gde[[i]]$avg_logFC,gde[[i]]$p_val),]
        gde[[i]] <- subset(x = gde[[i]], subset = p_val < return.thresh)
        gde[[i]]$cluster1.vs.cluster2 <- ident.1vs2
        gde[[i]]$gene <- rownames(x = gde[[i]])
        if(save.files){
            if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
            write.csv( gde[[i]], paste0(save.path,ident.1[i],"_vs_",ident.2[i],".csv"))
        }
    }
    return(bind_rows(gde))
}


# FilterGenes
#' filter gene names according to seurat object, produce uniformed gene format
#' pass non-existing genes or ill-formaed genes names to downstream analysis
#' will generate error
#' @param object Seurat object version 2
#' @param marker.genes gene names, marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
#' @param unique return unique name or not
#' @export marker.genes filtered, well-format gene names
#' @example FilterGenes(object, c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr"))
FilterGenes <- function(object, marker.genes, unique= TRUE){
        if(missing(object)) 
                stop("A seurat object must be provided first")
        if(class(object) != "seurat") 
                stop("A seurat object must be provided first")
        if(missing(marker.genes)) 
                stop("A list of marker genes must be provided")

                marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        
        species = CheckSpecies(object)
        if(species == "Human") marker.genes <- toupper(marker.genes)
        if(species == "Mouse") marker.genes <- Hmisc::capitalize(tolower(marker.genes)) 

        print(paste("Before filtration:",length(marker.genes)))
        marker.genes <- CaseMatch(search = marker.genes, match = rownames(x = object@raw.data))
        if(unique) marker.genes <- unique(marker.genes)
        print(paste("After filtration:",length(marker.genes)))
        return(as.character(marker.genes))
}


FPKM <- function(counts, lengths) {
    rownames(counts) = tolower(rownames(counts))
    names(lengths) = tolower(names(lengths))
    A = intersect(rownames(counts), names(lengths))
    counts = counts[A, ]
    lengths = lengths[A]
    rate = counts/lengths
    sweep(rate,2,FUN="/",STATS=colSums(counts)) * 1e+06
}

#=====Clean memory======================
GC <- function()
{
    while (gc()[2, 4] != gc()[2, 4] | gc()[1, 4] != gc()[1,
                                                         4]) {
    }
}

#' Extract ColorHexa from Seurat TSNE plot
#' @param object aligned seurat object with ident
#' @param return.vector TRUE/return full color vector, FALSE/return color levels
#' @param cells.use only return ColorHexa of selected cells
#' @param ... other TSNEPlot inputs 
#' @export colors: color vector named by cell ID
gg_colors <- function(object = object, return.vector=FALSE, cells.use = NULL,
                      print.plot = T, no.legend = TRUE, do.label = TRUE,colors.use =NULL,
                      do.return = TRUE, label.size = 6, gg_title="", ...){
    
    g1 <- Seurat::TSNEPlot(object = object, no.legend = no.legend,
                           do.label = do.label,do.return = do.return,
                           label.size = label.size, colors.use= colors.use,...)
    if(print.plot) print(g1)
    g <- ggplot2::ggplot_build(g1)
    #        print(unique(g$data[[1]]["colour"]))
    colors_df <- g$data[[1]][,c("colour","group")]
    
    #select color by cells ID
    cells <- Seurat::WhichCells(object)
    colors_df$cell.names <- cells
    if(!is.null(cells.use)) {
        colors_df <- colors_df[(colors_df$cell.names %in% cells.use),]
    }
    colors_df = colors_df[order(colors_df$group),]
    if(return.vector) {
        print(head(colors));print(length(colors))
        return(colors)
    } else {
        colors <- unique(colors_df$colour)
        print(colors)
        return(as.character(colors))
    }
}



# select ggplot color
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
#gg_color_hue(4)


#' group_by + top_n + mutate + re-arrange data frame
#' @param df  a data frame from FindAllMarkers
#' @param ... Name-value pairs of expressions, same as ... in dplyr::mutate
#' @param Top_n number of rows to return, same as n in dplyr::top_n. Default is NULL, return all rows.
#' @export top a re-arranged data frame, sorted by avg_logFC, then arrange by ...
#' @example 
#' major_cells <- c("Spermatogonia","Early Spermatocytes","Spermatocytes","Spermatids")
#' top <- group_top(df = All_markers, major_cells)
group_top_mutate <- function(df, ..., Top_n = 500){
    rownames(df) = NULL
    df <- df %>% dplyr::select("gene", everything()) # Moving the last column to the start
    new.col = deparse(substitute(...))
    new.order = assign(new.col,...)
    if(class(df) != "data.frame") df = as.data.frame(df)
    top <-  df %>% 
        dplyr::select("gene", everything()) %>%
        group_by(cluster) %>% 
        top_n(Top_n, wt = avg_logFC) %>%
        mutate(new.col = factor(cluster, levels = new.order)) %>%
        arrange(new.col)
    colnames(top)[which(colnames(top) == "new.col")] = new.col
    return(as.data.frame(top))
}


# modified GenePlot
# GenePlot.1(do.hover = TRUE) will return ggplot 
# including title
GenePlot.1 <- function (object, gene1, gene2, cell.ids = NULL, col.use = NULL, 
                        pch.use = 16, pt.size = 2, use.imputed = FALSE, use.scaled = FALSE, 
                        use.raw = FALSE, do.hover = TRUE, data.hover = "ident",  no.legend = TRUE,
                        do.identify = FALSE, dark.theme = FALSE, do.spline = FALSE, 
                        spline.span = 0.75, title = "",...) 
{
    cell.ids <- Seurat:::SetIfNull(x = cell.ids, default = object@cell.names)
    data.use <- as.data.frame(x = FetchData(object = object, 
                                            vars.all = c(gene1, gene2), cells.use = cell.ids, use.imputed = use.imputed, 
                                            use.scaled = use.scaled, use.raw = use.raw))
    data.plot <- data.use[cell.ids, c(gene1, gene2)]
    names(x = data.plot) <- c("x", "y")
    ident.use <- as.factor(x = object@ident[cell.ids])
    if (length(x = col.use) > 1) {
        col.use <- col.use[as.numeric(x = ident.use)]
    }
    else {
        col.use <- Seurat:::SetIfNull(x = col.use, default = as.numeric(x = ident.use))
    }
    gene.cor <- round(x = cor(x = data.plot$x, y = data.plot$y), 
                      digits = 2)
    if (dark.theme) {
        par(bg = "black")
        col.use <- sapply(X = col.use, FUN = function(color) ifelse(test = all(col2rgb(color) == 
                                                                                   0), yes = "white", no = color))
        axes = FALSE
        col.lab = "white"
    }
    else {
        axes = TRUE
        col.lab = "black"
    }
    
    if (dark.theme) {
        axis(side = 1, at = NULL, labels = TRUE, col.axis = col.lab, 
             col = col.lab)
        axis(side = 2, at = NULL, labels = TRUE, col.axis = col.lab, 
             col = col.lab)
    }
    if (do.spline) {
        spline.fit <- smooth.spline(x = data.plot$x, y = data.plot$y, 
                                    df = 4)
        loess.fit <- loess(formula = y ~ x, data = data.plot, 
                           span = spline.span)
        points(x = data.plot$x, y = loess.fit$fitted, col = "darkblue")
    }
    if (do.identify | do.hover) {
        p <- ggplot2::ggplot(data = data.plot, mapping = aes(x = x, 
                                                             y = y))
        p <- p + geom_point(mapping = aes(x = x,y=y, color = col.use), size = pt.size, 
                            shape = pch.use)
        p <- p + labs(title = paste0(title,"\n",gene.cor), x = gene1, y = gene2)
        if (no.legend) {
            p <- p + theme(legend.position = "none")
        }
        return(p)
    }
}


# HumanGenes
# turn list of gene character to uniformed Human gene list format
HumanGenes <- function(object, marker.genes, unique= TRUE){
    # marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
    if(missing(object)) 
        stop("A seurat object must be provided first")
    if(class(object) != "seurat") 
        stop("A seurat object must be provided first")
    if(missing(marker.genes)) 
        stop("A list of marker genes must be provided")
    if(rownames(object@raw.data)[1] == Hmisc::capitalize(tolower(rownames(object@raw.data)[1])))
        stop("This is mouse genome, use MouseGenes() instead!")
    
    marker.genes <- as.character(marker.genes)
    marker.genes <- unlist(strsplit(marker.genes,","))
    marker.genes <- toupper(marker.genes)
    print(paste("Before filtration:",length(marker.genes)))
    marker.genes <- CaseMatch(search = marker.genes, match = rownames(x = object@raw.data))
    if(unique) marker.genes <- unique(marker.genes)
    print(paste("After filtration:",length(marker.genes)))
    return(as.character(marker.genes))
}


#' Convert list to data frame
#'
#' This function will convert a list to a data frame, even if they are unequal length
#'
#' @param fname list
#' @export
#' @examples
#' library(GSVAdata)
#' data(brainTxDbSets)
#' brainTxDbSets_df <- list2df(brainTxDbSets)
list2df <- function(list){
    df <- do.call(rowr::cbind.fill, c(list, fill = NA))
    names(df) = names(list)
    return(df)
}


# make corlor bar for DoHeatmap
#' @param object Seurat object
#' @param marker_df FindAllMarkers results
#' @param Top_n top_n(Top_n, avg_logFC)
#' @param add.genes extra genes to add beyound FindAllMarkers results 
#' @param remove.legend TRUE/FALSE
#' @param color color scheme
#' @export g vertical ggplot
#' @example MakeCorlorBar(EC, top, Top_n=40)
MakeCorlorBar <- function(object, marker_df, Top_n = NULL, add.genes = NULL, color =NULL,
                          remove.legend =F, do.print = TRUE,do.return=FALSE){
    
        if(!is.null(Top_n)){
        colnames(marker_df)[8] ="cluster"
        marker_df <- marker_df %>% 
        group_by(cluster) %>% 
        top_n(Top_n, avg_logFC)
        }
        marker_df <- marker_df[,c("gene","cluster")]
        
        if(!is.null(add.genes)){
        marker_bar <- data.frame("gene" = add.genes,
                                 "cluster" = "marker.genes")
        marker_df = rbind.data.frame(marker_df[,c("gene","cluster")],
                                     marker_bar)
        }

    gene.use = rownames(object@scale.data)
    marker_df = marker_df[!duplicated(marker_df$gene),]
    marker_df = marker_df[marker_df$gene %in% gene.use,]
    
    marker_df$x <- 1:nrow(marker_df)
    
    wide <- table(marker_df$cluster) %>% as.data.frame
    colnames(wide) = c("cluster","w")
    
    marker_df %<>%   dplyr::left_join(wide, by = "cluster") %>%
            group_by(cluster) %>%
            dplyr::summarise(median = median(x, na.rm = TRUE)) %>%
            dplyr::inner_join(wide, by = "cluster")
    
    g <- ggplot(data = marker_df, aes(xmin = median - w / 2, xmax = median + w / 2,
                                      ymin = 0, ymax = 1, fill = cluster)) +
        geom_rect(colour = "white")+
        theme_bw() + 
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_blank(), 
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())+
        coord_flip() + scale_x_reverse()
    if(remove.legend) g = g + theme(legend.position="none")
    if(!is.null(color)) {
        colors_fill = color_scheme(color = color,
                                   names = unique(marker_df$cluster))
        g = g + scale_fill_manual(values = colors_fill)
    }
    if(do.print) {
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,deparse(substitute(object)),"vertical_bar.jpeg"), units="in", width=10, height=7,res=600)
        print(g)
        dev.off()
    } else if(do.return) {
        return(g)
    }
}


# Make color scheme vector
# names must be ordered vector or factor
color_scheme <- function(color,names){
    df_names <- as.data.frame(table(names))
    df_names_Var <- df_names$Freq
    color_code <- as.character(unlist(mapply(rep, color, df_names_Var)))
    names(color_code) = names
    return(color_code)
}

# clean up the gene names for downstream analysis
# turn list of gene character to uniform mouse gene list format
#' @param object Seurat object
#' @param marker.genes gene names, can be one gene or vector. Must be character
#' @param unique TRUE/FALSE, output unique gene name or not
#' @export marker.genes uniform mouse gene that exsit in object@data
#' @example MouseGenes(mouse_eyes,c("Cdh5","Pecam1","Flt1"))
MouseGenes <- function(object, marker.genes, unique = TRUE){
    # marker.genes =c("Cdh5,Pecam1,Flt1,Vwf,Plvap,Kdr") for example
    if(missing(object))
        stop("A seurat object must be provided first!")
    if(class(object) != "seurat")
        stop("A seurat object must be provided first!")
    if(missing(marker.genes))
        stop("A list of marker genes must be provided!")
    if(rownames(object@raw.data)[1] == toupper(rownames(object@raw.data)[1]))
        stop("This is human genome, use HumanGenes() instead!")
    
    marker.genes <- as.character(marker.genes)
    marker.genes <- Hmisc::capitalize(tolower(marker.genes))
    print(paste("Before filtration:",length(marker.genes)))
    marker.genes <- CaseMatch(search = marker.genes, match = rownames(x = object@raw.data))
    if(unique) marker.genes <- unique(marker.genes)
    print(paste("After filtration:",length(marker.genes)))
    return(as.character(marker.genes))
}


randomStrings <- function(n = 5000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


#' rename ident for certain cells only
#' One drawback, will not double check the duplicated cell names
#' @param object Seurat object
#' @param new.ident.name one character
#' @param cells.use a vector of cell names
#' @export object Seurat object with new ident
#' @example 
#' for(i in 1:length(SSC_labels_id)){
#'         SSCs <- RenameIdent.1(SSCs, new.ident.name = names(SSC_labels_id)[i],
#'                               cells.use = SSC_labels_id[[i]])
#' }
RenameIdent.1 <- function (object, new.ident.name, cells.use) 
{
    new.levels <- old.levels <- levels(x = object@ident)
    if (!(new.ident.name %in% old.levels)) {
        new.levels <- c(new.levels, new.ident.name)
    }
    
    ident.vector <- as.character(x = object@ident)
    names(x = ident.vector) <- names(object@ident)
    ident.vector[cells.use] <- new.ident.name
    object@ident <- factor(x = ident.vector, levels = new.levels)
    
    return(object)
}


# generate expression txt file for GSEA analysis
#' @param object Seurat object
#' @param k an integer for the number of folds. createFolds argment
#' @param do.return TRUE/FALSE
#' @param continuous.label NULL/continuous label #http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_Phenotype_Labels
#' @example PrepareGSEA(object, k = 50, continuous.label = major_cells)
PrepareGSEA <- function(object, k = 1, do.return = FALSE, continuous.label = NULL){
    
    set.seed(201)
    object@scale.data = NULL
    ident = FindIdentLabel(object)[1]
    
    if(k > 1){
        split_object <- SplitSeurat(object, split.by = ident)
        #Split object@ident to k group =====
        for(i in 1:length(split_object)){
            meta_index <- caret::createFolds(split_object[[i]]@meta.data[,ident],
                                         k = k, list = TRUE, 
                                         returnTrain = FALSE)
            for(n in 1:length(meta_index)){
                split_object[[i]]@meta.data[meta_index[[n]],"GSEA"] = 
                    paste(split_object[[i]]@meta.data[meta_index[[n]],ident],
                          n, sep = "_")
                }
        }
        new_object <-  Reduce(MergeSeurat,split_object) %>%
            SetAllIdent(id = "GSEA")
        if(!is.null(continuous.label))
            new_object@ident <- factor(x = new_object@ident,
                                    levels = paste0(rep(continuous.label,each = k),
                                                    "_",rep(1:k)))
        } else new_object <- object
    
        if(k == 1 & !is.null(continuous.label)) {
            new_object@ident <- factor(x = new_object@ident,
                                       levels = continuous.label)
        }

                
    print("#====Calculate Average Expression======")
    GSEA_expr <- AverageExpression(new_object)
    GSEA_name <- data.frame("NAME" = rownames(GSEA_expr),
                            "DESCRIPTION" = rep(NA,nrow(GSEA_expr)),
                            stringsAsFactors = F)
    GSEA_expr <- cbind.data.frame(GSEA_name,GSEA_expr)
    
    samples <- gsub("_([0-9]+).*$", "", colnames(GSEA_expr)[-c(1,2)])
    
    if(is.null(continuous.label)) {
        cls_list <- list(c(length(samples),length(unique(samples)), 1),
                     paste(c("#",unique(samples)), collapse = " "),
                     paste(samples, collapse = " "))
    } else 
        if(all(continuous.label %in% unique(object@ident))){
            numeric <- match(samples,continuous.label)
            cls_list <- list("#numeric",
                             paste(c("#",unique(samples)), collapse = "."),
                             paste(numeric, collapse = " "))
    } else stop("Incorrect continuous.label!")
    
    if(rownames(GSEA_expr)[1] == 
        Hmisc::capitalize(tolower(rownames(GSEA_expr)[1]))){
        print("#====Replace gene names ======")
        rownames.GSEA_expr = rownames(GSEA_expr)
        human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
        mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        
        genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), #filters = "mgi_symbol",
                         values = rownames(GSEA_expr) , mart = mouse,
                         attributesL = c("hgnc_symbol"), 
                         martL = human, uniqueRows=T)
        rm = duplicated(genesV2[,1])
        genesV2 = genesV2[!rm,]
        colnames(genesV2) = c("gene","NAME")
        colnames(GSEA_expr)[1] = "gene"
        GSEA_expr <- merge(genesV2,GSEA_expr,by = "gene")
        GSEA_expr = GSEA_expr[,-1]
    }
    file.name = paste(unique(object@ident), collapse = "_")
    if(!is.null(continuous.label)) file.name <- paste(continuous.label,collapse = "_")
    path <- paste0("output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    write.table(GSEA_expr, file = paste0(path, file.name,
                                         "_",k,".txt"),
                sep = "\t", quote = FALSE,row.names = FALSE)
    
    fn = paste0(path, file.name,"_",k,".cls")
    if (file.exists(fn)) file.remove(fn)
    lapply(cls_list, cat, "\n", file=fn,append=TRUE)
    
    if(do.return) return(GSEA_expr)
    
}


# scale with median, not mean
.scale <- function(df, fun1 = sd, fun2 = median){
    df1 <- sweep(df, 2, apply(df,2, fun1),"/")
    return(sweep(df1, 2, apply(df1,2,fun2),"-"))
}


# scale with median, not mean
ScaleData.1 <- function (object, genes.use = NULL, data.use = NULL, vars.to.regress, 
                         model.use = "linear", use.umi = FALSE, do.scale = TRUE, 
                         do.center = TRUE, scale.max = 10, block.size = 1000, min.cells.to.block = 3000, 
                         display.progress = TRUE, assay.type = "RNA", do.cpp = TRUE, 
                         check.for.norm = TRUE, do.par = FALSE, num.cores = 1) 
{
    data.use <- Seurat:::SetIfNull(x = data.use, default = GetAssayData(object = object, 
                                                               assay.type = assay.type, slot = "data"))
    if (check.for.norm) {
        if (!("NormalizeData" %in% names(object@calc.params))) {
            cat("NormalizeData has not been run, therefore ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.\n")
        }
        if (is.null(object@calc.params$NormalizeData$normalization.method)) {
            cat("ScaleData is running on non-normalized values. Recommended workflow is to run NormalizeData first.\n")
        }
    }
    genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = data.use))
    genes.use <- as.vector(x = intersect(x = genes.use, y = rownames(x = data.use)))
    data.use <- data.use[genes.use, ]
    if (!missing(x = vars.to.regress) && !is.null(x = vars.to.regress)) {
        data.use <- Seurat:::RegressOutResid(object = object, vars.to.regress = vars.to.regress, 
                                    genes.regress = genes.use, use.umi = use.umi, model.use = model.use, 
                                    display.progress = display.progress, do.par = do.par, 
                                    num.cores = num.cores)
        if (model.use != "linear") {
            use.umi <- TRUE
        }
        if (use.umi && missing(scale.max)) {
            scale.max <- 50
        }
    }
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("ScaleData"))]
    parameters.to.store$data.use <- NULL
    object <- Seurat:::SetCalcParams(object = object, calculation = "ScaleData", 
                            ... = parameters.to.store)
    if (!do.cpp) {
        return(.ScaleDataR(object = object, data.use = data.use, 
                          do.scale = do.scale, do.center = do.center, scale.max = scale.max, 
                          genes.use = genes.use))
    }
    scaled.data <- matrix(data = NA, nrow = length(x = genes.use), 
                          ncol = ncol(x = object@data))
    rownames(scaled.data) <- genes.use
    if (length(object@cell.names) <= min.cells.to.block) {
        block.size <- length(genes.use)
    }
    gc()
    colnames(scaled.data) <- colnames(object@data)
    max.block <- ceiling(x = length(x = genes.use)/block.size)
    gc()
    if (display.progress) {
        message("Scaling data matrix")
        pb <- txtProgressBar(min = 0, max = max.block, style = 3, 
                             file = stderr())
    }
    for (i in 1:max.block) {
        my.inds <- ((block.size * (i - 1)):(block.size * i - 
                                                1)) + 1
        my.inds <- my.inds[my.inds <= length(x = genes.use)]
        if (class(x = data.use) == "dgCMatrix" | class(x = data.use) == 
            "dgTMatrix") {
            data.scale <- Seurat:::FastSparseRowScale(mat = data.use[genes.use[my.inds], 
                                                            , drop = F], scale = do.scale, center = do.center, 
                                             scale_max = scale.max, display_progress = FALSE)
        }
        else {
            data.scale <- Seurat:::FastRowScale(mat = as.matrix(x = data.use[genes.use[my.inds], 
                                                                    , drop = F]), scale = do.scale, center = do.center, 
                                       scale_max = scale.max, display_progress = FALSE)
        }
        dimnames(x = data.scale) <- dimnames(x = data.use[genes.use[my.inds], 
                                                          ])
        scaled.data[genes.use[my.inds], ] <- data.scale
        rm(data.scale)
        gc()
        if (display.progress) {
            setTxtProgressBar(pb, i)
        }
    }
    if (display.progress) {
        close(pb)
    }
    object <- SetAssayData(object = object, assay.type = assay.type, 
                           slot = "scale.data", new.data = scaled.data)
    gc()
    object@scale.data[is.na(object@scale.data)] <- 0
    return(object)
}


# scale with median, not mean
.ScaleDataR <- function (object, genes.use = NULL, data.use = NULL, do.scale = TRUE, 
                         do.center = TRUE, scale.max = 10) 
{
    genes.use <- Seurat:::SetIfNull(x = genes.use, default = rownames(x = object@data))
    genes.use <- intersect(x = genes.use, y = rownames(x = object@data))
    data.use <- Seurat:::SSetIfNull(x = data.use, default = object@data[genes.use, 
                                                              ])
    object@scale.data <- matrix(data = NA, nrow = length(x = genes.use), 
                                ncol = ncol(x = object@data))
    dimnames(x = object@scale.data) <- dimnames(x = data.use)
    if (do.scale | do.center) {
        bin.size <- 1000
        max.bin <- floor(length(genes.use)/bin.size) + 1
        message("Scaling data matrix")
        pb <- txtProgressBar(min = 0, max = max.bin, style = 3, 
                             file = stderr())
        for (i in 1:max.bin) {
            my.inds <- ((bin.size * (i - 1)):(bin.size * i - 
                                                  1)) + 1
            my.inds <- my.inds[my.inds <= length(x = genes.use)]
            new.data <- t(x = .scale(x = t(x = as.matrix(x = data.use[genes.use[my.inds], 
                                                                     ])), fun1 = sd, fun2 = mean))
            new.data[new.data > scale.max] <- scale.max
            object@scale.data[genes.use[my.inds], ] <- new.data
            setTxtProgressBar(pb, i)
        }
        close(pb)
    }
    object@scale.data[is.na(object@scale.data)] <- 0
    return(object)
}


SetAllIdent.1 <- function (object, id = NULL) 
{
    id <- Seurat:::SetIfNull(x = id, default = "orig.ident")
    if (id %in% colnames(x = object@meta.data)) {
        cells.use <- rownames(x = object@meta.data)
        ident.use <- factor(x = object@meta.data[, id])
        names(ident.use) = cells.use
        object@ident <- ident.use
    }
    return(object)
}


# FeaturePlot doesn't return ggplot
# SingleFeaturePlot doesn't take seurat object as input
# modified SingleFeaturePlot, take seurat object and return ggplot
SingleFeaturePlot.1 <- function (object = object, feature = feature, pt.size = 1.0,
                                 dim.1 = 1, dim.2 = 2,pch.use = 16, cols.use  = c("lightgrey","blue"),
                                 gradient.use = c("orangered", "red4"),threshold= NA,text.size=15,
                                 cells.use = NULL,dim.codes, min.cutoff = 0, max.cutoff = Inf,
                                 use.imputed = FALSE, reduction.use = "tsne",no.axes = FALSE, no.legend = T, 
                                 dark.theme = FALSE,title=feature, do.return = TRUE,do.print=FALSE,x.lim = NULL,
                                 y.lim = NULL,legend.title = "log(UMI)", ...) 
{
    dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                slot = "key")
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(GetCellEmbeddings(object = object,
                                                 reduction.type = reduction.use, 
                                                 dims.use = c(dim.1,dim.2), 
                                                 cells.use = cells.use))
    x1 <- paste0(dim.code, dim.1)
    x2 <- paste0(dim.code, dim.2)
    data.plot$x <- data.plot[, grep(x1,colnames(data.plot),value = T)]
    data.plot$y <- data.plot[, grep(x2,colnames(data.plot),value = T)]
    data.plot$pt.size <- pt.size
    names(x = data.plot) <- c("x", "y")
    data.use <- t(x = FetchData(object = object, vars.all = feature, 
                                cells.use = cells.use, use.imputed = use.imputed,...))
    data.gene <- na.omit(object = data.frame(data.use[1,])) # Error in data.use[feature, ] : subscript out of bounds
    min.cutoff <- Seurat:::SetQuantile(cutoff = min.cutoff, data = data.gene)
    max.cutoff <- Seurat:::SetQuantile(cutoff = max.cutoff, data = data.gene)
    data.gene <- sapply(X = data.gene, FUN = function(x) {
        return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                      no = x))
    })
    data.gene <- sapply(X = data.gene, FUN = function(x) {
        return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                      no = x))
    })
    data.plot$gene <- data.gene
    
    if (length(x = cols.use) == 1) {
        brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
    }
    else {
        brewer.gran <- length(x = cols.use)
    }
    if (all(data.gene == 0)) {
        data.cut <- 0
    }
    else {
        data.cut <- as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.gene), 
                                                     breaks = brewer.gran)))
    }
    data.plot$col <- as.factor(x = data.cut)
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
    if (brewer.gran != 2) {
        if (length(x = cols.use) == 1) {
            p <- p + geom_point(mapping = aes(color = col), 
                                size = pt.size, shape = pch.use)
        }
        else {
            p <- p + geom_point(mapping = aes(color = col), 
                                size = pt.size, shape = pch.use)
        }
    }
    else {
        if (all(data.plot$gene == data.plot$gene[1])) {
            warning(paste0("All cells have the same value of ", 
                           feature, "."))
            p <- p + geom_point(color = cols.use[1], size = pt.size, 
                                shape = pch.use)
        }
        else {
            p <- p + geom_point(mapping = aes(color = gene), 
                                size = pt.size, shape = pch.use)
        }
    }
    if (no.axes) {
        p <- p + labs(title = feature, x = "", y = "") + theme(axis.line = element_blank(), 
                                                               axis.text.x = element_blank(), 
                                                               axis.text.y = element_blank(), 
                                                               axis.ticks = element_blank(), 
                                                               axis.title.x = element_blank(), 
                                                               axis.title.y = element_blank())
    }
    else {
        p <- p + labs(x = dim.codes[1], y = dim.codes[2])
    }
    if (no.legend) {
        p <- p + theme(legend.position = "none")
    }
    if (dark.theme) {
        p <- p + DarkTheme()
    }
    if(is.null(threshold)) {
        x = object@data[feature,]
        threshold <- histPeak(x)
    }
    p1 <- ChangeColorScale(p, alpha.use = 1,legend.title=legend.title,
                           scaled.expression.threshold = threshold,
                           gradient.use = gradient.use)
    if(!is.null(x.lim)) p1 = p1 + xlim(x.lim)
    if(!is.null(y.lim)) p1 = p1 + ylim(y.lim)
    p1 <- p1 +ggtitle(paste0(title))+
        theme(text = element_text(size=text.size),						
              axis.text.x = element_text(size=text.size*0.8),
              axis.text.y = element_text(size=text.size*0.8),
              plot.title = element_text(hjust = 0.5,size=text.size*1.5),
              legend.key.size = unit(text.size/4, "mm"),
              legend.key.width = unit(text.size/4, "mm"),
              legend.key.height = unit(text.size/4, "mm"))
    if(do.print) {
      path <- paste0("output/",gsub("-","",Sys.Date()),"/")
      if(!dir.exists(path)) dir.create(path, recursive = T)
      jpeg(paste0(path,"SingleFeaturePlot_",feature,".jpeg"), units="in", width=10, height=7,res=600)
      print(p1)
      dev.off()
    } else if(do.return) {
      return(p1)
    }
}


#' Split seurat cell names by certein criteria
#' A supporting funtion to SplitSeurat
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @export cell.subsets list of subseted cell names by certein condition
#' @example SplitCells(mouse_eyes, split.by = "orig.ident")
SplitCells <- function(object = object, split.by = "orig.ident"){
    
    cell.all <- FetchData(object = object, vars.all = split.by)
    if(class(cell.all[,1]) == "numeric"){
        cell.subsets <- list()
        cell.subsets[[1]] <- rownames(cell.all)[cell.all[,split.by] == 0]
        cell.subsets[[2]] <- rownames(cell.all)[cell.all[,split.by] > 0]
    }
    if(class(cell.all[,1]) == "factor"){
        df_table <- as.data.frame(table(cell.all[,1]))
        Levels <- df_table[df_table$Freq >0,"Var1"]
        cell.subsets <- lapply(Levels, function(x)
            rownames(cell.all)[cell.all[,split.by] == x])
    }
    return(cell.subsets)
}


#' Split Seurat by certein criteria and make tsne plot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @param select.plots output order, default to NULL. If want to change,use c(2,1) for example
#' @param return.data TRUE/FASLE, return splited ojbect or not.
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p ggplot object from TSNEplot
#' @example SplitDimPlot(mouse_eyes, split.by = "orig.ident")
#' @example SplitDimPlot(mouse_eyes, split.by = "orig.ident", select.plots = c(2,1))
SplitDimPlot <- function(object = object, split.by = "orig.ident",
                         reduction.use= "TSNE",
                         select.plots = NULL, return.plots = FALSE,
                         do.label = T, group.by = "ident", no.legend = TRUE,
                         pt.size = 1,label.size = 5,... ){
        
        subset.object <- SplitSeurat(object = object, split.by = split.by)
        levels <- object@meta.data[,split.by] %>% unique %>% sort
        
        p <- list()
        if(is.null(select.plots)) select.plots <- 1:length(subset.object)
        for(i in 1:length(select.plots)){
                p[[i]] <- DimPlot.1(object = subset.object[[select.plots[i]]],
                                    reduction.use = reduction.use,
                                    colors.use = ExtractMetaColor(subset.object[[select.plots[i]]]),
                                    do.label = do.label, group.by = group.by,
                                    do.return = T, no.legend = no.legend,
                                    pt.size = pt.size,label.size = label.size,...)+
                        ggtitle(levels[select.plots[i]])+
                        theme(text = element_text(size=20),     #larger text including legend title
                              plot.title = element_text(hjust = 0.5)) #title in middle
        }
        p <- p[lapply(p,length)>0] # remove NULL element
        if(return.plots) return(p) else print(do.call(cowplot::plot_grid, p))
}


#' Split Seurat object by certein criteria
#' A supporting funtion to SplitTSNEPlot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any colname in meta.data
#' @param select.seurat export subset.object in specific numeric order. Default is NULL
#' @export subset.object list of subseted object by certein conditions,
#' @example SplitCells(mouse_eyes, split.by = "orig.ident")
SplitSeurat <- function(object = object, split.by = "orig.ident", 
                        select.seurat = NULL){
    
    cell.subsets <- SplitCells(object = object, split.by = split.by)
    
    subset.object <- list()
    if(is.null(select.seurat)) select.seurat <- 1:length(cell.subsets)
    if(!is.null(select.seurat) & length(select.seurat) > length(cell.subsets))
        stop("Not enough sub seurat to be selected!")
    for(i in 1:length(select.seurat)){
        subset.object[[i]] <- SubsetData(object, cells.use =cell.subsets[[select.seurat[i]]])
    }
    if(class(FetchData(object = object, vars.all = split.by)[1]) == "numeric"){
            names(subset.object) =c("No","Yes")
            } else names(subset.object) <- sort(unique(object@meta.data[,split.by]))
    return(subset.object)
}


#' Split Seurat by certein criteria and make tsne plot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @param select.plots output order, default to NULL. If want to change,use c(2,1) for example
#' @param return.data TRUE/FASLE, return splited ojbect or not.
#' @param colors.use if colors.use is provided, must consider sychronize the color code
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p ggplot object from TSNEplot
#' @example SplitTSNEPlot(mouse_eyes, split.by = "orig.ident")
#' @example SplitTSNEPlot(mouse_eyes, split.by = "Rlbp1", select.plots = c(2,1))
SplitTSNEPlot <- function(object, split.by = "orig.ident",select.plots = NULL, 
                          do.return = TRUE, do.print = FALSE,
                          do.label = T, group.by = "ident", no.legend = TRUE,
                          pt.size = 1,label.size = 5,size=20,...){
    
    
    subset.object <- SplitSeurat(object = object, split.by = split.by)
    levels <- object@meta.data[,split.by] %>% unique %>% sort
    
    p <- list()
    if(is.null(select.plots)) select.plots <- 1:length(subset.object)
    for(i in 1:length(select.plots)){
        p[[i]] <- TSNEPlot.1(object = subset.object[[select.plots[i]]],
                             colors.use = ExtractMetaColor(subset.object[[select.plots[i]]]),
                             do.label = do.label, group.by = group.by,
                             do.return = T, no.legend = no.legend, 
                             pt.size = pt.size,label.size = label.size,...)+
            ggtitle(levels[select.plots[i]])+
            theme(text = element_text(size=size),     #larger text including legend title
                  plot.title = element_text(hjust = 0.5)) #title in middle
        
    }
    p <- p[lapply(p,length)>0] # remove NULL element
    if(do.print) {
        path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"SplitTSNEPlot_",paste(levels,collapse = "_"),"_",
                    FindIdentLabel(object),".jpeg"), units="in", width=10, height=7,res=600)
        print(do.call(cowplot::plot_grid, p))
        dev.off()
    } else if(do.return) {
            return(p)
    }
}


# find and print differentially expressed genes within all major cell types
# combine SubsetData, FindAllMarkers, write.csv
#' @param object Seurat object
#' @param split.by compatible to vars.all in Seurat::FetchData. If split.by is gene name,
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' @export gde.all data frame
#' @example FindAllMarkers.UMI(mouse_eyes)
SplitFindAllMarkers <- function(object, split.by = "orig.ident",test.use = "bimod",
                                write.csv = TRUE,...){
    
    subset.object <- SplitCells(object = object, split.by = split.by)
    conditions <- subset.object[[length(subset.object)]] # levels of conditions
    
    object.markers <- list()
    
    for(i in 1:length(conditions)){
        object.markers[[i]] <- FindAllMarkers.UMI(object = subset.object[[i]],
                                                  test.use = test.use,...)
        if(write.csv) write.csv(object.markers[[i]],
                                file=paste0("./output/",
                                            deparse(substitute(object)),
                                            "_",conditions[i],
                                            ".csv"))
    }
    return(object.markers)
}


# Split Seurat by certein criteria and make tsne plot
#' @param object Seurat object
#' @param split.by the criteria to split, can be gene name, or any variable in meta.data
#' @param select.plots output order, default to NULL. If want to change,use c(2,1) for example
#' @param markers marker vector
#' @param alias NX2 data frame c("gene", "alias")
#' @param threshold NULL or number. if NULL, use soft mean threshold 
#' @param return.data TRUE/FASLE, return splited ojbect or not.
#' @export p ggplot object from barchart
#' @example     SplitSingleFeaturePlot(subset.Glioma,group.by = "ident",split.by = "orig.ident",
#               no.legend = T,label.size=3,do.print =T,markers = markers,
#               threshold = 0.1)
SplitSingleFeaturePlot<- function(object, split.by = "orig.ident",select.plots = NULL, 
                                  markers, alias = NULL, do.return = TRUE, do.print = TRUE,
                                  do.label = T, group.by = "ident", threshold= NA,
                                  pt.size = 1,label.size = 5,size=20,nrow = NULL,
                                  ncol = NULL, x.lim = NULL, y.lim = NULL,...){
    
    
    subset.object <- SplitSeurat(object = object, split.by = split.by)
    levels <- object@meta.data[,split.by] %>% unique %>% sort
    
    if(is.null(select.plots)) select.plots <- 1:length(subset.object)
    for(marker in markers){
        g <- list()
        if(is.null(threshold)) {
                x = object@data[marker,]
                new.threshold <- histPeak(x)
                message(paste0("threshold for ",marker,":", new.threshold))
        } else new.threshold = threshold
        
        if(!is.null(alias)) {
                marker.alias = Alias(alias, marker)
                } else marker.alias = NULL
        
        
        for(i in 1:length(select.plots)){
                g[[i]] <- SingleFeaturePlot.1(object = subset.object[[select.plots[i]]],
                                          threshold= new.threshold,
                                          feature = marker,
                                          title = paste0(levels[select.plots[i]]),
                                          x.lim = x.lim, y.lim = y.lim,...)
                }
        g <- g[lapply(g ,length)>0] # remove NULL element
        if(do.print) {
            path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
            if(!dir.exists(path)) dir.create(path, recursive = T)
            jpeg(paste0(path,"Split_",paste(levels,collapse = "_"),
                        "_",marker,"_",new.threshold,".jpeg"), units="in", width=10, height=7,res=600)
            print(do.call(cowplot::plot_grid, c(g, align = "hv",
                                                nrow = nrow, ncol = ncol))+
                            ggtitle(paste0(marker, marker.alias))+
                            theme(plot.title = element_text(hjust = 0.5,size = 25, 
                                                            face = "bold")))
            print(paste0(which(markers == marker),":",length(markers)))
            dev.off()
        }
    }
}



# find histgram local maximam
histPeak <- function(vector){
        y <- hist(vector)
        df <- data.frame("counts" = y$counts, "breaks" = y$breaks[-length(y$breaks)])
        local.max = df$breaks[(which(diff(sign(diff(df$counts)))==-2))[1]] # local.max
        if(is.na(local.max)) {
                none_zeor_counts <- (df$counts[(df$breaks != 0)])
                local.max = df$breaks[which(none_zeor_counts == max(none_zeor_counts))[1]]
        }
        return(local.max)
}



SplitDotPlotGG.1 <- function (object, grouping.var, genes.plot,
                              gene.groups, cols.use = c("blue","red"),
                              col.min = -2.5, col.max = 2.5, 
                              dot.min = 0, dot.scale = 6,
                              group.by, plot.legend = FALSE,
                              do.return = FALSE, x.lab.rot = FALSE) 
{
    "Fix the mutate_impl error in orignal function"
    if (!missing(x = group.by)) {
        object <- SetAllIdent(object = object, id = group.by)
    }
    grouping.data <- FetchData(object = object, vars.all = grouping.var)[names(x = object@ident), 
                                                                         1]
    ncolor <- length(x = cols.use)
    ngroups <- length(x = unique(x = grouping.data))
    if (ncolor < ngroups) {
        stop(paste("Not enough colors supplied for number of grouping variables. Need", 
                   ngroups, "got", ncolor, "colors"))
    }
    else if (ncolor > ngroups) {
        cols.use <- cols.use[1:ngroups]
    }
    idents.old <- levels(x = object@ident)
    idents.new <- paste(object@ident, grouping.data, sep = "_")
    colorlist <- cols.use
    names(x = colorlist) <- levels(x = grouping.data)
    object@ident <- factor(x = idents.new, 
                           levels = unlist(x = lapply(X = idents.old,
                                                      FUN = function(x) {
                                                          lvls <- list()
                                                          for (i in seq_along(along.with = levels(x = grouping.data))) {
                                                              lvls[[i]] <- paste(x, levels(x = grouping.data)[i], 
                                                                                 sep = "_")
                                                          }
                                                          return(unlist(x = lvls))
                                                      })), ordered = TRUE)
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident
    data.to.plot <- data.to.plot %>% tidyr::gather(key = genes.plot, 
                                                   value = expression, -c(cell, id))
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
        summarize(avg.exp = ExpMean(x = expression), pct.exp = Seurat:::PercentAbove(x = expression, 
                                                                                     threshold = 0))
    data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
        mutate(avg.exp = scale(x = avg.exp)) %>% 
        mutate(avg.exp.scale = as.numeric(x = cut(x = MinMax(data = avg.exp,
                                                             max = col.max,
                                                             min = col.min),
                                                  breaks = 20)))
    data.to.plot <- data.to.plot %>% tidyr::separate(col = id, into = c("ident1", 
                                                                        "ident2"), sep = "_") %>% rowwise() %>% mutate(palette.use = colorlist[[ident2]],
                                                                                                                       ptcolor = colorRampPalette(colors = c("grey", palette.use))(20)[avg.exp.scale]) %>% 
        tidyr::unite("id", c("ident1", "ident2"), sep = "_")
    data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                      levels = rev(x = sub(pattern = "-", replacement = ".", 
                                                           x = genes.plot)))
    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    data.to.plot$id <- factor(x = data.to.plot$id, levels = levels(object@ident))
    palette.use <- unique(x = data.to.plot$palette.use)
    if (!missing(x = gene.groups)) {
        names(x = gene.groups) <- genes.plot
        data.to.plot <- data.to.plot %>% mutate(gene.groups = gene.groups[genes.plot])
    }
    p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                   y = id)) + geom_point(mapping = aes(size = pct.exp, 
                                                                                       color = ptcolor)) + scale_radius(range = c(0, dot.scale)) + 
        scale_color_identity() + theme(axis.title.x = element_blank(), 
                                       axis.title.y = element_blank())
    if (!missing(x = gene.groups)) {
        p <- p + facet_grid(facets = ~gene.groups, scales = "free_x", 
                            space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                         units = "lines"), strip.background = element_blank(), 
                                                                    strip.placement = "outside")
    }
    if (x.lab.rot) {
        p <- p + theme(axis.text.x = element_text(angle = 90, 
                                                  vjust = 0.5))
    }
    if (!plot.legend) {
        p <- p + theme(legend.position = "none")
    }
    else if (plot.legend) {
        plot.legend <- cowplot::get_legend(plot = p)
        palettes <- list()
        for (i in seq_along(along.with = colorlist)) {
            palettes[[names(colorlist[i])]] <- colorRampPalette(colors = c("grey", 
                                                                           colorlist[[i]]))(20)
        }
        gradient.legends <- mapply(FUN = Seurat:::GetGradientLegend, 
                                   palette = palettes, group = names(x = palettes), 
                                   SIMPLIFY = FALSE, USE.NAMES = FALSE)
        p <- p + theme(legend.position = "none")
        legends <- cowplot::plot_grid(plotlist = gradient.legends, 
                                      plot.legend, ncol = 1, rel_heights = c(1, rep.int(x = 0.5, 
                                                                                        times = length(x = gradient.legends))), scale = rep(0.5, 
                                                                                                                                            length(gradient.legends)), align = "hv")
        p <- cowplot::plot_grid(p, legends, ncol = 2, rel_widths = c(1, 
                                                                     0.3), scale = c(1, 0.8))
    }
    suppressWarnings(print(p))
    if (do.return) {
        return(p)
    }
}


# add x.lab
SingleVlnPlot.1 <- function (feature, data, cell.ident, do.sort, y.max, size.x.use,
                        size.y.use, size.title.use, adjust.use, point.size.use, x.lab,
                        cols.use, gene.names, y.log, x.lab.rot, y.lab.rot, legend.position,
                        remove.legend)
{
    feature.name <- colnames(data)
    colnames(data) <- "feature"
    feature <- "feature"
    set.seed(seed = 42)
    data$ident <- cell.ident
    if (do.sort) {
        data$ident <- factor(x = data$ident,
        levels = names(x = rev(x = sort(x = tapply(X = data[,feature],
        INDEX = data$ident,
        FUN = mean)))))
    }
    if (y.log) {
        noise <- rnorm(n = length(x = data[, feature]))/200
        data[, feature] <- data[, feature] + 1
    }
    else {
        noise <- rnorm(n = length(x = data[, feature]))/1e+05
    }
    if (all(data[, feature] == data[, feature][1])) {
        warning(paste0("All cells have the same value of ",
        feature, "."))
    }
    else {
        data[, feature] <- data[, feature] + noise
    }
    y.max <- Seurat:::SetIfNull(x = y.max, default = max(data[, feature]))
    plot <- ggplot(data = data, mapping = aes(x = factor(x = ident),
    y = feature)) +
    geom_violin(scale = "width", adjust = adjust.use,
    trim = TRUE, mapping = aes(fill = factor(x = ident))) +
    guides(fill = guide_legend(title = NULL)) + xlab(x.lab) +
    Seurat:::NoGrid() + ggtitle(feature) + theme(plot.title = element_text(size = size.title.use,
    face = "bold"), legend.position = legend.position, axis.title.x = element_text(face = "bold",
    colour = "#990000", size = size.x.use), axis.title.y = element_text(face = "bold",
    colour = "#990000", size = size.y.use))
    if (point.size.use != 0) {
        plot <- plot + geom_jitter(height = 0, size = point.size.use)
    }
    plot <- plot + ggtitle(feature.name)
    if (y.log) {
        plot <- plot + scale_y_log10()
    }
    else {
        plot <- plot + ylim(min(data[, feature]), y.max)
    }
    if (feature %in% gene.names) {
        if (y.log) {
            plot <- plot + ylab(label = "Log Expression level")
        }
        else {
            plot <- plot + ylab(label = "Expression level")
        }
    }
    else {
        plot <- plot + ylab(label = "")
    }
    if (!is.null(x = cols.use)) {
        plot <- plot + scale_fill_manual(values = cols.use)
    }
    if (x.lab.rot) {
        plot <- plot + theme(axis.text.x = element_text(angle = 45,
        hjust = 1, size = size.x.use))
    }
    if (y.lab.rot) {
        plot <- plot + theme(axis.text.x = element_text(angle = 90,
        size = size.y.use))
    }
    if (remove.legend) {
        plot <- plot + theme(legend.position = "none")
    }
    return(plot)
}


# add x.lab
VlnPlot.1 <- function (object, features.plot, ident.include = NULL, nCol = NULL,
do.sort = FALSE, y.max = NULL, same.y.lims = FALSE, size.x.use = 16,
size.y.use = 16, size.title.use = 20, adjust.use = 1, point.size.use = 1,
cols.use = NULL, group.by = NULL, y.log = FALSE, x.lab.rot = FALSE,
y.lab.rot = FALSE, legend.position = "right", single.legend = TRUE,
remove.legend = FALSE, do.return = FALSE, return.plotlist = FALSE,
x.lab = "Identity",...)
{
    if (is.null(x = nCol)) {
        if (length(x = features.plot) > 9) {
            nCol <- 4
        }
        else {
            nCol <- min(length(x = features.plot), 3)
        }
    }
    data.use <- data.frame(FetchData(object = object, vars.all = features.plot,
    ...), check.names = F)
    if (is.null(x = ident.include)) {
        cells.to.include <- object@cell.names
    }
    else {
        cells.to.include <- WhichCells(object = object, ident = ident.include)
    }
    data.use <- data.use[cells.to.include, , drop = FALSE]
    if (!is.null(x = group.by)) {
        ident.use <- as.factor(x = FetchData(object = object,
        vars.all = group.by)[cells.to.include, 1])
    }
    else {
        ident.use <- object@ident[cells.to.include]
    }
    gene.names <- colnames(x = data.use)[colnames(x = data.use) %in%
    rownames(x = object@data)]
    if (single.legend) {
        remove.legend <- TRUE
    }
    if (same.y.lims && is.null(x = y.max)) {
        y.max <- max(data.use)
    }
    plots <- lapply(X = features.plot, FUN = function(x) {
        return(SingleVlnPlot.1(feature = x, data = data.use[,x, drop = FALSE], cell.ident = ident.use, do.sort = do.sort,
        y.max = y.max, size.x.use = size.x.use, size.y.use = size.y.use,
        size.title.use = size.title.use, adjust.use = adjust.use,
        point.size.use = point.size.use, cols.use = cols.use,
        gene.names = gene.names, y.log = y.log, x.lab.rot = x.lab.rot,
        y.lab.rot = y.lab.rot, legend.position = legend.position,
        remove.legend = remove.legend, x.lab = x.lab))
    })
    if (length(x = features.plot) > 1) {
        plots.combined <- plot_grid(plotlist = plots, ncol = nCol)
        if (single.legend && !remove.legend) {
            legend <- get_legend(plot = plots[[1]] + theme(legend.position = legend.position))
            if (legend.position == "bottom") {
                plots.combined <- plot_grid(plots.combined,
                legend, ncol = 1, rel_heights = c(1, 0.2))
            }
            else if (legend.position == "right") {
                plots.combined <- plot_grid(plots.combined,
                legend, rel_widths = c(3, 0.3))
            }
            else {
                warning("Shared legends must be at the bottom or right of the plot")
            }
        }
    }
    else {
        plots.combined <- plots[[1]]
    }
    if (do.return) {
        if (return.plotlist) {
            return(plots)
        }
        else {
            return(plots.combined)
        }
    }
    else {
        if (length(x = plots.combined) > 1) {
            plots.combined
        }
        else {
            invisible(x = lapply(X = plots.combined, FUN = print))
        }
    }
}

# define topGOterms funtion ===
#' topGoterm incorporate the whole topgo analysis pipeline
#' @param int.genes genes of interest, extacted based DE analysis
#' @param bg.genes background genes
#' @param organism Select Human or mouse genes
#' @param getBM biomaRt::getBM result. A table with gene name and go id.
#' If no getBM is provided, function will find it automatically.
#' @param ontology Specify ontology database, select within in c("BP","CC","MF"). Default is "BP".
#' @param nodeSize an integer larger or equal to 1. Inherited from topGOdata ojbect
#' @param topNodes an integer larger or equal to 1. Show many rows in final GenTable results.
#' @export results GenTable result.
#' @example 
#' Spermatogonia.Go <-  topGOterms(int.genes = unique(Spermatogonia_markers),
#                                  bg.genes = rownames(SSCs@data),
#'                                 organism =  "Mouse",
#'                                 getBM = getBM)
topGOterms <- function(int.genes, bg.genes,organism =  "Mouse",getBM = NULL,
                       ontology = "BP",
                       nodeSize = 5,topNodes = 20){
    library(topGO)  
    if(is.null(getBM)){
        library(biomaRt)
        Database = useMart("ensembl") %>% listDatasets() %>% as.data.frame()
        # mmusculus_gene_ensembl is at 2nd
        dataset = tail(Database[grep(organism,Database$description),"dataset"],1) 
        print(paste("Use biomaRt",dataset,"dataset"))
        # collect gene names from biomart
        mart <- biomaRt::useMart("ensembl",dataset=dataset)
        # Get ensembl gene ids and GO terms
        getBM <- biomaRt::getBM(attributes = c("external_gene_name", "go_id"),
                                filters= "external_gene_name",
                                values = bg.genes,
                                mart = mart)
    }
    getBM <- getBM[getBM$go_id != '',]
    geneID2GO <- by(getBM$go_id, getBM$external_gene_name, function(x) as.character(x))
    all.genes <- sort(unique(as.character(getBM$external_gene_name)))
    print(paste("Total",length(all.genes),"genes.",length(int.genes), "gene of interest"))
    int.genes <- factor(as.integer(all.genes %in% int.genes))
    names(int.genes) = all.genes
    
    tab = as.list(ontology)
    names(tab) = ontology
    for(ont in ontology){
        go.obj <- new("topGOdata", ontology = ont,
                      allGenes = int.genes,
                      annot = annFUN.gene2GO,
                      gene2GO = geneID2GO,
                      nodeSize = nodeSize
        )
        resultFisher <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
        resultKS <- runTest(go.obj, algorithm = "classic", statistic = "ks")
        resultKS.elim <- runTest(go.obj, algorithm = "elim", statistic = "ks")
        
        tab$ont <- GenTable(go.obj, classicFisher = resultFisher, 
                            classicKS = resultKS, elimKS = resultKS.elim,topNodes = topNodes,
                            orderBy = "elimKS", ranksOf = "classicFisher")
    }
    topGOResultsSxT <- as.data.frame(do.call(rbind,tab))
    topGOResultsSxT <- topGOResultsSxT[(length(ontology)+1):nrow(topGOResultsSxT),]
    rownames(topGOResultsSxT) = NULL
    path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    jpeg(paste0(path,organism,"_go_analysis.jpeg"), units="in", width=10, height=7,
         res=600)
    showSigOfNodes(go.obj, score(resultFisher), firstSigNodes = 5, useInfo = 'def')
    dev.off()
    return(topGOResultsSxT)
}


Types2Markers <- function(df, by = "Cell_Type"){
    "Convert a list of N gene vector into a NX2 dataframe.
        Input
        df: N by 2 dataframe with gene at column1 and Cell_Type at column2.
        by: certerial to group by, chose between Cell_Type, or cluster.
        
        Output 
        List: list of N gene vector
        ---------------------------"
    # preparation
    df <- df[,c("gene","Cell_Type","cluster")]
    df[,"Cell_Type"] <- gsub(" ","_",df[,"Cell_Type"])
    df[,"cluster"] <- gsub(" ","_",df[,"cluster"])
    df <- df[gtools::mixedorder(df[,by]),]
    
    # Create a list contains marker gene vector
    x <- unique(df[,by])
    List <- list()
    List <- lapply(seq_along(x), function(y, df, i) {
        List[[i]]=unique(df[df[,by]==y[i],"gene"])},
        y=x,df=df)
    names(List) <- x
    return(List)
}


TSNEPlot.1 <- function(object, do.label = FALSE, pt.size = 1, label.size = 4, 
                       cells.use = NULL, colors.use = NULL, alpha =1,no.legend = TRUE,
                       text.repel = FALSE, label.repel = TRUE,force= 1,do.print=F,...)
{
    return(DimPlot.1(object = object, reduction.use = "tsne", no.legend = no.legend,
                     cells.use = cells.use, pt.size = pt.size, do.label = do.label, 
                     label.size = label.size, cols.use = colors.use,alpha = alpha,
                     text.repel = text.repel, label.repel = label.repel,force= force,
                     do.print=do.print,...))
}

#' Generate 3D TSNEplot
#' @param object Seurat object after performing RunTSNE(dim.embed = 3)
#' @param ... all other parameters are inherited from Seurat::TSNEPlot()
#' @export p/p3 ggplot object from TSNEplot
TSNEPlot.3D <- function (object, reduction.use = "tsne", dim.1 = 1, dim.2 = 2, dim.3 = 3,
                         cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE, 
                         cols.use = NULL, group.by = "ident", pt.shape = NULL, do.hover = FALSE, 
                         data.hover = "ident", do.identify = FALSE, do.label = FALSE, 
                         label.size = 4, no.legend = FALSE, no.axes = FALSE, dark.theme = FALSE, 
                         plot.order = NULL, plot.title = NULL, ...) 
{
    embeddings.use = GetDimReduction(object = object, reduction.type = reduction.use, 
                                     slot = "cell.embeddings")
    if (length(x = embeddings.use) == 0) {
        stop(paste(reduction.use, "has not been run for this object yet."))
    }
    if (ncol(x = embeddings.use) < 3) {
        stop(paste(reduction.use, "doesn't have the third dimension.
                           Suggest performing RunTSNE(dim.embed = 3)"))
    }
    cells.use <- Seurat:::SetIfNull(x = cells.use, default = colnames(x = object@data))
    dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                slot = "key")
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(x = embeddings.use)
    cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
    data.plot <- data.plot[cells.use, dim.codes]
    ident.use <- as.factor(x = object@ident[cells.use])
    if (group.by != "ident") {
        ident.use <- as.factor(x = FetchData(object = object, 
                                             vars.all = group.by)[cells.use, 1])
    }
    data.plot$ident <- ident.use
    data.plot$x <- data.plot[, dim.codes[1]]
    data.plot$y <- data.plot[, dim.codes[2]]
    data.plot$z <- data.plot[, dim.codes[3]]
    data.plot$pt.size <- pt.size
    if (!is.null(plot.order)) {
        if (any(!plot.order %in% data.plot$ident)) {
            stop("invalid ident in plot.order")
        }
        plot.order <- rev(c(plot.order, setdiff(unique(data.plot$ident), 
                                                plot.order)))
        data.plot$ident <- factor(data.plot$ident, levels = plot.order)
        data.plot <- data.plot[order(data.plot$ident), ]
    }
    rgl::open3d()
    
    car::scatter3d(x = data.plot$x, y = data.plot$y, z = data.plot$z)
    
    geom_point(mapping = aes(colour = factor(x = ident)), 
               size = pt.size)
    if (!is.null(x = pt.shape)) {
        shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use, 
                                                                     1]
        if (is.numeric(shape.val)) {
            shape.val <- cut(x = shape.val, breaks = 5)
        }
        data.plot[, "pt.shape"] <- shape.val
        p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
            geom_point(mapping = aes(colour = factor(x = ident), 
                                     shape = factor(x = pt.shape)), size = pt.size)
    }
    if (!is.null(x = cols.use)) {
        p <- p + scale_colour_manual(values = cols.use)
    }
    p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) + 
        zlab(label = dim.codes[[3]]) +
        scale_size(range = c(pt.size, pt.size))
    p3 <- p2 + Seurat:::SetXAxisGG() + Seurat:::SetYAxisGG()  + 
        Seurat:::SetLegendPointsGG(x = 6) + 
        Seurat:::SetLegendTextGG(x = 12) + no.legend.title + theme_bw() + 
        Seurat:::NoGrid()
    p3 <- p3 + theme(legend.title = element_blank())
    if (!is.null(plot.title)) {
        p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
    }
    if (do.label) {
        centers <- data.plot %>% dplyr::group_by(ident) %>% 
            summarize(x = median(x = x), y = median(x = y))
        p3 <- p3 + geom_point(data = centers, mapping = aes(x = x, 
                                                            y = y), size = 0, alpha = 0) + 
            geom_text(data = centers,
                      mapping = aes(label = ident), size = label.size)
    }
    if (dark.theme) {
        p <- p + DarkTheme()
        p3 <- p3 + DarkTheme()
    }
    if (no.legend) {
        p3 <- p3 + theme(legend.position = "none")
    }
    if (no.axes) {
        p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                         axis.text.y = element_blank(), axis.ticks = element_blank(), 
                         axis.title.x = element_blank(), axis.title.y = element_blank(), 
                         panel.background = element_blank(), panel.border = element_blank(), 
                         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                         plot.background = element_blank())
    }
    if (do.identify || do.hover) {
        if (do.bare) {
            plot.use <- p
        }
        else {
            plot.use <- p3
        }
        if (do.hover) {
            if (is.null(x = data.hover)) {
                features.info <- NULL
            }
            else {
                features.info <- FetchData(object = object, 
                                           vars.all = data.hover)
            }
            return(HoverLocator(plot = plot.use, data.plot = data.plot, 
                                features.info = features.info, dark.theme = dark.theme))
        }
        else if (do.identify) {
            return(FeatureLocator(plot = plot.use, data.plot = data.plot, 
                                  dark.theme = dark.theme, ...))
        }
    }
    if (do.return) {
        if (do.bare) {
            return(p)
        }
        else {
            return(p3)
        }
    }
}


# FeaturePlot in Seurat 3 doesn't support two gene overlay.
# modify FeaturePlot in Seurat 2 and support Seurat 3 object
FeaturePlot.1 <- function (object, features.plot, min.cutoff = NA, max.cutoff = NA,
dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1,
cols.use = c("yellow",  "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE,
data.hover = "ident", do.identify = FALSE, reduction = "tsne",
use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE,
dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
{
    cells.use <- cells.use %||% colnames(x = object@data)
    if (is.null(x = nCol)) {
        nCol <- 2
        if (length(x = features.plot) == 1) {
            nCol <- 1
        }
        if (length(x = features.plot) > 6) {
            nCol <- 3
        }
        if (length(x = features.plot) > 9) {
            nCol <- 4
        }
    }
    num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) +
    1
    if (overlay | do.hover) {
        num.row <- 1
        nCol <- 1
    }
    par(mfrow = c(num.row, nCol))
    dim.code <- Embeddings(object = object[[reduction]])
    dim.codes <- paste0(dim.code, c(dim.1, dim.2))
    data.plot <- as.data.frame(GetCellEmbeddings(object = object,
    reduction.type = reduction.use, dims.use = c(dim.1,
    dim.2), cells.use = cells.use))
    x1 <- paste0(dim.code, dim.1)
    x2 <- paste0(dim.code, dim.2)
    data.plot$x <- data.plot[, x1]
    data.plot$y <- data.plot[, x2]
    data.plot$pt.size <- pt.size
    names(x = data.plot) <- c("x", "y")
    data.use <- t(x = FetchData(object = object, vars.all = features.plot,
    cells.use = cells.use, use.imputed = use.imputed))
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = min(data.use[feature,
        ]), no = cutoff)
    }, cutoff = min.cutoff, feature = features.plot)
    
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        ifelse(test = is.na(x = cutoff), yes = max(data.use[feature,
        ]), no = cutoff)
    }, cutoff = max.cutoff, feature = features.plot)
    check_lengths = unique(x = vapply(X = list(features.plot,
    min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check_lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    if (overlay) {
        pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot,
        data.plot = data.plot, pt.size = pt.size, pch.use = pch.use,
        cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff,
        max.cutoff = max.cutoff, no.axes = no.axes, no.legend = no.legend,
        dark.theme = dark.theme))
    }
    else {
        pList <- mapply(FUN = SingleFeaturePlot, feature = features.plot,
        min.cutoff = min.cutoff, max.cutoff = max.cutoff,
        MoreArgs = list(data.use = data.use, data.plot = data.plot,
        pt.size = pt.size, pch.use = pch.use, cols.use = cols.use,
        dim.codes = dim.codes, no.axes = no.axes, no.legend = no.legend,
        dark.theme = dark.theme, vector.friendly = vector.friendly),
        SIMPLIFY = FALSE)
    }
    if (do.hover) {
        if (length(x = pList) != 1) {
            stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
        }
        if (is.null(x = data.hover)) {
            features.info <- NULL
        }
        else {
            features.info <- FetchData(object = object, vars.all = data.hover)
        }
        return(HoverLocator(plot = pList[[1]], data.plot = data.plot,
        features.info = features.info, dark.theme = dark.theme,
        title = features.plot))
    }
    else if (do.identify) {
        if (length(x = pList) != 1) {
            stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
        }
        return(FeatureLocator(plot = pList[[1]], data.plot = data.plot,
        dark.theme = dark.theme))
    }
    else {
        print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
    }
    par(mfrow = c(1, 1))
    if (do.return) {
        return(pList)
    }
}



# Support FeaturePlot.1
BlendPlot <- function (data.use, features.plot, data.plot, pt.size, pch.use,
cols.use, dim.codes, min.cutoff, max.cutoff, no.axes, no.legend,
dark.theme)
{
    num.cols <- length(x = cols.use)
    cols.not.provided <- colors(distinct = TRUE)
    cols.not.provided <- cols.not.provided[!(grepl(pattern = paste(cols.use,
    collapse = "|"), x = cols.not.provided, ignore.case = TRUE))]
    if (num.cols > 4) {
        cols.use <- cols.use[c(1:4)]
    } else if ((num.cols == 2) || (num.cols == 3)) {
        blend <- BlendColors(cols.use[c(num.cols - 1, num.cols)])
        cols.use <- c(cols.use, blend)
        if (num.cols == 2) {
            cols.use <- c(sample(x = cols.not.provided, size = 1),
            cols.use)
        }
    } else if ((num.cols == 1)) {
        if (cols.use %in% rownames(x = brewer.pal.info)) {
            palette <- brewer.pal(n = 3, name = cols.use)
            cols.use <- c(palette, BlendColors(palette[c(2,
            3)]))
        }
        else {
            cols.high <- sample(x = cols.not.provided, size = 2,
            replace = FALSE)
            cols.use <- c(cols.use, cols.high, BlendColors(cols.high))
        }
    } else if (num.cols <= 0) {
        cols.use <- c("yellow", "red", "blue", BlendColors("red",
        "blue"))
    }
    names(x = cols.use) <- c("low", "high1", "high2", "highboth")
    length.check <- vapply(X = list(features.plot, min.cutoff,
    max.cutoff), FUN = function(x) {
        return(length(x = x) != 2)
    }, FUN.VALUE = logical(length = 1))
    if (any(length.check)) {
        stop("An overlayed FeaturePlot only works with two features and requires two minimum and maximum cutoffs")
    }
    min.cutoff <- c(SetQuantile(cutoff = min.cutoff[1], data = data.gene[features.plot[1],
    ]), SetQuantile(cutoff = min.cutoff[2], data = data.gene[features.plot[2],
    ]))
    max.cutoff <- c(SetQuantile(cutoff = max.cutoff[1], data = data.gene[features.plot[1],
    ]), SetQuantile(cutoff = max.cutoff[2], data = data.gene[features.plot[2],
    ]))
    data.gene <- stats::na.omit(object = data.frame(data.use[features.plot,
    ]))
    cell.names <- colnames(x = data.gene)
    data.gene <- matrix(data = vapply(X = data.gene, FUN = function(x) ifelse(test = x <
    min.cutoff, yes = min.cutoff, no = x), FUN.VALUE = c(1,
    1)), nrow = 2)
    data.gene <- matrix(data = vapply(X = as.data.frame(x = data.gene),
    FUN = function(x) ifelse(test = x > max.cutoff, yes = max.cutoff,
    no = x), FUN.VALUE = c(1, 1)), nrow = 2)
    data.gene <- as.data.frame(x = data.gene)
    rownames(x = data.gene) <- features.plot
    colnames(x = data.gene) <- cell.names
    if (all(data.gene == 0)) {
        data.cut <- 0
    } else {
        cuts <- apply(X = data.gene, MARGIN = 1, FUN = cut,
        breaks = 2, labels = FALSE)
        cuts.dim <- dim(x = cuts)
        if (cuts.dim[1] > cuts.dim[2]) {
            cuts <- t(x = cuts)
        }
        data.cut = apply(X = cuts, MARGIN = 2, FUN = function(x) {
            return(if ((x[1] == 1) && (x[2] == 2)) {
                "high2"
            } else if ((x[1] == 2) && (x[2] == 1)) {
                "high1"
            } else if ((x[1] == 2) && (x[2] == 2)) {
                "highboth"
            } else {
                "low"
            })
        })
        data.cut <- as.factor(x = data.cut)
    }
    data.plot$colors <- data.cut
    legend.names <- c(high1 = paste("High", features.plot[1]),
    high2 = paste("High", features.plot[2]), highboth = "High both")
    title <- paste0(features.plot, collapse = " + ")
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
    p <- p + geom_point(mapping = aes(color = colors), size = pt.size,
    shape = pch.use)
    p <- p + scale_color_manual(values = cols.use, limits = c("high1",
    "high2", "highboth"), labels = legend.names, guide = guide_legend(title = NULL,
    override.aes = list(size = 2)))
    if (no.axes) {
        p <- p + labs(title = title, x = "", y = "") + theme(axis.line = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank())
    } else {
        p <- p + labs(title = title, x = dim.codes[1], y = dim.codes[2])
    }
    if (no.legend) {
        p <- p + theme(legend.position = "none")
    }
    if (dark.theme) {
        p <- p + DarkTheme()
    }
    return(p)
}


SetQuantile <- function (cutoff, data) {
    if (grepl(pattern = "^q[0-9]{1,2}$", x = as.character(x = cutoff),
    perl = TRUE)) {
        this.quantile <- as.numeric(x = sub(pattern = "q", replacement = "",
        x = as.character(x = cutoff)))/100
        data <- unlist(x = data)
        data <- data[data > 0]
        cutoff <- quantile(x = data, probs = this.quantile)
    }
    return(as.numeric(x = cutoff))
}

