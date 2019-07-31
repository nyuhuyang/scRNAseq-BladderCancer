library(dplyr)
library(magrittr)
library(Seurat)
library(ggplot2)

# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
        if (!is.null(x = lhs)) {
                return(lhs)
        } else {
                return(rhs)
        }
}

# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
        return(rhs)
    } else {
        return(lhs)
    }
}

# Generic GetAssayData functions will generate Error in t.default(x = GetAssayData(object = object, slot = slot)[features,  : 
# argument is not a matrix
GetAssayData <- function(object, slot = "data", assay = NULL){
    assay <- assay %||% DefaultAssay(object = object)
    assay.data <- GetAssay(object = object, assay = assay)
    mtrix <- slot(object = assay.data, name = slot)
    if(class(mtrix) =="dgCMatrix") mtrix = as.matrix(mtrix)
    return(mtrix)
}


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
                     row.names = colnames(object))
    colnames(mat) = get("label")
    newMetaData = .AddMetaColor(mat = mat, colors = colors)
    object <- AddMetaData(object = object, metadata = newMetaData)
    
    return(object)
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
        if(is.null(gene.alias)) return(NULL)
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
        if(class(subject)=="Seurat"){
                subject = rownames(subject)[1]
        }
        # Human gene
        if(subject == toupper(subject))
                return("Human")
        # Mouse gene
        if(subject == Hmisc::capitalize(tolower(subject)))
                return("Mouse") 

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


# Basic function to convert human to mouse gene names
#' @example genes <- Human2Mouse(humGenes)
Human2Mouse <- function(x){
    
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                     values = x , mart = human, attributesL = c("mgi_symbol"), 
                     martL = mouse, uniqueRows=T)
    
    humanx <- unique(genesV2[, 2])
    
    # Print the first 6 genes found to the screen
    print(head(humanx))
    return(humanx)
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


#' DoHeatmap.1, automatically group top DE genes from FindMakers output
#example   DoHeatmap.1(SSCs,top,Top_n = 15, 
#                   group.order = major_cells,ident.use = "all cell types",
#                   group.label.rot = T,cex.row = 5,remove.key =T)
#
DoHeatmap.1 <- function(object, marker_df,features = NULL, no.legend =F,unique.name= F,
                        Top_n = 10, cells = NULL, group.by = "ident", 
                         group.bar = TRUE, disp.min = -2.5, disp.max = NULL, slot = "scale.data", 
                         assay = NULL, label = TRUE, cols=NULL, size = 5.5, hjust = 0, angle = 45, 
                         raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02, 
                         combine = TRUE,title = NULL,title.size = 14,do.print = FALSE,
                        pal_gsea = TRUE,position = "right",
                        cex.row=12,legend.size = NULL,units="in", width=10, height=7,res=600,...){
    if(unique.name) {
        v <- paste(unique(object$orig.ident),collapse = "_")
    } else v <- deparse(substitute(object))
    v = paste0(v,"_",FindIdentLabel(object))
    if(!no.legend) v = paste0(v, "_Legend")
    if (!missing(x = marker_df)) {
        
        colnames(marker_df)[grep("cluster*.",colnames(marker_df))[1]] = "cluster"
        top <-  marker_df %>% 
            group_by(cluster) %>% 
            top_n(Top_n, avg_logFC)
        features = c(as.character(top$gene),features)
        }
    heatmap <- DoHeatmap(object = object, features = features, cells = cells, group.by = group.by,
                         group.bar = group.bar, disp.min = disp.min, disp.max = disp.max,
                         slot = slot, assay = assay, label = label, size = size,
                         hjust = hjust, angle = angle, raster = raster, draw.lines = draw.lines,
                         lines.width = lines.width, group.bar.height = group.bar.height, 
                         combine = combine)+
        scale_y_discrete(position = position)
    if(pal_gsea) heatmap = heatmap + scale_fill_gradientn(colors = ggsci::pal_gsea()(12))
    if(!is.null(title)) {
        heatmap = heatmap+ ggtitle(title)+ 
            theme(plot.title = element_text(size=title.size, hjust = 0.5,face="plain"))
        }
    heatmap = heatmap + theme(axis.text.y = element_text(size = cex.row))
    if(!is.null(legend.size)) {
        heatmap = heatmap + theme(legend.text = element_text(size = legend.size),
                                  legend.title = element_text(size = legend.size*2),
                                  legend.key.size = unit(legend.size/10,"cm"))
    }
    if(no.legend) heatmap = heatmap + NoLegend()
    if(do.print){
            path <- paste0("output/",gsub("-","",Sys.Date()),"/")
            if(!dir.exists(path)) dir.create(path, recursive = T)
            jpeg(paste0(path,"Doheatmap_top",Top_n,"_",v,".jpeg"),
                 units=units, width=width, height=height,res=res)
            print(heatmap)
            dev.off()
    } else return(heatmap)
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
    
    return(base::as.character(df_colors[,color_index]))
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

# Support the manual color scheme in blend figure
TSNEPlot.1 <- function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
    c("lightgrey", "#ff0000", "#00ff00") } else {
    c("lightgrey", "blue") }, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, 
    reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data", 
    blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4, 
    repel = FALSE, ncol = NULL, combine = TRUE, coord.fixed = FALSE, 
    by.col = TRUE, sort.cell = FALSE)
{
    no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                      axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
                                                                                             size = 14, 
                                                                                             margin = margin(r = 7)))
    reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
    if (length(x = dims) != 2 || !is.numeric(x = dims)) {
        stop("'dims' must be a two-length integer vector")
    }
    if (blend && length(x = features) != 2) {
        stop("Blending feature plots only works with two features")
    }
    if (blend) {
        default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
        cols <- switch(EXPR = as.character(x = length(x = cols)), 
                       `0` = {
                           warning("No colors provided, using default colors", 
                                   call. = FALSE, immediate. = TRUE)
                           default.colors
                       }, `1` = {
                           warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                                   call. = FALSE, immediate. = TRUE)
                           c(cols, default.colors[2:3])
                       }, `2` = {
                           warning("Only two colors provided, assuming specified are for features and agumenting with '", 
                                   default.colors[1], "' for double-negatives", 
                                   call. = FALSE, immediate. = TRUE)
                           c(default.colors[1], cols)
                       }, `3` = cols,
                       `4` = cols,
                       {
                           warning("More than four colors provided, using only first four", 
                                   call. = FALSE, immediate. = TRUE)
                           cols[1:4]
                       })
    }
    if (blend && length(x = cols) < 3) {
        stop("Blending feature plots only works with three colors; first one for negative cells")
    }
    dims <- paste0(Key(object = object[[reduction]]), dims)
    cells <- cells %||% colnames(x = object)
    data <- FetchData(object = object, vars = c(dims, "ident", 
                                                features), cells = cells, slot = slot)
    if (ncol(x = data) < 4) {
        stop("None of the requested features were found: ", 
             paste(features, collapse = ", "), " in slot ", slot, 
             call. = FALSE)
    } else if (!all(dims %in% colnames(x = data))) {
        stop("The dimensions requested were not found", call. = FALSE)
    }
    features <- colnames(x = data)[4:ncol(x = data)]
    min.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = min(data[, feature]), no = cutoff))
    }, cutoff = min.cutoff, feature = features)
    max.cutoff <- mapply(FUN = function(cutoff, feature) {
        return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                               feature]), no = cutoff))
    }, cutoff = max.cutoff, feature = features)
    check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                                max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
    if (length(x = check.lengths) != 1) {
        stop("There must be the same number of minimum and maximum cuttoffs as there are features")
    }
    brewer.gran <- ifelse(test = length(x = cols) == 1, 
                          yes = brewer.pal.info[cols, ]$maxcolors, 
                          no = length(x = cols))
    data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
                                       FUN = function(index) {
                                           data.feature <- as.vector(x = data[, index])
                                           min.use <- Seurat:::SetQuantile(cutoff = min.cutoff[index - 
                                                                                          3], data.feature)
                                           max.use <- Seurat:::SetQuantile(cutoff = max.cutoff[index - 
                                                                                          3], data.feature)
                                           data.feature[data.feature < min.use] <- min.use
                                           data.feature[data.feature > max.use] <- max.use
                                           if (brewer.gran == 2) {
                                               return(data.feature)
                                           }
                                           data.cut <- if (all(data.feature == 0)) {
                                               0
                                           }
                                           else {
                                               as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                                                                breaks = brewer.gran)))
                                           }
                                           return(data.cut)
                                       })
    colnames(x = data)[4:ncol(x = data)] <- features
    rownames(x = data) <- cells
    data$split <- if (is.null(x = split.by)) {
        Seurat:::RandomName()
    } else {
        switch(EXPR = split.by, ident = Idents(object = object)[cells], 
               object[[split.by, drop = TRUE]][cells])
    }
    if (!is.factor(x = data$split)) {
        data$split <- factor(x = data$split)
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    plots <- vector(mode = "list", length = ifelse(test = blend, 
                                                   yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
    xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
    ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
    if (blend) {
        ncol <- 4
        color.matrix <- Seurat:::BlendMatrix(two.colors = cols[2:3], 
                                    col.threshold = blend.threshold, negative.color = cols[1])
        if(length(cols)>3) col.both <- cols[4]
        cols <- cols[2:3]
        colors <- list(color.matrix[, 1], color.matrix[1, ], 
                       as.vector(x = color.matrix))
    }
    for (i in 1:length(x = levels(x = data$split))) {
        ident <- levels(x = data$split)[i]
        data.plot <- data[as.character(x = data$split) == ident, 
                          , drop = FALSE]
        if (blend) {
            features <- features[1:2]
            no.expression <- features[colMeans(x = data.plot[, features]) == 0]
            if (length(x = no.expression) != 0) {
                stop("The following features have no value: ", 
                     paste(no.expression, collapse = ", "), call. = FALSE)
            }
            data.plot <- cbind(data.plot[, c(dims, "ident")], 
                               Seurat:::BlendExpression(data = data.plot[, features[1:2]]))
            features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
        }
        for (j in 1:length(x = features)) {
            feature <- features[j]
            if (blend) {
                cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
                cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
            } else {
                cols.use <- NULL
            }
            data.single <- data.plot[, c(dims, "ident", feature, 
                                         shape.by)]
            if (sort.cell) {
                data.single <- data.single[order(data.single[, 
                                                             feature]), ]
            }
            plot <- Seurat:::SingleDimPlot(data = data.single, dims = dims, 
                                  col.by = feature, order = order, pt.size = pt.size, 
                                  cols = cols.use, shape.by = shape.by, label = FALSE) + 
                scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + 
                theme_cowplot()
            if (label) {
                plot <- LabelClusters(plot = plot, id = "ident", 
                                      repel = repel, size = label.size)
            }
            if (length(x = levels(x = data$split)) > 1) {
                plot <- plot + theme(panel.border = element_rect(fill = NA, 
                                                                 colour = "black"))
                plot <- plot + if (i == 1) {
                    labs(title = feature)
                }
                else {
                    labs(title = NULL)
                }
                if (j == length(x = features) && !blend) {
                    suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident)) + 
                                         no.right)
                }
                if (j != 1) {
                    plot <- plot + theme(axis.line.y = element_blank(), 
                                         axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                                         axis.title.y.left = element_blank())
                }
                if (i != length(x = levels(x = data$split))) {
                    plot <- plot + theme(axis.line.x = element_blank(), 
                                         axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                                         axis.title.x = element_blank())
                }
            }
            else {
                plot <- plot + labs(title = feature)
            }
            if (!blend) {
                plot <- plot + guides(color = NULL)
                cols.grad <- cols
                if (length(x = cols) == 1) {
                    plot <- plot + scale_color_brewer(palette = cols)
                }
                else if (length(x = cols) > 1) {
                    unique.feature.exp <- unique(data.plot[, feature])
                    if (length(unique.feature.exp) == 1) {
                        warning("All cells have the same value (", 
                                unique.feature.exp, ") of ", feature, 
                                ".")
                        if (unique.feature.exp == 0) {
                            cols.grad <- cols[1]
                        }
                        else {
                            cols.grad <- cols
                        }
                    }
                    plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                                                                                 guide = "colorbar"))
                }
            }
            if (coord.fixed) {
                plot <- plot + coord_fixed()
            }
            plot <- plot
            plots[[(length(x = features) * (i - 1)) + j]] <- plot
        }
    }
    if (blend) {
        blend.legend <- Seurat:::BlendMap(color.matrix = color.matrix)
        for (ii in 1:length(x = levels(x = data$split))) {
            suppressMessages(expr = plots <- append(x = plots, 
                values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 1, 
                                                                                                  yes = levels(x = data$split)[ii], 
                                                                                                  no = "")), 
                                                                expand = c(0, 0)) + labs(x = features[1],y = features[2],
                                                                                         title = if (ii == 1) { paste("Color threshold:", blend.threshold) 
                                                                                             } else {NULL}) + 
                                  no.right), after = 4 * ii - 1))
        }
    }
    plots <- Filter(f = Negate(f = is.null), x = plots)
    if (combine) {
        if (is.null(x = ncol)) {
            ncol <- 2
            if (length(x = features) == 1) {
                ncol <- 1
            }
            if (length(x = features) > 6) {
                ncol <- 3
            }
            if (length(x = features) > 9) {
                ncol <- 4
            }
        }
        ncol <- ifelse(test = is.null(x = split.by) || blend, 
                       yes = ncol, no = length(x = features))
        legend <- if (blend) {
            "none"
        }
        else {
            split.by %iff% "none"
        }
        if (by.col && !is.null(x = split.by) && !blend) {
            plots <- lapply(X = plots, FUN = function(x) {
                return(suppressMessages(expr = x + theme_cowplot() + 
                                            ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = "")) + 
                                            no.right))
            })
            nsplits <- length(x = levels(x = data$split))
            idx <- 1
            for (i in (length(x = features) * (nsplits - 1) + 
                       1):(length(x = features) * nsplits)) {
                plots[[i]] <- suppressMessages(plots[[i]] + 
                                                   scale_y_continuous(sec.axis = dup_axis(name = features[[idx]])) + 
                                                   no.right)
                idx <- idx + 1
            }
            idx <- 1
            for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                            1)) {
                plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]])
                idx <- idx + 1
            }
            idx <- 1
            if (length(x = features) == 1) {
                for (i in 1:length(x = plots)) {
                    plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]])
                    idx <- idx + 1
                }
            }
            plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                                                                f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
            plots <- CombinePlots(plots = plots, ncol = nsplits, 
                                  legend = legend)
        } else {
            plots <- CombinePlots(plots = plots, ncol = ncol, 
                                  legend = legend, nrow = split.by %iff% length(x = levels(x = data$split)))
        }
    }
    return(plots)
}


FgseaBarplot <- function(pathways=hallmark, stats=res, nperm=1000,cluster = 1,
                         sample="",pathway.name = "Hallmark", hjust=0.5, 
                         width=10, height = 7, no.legend = FALSE,
                         cut.off = c("pval","padj"),cut.off.value = 0.25,
                         do.print = TRUE, do.return = FALSE){
    
    res = stats[order(stats["avg_logFC"]),]
    res1 = res[res$cluster == cluster,c("gene","avg_logFC")] %>% deframe
    
    fgseaRes <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
    print(dim(fgseaRes))
    
    if(cut.off == "pval"){
        (topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=25), pathway])
        (topPathwaysDown <- fgseaRes[NES < 0][head(order(pval), n=25), pathway])
    }
    if(cut.off == "padj"){
        (topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=25), pathway])
        (topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=25), pathway])
    }
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    out_write=fgseaRes[match(topPathways,pathway)]
    colnames(out_write)[1] <- 'Pathways'
    out_write$To_Plot_value <- -log10(out_write$pval)
    out_write$sign <- ifelse(out_write$NES >0,1,-1)
    out_write$To_Plot_value <- out_write$To_Plot_value*out_write$sign
    out_write$sign <- ifelse(out_write$NES >0,"Upregulated", "Downregulated")
    
    legend_title = ifelse(cut.off == "padj",yes = "Adjusted p value", no = "P value")
    
    p<-ggbarplot(out_write,
                 x = "Pathways",
                 y = "NES",
                 #fill = "sign",           # change fill color by mpg_level
                 color = "white",            # Set bar border colors to white
                 palette = "jco",            # jco journal color palett. see ?ggpar
                 sort.val = "asc",          # Sort the value in descending order
                 sort.by.groups = FALSE,     # Don't sort inside each group
                 ylab = 'Normalized Enrichment Score',
                 legend.title = paste(legend_title,"<",cut.off.value),
                 rotate = TRUE,
                 title = paste(pathway.name,"pathways in",sample,cluster),
                 ggtheme = theme_minimal(base_size = 15))+
        guides(fill = guide_legend(reverse = TRUE))+
        theme(text = element_text(size=12),
              plot.title = element_text(size = 16,hjust = hjust))
    if(cut.off == "pval") p = p + geom_col(aes(fill = pval < cut.off.value))
    if(cut.off == "padj") p = p + geom_col(aes(fill = padj < cut.off.value))
    if(no.legend) p = p + NoLegend()
    path <- paste0("output/",gsub("-","",Sys.Date()),"/")
    if(do.print){
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,sample,"_",cluster,"-",pathway.name,".jpeg"), units="in", width=width, height=height,res=600)
        print(p)
        dev.off()
    }
    if(do.return) return(p)

}

#' FgseaDotPlot generate Dot plot using findmarker results based on FGSEA
#' @param stats findmarker results
#' @param pathways pathway list
#' @param sample add to title names
#' @param order.by c(1,"pval") means order by pval in cluster 1
#' @param order specify order of y axis
#' @param nperm fgsea param, default 1000
#' @param ... ggballoonplot param
#' @example FgseaDotPlot(stats=res, pathways=hallmark,sample = "each B_MCL clusters")
FgseaDotPlot <- function(stats=results, pathways=hallmark, x = "cluster", y = "pathway",
                         nperm=1000,size = "-log10(pval)", order = NULL,
                         fill = "NES",sample="",order.by = c(1,"pval"),decreasing = T,
                         pathway.name = "Hallmark", padj = 0.25, pval=0.05,verbose=T,...){
    
    #clusters = ifelse(class(stats$cluster)=="factor", 
    #                  as.character(levels(stats$cluster)),
    #                  as.character(unique(stats$cluster)))
    clusters = unique(as.character(stats$cluster))
    fgseaRes <- list()
    for(i in 1:length(clusters)){
        res1 = stats[stats$cluster == clusters[i],]
        res1 = res1[order(res1["p_val_adj"]),c("gene","avg_logFC")]  %>% deframe
        fgseaRes[[i]] <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
        fgseaRes[[i]] = as.data.frame(fgseaRes[[i]])
        fgseaRes[[i]] = fgseaRes[[i]][,c("pathway","pval","padj","NES")]
        if(i == order.by[1]) {
            order = fgseaRes[[i]][order(fgseaRes[[i]][,order.by[2]],
                                      decreasing = decreasing),y]
            
        }
        if(!is.null(pval)) fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$pval < pval,]
        if(!is.null(padj)) fgseaRes[[i]] = fgseaRes[[i]][fgseaRes[[i]]$padj < padj,]
        fgseaRes[[i]]$cluster = clusters[i]
    }
    df_fgseaRes <- data.table::rbindlist(fgseaRes) %>% as.data.frame()
    df_fgseaRes[,"-log10(pval)"] = -log10(df_fgseaRes$pval)
    df_fgseaRes[,"-log10(padj)"] = -log10(df_fgseaRes$padj)
    if(!is.null(order)) df_fgseaRes[,y] = factor(df_fgseaRes[,y],order)
    if(verbose) print(round(dim(df_fgseaRes)/length(clusters)))
    plot<- ggballoonplot(df_fgseaRes, x = x, y = y,
                         size = size, fill = fill,
                         size.range = c(1, 5),
                         font.xtickslab=14, 
                         font.ytickslab=round(50/dim(df_fgseaRes)[1]*14),
                         title = paste(pathway.name,"pathways enriched in",sample),
                         legend.title = ifelse(fill =="NES",
                                               "Normalized\nenrichment\nscore",
                                               NULL),
                         xlab = "", ylab = "",
                         font.x =14,font.y= 14,font.main=16,...) +
        scale_fill_gradientn(colors = pal_gsea()(10))+
        theme(plot.title = element_text(hjust = 0.5))
    if(size == "padj") plot = plot + scale_size(breaks=c(0,0.05,0.10,0.15,0.2,0.25),
                                               labels=rev(c(0,0.05,0.10,0.15,0.2,0.25)))
    path <- paste0("output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    jpeg(paste0(path,"Dotplot_",sample,"_",pathway.name,
                "_",padj,"_",pval,".jpeg"), units="in", width=10, height=7,res=600)
    print(plot)
    dev.off()
}


#' Combine FindAllMarkers and calculate average UMI
#' Modified Seurat::FindAllMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindAllMarkers
#' #' @export gde.all data frame
#' @example FindAllMarkers.UMI(object)
FindAllMarkers.UMI <- function (object, assay = NULL, features = NULL, logfc.threshold = 0.25, 
                              test.use = "wilcox", slot = "data", min.pct = 0.1, min.diff.pct = -Inf, 
                              node = NULL, verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf, 
                              random.seed = 1, latent.vars = NULL, min.cells.feature = 3, 
                              min.cells.group = 3, pseudocount.use = 1, return.thresh = 0.01, 
                              ...) 
    {
        MapVals <- function(vec, from, to) {
            vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
            vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
            return(unname(obj = vec2))
        }
        if ((test.use == "roc") && (return.thresh == 0.01)) {
            return.thresh <- 0.7
        }
        if (is.null(x = node)) {
            idents.all <- sort(x = unique(x = Idents(object = object)))
        }
        else {
            tree <- Tool(object = object, slot = "BuildClusterTree")
            if (is.null(x = tree)) {
                stop("Please run 'BuildClusterTree' before finding markers on nodes")
            }
            descendants <- DFT(tree = tree, node = node, include.children = TRUE)
            all.children <- sort(x = tree$edge[, 2][!tree$edge[, 
                                                               2] %in% tree$edge[, 1]])
            descendants <- MapVals(vec = descendants, from = all.children, 
                                   to = tree$tip.label)
            drop.children <- setdiff(x = tree$tip.label, y = descendants)
            keep.children <- setdiff(x = tree$tip.label, y = drop.children)
            orig.nodes <- c(node, as.numeric(x = setdiff(x = descendants, 
                                                         y = keep.children)))
            tree <- drop.tip(phy = tree, tip = drop.children)
            new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
            idents.all <- (tree$Nnode + 2):max(tree$edge)
        }
        genes.de <- list()
        messages <- list()
        for (i in 1:length(x = idents.all)) {
            if (verbose) {
                message("Calculating cluster ", idents.all[i])
            }
            genes.de[[i]] <- tryCatch(expr = {
                FindMarkers.UMI(object = object, assay = assay, ident.1 = if (is.null(x = node)) {
                    idents.all[i]
                }
                else {
                    tree
                }, ident.2 = if (is.null(x = node)) {
                    NULL
                }
                else {
                    idents.all[i]
                }, features = features, logfc.threshold = logfc.threshold, 
                test.use = test.use, slot = slot, min.pct = min.pct, 
                min.diff.pct = min.diff.pct, verbose = verbose, 
                only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
                random.seed = random.seed, latent.vars = latent.vars, 
                min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
                pseudocount.use = pseudocount.use, ...)
            }, error = function(cond) {
                return(cond$message)
            })
            if (class(x = genes.de[[i]]) == "character") {
                messages[[i]] <- genes.de[[i]]
                genes.de[[i]] <- NULL
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
                else if (is.null(x = node) || test.use %in% c("bimod", 
                                                              "t")) {
                    gde <- gde[order(gde$p_val, -gde[, 2]), ]
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
        if ((only.pos) && nrow(x = gde.all) > 0) {
            return(subset(x = gde.all, subset = gde.all[, 2] > 0))
        }
        rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
        if (nrow(x = gde.all) == 0) {
            warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
        }
        if (length(x = messages) > 0) {
            warning("The following tests were not performed: ", 
                    call. = FALSE, immediate. = TRUE)
            for (i in 1:length(x = messages)) {
                if (!is.null(x = messages[[i]])) {
                    warning("When testing ", idents.all[i], " versus all:\n\t", 
                            messages[[i]], call. = FALSE, immediate. = TRUE)
                }
            }
        }
        if (!is.null(x = node)) {
            gde.all$cluster <- MapVals(vec = gde.all$cluster, from = new.nodes, 
                                       to = orig.nodes)
        }
        return(gde.all)
    }


#' FindIdentLabel: Find identical label between ident and metadata
#' @object seurat object
#' @label colname in metadata
FindIdentLabel <- function(object){
    ident.label <- as.vector(Idents(object))
    labels <- sapply(object@meta.data,
                     function(x) all(ident.label == x)) %>% .[.] %>% .[!is.na(.)]
    label <- names(labels[labels])
    label =  label[!(label %in% c("seurat_clusters","ident"))][1]
    return(label)
}


#' Seurat 3
#' Calculate average UMI and attach to FindMarkers results 
#' Modified Seurat::FindMarkers function to add average UMI for group1 (UMI.1) and group 2 (UMI.2)
#' @param ... all paramethers are the same as Seurat::FindMarkers
#' @export gde.all data frame
#' @example FindMarkers.UMI(object,ident.1 = "1")
FindMarkers.UMI <- function (object, ident.1 = NULL, ident.2 = NULL, group.by = NULL, 
                              subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL, 
                              features = NULL, logfc.threshold = 0.25, test.use = "wilcox", 
                              min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE, 
                              max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL, 
                              min.cells.feature = 3, min.cells.group = 3, pseudocount.use = 1, 
                              ...){
    if (!is.null(x = group.by)) {
        if (!is.null(x = subset.ident)) {
            object <- subset(x = object, idents = subset.ident)
        }
        Idents(object = object) <- group.by
    }
    if (!is.null(x = assay) && !is.null(x = reduction)) {
        stop("Please only specify either assay or reduction.")
    }
    data.slot <- ifelse(test = test.use %in% c("negbinom", "poisson", 
                                               "DESeq2"), yes = "counts", no = slot)
    if (is.null(x = reduction)) {
        assay <- assay %||% DefaultAssay(object = object)
        data.use <- GetAssayData(object = object, slot = data.slot,
                                 assay = assay)
    }
    else {
        if (data.slot == "counts") {
            stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
        }
        data.use <- t(x = Embeddings(object = object, reduction = reduction))
    }
    if (is.null(x = ident.1)) {
        stop("Please provide ident.1")
    }
    else if ((length(x = ident.1) == 1 && ident.1[1] == "clustertree") || 
             is(object = ident.1, class2 = "phylo")) {
        if (is.null(x = ident.2)) {
            stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
        }
        tree <- if (is(object = ident.1, class2 = "phylo")) {
            ident.1
        }
        else {
            Tool(object = object, slot = "BuildClusterTree")
        }
        if (is.null(x = tree)) {
            stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
        }
        ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, 
                                                     node = ident.2)]
        ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, 
                                                      node = ident.2)]
    }
    if (length(x = as.vector(x = ident.1)) > 1 && any(as.character(x = ident.1) %in% 
                                                      colnames(x = data.use))) {
        bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% 
                                                      colnames(x = data.use))]
        if (length(x = bad.cells) > 0) {
            stop(paste0("The following cell names provided to ident.1 are not present in the object: ", 
                        paste(bad.cells, collapse = ", ")))
        }
    }
    else {
        ident.1 <- WhichCells(object = object, idents = ident.1)
    }
    if (length(x = as.vector(x = ident.2)) > 1 && any(as.character(x = ident.2) %in% 
                                                      colnames(x = data.use))) {
        bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% 
                                                      colnames(x = data.use))]
        if (length(x = bad.cells) > 0) {
            stop(paste0("The following cell names provided to ident.2 are not present in the object: ", 
                        paste(bad.cells, collapse = ", ")))
        }
    }
    else {
        if (is.null(x = ident.2)) {
            ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
        }
        else {
            ident.2 <- WhichCells(object = object, idents = ident.2)
        }
    }
    if (!is.null(x = latent.vars)) {
        latent.vars <- FetchData(object = object, vars = latent.vars, 
                                 cells = c(ident.1, ident.2))
    }
    counts <- switch(EXPR = data.slot, scale.data = GetAssayData(object = object, 
                                                                 slot = "counts", assay = assay), numeric())
    de.results <- FindMarkers(object = data.use, slot = data.slot, 
                              counts = counts, cells.1 = ident.1, cells.2 = ident.2, 
                              features = features, reduction = reduction, logfc.threshold = logfc.threshold/log2(exp(1)), 
                              test.use = test.use, min.pct = min.pct, min.diff.pct = min.diff.pct, 
                              verbose = verbose, only.pos = only.pos, max.cells.per.ident = max.cells.per.ident, 
                              random.seed = random.seed, latent.vars = latent.vars, 
                              min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
                              pseudocount.use = pseudocount.use, ...)
    
    de.results$avg_logFC = log2(exp(1)) * de.results$avg_logFC
    avg_UMI.1 <- Matrix::rowMeans(expm1(x = data.use[, ident.1]))
    avg_UMI.2 <- Matrix::rowMeans(expm1(x = data.use[, ident.2]))
    avg_UMI <-data.frame(avg_UMI.1, avg_UMI.2)
    de.results <- cbind(de.results,avg_UMI[match(rownames(de.results),rownames(avg_UMI)),])
    
    return(de.results)
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
                            max.cells.per.ident = Inf, random.seed = 1, latent.vars = "nCount_RNA",
                            min.cells.gene = 3,min.cells.group=3, pseudocount.use = 1, 
                            assay.type = "RNA",save.path = NULL,save.files = TRUE,...){
    #prepare save folder name
    if(save.files){
        if(class(ident.1) == "numeric" & class(ident.2) == "numeric") {
            ident1="";ident2=""
        } else {
            ident1 <- unique(gsub('\\_.*', '', ident.1))
            ident2 <- unique(gsub('\\_.*', '', ident.2))
        }
        if(is.null(save.path)){
            path <- paste0("output/",gsub("-","",Sys.Date()),"/")
            save.path <- paste0(path,ident1,"_vs_",ident2,"/")
        }
    }
    gde <- list()
    for(i in 1:length(ident.1)) {
        ident.1vs2 <- paste(ident.1[i], ident.2[i], sep = ".vs.")
        print(ident.1vs2)
        ifelse(class(ident.1[i])=="list", ident1 <- ident.1[[i]], ident1 <- ident.1[i])
        ifelse(class(ident.2[i])=="list", ident2 <- ident.2[[i]], ident2 <- ident.2[i])
        gde[[i]] <- FindMarkers.UMI(object = object, 
                                    ident.1 = ident1,
                                    ident.2 = ident2, 
                                    assay.type = assay.type, 
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
            write.csv( gde[[i]], paste0(save.path,ident1,"_vs_",ident2,".csv"))
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
#' @example FilterGenes(object, c("cdh5,PECAM1,flt1,Vwf,Plvap,Kdr"))
FilterGenes <- function(object, marker.genes, unique= TRUE){
        if(missing(object)) 
                stop("A seurat object must be provided first")
        if(class(object) != "Seurat")
                stop("A seurat object must be provided first")
        if(missing(marker.genes)) 
                stop("A list of marker genes must be provided")

                marker.genes <- as.character(marker.genes)
        marker.genes <- unlist(strsplit(marker.genes,","))
        marker.genes <- gsub(" ","",marker.genes)
        species = CheckSpecies(object)
        if(species == "Human") marker.genes <- toupper(marker.genes)
        if(species == "Mouse") marker.genes <- Hmisc::capitalize(tolower(marker.genes)) 

        print(paste("Before filtration:",length(marker.genes)))
        marker.genes <- CaseMatch(search = marker.genes, match = rownames(object))
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
#' @param features extra genes to add beyound FindAllMarkers results 
#' @param remove.legend TRUE/FALSE
#' @param color color scheme
#' @export g vertical ggplot
#' @example MakeCorlorBar(EC, top, Top_n=40)
MakeCorlorBar <- function(object, marker_df, Top_n = NULL, features = NULL, color =NULL,unique.name = F,
                          no.legend =F, legend.size = 10,do.print = TRUE,do.return=FALSE,width=10, height=7){
    if(unique.name) {
        v <- paste(unique(object$orig.ident),collapse = "_")
    } else v <- deparse(substitute(object))
    v = paste0(v,"_",FindIdentLabel(object))
    if(!no.legend) v = paste0(v, "_Legend")
    
    if(!is.null(features)){
        marker_bar <- data.frame("gene" = features,
                                 "cluster" = "marker.genes")
        rownames(marker_bar) = marker_bar$gene 
    } else marker_bar <- NULL
    
    if(!missing(marker_df)){
        colnames(marker_df)[grep("cluster",colnames(marker_df))] ="cluster"
        marker_df %<>% group_by(cluster)
        if(!is.null(Top_n)) marker_df %<>% top_n(Top_n, avg_logFC)
        marker_df = rbind(marker_df[,c("gene","cluster")],
                                     marker_bar)
    } else marker_df <- marker_bar

    gene.use = rownames(object@assays$RNA@data)
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
        
    if(!no.legend) g = g + theme(legend.title = element_text(size=legend.size*1.2),
                                 legend.text = element_text(size=legend.size),
                                 legend.key.size = unit(legend.size/5,"line"))
    if(no.legend) g = g + theme(legend.position="none")
    if(!is.null(color)) {
        colors_fill = color_scheme(color = color,
                                   names = unique(marker_df$cluster))
        g = g + scale_fill_manual(values = colors_fill)
    }
    if(do.print) {
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"Doheatmap_vcolorbar_",v,".jpeg"),
             units="in", width=width, height=height,res=600)
        print(g)
        dev.off()
    } else if(do.return) {
        return(g)
    }
}


# Support function for TSNEPlot.1, modified from Seurate:::LabelClusters function
LabelRepel <- function(plot, id, clusters = NULL, labels = NULL, split.by = NULL, 
                       repel = TRUE,color = NULL,alpha= NULL, ...){
    xynames <- unlist(x = Seurat:::GetXYAesthetics(plot = plot), use.names = TRUE)
    if (!id %in% colnames(x = plot$data)) {
        stop("Cannot find variable ", id, " in plotting data")
    }
    if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
        warning("Cannot find splitting variable ", id, " in plotting data")
        split.by <- NULL
    }
    data <- plot$data[, c(xynames, id, split.by)]
    g <- ggplot_build(plot)
    data_color <- g$data[[1]]
    data = cbind(data, color = data_color$colour[match(data$tSNE_1, data_color$x)])
    possible.clusters <- as.character(x = na.omit(object = unique(x = data[,id])))
    groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[,id])))
    if (any(!groups %in% possible.clusters)) {
        stop("The following clusters were not found: ", 
             paste(groups[!groups %in% 
                              possible.clusters], collapse = ","))
    }
    labels.loc <- lapply(X = groups, FUN = function(group) {
        data.use <- data[data[, id] == group, , drop = FALSE]
        data.medians <- if (!is.null(x = split.by)) {
            do.call(what = "rbind", 
                    args = lapply(X = unique(x = data.use[, split.by]), 
                                  FUN = function(split) {
                                      medians <- apply(X = data.use[data.use[, split.by] == 
                                                                        split, xynames, drop = FALSE], 
                                                       MARGIN = 2, 
                                                       FUN = median, na.rm = TRUE)
                                      medians <- as.data.frame(x = t(x = medians))
                                      medians[, split.by] <- split
                                      return(medians)
                                      }))
        } else {
            as.data.frame(x = t(x = apply(X = data.use[, xynames, 
                                                       drop = FALSE], MARGIN = 2, 
                                          FUN = median, na.rm = TRUE)))
        }
        data.medians[, id] <- group
        return(data.medians)
    })
    labels.loc <- do.call(what = "rbind", args = labels.loc)
    labels <- labels %||% groups
    if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
        stop("Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (", 
             nrow(x = labels.loc), ").")
    }
    names(x = labels) <- groups
    for (group in groups) {
        labels.loc[labels.loc[, id] == group, id] <- labels[group]
    }
    labels.loc = labels.loc[order(labels.loc[,id]),]
    labels.loc = cbind.data.frame(labels.loc, 
                                 color = data$color[match(labels.loc[,id], data[,id])])
    geom.use <- ifelse(test = repel, yes = ggrepel::geom_label_repel, 
                       no = geom_label)
    plot = plot + geom.use(data = labels.loc, mapping = aes_string(x = xynames["x"], 
                                                                   y = xynames["y"], 
                                                                   label = id),
                           colour = labels.loc$color,
                           alpha = alpha)
    return(plot)
    
}


# modify it for Seurat v3
RunHarmony.1 <- function (object, group.by.vars, dims.use, theta = NULL, lambda = NULL, sigma = 0.1, 
                          nclust = 100, tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                          max.iter.cluster = 20, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                          plot_convergence = FALSE, verbose = TRUE, reference_values = NULL)
{
        if (!"Seurat" %in% class(object)) {
                stop("This Function is meant to be run on a Seurat object!")
        }
        if (!"pca" %in% names(object@reductions)) {
                stop("PCA must be computed before running Harmony.")
        }
        if (missing(dims.use)) {
                dims.use <- 1:ncol(object@reductions$pca@cell.embeddings)
        }
        else if (!all(dims.use %in% 1:ncol(object@reductions$pca@cell.embeddings))) {
                stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
        }
        if (length(dims.use) == 1) {
                stop("only specified one dimension in dims.use")
        }
        if (!group.by.vars %in% colnames(object@meta.data)) {
                stop(sprintf("ERROR: Primary integration variable [%s] is not in meta.data"))
        }
        missing.vars <- setdiff(group.by.vars, colnames(object@meta.data))
        if (length(missing.vars) > 0) {
                msg <- gettextf(ngettext(length(missing.vars),
                                         "trying to integrate over missing variable: %s",
                                         "trying to integrate over missing variables: %s",
                                         domain = "R-base"),
                                paste(missing.vars, collapse = ", "))
                stop(msg)
        }
        if (verbose) {
                message(gettextf("running Harmony using %d PCs", length(dims.use)))    
        }
        harmonyEmbed <- HarmonyMatrix(object@reductions$pca@cell.embeddings,
                                      object@meta.data, group.by.vars, FALSE, 0, 
                                      theta, lambda, sigma, nclust, tau, block.size, max.iter.harmony, 
                                      max.iter.cluster, epsilon.cluster, epsilon.harmony,
                                      plot_convergence, FALSE, verbose, reference_values)
        rownames(harmonyEmbed) <- row.names(object@meta.data)
        colnames(harmonyEmbed) <- paste0("harmony_", 1:ncol(harmonyEmbed))
        
        #data.use <- PrepDR(object = object,features = NULL,verbose = F)
        #feature.loadings <- (as.matrix(x = data.use) %*% as.matrix(x = harmonyEmbed))
        
        object[["harmony"]] <- CreateDimReducObject(
        embeddings = harmonyEmbed,
        #loadings = feature.loadings,
        key = "harmony_",
        assay = DefaultAssay(object = object)
        )
        return(object)
}


# modify it for monocle v3
RunHarmony.2 <- function (cds, group.by.vars, dims.use, theta = NULL, lambda = NULL, sigma = 0.1, 
                          nclust = 100, tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                          max.iter.cluster = 20, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                          plot_convergence = FALSE, verbose = TRUE, reference_values = NULL)
{
    if (!"CellDataSet" %in% class(cds)) {
        stop("This Function is meant to be run on a Seurat object!")
    }
    if (ncol(cds@normalized_data_projection) == 0) {
        stop("PCA must be computed before running Harmony.")
    }
    if (missing(dims.use)) {
        dims.use <- 1:ncol(cds@normalized_data_projection)
    }
    else if (!all(dims.use %in% 1:ncol(cds@normalized_data_projection))) {
        stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
    }
    if (length(dims.use) == 1) {
        stop("only specified one dimension in dims.use")
    }
    if (!group.by.vars %in% colnames(pData(cds))) {
        stop(sprintf("ERROR: Primary integration variable [%s] is not in meta.data",group.by.vars))
    }
    if (verbose) {
        message(gettextf("running Harmony using %d PCs", length(dims.use)))    
    }
    harmonyEmbed <- HarmonyMatrix(cds@normalized_data_projection,
                                  pData(cds), group.by.vars, FALSE, 0, 
                                  theta, lambda, sigma, nclust, tau, block.size, max.iter.harmony, 
                                  max.iter.cluster, epsilon.cluster, epsilon.harmony,
                                  plot_convergence, FALSE, verbose, reference_values)

    cds@normalized_data_projection = harmonyEmbed  
    return(cds)
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


# Run GSEA and generate reports
#'@example ReportGSEA(file = "c5.all.827-D14.827-D16.827-D28.827-D0.Gsea.1552707417126",pos=T,ncol=3)
ReportGSEA <- function(file, pos=T,ncol = 3){
        (gsea_path <- paste("~/gsea_home/output",tolower(format(Sys.Date(), "%b%d")),
                            file,sep ='/'))
        #(pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*pos.*xls"))
        p<-1; if(pos) p <-2
        (pos.xls.path <- list.files(gsea_path,pattern="gsea_report_for_.*xls"))
        GSEA_output <- readr::read_delim(paste0(gsea_path,"/",pos.xls.path[p]),"\t", 
                                         escape_double = FALSE, trim_ws = TRUE)
        print(GSEA_output %>% .[-c(2,3,12)] %>% head(50) %>% kable() %>% kable_styling())
        
        (GSEA.plots <- sapply(GSEA_output$NAME[1:9], function(name) {
                paste0("enplot_",name, "_([0-9]+)*\\.png$")}) %>%
                        sapply(function(x) list.files(path = gsea_path, pattern =x)) %>%
                        .[sapply(.,length)>0] %>% #Remove empty elements from list with character(0)
                        paste(gsea_path, ., sep = "/")) 
        CombPngs(GSEA.plots, ncol = ncol)
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        GSEA_path <- paste0(path, "GSEA_xls")
        if(!dir.exists(GSEA_path)) dir.create(GSEA_path, recursive = T)
        file.copy(paste0(gsea_path,"/",pos.xls.path),GSEA_path)
        file.rename(paste0(path,list.files(path,pattern=paste0("GSEA.plots_CombPngs"))), 
                    paste0(path,file,"-",p,".jpeg"))
        
}


#' re-order seurat idents factors
sortIdent <- function(object,numeric=F){
    
    levels  <- Idents(object) %>% unique %>% as.character
    if(numeric) levels %<>% as.numeric 
    if(!is.null(FindIdentLabel(object))) {
        object@meta.data[,FindIdentLabel(object)] %<>% factor(levels = sort(levels))
    }
    Idents(object) %<>% factor(levels = sort(levels))
    
    return(object)
}


#' Modified TSNEPlot
#' @param label.repel 
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
TSNEPlot.1 <- function(object,dims = c(1, 2),cells = NULL,cols = NULL,pt.size = NULL,
                       reduction = "tsne",group.by = NULL,split.by = NULL,shape.by = NULL,
                       order = NULL,label = FALSE,label.repel= TRUE, label.size = 4,
                       repel = TRUE,alpha = 0.85, 
                       cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
                       na.value = 'grey50',combine = TRUE,ncol = NULL,title = NULL,
                       no.legend = F,do.print = F,do.return = T,unique.name=F,
                       width=10, height=7,border = FALSE, ...) {
    if(unique.name) {
        VarName <- unique(object$orig.ident)
        VarName <- paste(VarName[1:min(5,length(VarName))],collapse = "_")
    } else VarName <- deparse(substitute(object))
    VarName = paste0(VarName,"_",group.by %||% FindIdentLabel(object),
                     ifelse(!is.null(split.by), yes = paste0("_",split.by), no =""))
    
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    reduction <- reduction %||% Seurat:::DefaultDimReduc(object = object)
    cells <- cells %||% colnames(x = object)
    data <- object@reductions$tsne@cell.embeddings[cells,dims]
    data <- as.data.frame(x = data)
    dims <- paste0(object@reductions$tsn@key, dims)
    object[["ident"]] <- Idents(object = object)
    group.by <- group.by %||% "ident"
    data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    plots <- lapply(X = group.by, FUN = function(x) {
        plot <- Seurat:::SingleDimPlot(data = data[, c(dims, x, split.by, shape.by)], 
                                       dims = dims, col.by = x, cols = cols,
        pt.size = pt.size, shape.by = shape.by, order = order,
        label = FALSE, cells.highlight = cells.highlight,
        cols.highlight = cols.highlight, sizes.highlight = sizes.highlight,
        na.value = na.value, ...)
        if (label & label.repel == F ) {
            plot <- LabelClusters(plot = plot, id = x, repel = repel,
            size = label.size, split.by = split.by)
        }
        if (label & label.repel) {
            plot <- LabelRepel(plot = plot, id = x, repel = repel,
                               size = label.size, split.by = split.by, color= cols,
                               alpha = alpha)
        }
        if (!is.null(x = split.by)) {
            plot <- plot + Seurat:::FacetTheme() + 
                facet_wrap(facets = vars(!!sym(x = split.by)),
                           ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                               length(x = unique(x = data[, split.by]))
                               } else ncol )+
                theme(strip.text = element_text(face="plain",size=15))
            
        }
        if(border == TRUE) plot = plot + theme(panel.border = element_rect(colour = "black"))
        return(plot)
    })
    if (combine) {
        plots <- CombinePlots(plots = plots, ncol = if (!is.null(x = split.by) &&
        length(x = group.by) > 1) {
            1
        }
        else {
            ncol
        }, ...)
    }
    if(!is.null(title)){
        plots = plots + ggtitle(title)+
            theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
    }
    if(no.legend) {
        plots = plots + NoLegend()
        L = ""
    } else L = "_Legend"
    if(do.print) {
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"TSNEPlot_",VarName,L,".jpeg"), 
             units="in", width=width, height=height,res=600)
        print(plots)
        dev.off()
    }
    if(do.return) return(plots)
}


#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
UMAPPlot.1 <- function(
    object,dims = c(1, 2),cells = NULL,cols = NULL,pt.size = NULL,
    reduction = NULL,group.by = NULL,split.by = NULL,shape.by = NULL,
    order = NULL,label = FALSE,label.size = 4,repel = FALSE,
    cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
    na.value = 'grey50',combine = TRUE,ncol = NULL,title = NULL,
    no.legend = F,do.print = F,do.return = T,unique.name=F,...) {
    if(unique.name) {
        VarName <- unique(object$orig.ident)
        VarName <- paste(VarName[1:min(5,length(VarName))],collapse = "_")
    } else VarName <- deparse(substitute(object))
    VarName = paste0(VarName,"_",group.by %||% FindIdentLabel(object))
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    cells <- cells %||% colnames(x = object)
    data <- object@reductions$umap@cell.embeddings[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0("UMAP_", dims)
    object[['ident']] <- Idents(object = object)
    group.by <- group.by %||% 'ident'
    data[, group.by] <- object[[group.by]][cells, , drop = FALSE]
    for (group in group.by) {
        if (!is.factor(x = data[, group])) {
            data[, group] <- factor(x = data[, group])
        }
    }
    if (!is.null(x = shape.by)) {
        data[, shape.by] <- object[[shape.by, drop = TRUE]]
    }
    if (!is.null(x = split.by)) {
        data[, split.by] <- object[[split.by, drop = TRUE]]
    }
    plots <- lapply(
        X = group.by,
        FUN = function(x) {
            plot <- SingleDimPlot(
                data = data[, c(dims, x, split.by, shape.by)],
                dims = dims,
                col.by = x,
                cols = cols,
                pt.size = pt.size,
                shape.by = shape.by,
                order = order,
                label = FALSE,
                cells.highlight = cells.highlight,
                cols.highlight = cols.highlight,
                sizes.highlight = sizes.highlight,
                na.value = na.value
            )
            if (label) {
                plot <- LabelClusters(
                    plot = plot,
                    id = x,
                    repel = repel,
                    size = label.size,
                    split.by = split.by
                )
            }
            if (!is.null(x = split.by)) {
                plot <- plot + Seurat:::FacetTheme() +
                    facet_wrap(
                        facets = vars(!!sym(x = split.by)),
                        ncol = if (length(x = group.by) > 1) {
                            length(x = unique(x = data[, split.by]))
                        } else {
                            NULL
                        }
                    )
            }
            return(plot)
        }
    )
    if (combine) {
        plots <- CombinePlots(
            plots = plots,align= align,
            ncol = if (!is.null(x = split.by) && length(x = group.by) > 1) {
                1
            } else {
                ncol
            },
            ...
        )
    }
    if(!is.null(title)){
        plots = plots + ggtitle(title)+
            theme(plot.title = element_text(hjust = 0.5,size=15,face = "plain"))
    }
    if(no.legend) {
        plots = plots + NoLegend()
        L = ""
    } else L = "_Legend"
    if(do.print) {
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"UMAPPlot_",VarName,L,".jpeg"), 
             units="in", width=10, height=7,res=600)
        print(plots)
        dev.off()
    }
    if(do.return) return(plots)
}
