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


#' DoHeatmap.1, automatically group top DE genes from FindMakers output
#example   DoHeatmap.1(SSCs,top,Top_n = 15, 
#                   group.order = major_cells,ident.use = "all cell types",
#                   group.label.rot = T,cex.row = 5,remove.key =T)
#
DoHeatmap.1 <- function(object, marker_df, add.genes = NULL, no.legend =F,unique.name= F,
                        Top_n = 10, features = NULL, cells = NULL, group.by = "ident", 
                         group.bar = TRUE, disp.min = -2.5, disp.max = NULL, slot = "scale.data", 
                         assay = NULL, label = TRUE, cols=NULL, size = 5.5, hjust = 0, angle = 45, 
                         raster = TRUE, draw.lines = TRUE, lines.width = NULL, group.bar.height = 0.02, 
                         combine = TRUE,title = NULL,title.size = 14,do.print = FALSE,
                        cex.row=12,legend.size = NULL,units="in", width=10, height=7,res=600,...){
    if(unique.name) {
        v <- paste(unique(object$orig.ident),collapse = "_")
    } else v <- deparse(substitute(object))
    v = paste0(v,"_",FindIdentLabel(object))
    if (!missing(x = marker_df)) {
        
        colnames(marker_df)[grep("cluster*.",colnames(marker_df))[1]] = "cluster"
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
    heatmap <- DoHeatmap(object = object, features = as.vector(top_gene), cells = cells, group.by = group.by,
                         group.bar = group.bar, disp.min = disp.min, disp.max = disp.max,
                         slot = slot, assay = assay, label = label, size = size,
                         hjust = hjust, angle = angle, raster = raster, draw.lines = draw.lines,
                         lines.width = lines.width, group.bar.height = group.bar.height, 
                         combine = combine)+
        scale_y_discrete(position = "right")
    heatmap = heatmap + scale_fill_gradientn(colors = ggsci::pal_gsea()(12))
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


FgseaBarplot <- function(pathways=hallmark, stats=res, nperm=1000,cluster = 1,
                         sample="",pathway.name = "Hallmark", hjust=0.5){
    
    res = stats[order(stats["p_val_adj"]),]
    res1 = res[res$cluster == cluster,c("gene","avg_logFC")] %>% deframe
    
    fgseaRes <- fgsea(pathways=pathways, stats=res1, nperm=nperm)
    print(dim(fgseaRes))
    topPathwaysUp <- fgseaRes[NES > 0][head(order(padj), n=25), pathway]
    topPathwaysDown <- fgseaRes[NES < 0][head(order(padj), n=25), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
    out_write=fgseaRes[match(topPathways,pathway)]
    colnames(out_write)[1] <- 'Pathways'
    out_write$To_Plot_value <- -log10(out_write$pval)
    out_write$sign <- ifelse(out_write$NES >0,1,-1)
    out_write$To_Plot_value <- out_write$To_Plot_value*out_write$sign
    out_write$sign <- ifelse(out_write$NES >0,"Upregulated", "Downregulated")
    
    p<-ggbarplot(out_write,
                 x = "Pathways",
                 y = "NES",
                 #fill = "sign",           # change fill color by mpg_level
                 color = "white",            # Set bar border colors to white
                 palette = "jco",            # jco journal color palett. see ?ggpar
                 sort.val = "asc",          # Sort the value in descending order
                 sort.by.groups = FALSE,     # Don't sort inside each group
                 x.text.angle = 90,          # Rotate vertically x axis texts
                 ylab = 'Normalized Enrichment Score',
                 legend.title = "padj < 0.25",
                 rotate = TRUE,
                 title = paste(pathway.name,"pathways in",sample,cluster),
                 ggtheme = theme_minimal(base_size = 15))+
        geom_col(aes(fill=padj<0.25))+
        guides(fill = guide_legend(reverse = TRUE))+
        theme(text = element_text(size=12),
              plot.title = element_text(size = 12,hjust = hjust))
    path <- paste0("output/",gsub("-","",Sys.Date()),"/")
    if(!dir.exists(path)) dir.create(path, recursive = T)
    jpeg(paste0(path,sample,"_",cluster,"-",pathway.name,".jpeg"), units="in", width=10, height=7,res=600)
    print(p)
    dev.off()
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
    ave_UMI.1 <- Matrix::rowMeans(expm1(x = data.use[, ident.1]))
    ave_UMI.2 <- Matrix::rowMeans(expm1(x = data.use[, ident.2]))
    avg_UMI <-data.frame(ave_UMI.1, ave_UMI.2)
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
    gde <- list()
    for(i in 1:length(ident.1)) {
        ident.1vs2 <- paste(ident.1[i], ident.2[i], sep = ".vs.")
        print(ident.1vs2)
        ifelse(class(ident.1[i])=="list", ident1 <- ident.1[[i]], ident1 <- ident.1[i])
        ifelse(class(ident.2[i])=="list", ident2 <- ident.2[[i]], ident2 <- ident.2[i])
        num1 <- as.numeric(unique(gsub('.*\\_', '', ident1)))
        num2 <- as.numeric(unique(gsub('.*\\_', '', ident2)))
        suppressWarnings(if(num1 == num2) num1 = "")
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
#' @param add.genes extra genes to add beyound FindAllMarkers results 
#' @param remove.legend TRUE/FALSE
#' @param color color scheme
#' @export g vertical ggplot
#' @example MakeCorlorBar(EC, top, Top_n=40)
MakeCorlorBar <- function(object, marker_df, Top_n = NULL, add.genes = NULL, color =NULL,
                          no.legend =F, legend.size = NULL,do.print = TRUE,do.return=FALSE,width=10, height=7){
    
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

    gene.use = rownames(object@assays$RNA@scale.data)
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
              axis.ticks.y=element_blank(),
              legend.title = element_text(size=legend.size*1.2),
              legend.text = element_text(size=legend.size),
              legend.key.size = unit(legend.size/5,"line"))+
        coord_flip() + scale_x_reverse()+
    if(no.legend) g = g + theme(legend.position="none")
    if(!is.null(color)) {
        colors_fill = color_scheme(color = color,
                                   names = unique(marker_df$cluster))
        g = g + scale_fill_manual(values = colors_fill)
    }
    if(do.print) {
        v <- deparse(substitute(object))
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,"Doheatmap_",v,"_top",Top_n,"_",FindIdentLabel(object),"_vbar.jpeg"),
             units="in", width=width, height=height,res=600)
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

# Seurat 3
SingleDimPlot <- function (data, dims, col.by = NULL, cols = NULL, pt.size = NULL, 
                           shape.by = NULL, order = NULL, label = FALSE, repel = FALSE, 
                           label.size = 4, cells.highlight = NULL, cols.highlight = "red", 
                           sizes.highlight = 1, na.value = "grey50") 
{
    pt.size <- pt.size %||% Seurat:::AutoPointSize(data = data)
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    if (!is.data.frame(x = data)) {
        data <- as.data.frame(x = data)
    }
    if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
        stop("Cannot find dimensions to plot in data")
    }
    else if (is.numeric(x = dims)) {
        dims <- colnames(x = data)[dims]
    }
    if (!is.null(x = cells.highlight)) {
        highlight.info <- SetHighlight(cells.highlight = cells.highlight, 
                                       cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||% 
                                           pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||% 
                                           "black", pt.size = pt.size)
        order <- highlight.info$plot.order
        data$highlight <- highlight.info$highlight
        col.by <- "highlight"
        pt.size <- highlight.info$size
        cols <- highlight.info$color
    }
    if (!is.null(x = order) && !is.null(x = col.by)) {
        if (typeof(x = order) == "logical") {
            if (order) {
                data <- data[order(data[, col.by]), ]
            }
        }
        else {
            order <- rev(x = c(order, setdiff(x = unique(x = data[, 
                                                                  col.by]), y = order)))
            data[, col.by] <- factor(x = data[, col.by], levels = order)
            new.order <- order(x = data[, col.by])
            data <- data[new.order, ]
            if (length(x = pt.size) == length(x = new.order)) {
                pt.size <- pt.size[new.order]
            }
        }
    }
    if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
        warning("Cannot find ", col.by, " in plotting data, not coloring plot")
        col.by <- NULL
    }
    else {
        col.index <- match(x = col.by, table = colnames(x = data))
        if (grepl(pattern = "^\\d", x = col.by)) {
            col.by <- paste0("x", col.by)
        }
        else if (grepl(pattern = "-", x = col.by)) {
            col.by <- gsub(pattern = "-", replacement = ".", 
                           x = col.by)
        }
        colnames(x = data)[col.index] <- col.by
    }
    if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
        warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
    }
    plot <- ggplot(data = data) + geom_point(mapping = aes_string(x = dims[1], 
                                                                  y = dims[2], color = paste0("`", col.by, "`"), shape = shape.by), 
                                             size = pt.size) + guides(color = guide_legend(override.aes = list(size = 3))) + 
        labs(color = NULL)
    if (label && !is.null(x = col.by)) {
        plot <- LabelClusters(plot = plot, id = col.by, repel = repel, 
                              size = label.size)
    }
    if (!is.null(x = cols)) {
        plot <- plot + if (length(x = cols) == 1) {
            scale_color_brewer(palette = cols, na.value = na.value)
        }
        else {
            scale_color_manual(values = cols, na.value = na.value)
        }
    }
    plot <- plot + cowplot::theme_cowplot()
    return(plot)
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


#' Add several args from Seurat 2
#' @param no.legend remove legend
#' @param title add ggplot title
#' @param do.print save jpeg file
#' @param unique.name save jpeg file with unique name
#' @param do.return return plot
TSNEPlot.1 <- function(
    object,dims = c(1, 2),cells = NULL,cols = NULL,pt.size = NULL,
    reduction = NULL,group.by = NULL,split.by = NULL,shape.by = NULL,
    order = NULL,label = FALSE,label.size = 4,repel = FALSE,
    cells.highlight = NULL,cols.highlight = 'red',sizes.highlight = 1,
    na.value = 'grey50',combine = TRUE,ncol = NULL,title = NULL,
    no.legend = F,do.print = F,do.return = T,unique.name=F,...) {
    if(unique.name) {
        v <- unique(object$orig.ident)
        v <- paste(v[1:min(5,length(v))],collapse = "_")
    } else v <- deparse(substitute(object))
    v = paste0(v,"_",FindIdentLabel(object))
    if (length(x = dims) != 2) {
        stop("'dims' must be a two-length vector")
    }
    cells <- cells %||% colnames(x = object)
    data <- object@reductions$tsne@cell.embeddings[cells, dims]
    data <- as.data.frame(x = data)
    dims <- paste0("tSNE_", dims)
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
        jpeg(paste0(path,"TSNEPlot_",v,L,".jpeg"), 
             units="in", width=10, height=7,res=600)
        print(plots)
        dev.off()
    }
    if(do.return) return(plots)
}


#' Add several args from Seurat 2
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
        v <- unique(object$orig.ident)
        v <- paste(v[1:min(5,length(v))],collapse = "_")
    } else v <- deparse(substitute(object))
    v = paste0(v,"_",FindIdentLabel(object))
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
        jpeg(paste0(path,"UMAPPlot_",v,L,".jpeg"), 
             units="in", width=10, height=7,res=600)
        print(plots)
        dev.off()
    }
    if(do.return) return(plots)
}
