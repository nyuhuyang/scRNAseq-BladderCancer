# FeaturePlot in Seurat 3 doesn't support two gene overlay.
# modify FeaturePlot in Seurat 2 and support Seurat 3 object
FeaturePlot.2 <- function (object, features, min.cutoff = NA, max.cutoff = NA,
                           dims = c(1,2), cells = NULL, pt.size = 1,slot = "data",
                           cols.use = c("yellow",  "red"), pch.use = 16, overlay = FALSE, do.hover = FALSE,
                           data.hover = "ident", do.identify = FALSE, reduction = "tsne",
                           use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE,
                           dark.theme = FALSE, do.return = FALSE, vector.friendly = FALSE)
{
        cells <- cells %||% colnames(x = object)
        if (is.null(x = nCol)) {
                nCol <- 2
                if (length(x = features) == 1) {
                        nCol <- 1
                }
                if (length(x = features) > 6) {
                        nCol <- 3
                }
                if (length(x = features) > 9) {
                        nCol <- 4
                }
        }
        num.row <- floor(x = length(x = features)/nCol - 1e-05) +
                1
        if (overlay | do.hover) {
                num.row <- 1
                nCol <- 1
        }
        par(mfrow = c(num.row, nCol))
        dims <- paste0(Key(object = object[[reduction]]), dims)
        data.plot <- FetchData(object = object, vars = dims, cells = cells, slot = slot)
        x1 <- dims[1]
        x2 <- dims[2]
        data.plot$x <- data.plot[, x1]
        data.plot$y <- data.plot[, x2]
        data.plot$pt.size <- pt.size
        names(x = data.plot) <- c("x", "y")
        data.use <- t(x = FetchData(object = object, vars = features, cells = cells, slot = slot))
        min.cutoff <- mapply(FUN = function(cutoff, feature) {
                ifelse(test = is.na(x = cutoff), yes = min(data.use[feature,
                                                                    ]), no = cutoff)
        }, cutoff = min.cutoff, feature = features)
        
        max.cutoff <- mapply(FUN = function(cutoff, feature) {
                ifelse(test = is.na(x = cutoff), yes = max(data.use[feature,
                                                                    ]), no = cutoff)
        }, cutoff = max.cutoff, feature = features)
        check_lengths = unique(x = vapply(X = list(features,
                                                   min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
        if (length(x = check_lengths) != 1) {
                stop("There must be the same number of minimum and maximum cuttoffs as there are features")
        }
        if (overlay) {
                pList <- list(BlendPlot(data.use = data.use, features = features,
                                        data.plot = data.plot, pt.size = pt.size, pch.use = pch.use,
                                        cols.use = cols.use, dim.codes = dims, min.cutoff = min.cutoff,
                                        max.cutoff = max.cutoff, no.axes = no.axes, no.legend = no.legend,
                                        dark.theme = dark.theme))
        }
        else {
                pList <- mapply(FUN = SingleFeaturePlot, feature = features,
                                min.cutoff = min.cutoff, max.cutoff = max.cutoff,
                                MoreArgs = list(data.use = data.use, data.plot = data.plot,
                                                pt.size = pt.size, pch.use = pch.use, cols.use = cols.use,
                                                dim.codes = dims, no.axes = no.axes, no.legend = no.legend,
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
                                    title = features))
        }
        else if (do.identify) {
                if (length(x = pList) != 1) {
                        stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
                }
                return(FeatureLocator(plot = pList[[1]], data.plot = data.plot,
                                      dark.theme = dark.theme))
        }
        #else {
        #        print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
        #}
        par(mfrow = c(1, 1))
        if (do.return) {
                return(pList[[1]])
        }
}



# Support FeaturePlot.1
BlendPlot <- function (data.use, features, data.plot, pt.size, pch.use,
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
        length.check <- vapply(X = list(features, min.cutoff,
                                        max.cutoff), FUN = function(x) {
                                                return(length(x = x) != 2)
                                        }, FUN.VALUE = logical(length = 1))
        if (any(length.check)) {
                stop("An overlayed FeaturePlot only works with two features and requires two minimum and maximum cutoffs")
        }
        min.cutoff <- c(SetQuantile(cutoff = min.cutoff[1], data = data.gene[features[1],
                                                                             ]), SetQuantile(cutoff = min.cutoff[2], data = data.gene[features[2],
                                                                                                                                      ]))
        max.cutoff <- c(SetQuantile(cutoff = max.cutoff[1], data = data.gene[features[1],
                                                                             ]), SetQuantile(cutoff = max.cutoff[2], data = data.gene[features[2],
                                                                                                                                      ]))
        data.gene <- stats::na.omit(object = data.frame(data.use[features,
                                                                 ]))
        cell.names <- colnames(x = data.gene)
        data.gene <- matrix(data = vapply(X = data.gene, FUN = function(x) ifelse(test = x <
                                                                                          min.cutoff, yes = min.cutoff, no = x), FUN.VALUE = c(1,
                                                                                                                                               1)), nrow = 2)
        data.gene <- matrix(data = vapply(X = as.data.frame(x = data.gene),
                                          FUN = function(x) ifelse(test = x > max.cutoff, yes = max.cutoff,
                                                                   no = x), FUN.VALUE = c(1, 1)), nrow = 2)
        data.gene <- as.data.frame(x = data.gene)
        rownames(x = data.gene) <- features
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
        legend.names <- c(high1 = paste("High", features[1]),
                          high2 = paste("High", features[2]), highboth = "High both")
        title <- paste0(features, collapse = " + ")
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
        p <- p + theme_bw() + theme(panel.border = element_blank(), 
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(), 
                                    axis.line = element_line(colour = "black")) 
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
