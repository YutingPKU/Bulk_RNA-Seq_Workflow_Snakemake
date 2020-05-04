######Some supplemental fns that are used by heatmap_plot*.R and cluster_plot.R
#BASED on heatmap_supp_funcs.R written by Henry Long
library(methods)

zscore = function(x){
    y=(x-mean(x))/sd(x)
    return(y)
}

cmap <- function(x, colorstart=NULL, use_viridis=FALSE) {
    colors = c("#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363", "#BD4931", "#74c476", "#6baed6", "#fd8d3c", "#9e9ac8", "#969696", "#D67D6B", "#9ecae1", "#fdae6b", "#bcbddc", "#bdbdbd", "#E0A89D", "#c6dbef", "#fdd0a2", "#c7e9c0", "#d9d9d9", "#F0CEC7","#006400", "#FFD700","#1874CD","#575757","#FF3030", "#a1d99b","#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363", "#BD4931", "#6baed6", "#fd8d3c", "#74c476", "#9e9ac8", "#969696", "#D67D6B", "#9ecae1", "#fdae6b", "#a1d99b", "#bcbddc", "#bdbdbd", "#E0A89D", "#c6dbef", "#fdd0a2", "#c7e9c0", "#dadaeb", "#d9d9d9", "#F0CEC7","#006400", "#FFD700","#1874CD","#575757","#FF3030")
    x <- sort(unique(na.omit(as.vector(x))))
    if(use_viridis) {
        col <- viridis(length(x))
    } else {
        col <- colors[(colorstart+1):(colorstart+length(x))]
    }    
    names(col) <- x
    return(col)
}

make_complexHeatmap_annotation <- function(annotation){
    MIN_UNIQUE <- 6
    global_gp = gpar(fontsize = 3)
    title_gp = gpar(fontsize = 3)

    colorlist <- list()
    colorcount = 0
    nn<-length(annotation)
    for (i in 1:nn) {
        ann <- as.matrix(annotation[,i])
        if(length(sort(unique(na.omit(as.vector(ann))))) < MIN_UNIQUE | is.numeric(ann)==FALSE) {
            colorlist[[i]] <- cmap(ann, colorstart=colorcount)
            colorcount = colorcount + length(unique(ann))
        } else {
            colorlist[[i]] <- colorRamp2(seq(min(ann, na.rm = TRUE), max(ann, na.rm = TRUE), length = 3), c("white","yellow", "red"))
        }
    }
    names(colorlist) <- c(colnames(annotation)[1:nn])
    
    ha1 = HeatmapAnnotation(df = annotation[,1:nn,drop=FALSE], gap=unit(0.5,"mm"), col = colorlist,annotation_name_gp = gpar(fontsize = 3),
                            annotation_height = unit(1,"mm"), height = unit(1,"mm"),simple_anno_size = unit(1.5, "mm"),
                            annotation_legend_param = list(title_gp=gpar(fontsize=3), grid_height = unit(3,"mm"), labels_gp=gpar(fontsize=3)))

    return(ha1)
}

#NOTE: LEN removed the threeD code
#make_pca_plots <- function(data_matrix, threeD = TRUE, labels = TRUE, pca_title = "Data Matrix", legend_title = "", ClassColorings) {
#
#    #Standard PCA analysis
#    pca_out <- prcomp(data_matrix, scale. = TRUE, tol = 0.05)
#    pc_var <- signif(100.0 * summary(pca_out)[[6]][2,], digits = 3)
#    #pc_var <- signif(100.0 * summary(pca_out)[[6]][2,1:3], digits = 3)
#        
#    #### NEW
#    #pc_var <- signif(100.0 * summary(pca_out)[[6]][2,1:3], digits = 3)
#    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#
#    plot(pca_out$x[,"PC1"], pca_out$x[,"PC2"],  col=ClassColorings, pch=16, xlab=paste0("PC1 (", pc_var[1], "% of variance)"), ylab=paste0("PC2 (", pc_var[2], "% of variance)"), main = paste0('PCA analysis of ',pca_title))
#    if(labels == TRUE) {text(pca_out$x[,"PC1"], pca_out$x[,"PC2"], labels=row.names(data_matrix), cex= 0.7, pos=3)}
#    if(legend_title != "") {
#        mycols = unique(ClassColorings)
#        mynames = unique(names(ClassColorings))
#        t = sort(mynames)
#        tt = ClassColorings[t]
#        legend("topright", inset=c(-0.23,0), legend = unique(t), col=unique(tt), pch = 16, title = legend_title)
#    }
#
#    if(threeD==TRUE){
#        #try 3D plot
#        library("rgl", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
#        pca3d<-cbind(pca_out$x[,1], pca_out$x[,2], pca_out$x[,3])
#        plot3d(pca3d, type="s",col=ClassColorings, size=1, scale=0.2)
#    }
#
#    return(pca_out)
#}


######### FROM https://github.com/vqv/ggbiplot/blob/master/R/ggscreeplot.r
ggscreeplot <- function(pcobj, type = c('pev', 'cev')) {
    type <- match.arg(type)
    d <- pcobj$sdev^2
    yvar <- switch(type,
                   pev = d / sum(d),
                   cev = cumsum(d) / sum(d))
    yvar.lab <- switch(type,
                       pev = 'proportion of explained variance',
                       cev = 'cumulative proportion of explained variance')
    df <- data.frame(PC = 1:length(d), yvar = yvar)
    fontsize = 4
    linesize = .5
    ggplot(data = df, aes(x = PC, y = yvar)) +
      xlab('principal component number') + ylab(yvar.lab) +
      geom_point(size = .5) + geom_path() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank(), 
            panel.background = element_rect(fill = "white", colour = NA),
            panel.border = element_rect(size = linesize, color = 'black',fill = NA  ))+
      theme(axis.text.x = element_text(size=fontsize),
            axis.text.y = element_text(size=fontsize),  
            axis.title.x = element_text(size=fontsize),
            axis.title.y = element_text(size=fontsize),
            axis.line = element_blank(), axis.ticks = element_line(size = .3),
            legend.position = "none", 
            #legend.text = element_text(size = fontsize),
            #legend.title  = element_text(size = fontsize),legend.key.size = unit(0.2, "cm"),
            plot.title = element_text(size=fontsize, hjust = 0.5))
    
}
######### FROM https://github.com/vqv/ggbiplot/blob/master/R/ggbiplot.r - with edits by Mahesh Vangala

ggbiplot <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE,
                      obs.scale = 1 - scale, var.scale = scale,
                      groups = NULL, ellipse = FALSE, ellipse.prob = 0.68,
                      labels = NULL, labels.size = 3, alpha = 1,
                      var.axes = TRUE,
                      circle = FALSE, circle.prob = 0.69,
                      varname.size = 3, varname.adjust = 1.5,
                      varname.abbrev = FALSE, ...)
{
#    library(ggplot2)
#    library(plyr)
#    library(scales)
#    library(grid)
#    library(ggrepel)

    stopifnot(length(choices) == 2)

    # Recover the SVD
    if(inherits(pcobj, 'prcomp')){
        nobs.factor <- sqrt(nrow(pcobj$x) - 1)
        d <- pcobj$sdev
        u <- sweep(pcobj$x, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- pcobj$rotation
    } else if(inherits(pcobj, 'princomp')) {
        nobs.factor <- sqrt(pcobj$n.obs)
        d <- pcobj$sdev
        u <- sweep(pcobj$scores, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- pcobj$loadings
    } else if(inherits(pcobj, 'PCA')) {
        nobs.factor <- sqrt(nrow(pcobj$call$X))
        d <- unlist(sqrt(pcobj$eig)[1])
        u <- sweep(pcobj$ind$coord, 2, 1 / (d * nobs.factor), FUN = '*')
        v <- sweep(pcobj$var$coord,2,sqrt(pcobj$eig[1:ncol(pcobj$var$coord),1]),FUN="/")
    } else if(inherits(pcobj, "lda")) {
        nobs.factor <- sqrt(pcobj$N)
        d <- pcobj$svd
        u <- predict(pcobj)$x/nobs.factor
        v <- pcobj$scaling
        d.total <- sum(d^2)
    } else {
         stop('Expected a object of class prcomp, princomp, PCA, or lda')
    }

    # Scores
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[,choices], 2, d[choices]^obs.scale, FUN='*'))

    # Directions
    v <- sweep(v, 2, d^var.scale, FUN='*')
    df.v <- as.data.frame(v[, choices])

    names(df.u) <- c('xvar', 'yvar')
    names(df.v) <- names(df.u)

    if(pc.biplot) {
        df.u <- df.u * nobs.factor
    }

    # Scale the radius of the correlation circle so that it corresponds to
    # a data ellipse for the standardized PC scores
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)

    # Scale directions
    v.scale <- rowSums(v^2)
    df.v <- r * df.v / sqrt(max(v.scale))

    # Change the labels for the axes
    if(obs.scale == 0) {
        u.axis.labs <- paste('standardized PC', choices, sep='')
    } else {
        u.axis.labs <- paste('PC', choices, sep='')
    }

    # Append the proportion of explained variance to the axis labels
    u.axis.labs <- paste(u.axis.labs,
                         sprintf('(%0.1f%% explained var.)',
                                 100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))

    # Score Labels
    if(!is.null(labels)) {
        df.u$labels <- labels
    }

    # Grouping variable
    if(!is.null(groups)) {
        df.u$groups <- groups
    }

    # Variable Names
    if(varname.abbrev) {
        df.v$varname <- abbreviate(rownames(v))
    } else {
        df.v$varname <- rownames(v)
    }

    # Variables for text label placement
    df.v$angle <- with(df.v, (180/pi) * atan(yvar / xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar)) / 2)

    # Base plot
    g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) +
                xlab(u.axis.labs[1]) + ylab(u.axis.labs[2]) + coord_equal()

    if(var.axes) {
        # Draw circle
        if(circle)
        {
            theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
            circle <- data.frame(xvar = r * cos(theta), yvar = r * sin(theta))
            g <- g + geom_path(data = circle, color = muted('white'),
                               size = 1/2, alpha = 1/3)
        }

        # Draw directions
        g <- g +
             geom_segment(data = df.v,
             aes(x = 0, y = 0, xend = xvar, yend = yvar),
             arrow = arrow(length = unit(1/2, 'picas')),
             color = muted('red'))
    }

    # Draw either labels or points
    if(!is.null(df.u$labels)) {
        if(!is.null(df.u$groups)) {
            g <- g + geom_point(aes(xvar, yvar), color = "black",
                                size = 2) + geom_text_repel(aes(xvar, yvar, color = groups,
                                label = labels), size = 3, fontface = "bold",
                                box.padding = unit(0.5, "lines"), point.padding = unit(1.6,
                                    "lines"), segment.color = "#555555", segment.size = 0.5,
                                arrow = arrow(length = unit(0.01, "npc")), force = 1,
                                max.iter = 2000) + geom_point(aes(color=groups))
        } else {
            g <- g + geom_text(aes(label = labels), size = labels.size)
        }
    } else {
            if(!is.null(df.u$groups)) {
                g <- g + geom_point(aes(color = groups), alpha = alpha)
            } else {
                g <- g + geom_point(alpha = alpha)
            }
    }

    # Overlay a concentration ellipse if there are groups
    if(!is.null(df.u$groups) && ellipse) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
        circle <- cbind(cos(theta), sin(theta))

        ell <- ddply(df.u, 'groups', function(x) {
            if(nrow(x) <= 2) {
                return(NULL)
            }
            sigma <- var(cbind(x$xvar, x$yvar))
            mu <- c(mean(x$xvar), mean(x$yvar))
            ed <- sqrt(qchisq(ellipse.prob, df = 2))
            data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'),
                       groups = x$groups[1])
            })
            names(ell)[1:2] <- c('xvar', 'yvar')
            g <- g + geom_path(data = ell, aes(color = groups, group = groups))
        }

        # Label the variable axes
        if(var.axes) {
            g <- g +
                 geom_text(data = df.v,
                           aes(label = varname, x = xvar, y = yvar,
                               angle = angle, hjust = hjust),
                           color = 'darkred', size = varname.size)
        }
        # Change the name of the legend for groups
        # if(!is.null(groups)) {
        #   g <- g + scale_color_brewer(name = deparse(substitute(groups)),
        #                               palette = 'Dark2')
        # }

        # TODO: Add a second set of axes

    return(g)
}






######END SUPPLEMENTAL Fns #######
