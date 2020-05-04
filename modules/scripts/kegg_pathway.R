suppressMessages(library("dplyr"))
suppressMessages(library("AnnotationDbi"))
suppressMessages(library("org.Mm.eg.db"))
suppressMessages(library("org.Hs.eg.db"))
suppressMessages(library("gage"))
suppressMessages(library("gageData"))
suppressMessages(library("pathview"))
suppressWarnings(suppressMessages(library("clusterProfiler")))
suppressMessages(library("XML"))
suppressMessages(library("ggplot2"))
suppressMessages(library("stringr"))

## The traceback is actually necessary to not break pipe at the stop step, so leave on
options(error = function() traceback(2))

## Example inputs
#deseq_file = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/test_all/diffexp/BCellvsTCell/BCellvsTCell.deseq.csv"
#keggpvalcutoff = 0.1
#numkeggpathways = 5
#reference = "hg19"
#kegg_dir = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/test_all/diffexp/BCellvsTCell/kegg_pathways/"
#temp_dir = "/mnt/cfce-stor1/home/mgc31/code/viperproject/analysis/test_all/diffexp/BCellvsTCell/temp/"
#kegg_table_up = "/mnt/cfce-stor1/home/mgc31/code/viperproject/uptest.csv"
#kegg_table_down = "/mnt/cfce-stor1/home/mgc31/code/viperproject/downtest.csv"
#keggsummary_pdf = "/mnt/cfce-stor1/home/mgc31/code/viperproject/test.pdf"
#up_kegg_png = "/mnt/cfce-stor1/home/mgc31/code/viperproject/upkeggtest.png"
#down_kegg_png = "/mnt/cfce-stor1/home/mgc31/code/viperproject/downkeggtest.png"
#gsea_table = "/mnt/cfce-stor1/home/mgc31/code/viperproject/gseatest.csv"
#gsea_pdf = "/mnt/cfce-stor1/home/mgc31/code/viperproject/gseatest.pdf"

kegg_pathway_f<- function(deseq_file, numkeggpathways,kegg_dir,reference,temp_dir, kegg_table_up,kegg_table_down,keggsumary_pdf,up_kegg_png,down_kegg_png,gsea_table,gsea_pdf) {

    ## These are here until we update snakemake
    numkeggpathways = as.numeric(numkeggpathways)
    #keggpvalcutoff = as.numeric(keggpvalcutoff)

    ## Will need this path stuff for later as kegg output is very messy
    mainDir = substr(kegg_dir, 1, nchar(kegg_dir)-14)
    dir.create(file.path(mainDir, "kegg_pathways/"), showWarnings = FALSE)

    ## Read in deseq table
    detable = read.table(deseq_file, header=TRUE, sep=",", fill=TRUE)
    rownames(detable) <- detable[,1]

    ## Append ENSEMBL and ENTREZ IDs from loaded in db
    if (reference == "hg19") {IDdb = org.Hs.eg.db}
    if (reference == "mm9") {IDdb = org.Mm.eg.db}
        
    detable$entrez = mapIds(IDdb,
                        keys=row.names(detable),
                        column="ENTREZID",
                        keytype="SYMBOL",
                        multiVals="first")

    ## Couple failsafes
    detable = na.omit(detable)
    detable = detable[is.finite(detable$log2FoldChange),]
    
    ## Setting up gage input, needs the log2fc with the entrez id
    gageinput = detable$log2FoldChange
    names(gageinput) = detable$entrez

    ## Run gage
    keggres = gage(gageinput, gsets = kss, same.dir=TRUE)

    ## Output kegg results with respect to  upregulation
    kegg_up = keggres$greater
    kegg_up = cbind(rownames(kegg_up), kegg_up)
    colnames(kegg_up)[1] = "Kegg_pathway"
    xx = gsub(",","", as.matrix(kegg_up[,1]))
    kegg_up[,1] = xx
    write.table(kegg_up, file = kegg_table_up, quote=F, col.names=TRUE, row.names=FALSE, sep=",")

    ## Output kegg results with respect to  downregulation
    kegg_down= keggres$less
    kegg_down= cbind(rownames(kegg_down), kegg_down)
    colnames(kegg_down)[1] = "Kegg_pathway"
    xx = gsub(",","", as.matrix(kegg_down[,1]))
    kegg_down[,1] = xx
    write.table(kegg_down, file = kegg_table_down, quote=F, col.names=TRUE, row.names=FALSE, sep=",")

    ## Concat tables for filtering and testing
    full = rbind(kegg_up, kegg_down)
    fullorder = full[order(full[,4]),]
    fullkegg = fullorder[!duplicated(fullorder[,1]),]
        
    ## Stop run if the params don't match for output
    #kegg_output_filter = subset(fullkegg, fullkegg[,4] < keggpvalcutoff)
    #if(nrow(kegg_output_filter) < numkeggpathways) {stop(paste("Only ",nrow(kegg_output_filter), " pathways pass the current keggpvalcutoff of ", keggpvalcutoff, ", please run again with increased pval. Check comp.kegg.txt for details", sep="")) }

    ## Get the pathways
    keggrespathways = keggres$stats
    keggrespathways = keggrespathways[order(-abs(keggrespathways[,1])),]
    keggrespathways = rownames(keggrespathways)[1:numkeggpathways]
    keggresids = substr(keggrespathways, start=1, stop=8)

    ## Plot using pathview
    if (reference == "hg19") {orgparam = "hsa"}
    if (reference == "mm9") {orgparam = "mmu"}
    
    normwd = getwd()
    setwd(temp_dir)
    
    suppressMessages(for ( i in 1:numkeggpathways) {
        pvout <- pathview(gene.data=gageinput,              ## Gene list
                          pathway.id=keggresids[i],         ## Which pathway
                          species = orgparam,                  ## Species
                          #limit = list(gene=max(abs(gageinput)),cpd=1),
                          #kegg.dir = temp_dir               ## Save directory
                          low = list(gene = "blue", cpd = "purple"),
                          mid = list(gene = "gray", spd = "gray"),
                          high = list(gene = "red", cpd = "gray")  ## Color scale
                     )
    })
    setwd(normwd)

    ## Renaming files
    # Create variable with keggrespathways sorted and pull out name
    sortkeggrespathways = sort(keggrespathways)
    if (reference == "hg19") {
      newnames = substr(sortkeggrespathways, 10, nchar(sortkeggrespathways))
      newnames = gsub(" ", "_", newnames)
      newnames = paste0(newnames, "_", match(sortkeggrespathways,keggrespathways))
    }
    
    # Read in the list of made png files
    png_files <- list.files(temp_dir, pattern=glob2rx("*.png"))
    if (reference == "hg19") { file.rename(paste0(temp_dir,png_files), paste0(kegg_dir, newnames, ".png")) } else { file.rename(paste0(temp_dir, png_files), paste0(kegg_dir, png_files)) }
    # Repeat for xml files
    xml_files <- list.files(temp_dir, pattern=glob2rx("*.xml"))
    if (reference == "hg19") { file.rename(paste0(temp_dir,xml_files), paste0(kegg_dir, newnames, ".xml")) } else { file.rename(paste0(temp_dir, xml_files), paste0(kegg_dir, xml_files)) }

    ## Create KEGG Summary Tables for up and downregulated pathways
    plot_list = list()
    
    ## Create Up KEGG Summary Table
    upvalues = sapply(kegg_up[,5], as.numeric)
    uplogqval = -log(upvalues)
    upkeggsummary = data.frame(Keggpathway=substr(names(uplogqval), 10, nchar(names(uplogqval))), logqval=uplogqval)
    
    ## Create title for plot
    temptitle = tail(unlist(strsplit(keggsummary_pdf, split="/")), n=1)
    temptitle = head(unlist(strsplit(temptitle, split="[.]")), n=1)
    uptitle = paste(temptitle, "_Top_", numkeggpathways, "_Kegg_Pathways_UP", sep="")
    downtitle = paste(temptitle, "_Top_", numkeggpathways, "_Kegg_Pathways_DOWN", sep="")
    
    ## Create Up Plot
    kegg.term.width = 20
    upkegg.df <- as.data.frame(upkeggsummary[numkeggpathways:1,], stringsAsFactors=FALSE)
    to_fold <- which(lapply(as.character((upkegg.df$Keggpathway)), nchar) > kegg.term.width)

    foldednames <- sapply(upkegg.df$Keggpathway[to_fold],function(x) paste(strwrap(x,width = kegg.term.width),collapse = "\n"))
    levels(upkegg.df$Keggpathway) <- c(levels(upkegg.df$Keggpathway), foldednames)
    upkegg.df[to_fold,1] <- foldednames
    
    ggup <- ggplot(upkegg.df, aes(y=reorder(Keggpathway, logqval), x=logqval)) +
        geom_segment(aes(y=reorder(Keggpathway, logqval), yend=reorder(Keggpathway, logqval), x=0, xend=logqval)) +
        geom_point(aes(y=reorder(Keggpathway, logqval), x=logqval), color = "steelblue", size=5) +
        labs(x="- Log(Q-value)", y=NULL, title=uptitle) +
        theme(panel.grid.major.y=element_blank()) +
        theme(panel.grid.minor=element_blank()) +
        theme(axis.text.y=element_text(size=10)) +
        theme(plot.title = element_text(size=8, margin = margin(10, 0, 20, 0)))
    plot_list[[1]] = ggup                                                   
    

    ## Create Down KEGG Summary Table
    downvalues = sapply(kegg_down[,5], as.numeric)
    downlogqval = -log(downvalues)
    downkeggsummary = data.frame(Keggpathway=substr(names(downlogqval), 10, nchar(names(downlogqval))), logqval=downlogqval)

    ## Create Down Plot
    kegg.term.width = 20
    downkegg.df <- as.data.frame(downkeggsummary[numkeggpathways:1,], stringsAsFactors=FALSE)
    to_fold <- which(lapply(as.character((downkegg.df$Keggpathway)), nchar) > kegg.term.width)

    foldednames <- sapply(downkegg.df$Keggpathway[to_fold],function(x) paste(strwrap(x,width = kegg.term.width),collapse = "\n"))
    levels(downkegg.df$Keggpathway) <- c(levels(downkegg.df$Keggpathway), foldednames)
    downkegg.df[to_fold,1] <- foldednames

    ggdown <- ggplot(downkegg.df, aes(y=reorder(Keggpathway, logqval), x=logqval)) +
        geom_segment(aes(y=reorder(Keggpathway, logqval), yend=reorder(Keggpathway, logqval), x=0, xend=logqval)) +
        geom_point(aes(y=reorder(Keggpathway, logqval), x=logqval), color = "steelblue", size=5) +
        labs(x="- Log(Q-value)", y=NULL, title=downtitle) +
        theme(panel.grid.major.y=element_blank()) +
        theme(panel.grid.minor=element_blank()) +
        theme(axis.text.y=element_text(size=10)) +
        theme(plot.title = element_text(size=8, margin = margin(10, 0, 20, 0)))
    plot_list[[2]] = ggdown

    
    pdf(keggsummary_pdf)
    for (i in 1:2) { print(plot_list[[i]]) }
    junk <- dev.off()

    png(up_kegg_png, width = 8, height = 8, unit="in",res=300)
    print(plot_list[[1]])
    junk <- dev.off()

    png(down_kegg_png, width = 8, height = 8, unit="in",res=300)
    print(plot_list[[2]])
    junk <- dev.off()
    
    
    ## GSEA Analysis
    gseainput = sort(gageinput, decreasing=TRUE)
    gseainput = gseainput[is.finite(gseainput)]
    
    if (reference == "hg19") {orgparam = "hsa"}
    if (reference == "mm9") {orgparam = "mmu"}

    fullgsea <- gseKEGG(geneList = gseainput,
                        organism     = orgparam,
                        nPerm        = 100,
                        minGSSize    = 1,
                        pvalueCutoff = 0.99,
                        verbose      = FALSE,
                        use_internal_data = FALSE)

    gsea_data = summary(fullgsea)
    gsea_data = gsea_data[order(-abs(gsea_data$NES)),]
    xx = gsub(",","", as.matrix(gsea_data[,2]))
    gsea_data[,2] = xx
    write.table(gsea_data, file = gsea_table, quote=FALSE, sep= ",", row.names=FALSE, col.names=TRUE)

    pdf(gsea_pdf)
    plot.new()
    mtext("GSEA_plots")
    for ( i in 1:10 ) {
        gseaplot(fullgsea, geneSetID = gsea_data[i,1])
        mtext(gsea_data[i,2])
    }
    junk <- dev.off()
    
    }

args <- commandArgs( trailingOnly = TRUE )
deseq_file = args[1]
numkeggpathways = args[2]
kegg_dir = args[3]
reference = args[4]
temp_dir = args[5]
kegg_table_up = args[6]
kegg_table_down = args[7]
keggsummary_pdf = args[8]
up_kegg_png = args[9]
down_kegg_png = args[10]
gsea_table = args[11]
gsea_pdf = args[12]

species <- ifelse(reference == 'hg19', 'hsa', 'mmu')
kegg.set = kegg.gsets(species = species)
ks = kegg.set$kg.sets
kss = kegg.set$kg.sets[kegg.set$sigmet.idx]
names(kss) = gsub("/","",names(kss))

## Removing pathways that I know don't load properly... no idea why
pathway_errors = c("hsa01200 Carbon metabolism", "hsa01230 Biosynthesis of amino acids", "hsa01212 Fatty acid metabolism", "hsa01210 2-Oxocarboxylic acid metabolism", "hsa01100 Metabolic pathways", "hsa00533 Glycosaminoglycan biosynthesis - keratan sulfate", "hsa00514 Other types of O-glycan biosynthesis", "hsa00511 Other glycan degradation")
kss[which(names(kss) %in% pathway_errors)] <- NULL


kegg_pathway_f(deseq_file, numkeggpathways,kegg_dir,reference,temp_dir, kegg_table_up,kegg_table_down,keggsumary_pdf,up_kegg_png,down_kegg_png,gsea_table,gsea_pdf)
