
#generates a boxplot of CPK col in cdr3/CPK.csv
cpk_plot <- function(inputFile, plot_out) {
    png(file=plot_out)
    cpk_table = read.table(inputFile, header=TRUE, sep=",", row.names=1, check.names=FALSE)
    #print(cpk_table)
    #print(cpk_table['CPK'])
    boxplot(cpk_table['CPK'], ylab="CPK of TCR")
    junk <- dev.off()

}

args <- commandArgs( trailingOnly = TRUE )
cpkInput=args[1]
pngFile=args[2]
cpk_plot(cpkInput, pngFile)


