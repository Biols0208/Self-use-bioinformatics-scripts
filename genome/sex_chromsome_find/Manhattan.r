#!/R-4.2/bin/Rscript

### choose libPaths ###
.libPaths("/R-4.2/lib/R/library/")

args <- commandArgs (T)
library(CMplot)

indata <- read.table(args[1], header=T, sep="\t", quote="")
#pdf(args[2], width = 15, height = 6)

CMplot(indata, plot.type="m", 
#        col=c("grey30","grey60"), 
        LOG10=T, 
        ylim=c(0,10), 
        cex=c(0.0001,0.0001),
        threshold=c(as.numeric(args[3]),as.numeric(args[4])), 
        threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c("black","grey"), 
        amplify=T, chr.den.col=NULL, 
        signal.col=c("red","green"), signal.cex=c(0.5,0.5), signal.pch=c(19,19),
        file="jpg", dpi=600, file.name=args[2], file.output=TRUE, verbose=F)

