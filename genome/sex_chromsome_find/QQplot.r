#!/R-4.2/bin/Rscript

### choose libPaths ###
.libPaths("/R-4.2/lib/R/library/")

args <- commandArgs (T)
library(CMplot)
library(qqman)

results_log <- read.table(args[1], header=T)
p_value=results_log$P
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda

#pdf(args[2], width = 6, height = 6)

CMplot(results_log, plot.type = "q", threshold = 0.05, signal.cex=0.5, conf.int.col="grey", file="jpg", dpi=600, file.name=args[2], file.output=TRUE, verbose=F,cex=c(0.3,0.3))

#dev.off()
