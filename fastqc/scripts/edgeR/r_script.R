
# zz <- file("/researchers/swati.dawar/Analysis/210506_NB501056_0631_AHGNWGBGXJ/DEG/DEG_RESULTS/Casp2KO.GFP-Cas9.GFP_plots/r_script.err", open="wt")
# sink(zz)
# sink(zz, type="message")
source("https://www.bioconductor.org/biocLite.R")
biocLite(c('limma', 'edgeR'), ask = FALSE)

library("edgeR")
library("limma")

# read in matrix and groups
toc <- read.table("/data/htseq/output/Casp2KO.GFP-Cas9.GFP_HTSeq_Expression_Matrix.txt", sep="\t", comment="", as.is=T)
groups <- sapply(toc[1, -1], strsplit, ":")
for(i in 1:length(groups)) { g <- make.names(groups[[i]][2]); names(groups)[i] <- g; groups[[i]] <- groups[[i]][-2] }
colnames(toc) <- make.names(toc[2,])
toc[,1] <- gsub(",", ".", toc[,1])
tagnames <- toc[-(1:2), 1]
rownames(toc) <- toc[,1]
toc <- toc[-(1:2), -1]
for(i in colnames(toc)) toc[, i] <- as.numeric(toc[,i])
norm_factors <- calcNormFactors(as.matrix(toc))

pw_tests <- list()
uniq_groups <- unique(names(groups))
for(i in 1:(length(uniq_groups)-1)) for(j in (i+1):length(uniq_groups)) pw_tests[[length(pw_tests)+1]] <- c(uniq_groups[i], uniq_groups[j])
DGE <- DGEList(toc, lib.size=norm_factors*colSums(toc), group=names(groups))
# pdf("/data/edgeR/output/MA_plots_normalisation.pdf", width=14)
# for(i in 1:length(pw_tests)) {
# 	j <- c(which(names(groups) == pw_tests[[i]][1])[1], which(names(groups) == pw_tests[[i]][2])[1])
# 	par(mfrow = c(1, 2))
# 	maPlot(toc[, j[1]], toc[, j[2]], normalize = TRUE, pch = 19, cex = 0.2, ylim = c(-10, 10), main=paste("MA Plot", colnames(toc)[j[1]], "vs", colnames(toc)[j[2]]))
# 	grid(col = "blue")
# 	abline(h = log2(norm_factors[j[2]]), col = "red", lwd = 4)
# 	maPlot(DGE$counts[, j[1]]/DGE$samples$lib.size[j[1]], DGE$counts[, j[2]]/DGE$samples$lib.size[j[2]], normalize = FALSE, pch = 19, cex = 0.2, ylim = c(-8, 8), main=paste("MA Plot", colnames(toc)[j[1]], "vs", colnames(toc)[j[2]], "Normalised"))
# 	grid(col = "blue")
# }
# dev.off()
pdf(file="/data/edgeR/output/MDSplot.pdf")
plotMDS(DGE, main="MDS Plot", col=as.numeric(factor(names(groups)))+1, xlim=c(-3,3))
dev.off()
tested <- list()

group_fact <- factor(names(groups))
design <- model.matrix(~ -1 + group_fact)
colnames(design) <- sub("group_fact", "", colnames(design))

isexpr <- rowSums(cpm(toc)>1) >= 2
toc <- toc[isexpr, ]
pdf(file="/data/edgeR/output/LIMMA_voom.pdf")
y <- voom(toc, design, plot=TRUE, lib.size=colSums(toc)*norm_factors)
dev.off()

pdf(file="/data/edgeR/output/LIMMA_MDS_plot.pdf")
plotMDS(y, labels=colnames(toc), col=as.numeric(factor(names(groups)))+1, gene.selection="common")
dev.off()
fit <- lmFit(y, design)

tab <- data.frame(Feature=rownames(y$E), y$E, stringsAsFactors=F)

write.table(tab, "/data/edgeR/output/Casp2KO.GFP-Cas9.GFP_log_dig_matrix.txt", quote=F, sep="\t", row.names=F)

cont <- c("G1-G2")
for(i in uniq_groups)  cont <- gsub(paste(groups[[i]], "([^0-9])", sep=""), paste(i, "\\1", sep=""), cont)
for(i in uniq_groups)  cont <- gsub(paste(groups[[i]], "$", sep=""), i, cont)

cont <- makeContrasts(contrasts=cont, levels=design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

tab <- NULL
options(digits = 6)
for(i in colnames(fit2)) {
	tab_tmp <- topTable(fit2, coef=i, n=Inf, sort.by="none", adjust.method="BH")
	colnames(tab_tmp) <- paste(i, colnames(tab_tmp), sep=":")
	if(is.null(tab)) {
		tab <- cbind(Feature=rownames(tab_tmp), tab_tmp)
	} else tab <- cbind(tab, tab_tmp[tab[,1],])
}

write.table(tab, "/data/edgeR/output/Casp2KO.GFP-Cas9.GFP_results.tsv", quote=F, sep="\t", row.names=F)
# sink(type="message")
# sink()
