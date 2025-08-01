library(VariantAnnotation)
library(data.table)
library(ggplot2)
library(ggVennDiagram)
library(patchwork)

make_table <- function(filename) {
    vcf <- readVcf(filename)

    gt <- as.data.table(geno(vcf)$GT)
    fixed.names <- colnames(gt)
    fixed.names <- gsub("[0-9]+_", "", fixed.names)
    colnames(gt) <- fixed.names
    x <- cbind(data.table(chrom = as.character(seqnames(vcf)),
                          start = start(vcf),
                          vid = rownames(vcf)),
                gt)

               
    x.long <- melt(x, id.vars=c("chrom", "start", "vid"))[value != "0/0" & value != "0|0"]

    x.long
}

make_venn <- function(filename, name) {
    x.long <- make_table(filename)

    x.ids <- split(x.long$vid, x.long$variable)

    w <- 12
    g <- ggVennDiagram(x.ids) + guides(fill = "none") +
         ggtitle(name) +
         coord_cartesian(clip = 'off') +
         theme(plot.margin=unit(c(w, w, w, w), "pt"))
    g
}

p <- make_venn("results/out.svelt.vcf", "Svelt") +
     make_venn("results/out.jasmine.vcf", "Jasmine")

pdf("results/merge-venn.pdf", width=10, height=5.8)
print(p)
dev.off()