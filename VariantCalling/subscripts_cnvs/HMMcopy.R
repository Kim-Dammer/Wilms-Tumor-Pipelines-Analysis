require(HMMcopy); require(karyoploteR)
if (!require("HMMcopy", quietly = TRUE)) {
    stop("HMMcopy is not installed.")
}
if (!require("karyoploteR", quietly = TRUE)) {
    stop("karyoploteR is not installed.")
}

check_file <- function(file, name) {
    if (length(file) == 1) {
        return(file[1])
    } else if (length(file) == 0) {
        stop(parse0("no ", name, " ref file found"))
    } else {
        stop(parse0("multiple ", name, " ref files found"))
    }
}

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
refgenome <- args[2]
refgenomedir = file.path("/omics/groups/OE0132/internal/SequenceAlignment/Reference_Genome",
                         refgenome)
#readcounts
rfile_tumor <- file.path(outdir, "tumor_readcounts.wig")
# gc ref
gfile = list.files(refgenomedir, pattern="*gc.wig$", full.names=TRUE)
gfile <- check_file(gfile, "gc")
# mapping ref
mfile = list.files(refgenomedir, pattern="*fa.map.ws_1000.wig$", full.names=TRUE)
mfile <- check_file(mfile, "mapping")
tumor_reads <- wigsToRangedData(rfile_tumor, gfile, mfile)
# Correction of Copy Number Profile by Mapping, GC Ref Data
tumor_copy <- correctReadcount(tumor_reads)
tumor_copy$chr <- as.factor(tumor_copy$chr)
write.csv(data.frame(tumor_copy), sprintf("%s/TumorCopy.csv", outdir))

# plot CN Profile
set.seed(42)
tumor_copy2 = tumor_copy[sample(nrow(tumor_copy), size = .03 * nrow(tumor_copy), replace = FALSE), ]
## tumor CN karyo
png(paste0(outdir, "/CNTumorKaryo.png"), width = 2000, height = 300, res=150)
kp = plotKaryotype(plot.type=4, genome = "mm10")
kpDataBackground(kp, data.panel=1)
kpPoints(kp,
         chr = paste0("chr", tumor_copy2$chr),
         x = tumor_copy2$start,
         y = tumor_copy2$cor.map,
         ymin = 0,
         ymax = 2,
         data.panel = 1,
         col = "blue",
         cex = .15)
kpAbline(kp,
         chr = NULL,
         h = c(0.5, 1, 1.5),
         ymin = 0,
         ymax = 2,
         v = NULL,
         data.panel = 1,
         r0 = NULL,
         r1 = NULL,
         clipping = TRUE,
         lwd = 1,
         lty = 3,
         col = "black")
dev.off()
