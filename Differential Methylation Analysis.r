# ============================================================
# Differential Methylation Analysis with Permutation-Derived
# Significance Thresholds for C9orf72
# Files needed: normalized.Rdata, mSetSqFlt.Rdata, cell-type-proportions.csv, C9_DNAmFTD_lengths.txt
# These files were not provided on GitHub due to concerns about sharing patient data.
# ============================================================

# -----------------------------
# Libraries
# -----------------------------
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)   # EPICv2 annotation
library(IlluminaHumanMethylationEPICv2manifest)         # EPICv2 manifest
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(ggplot2)
library(DMRcate)
library(Gviz)
library(EpiSmokEr)                                      # Smoking score (SSc)
library(qqman)                                          # Manhattan/QQ
library("CONFINED")
library(ggrepel)                                        # Non-overlapping text labels
library(topr)                                           # Manhattan helper

# -----------------------------
# Load normalized objects
# -----------------------------
load("normalized.Rdata")
load("mSetSqFlt.Rdata")

# -----------------------------
# Compute M-values / Beta-values
# -----------------------------
mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

# Harmonize probe IDs (drop last 5 characters) to match annotation
rownames(mVals) <- substr(rownames(mVals), 1, nchar(rownames(mVals)) - 5)
rownames(bVals) <- substr(rownames(bVals), 1, nchar(rownames(bVals)) - 5)

# Fix diagnosis label per original code
targets$diagnosis[targets$diagnosis == "Subjectieve klachten"] <- "SC"

# -----------------------------
# Probe annotation aligned to mVals
# -----------------------------
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno$Name <- substr(anno$Name, 1, nchar(anno$Name) - 5)
anno <- anno[match(rownames(mVals), anno$Name), ]

# Quick peek (sanity check)
head(mVals[1:10, 1:10])

# -----------------------------
# Figure sizing helper (for notebooks)
# -----------------------------
fig <- function(width, heigth) {
  options(repr.plot.width = width, repr.plot.height = heigth)
}

# ============================================================
# Cell-type proportions and target data merge
# ============================================================
proportions <- read.csv("cell-type-proportions.csv", header = TRUE, row.names = 1)
head(proportions)

# Match to targets via Basename
proportions <- proportions[targets$Basename, ]
targets <- cbind(targets, proportions)
rownames(targets) <- targets$Basename
head(targets)

# ============================================================
# Identify individuals with missing C9orf72 repeat lengths
# ============================================================
repeat.lengths <- read.table("C9_DNAmFTD_lengths.txt", header = TRUE, sep = "\t")

indices.missing.lengths <- is.na(repeat.lengths$C9ORF1_harmonized) &
  is.na(repeat.lengths$C9ORF2_harmonized)
missing.lengths <- paste(repeat.lengths$Sentrix_ID[indices.missing.lengths],
                         repeat.lengths$Sentrix_Position[indices.missing.lengths],
                         sep = "_")

# Restrict matrices and targets to samples with known lengths
c9orf.mVals   <- mVals[, !(colnames(mVals) %in% missing.lengths)]
c9orf.bVals   <- bVals[, !(colnames(mVals) %in% missing.lengths)]
c9orf.targets <- targets[!(targets$Basename %in% missing.lengths), ]
c9orf <- factor(c9orf.targets$interpretation)

# Covariates used in the design matrix
sex   <- as.numeric(factor(c9orf.targets$sex))
plate <- as.numeric(factor(c9orf.targets$Sample_Plate))

# EpiSmokEr smoking score (SSc) from beta-values
result_SSc <- epismoker(dataset = c9orf.bVals, method = "SSc")
smoking <- result_SSc$smokingScore

# -----------------------------
# Design matrix 
#   Model: methylation ~ C9orf72 group + covariates (cell fractions)
# -----------------------------
design <- model.matrix(
  ~0 + c9orf + sex + age + smoking + plate + CD8T + CD4T + NK + Bcell + Mono + Neu,
  data = c9orf.targets
)

# ============================================================
# Differential methylation with limma
# ============================================================
fit <- lmFit(c9orf.mVals, design)
contMatrix <- makeContrasts(c9orfexpanded - c9orfno, levels = design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# Summary of DE status across probes
summary(decideTests(fit2))

# Top table of DMPs (no annotation)
results <- topTable(fit2, number = Inf, coef = 1)

# Re-attach annotation (as in original)
anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno$Name <- substr(anno$Name, 1, nchar(anno$Name) - 5)
anno <- anno[match(rownames(c9orf.mVals), anno$Name), ]
DMPs <- topTable(fit2, num = Inf, coef = 1, genelist = anno)

# ------------------------------------------------------------
# Genomic inflation factor (lambda) and QQ
# ------------------------------------------------------------
DMPs$CHR <- as.numeric(gsub(pattern = 'chr', replacement = '', x = DMPs$chr))
chisq <- qchisq(1 - DMPs$P.Value, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

# QQ plot of P-values
fig(6, 5)
qq(DMPs$P.Value)

# Density of observed -log(p) (natural log per original)
fig(5, 3)
ggplot(DMPs, aes(x = -log(P.Value))) +
  geom_density(fill = "gray") +
  labs(title = "", x = "Observed -log(p)", y = "Density") +
  theme_classic()

# ============================================================
# Manhattan plot (qqman/topr formatting)
# ============================================================
DMPs.rearranged <- cbind(
  chrom = substr(DMPs$chr, 4, nchar(DMPs$chr)),  # drop "chr"
  pos   = DMPs$pos,
  rsid  = DMPs$Name,
  p     = DMPs$P.Val,
  B     = DMPs$B
)
DMPs.rearranged <- as.data.frame(DMPs.rearranged)
DMPs.rearranged$chrom <- as.numeric(DMPs.rearranged$chrom)
DMPs.rearranged$pos   <- as.numeric(DMPs.rearranged$pos)
DMPs.rearranged$p     <- as.numeric(DMPs.rearranged$p)
DMPs.rearranged$B     <- as.numeric(DMPs.rearranged$B)

head(DMPs.rearranged)

fig(10, 5)
manhattan(DMPs.rearranged, annotate = 5e-9, sign_thresh = 4.4e-10)

# ============================================================
# Locus zoom around C9orf72 (if EnsDb available)
# ============================================================
fig(6, 6)
library(locuszoomr)
if (require(EnsDb.Hsapiens.v86)) {
  loc <- locus(
    data   = DMPs.rearranged,
    gene   = 'C9orf72',
    flank  = 5e4,
    ens_db = "EnsDb.Hsapiens.v86"
  )
  summary(loc)
  locus_plot(loc)
}

# ============================================================
# Volcano plot with labels (threshold 4.4e-10; natural -log)
# ============================================================
DMPs$Point.Color <- ifelse(DMPs$P.Value > 4.4e-10, 'gray',
                           ifelse(DMPs$logFC > 0, "blue", "red"))
DMPs$Label <- NA

# Gene label: first symbol before ';'; default "CHID1" (as in original)
split_vec <- strsplit(DMPs$UCSC_RefGene_Name[DMPs$P.Value <= 4.4e-10], ";")
gene <- sapply(split_vec, `[`, 1)
gene[is.na(gene)] <- "CHID1"
DMPs$Label[DMPs$P.Value <= 4.4e-10] <-
  paste(gene, DMPs$Name[DMPs$P.Value <= 4.4e-10], sep = "\n")

fig(10, 5)
ggplot(data = DMPs, aes(x = logFC, y = -log(P.Value), label = Label)) +
  geom_point(color = DMPs$Point.Color, size = 2) +
  geom_hline(yintercept = -log(4.4e-10), color = "red", linetype = "dashed") +
  geom_label_repel(box.padding = 3, max.overlaps = Inf,
                   min.segment.length = 0, size = 3.5, point.padding = NA) +
  annotate("text",
           x = min(DMPs$logFC),
           y = -log(4.4e-10) + 0.6,
           label = "4.4e-10", hjust = 0.8, vjust = 0, color = "red") +
  theme_classic()

# ============================================================
# Permutation loop (100 permutations) to derive thresholds
# ============================================================
orig.design   <- design
permuted.DMPs <- list()

for (i in 1:100) {
  set.seed(i)
  
  # Randomly select 27 samples to be the "expanded" group
  sampled.expanded <- sample(x = rownames(orig.design), size = 27, replace = FALSE)
  
  # Copy the original design and flip group indicators
  this.design <- orig.design
  this.design[which(rownames(this.design) %in% sampled.expanded), 1] <- 1
  this.design[which(rownames(this.design) %in% sampled.expanded), 2] <- 0
  
  this.design[which(!(rownames(this.design) %in% sampled.expanded)), 1] <- 0
  this.design[which(!(rownames(this.design) %in% sampled.expanded)), 2] <- 1
  
  # Refit limma with permuted labels
  fit <- lmFit(c9orf.mVals, this.design)
  contMatrix <- makeContrasts(c9orfexpanded - c9orfno, levels = this.design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  
  summary(decideTests(fit2))
  
  results <- topTable(fit2, number = Inf, coef = 1)
  
  # Attach annotation as done above
  anno <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  anno$Name <- substr(anno$Name, 1, nchar(anno$Name) - 5)
  anno <- anno[match(rownames(c9orf.mVals), anno$Name), ]
  permuted.DMPs[[i]] <- topTable(fit2, num = Inf, coef = 1, genelist = anno)
  
  cat("Permutation ", i, "\n")
}

# Save permutation results object
save(permuted.DMPs, file = "permuted.ewas.results.Rdata")

# ============================================================
# Post-permutation threshold assessment (three approaches)
# ============================================================

# 1) Count permutations with hits at threshold = 3.984926e-08
threshold <- 3.984926e-08
permutations <- as.data.frame(matrix(nrow = 100, ncol = 2))
permutations[, 1] <- 1:100
for (i in 1:100) {
  this.DMP <- permuted.DMPs[[i]]
  permutations[i, 2] <- sum(this.DMP$P.Value <= threshold)
}
sum(permutations[, 2] > 0)  # number of permutations with ≥1 hit

# 2) Repeat for threshold = 3.397500e-10
threshold <- 3.397500e-10
permutations <- as.data.frame(matrix(nrow = 100, ncol = 2))
permutations[, 1] <- 1:100
for (i in 1:100) {
  this.DMP <- permuted.DMPs[[i]]
  permutations[i, 2] <- sum(this.DMP$P.Value <= threshold)
}
sum(permutations[, 2] > 0)

# 3) Empirical minimum P across permutations
threshold <- 1
for (i in 1:100) {
  this.DMP <- permuted.DMPs[[i]]
  if (min(this.DMP$P.Value) < threshold) {
    threshold <- min(this.DMP$P.Value)
  }
}
threshold  # empirical min P across permutations

# ============================================================
# Manhattan with extra reference lines
# ============================================================
fig(6, 6)
p <- manhattan(DMPs.rearranged, annotate = 5e-9, sign_thresh = 9e-8)
p + geom_hline(yintercept = 9.35, color = "blue", linetype = "dashed") +
  annotate("text", x = 1, y = 9.35 + 0.6, label = "4.4e-10",
           hjust = 0, vjust = 0, color = "blue")

# ============================================================
# Volcano with revised thresholds (base-10 scale here)
# ============================================================
DMPs$Point.Color <- ifelse(DMPs$P.Value > 4.37039646887794e-10, 'gray',
                           ifelse(DMPs$logFC > 0, "blue", "red"))
DMPs$Label <- NA

split_vec <- strsplit(DMPs$UCSC_RefGene_Name[DMPs$P.Value <= 4.37039646887794e-10], ";")
gene <- sapply(split_vec, `[`, 1)
gene[is.na(gene)] <- "CHID1"

DMPs$Label[DMPs$P.Value <= 4.37039646887794e-10] <-
  paste(gene, DMPs$Name[DMPs$P.Value <= 4.37039646887794e-10], sep = "\n")

fig(6, 6)
ggplot(data = DMPs, aes(x = logFC, y = -log10(P.Value), label = Label)) +
  geom_point(color = DMPs$Point.Color, size = 2) +
  geom_hline(yintercept = -log10(4.37039646887794e-10), color = "blue", linetype = "dashed") +
  annotate("text",
           x = min(DMPs$logFC),
           y = -log10(4.37039646887794e-10) + 0.6,
           label = "4.4e-10", hjust = 0.3, vjust = 0, color = "blue") +
  geom_hline(yintercept = -log10(9e-8), color = "red", linetype = "dashed") +
  annotate("text",
           x = min(DMPs$logFC),
           y = -log10(9e-8) + 0.6,
           label = "9e-08", hjust = 0.4, vjust = 0, color = "red") +
  geom_label_repel(box.padding = 3, max.overlaps = Inf,
                   min.segment.length = 0, size = 3.5, point.padding = NA) +
  theme_classic()