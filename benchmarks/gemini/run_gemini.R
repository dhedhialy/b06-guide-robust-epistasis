## GEMINI (sellerslab/gemini, Bioc 1.24.0) on the COLO1 dual-guide screen.
## Input contract: counts.matrix rows = guide pairs (rownames "sg1;sg2"),
## columns = samples; guide.annotation carries the two per-row genes;
## sample.replicate.annotation maps columns to samples. Output: per gene
## pair the median "strong" combination score across the 3 endpoint samples
## and the minimum FDR, exported as gemini_scores.csv. Rerunnable, commit-free,
## no git metadata.

suppressMessages(library(gemini))

args <- commandArgs(trailingOnly = TRUE)
annotated <- args[1]
if (is.na(annotated)) stop("usage: run_gemini.R <COLO1_RUNMERGED_EXACT_ANNOTATED.txt>")
outdir <- ifelse(length(args) >= 2, args[2], "benchmarks/gemini")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

raw <- read.delim(annotated, stringsAsFactors = FALSE, check.names = FALSE)
raw <- raw[!duplicated(raw[c("sgRNA1_ID", "sgRNA2_ID")]), ]
raw$sgRNA1_ID <- trimws(raw$sgRNA1_ID)
raw$sgRNA2_ID <- trimws(raw$sgRNA2_ID)

# counts by guide pair; lib.COLO.1 = ETP (library) reference, the three
# CPID1020/23/26 columns = endpoint (HT29-selected) replicates.
pairid <- paste(raw$sgRNA1_ID, raw$sgRNA2_ID, sep = ";")
count_cn <- c("lib.COLO.1", "SIDM00136_CPID1020",
              "SIDM00136_CPID1023", "SIDM00136_CPID1026")
counts.matrix <- data.frame(lapply(raw[count_cn], as.numeric),
                            row.names = pairid, check.names = FALSE)
colnames(counts.matrix) <- count_cn
counts.matrix <- as.matrix(counts.matrix)
counts.matrix <- counts.matrix[rowSums(counts.matrix) > 0, , drop = FALSE]

guide.annotation <- data.frame(
  rowname = rownames(counts.matrix),
  U6.guide = sub(";.*", "", rownames(counts.matrix)),
  H1.guide = sub(".*;", "", rownames(counts.matrix)),
  U6.gene = raw$Gene1[match(rownames(counts.matrix), pairid)],
  H1.gene = raw$Gene2[match(rownames(counts.matrix), pairid)],
  stringsAsFactors = FALSE)

sample.replicate.annotation <- data.frame(
  colname = count_cn,
  samplename = c("lib", "HT29", "HT29", "HT29"),
  replicate = c(NA, "RepA", "RepB", "RepC"),
  stringsAsFactors = FALSE)

cat("counts dim:", dim(counts.matrix), "\n")
cat("NONTARGET00741 pairs:", sum(raw$Gene1 == "NONTARGET00741" |
                                   raw$Gene2 == "NONTARGET00741"), "\n")

Input <- gemini_create_input(
  counts.matrix = counts.matrix,
  sample.replicate.annotation = sample.replicate.annotation,
  guide.annotation = guide.annotation,
  sample.column.name = "samplename",
  gene.column.names = c("U6.gene", "H1.gene"),
  ETP.column = "lib.COLO.1",
  verbose = FALSE)

Input <- gemini_calculate_lfc(Input, normalize = TRUE, CONSTANT = 32)

Model <- gemini_initialize(
  Input = Input,
  nc_gene = "NONTARGET00741",
  pattern_join = ";",
  pattern_split = ";",
  cores = 15,
  verbose = FALSE)

Model <- gemini_inference(Model, n_iterations = 50,
                          force_results = TRUE, cores = 15, verbose = FALSE)

# FDR needs non-interacting pairs; use the NONTARGET-containing gene pairs.
srows <- rownames(Model$s)
nc_pairs <- srows[grepl("NONTARGET", srows, fixed = TRUE)]
cat("nc_pairs:", length(nc_pairs), "\n")
Score <- gemini_score(Model, pc_gene = NA, nc_pairs = nc_pairs)

strong <- Score$strong        # genes x samples; higher = stronger interaction
fdr <- Score$fdr_strong
strong_med <- apply(strong, 1, median, na.rm = TRUE)
fdr_min <- apply(fdr, 1, min, na.rm = TRUE)
out <- data.frame(gene_pair = names(strong_med),
                  gemini_strong_median = unname(strong_med),
                  gemini_fdr = unname(fdr_min),
                  stringsAsFactors = FALSE)
write.csv(out, file.path(outdir, "gemini_scores.csv"), row.names = FALSE)
cat("rows exported:", nrow(out), "\n")
cat("sessionInfo():\n")
print(sessionInfo())