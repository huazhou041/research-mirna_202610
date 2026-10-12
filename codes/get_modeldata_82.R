rm(list=ls())
library(readxl)
setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

outdir = 'data/model_data'
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

# ---- Read normalized expression (82 screening samples) ----
# Row 1 = sample IDs, row 2 = labels (Ctrl/Lung Ca), rows 3+ = data
# Cols 1-7 = annotation, col 1 = Probe Name, cols 8+ = samples
raw = read_excel('data/normalized data_screening set 82.xlsx', sheet=1, col_names=FALSE)

sample_ids  = unlist(raw[1, 8:ncol(raw)])   # full RCC filenames
labels      = unlist(raw[2, 8:ncol(raw)])   # "Ctrl" or "Lung Ca"
# Col 5 = Class; keep only Endogenous1 probes (excludes ligation, housekeeping,
# spike-in, positive, negative controls)
probe_class = unlist(raw[3:nrow(raw), 5])
keep        = !is.na(probe_class) & probe_class == 'Endogenous1'
probe_names = unlist(raw[3:nrow(raw), 1])[keep]
expr_mat    = raw[3:nrow(raw), 8:ncol(raw)][keep, ]
cat('Probes kept (Endogenous1):', sum(keep), 'of', length(keep), '\n')

# Numeric expression matrix: expr_mat is 828 miRNA rows x 82 sample cols; transpose to 82 x 828
expr_num = matrix(as.numeric(unlist(expr_mat)), nrow=nrow(expr_mat), ncol=ncol(expr_mat))
expr = data.frame(t(expr_num), stringsAsFactors=FALSE)

# Replace + with _ in probe names to match old EDRN naming convention
probe_names = gsub('[+]', '_', trimws(probe_names))
colnames(expr) = probe_names

# Extract numeric sample key: e.g. "B3_18.Ctrl_01.RCC" -> "18", "B3_189.LC_09.RCC" -> "189"
extract_key = function(s) regmatches(s, regexpr('[0-9]+(?=[.](Ctrl|LC)_)', s, perl=TRUE))
sample_keys = sapply(sample_ids, extract_key)
rownames(expr) = sample_keys

# Outcome: Ctrl=0, Lung Ca=1
recurrence = ifelse(trimws(labels) == 'Lung Ca', 1L, 0L)
names(recurrence) = sample_keys

cat('Expression matrix:', nrow(expr), 'samples x', ncol(expr), 'miRNAs\n')
cat('Outcome counts:\n'); print(table(recurrence))

# ---- miRNA-only RDS ----
data_mirna = cbind.data.frame(recurrence=recurrence, expr)

saveRDS(data_mirna, file.path(outdir, 'nyu_82_miRNA.RDS'))

cat('\nSaved:\n')
cat(' nyu_82_miRNA.RDS:', nrow(data_mirna), 'x', ncol(data_mirna), '\n')
