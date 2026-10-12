rm(list=ls())
library(readxl)
setwd('/Volumes/hua_mac/research/aris/harvey/mirna_202605')

outdir = 'data/model_data'
if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)

# ---- Read normalized expression (166 samples) ----
# Row 1 = sample IDs, row 2 = labels (Benign/Lu Ca), rows 3+ = data
# Cols 1-7 = annotation, col 1 = Probe Name, cols 8+ = samples
raw = read_excel('data/normalized data_validation set 166.xlsx', sheet=1, col_names=FALSE)

sample_ids = unlist(raw[1, 8:ncol(raw)])   # full RCC filenames
labels     = unlist(raw[2, 8:ncol(raw)])   # "Benign" or "Lu Ca"
# Col 5 = Class; keep only Endogenous1 probes (excludes ligation, housekeeping,
# spike-in, positive, negative controls)
probe_class = unlist(raw[3:nrow(raw), 5])
keep        = !is.na(probe_class) & probe_class == 'Endogenous1'
probe_names = unlist(raw[3:nrow(raw), 1])[keep]
expr_mat    = raw[3:nrow(raw), 8:ncol(raw)][keep, ]
cat('Probes kept (Endogenous1):', sum(keep), 'of', length(keep), '\n')

# Numeric expression matrix: expr_mat is 828 miRNA rows x 166 sample cols; transpose to 166 x 828
expr_num = matrix(as.numeric(unlist(expr_mat)), nrow=nrow(expr_mat), ncol=ncol(expr_mat))
expr = data.frame(t(expr_num), stringsAsFactors=FALSE)

# Replace + with _ in probe names to match old EDRN naming convention
probe_names = gsub('[+]', '_', trimws(probe_names))
colnames(expr) = probe_names

# Extract numeric sample key from RCC filename: e.g. "NYU2239" -> "2239"
extract_key = function(s) gsub('NYU', '', regmatches(s, regexpr('NYU[0-9]+', s, perl=TRUE)))
sample_keys = sapply(sample_ids, extract_key)

# NYU2208 and NYU3044 each have two time-point collections — keep both as separate rows
# make.unique appends .1 to the second occurrence: "2208", "2208.1"
rownames(expr) = make.unique(sample_keys)
sample_keys    = rownames(expr)

# Outcome: Benign=0, Lu Ca=1
recurrence = ifelse(trimws(labels) == 'Lu Ca', 1L, 0L)
names(recurrence) = sample_keys

cat('Expression matrix:', nrow(expr), 'samples x', ncol(expr), 'miRNAs\n')
cat('Outcome counts:\n'); print(table(recurrence))

# ---- Read clinical metadata ----
# miRNA_mapping_surgical_166-2.xlsx
# Col 1=Sample_ID, 13=Gender, 14=Race, 15=Smoking, 16=PackYears, 19=Age
map = read_excel('data/miRNA_mapping_surgical_166-2.xlsx', sheet=1, col_names=FALSE)
map_ids  = gsub('NYU', '', unlist(map[2:nrow(map), 1]))
map_data = data.frame(
  key     = map_ids,
  Gender  = unlist(map[2:nrow(map), 13]),
  Race    = unlist(map[2:nrow(map), 14]),
  Smoking = unlist(map[2:nrow(map), 15]),
  PackYrs = as.numeric(unlist(map[2:nrow(map), 16])),
  Age     = as.numeric(unlist(map[2:nrow(map), 19])),
  stringsAsFactors = FALSE
)

# Match metadata preserving timepoint order:
# make.unique gives "2208" for 1st collection and "2208.1" for 2nd.
# Use within-group occurrence rank so each timepoint gets its own metadata row.
lookup_keys = sub('[.]\\d+$', '', sample_keys)
occ_rank    = as.integer(ave(seq_along(lookup_keys), lookup_keys, FUN=seq_along))
idx = mapply(function(key, rank) {
  rows = which(map_data$key == key)
  if (length(rows) == 0L) return(NA_integer_)
  rows[min(rank, length(rows))]
}, lookup_keys, occ_rank)
meta = map_data[idx, ]
rownames(meta) = sample_keys

# Encode clinical variables (same scheme as previous EDRN pipeline)
# Gender: Female=0, Male=1, NA or unmatched=-9
sex = rep(-9L, nrow(meta))
sex[!is.na(meta$Gender) & meta$Gender == 'Female'] = 0L
sex[!is.na(meta$Gender) & meta$Gender == 'Male']   = 1L

# Race: White=5, Asian=4, Black=3, Hispanic=2, Other/Unknown/NA=1, unmatched=-9
race = rep(-9L, nrow(meta))
race[!is.na(meta$Race)] = 1L
race[!is.na(meta$Race) & grepl('Hispanic', meta$Race, ignore.case=TRUE)] = 2L
race[!is.na(meta$Race) & grepl('Black',    meta$Race, ignore.case=TRUE)] = 3L
race[!is.na(meta$Race) & grepl('Asian',    meta$Race, ignore.case=TRUE)] = 4L
race[!is.na(meta$Race) & grepl('White',    meta$Race, ignore.case=TRUE)] = 5L

# Smoking: Current=1, Former/Never=0, NA=-9
smoking = rep(-9L, nrow(meta))
smoking[!is.na(meta$Smoking) & meta$Smoking %in% c('Former','Never')] = 0L
smoking[!is.na(meta$Smoking) & meta$Smoking == 'Current']             = 1L

# Pack years: numeric, NA=-9
pack_yrs = ifelse(is.na(meta$PackYrs), -9, meta$PackYrs)

# Age: numeric, NA=-9
age = ifelse(is.na(meta$Age), -9, meta$Age)

clincial = data.frame(
  recurrence        = recurrence,
  Tissue_Sex        = sex,
  Tissue_Race       = race,
  Tissue_Smoking    = smoking,
  Tissue_Pack_years = pack_yrs,
  Tissue_Age        = age,
  row.names         = sample_keys,
  stringsAsFactors  = FALSE
)

cat('\nClinical data summary:\n')
cat('Samples matched to metadata:', sum(!is.na(idx)), 'of', length(sample_keys), '\n')
cat('Sex counts:\n');     print(table(clincial$Tissue_Sex))
cat('Race counts:\n');    print(table(clincial$Tissue_Race))
cat('Smoking counts:\n'); print(table(clincial$Tissue_Smoking))

# ---- Build and save RDS files ----
data_clinical = clincial
data_mirna    = cbind.data.frame(recurrence=recurrence, expr)
data_both     = cbind.data.frame(clincial, expr)

saveRDS(data_clinical, file.path(outdir, 'nyu_166_Clinical.RDS'))
saveRDS(data_mirna,    file.path(outdir, 'nyu_166_miRNA.RDS'))
saveRDS(data_both,     file.path(outdir, 'nyu_166_Both.RDS'))

cat('\nSaved:\n')
cat(' nyu_166_Clinical.RDS:', nrow(data_clinical), 'x', ncol(data_clinical), '\n')
cat(' nyu_166_miRNA.RDS:   ', nrow(data_mirna),    'x', ncol(data_mirna),    '\n')
cat(' nyu_166_Both.RDS:    ', nrow(data_both),     'x', ncol(data_both),     '\n')
