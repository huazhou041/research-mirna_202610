#!/usr/bin/env Rscript
# Extended modeling script: GLM, Random Forest, XGBoost with unified bootstrap CV
# Runs all three learners on the same data split (no learner parameter)
library(argparser)

p = arg_parser("Run training pipeline: GLM, RF, and XGBoost")
p = add_argument(p, "drug", help = "Drug to learn (outcome column)")
p = add_argument(p, "matrix", help = "Matrix to train on (RDS file)")
p = add_argument(p, "hyperparams", help = "json file with hyperparameters")
p = add_argument(p, "prefix", help = "prefix to write files (base; _GLM/_RF/_XGB appended)")
p = add_argument(p, "--varcut", default = 100,
                 help = "Number of genes to keep based on variance.")
p = add_argument(p, "--nboot", default = 50,
                 help = "Number of bootstrap iterations.")
p = add_argument(p, "--kfold", default = 10L, help = "Number of CV folds")
p = add_argument(p, "--cores", default = "Slurm",
                 help = "Method to determine cores")

argv = parse_args(p)

suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(caret)
  library(pROC)
})

drug = argv$drug
mat_file = argv$matrix
hyperparam_file = argv$hyperparams
prefix = argv$prefix
varcut = argv$varcut
learners = c("GLM", "RF", "XGB")
kfold = argv$kfold
kboot = argv$nboot

Ncores = future::availableCores(method = argv$cores)
options(mc.cores = Ncores)

cat('*' %+% rep('*', 40) %+% '*\n')
cat('Running all learners on:', mat_file, '\n')
cat('Varcut:', varcut, '| Bootstrap:', kboot, '| Folds:', kfold, '\n')
cat('*' %+% rep('*', 40) %+% '*\n')

# ============================================================
# DATA PREP (same as tune_drug_hua.R)
# ============================================================
message(paste0('Reading ', mat_file, '...'))
mat = readRDS(mat_file)
mat[[drug]] = as.factor(mat[[drug]])

tissues = grep("^Tissue_", colnames(mat), value = TRUE)
genes = setdiff(colnames(mat), c(drug, tissues))

message('Subsetting matrix...')
mat = mat[!is.na(mat[[drug]]), ]
origmat = mat

set.seed(123)

# ============================================================
# RUN ALL LEARNERS ON SAME DATA SPLIT
# ============================================================
for (learner in learners) {
  cat('\n' %+% rep('=', 50) %+% '\n')
  cat('LEARNER:', learner, '\n')
  cat(rep('=', 50) %+% '\n\n')

  cvtotal = vector("list", kboot)
  cvtotal2 = vector("list", kboot)

  # ============================================================
  # BOOTSTRAP LOOP (100 iterations)
  # ============================================================
  for (tt in 1:kboot) {
    mat = origmat
    fold_id = rep(NA, nrow(mat))

    flds = createFolds(mat[[drug]], k = kfold, list = TRUE, returnTrain = FALSE)
    for (i in 1:kfold) {
      fold_id[flds[[i]]] = i
    }

    # Feature selection: top-N by MAD per fold
    if (tt == 1) message(sprintf('Selecting top %d features per fold', varcut))

    FS_top = function(mat_fold, top_n) {
      mad_vals = apply(mat_fold[, genes], 2, mad, na.rm = TRUE)
      top_idx = order(mad_vals, decreasing = TRUE)[1:min(top_n, length(genes))]
      selected_genes = genes[top_idx]
      c(tissues, selected_genes)
    }

    features = mclapply(1:kfold, function(i) {
      FS_top(mat[fold_id != i, ], varcut)
    })
    features[[kfold + 1L]] = FS_top(mat, varcut)

    # ============================================================
    # CV LOOP
    # ============================================================
    cvSummary = vector("list", kfold)
    cvSummary2 = data.frame(matrix(NaN, 0, 4))
    names(cvSummary2) = c('subj', 'foldid', 'bootid', 'prob')

    for (i in 1:kfold) {
      kmat = mat[, c(drug, features[[i]]), drop = FALSE]
      ix = which(fold_id == i)

      train_data = data.frame(kmat[-ix, ])
      valid_data = data.frame(kmat[ix, ])

      if (learner == "GLM") {
        # ========== GLM: elastic net ==========
        library(glmnet)

        hyper_params = jsonlite::fromJSON(hyperparam_file,
                                          simplifyVector = TRUE,
                                          simplifyDataFrame = FALSE,
                                          simplifyMatrix = FALSE)

        X_train = as.matrix(train_data[, -which(colnames(train_data) == drug)])
        y_train = as.numeric(train_data[[drug]]) - 1  # 0/1

        X_valid = as.matrix(valid_data[, -which(colnames(valid_data) == drug)])
        y_valid = as.numeric(valid_data[[drug]]) - 1

        kmodel = lapply(hyper_params, function(alpha) {
          tryCatch(
            cv.glmnet(X_train, y_train, family = "binomial", alpha = alpha, nfolds = 5),
            error = function(e) {
              cv.glmnet(X_train, y_train, family = "binomial", alpha = alpha,
                        lambda = exp(seq(log(0.001), log(5), length.out = 100)))
            }
          )
        })

        cvms = sapply(kmodel, function(m) m$cvm[m$lambda == m$lambda.min])
        imod = which.min(cvms)
        kmodel = kmodel[[imod]]

        pred_t = predict(kmodel, X_valid, s = 'lambda.min', type = 'response')
        pred_t = as.vector(pred_t)

        cauc = roc(y_valid, pred_t, quiet = TRUE)$auc
        kvalid = data.frame(AUC = cauc)

      } else if (learner == "RF") {
        # ========== Random Forest ==========
        library(randomForest)

        kmodel = randomForest(
          as.formula(paste0(drug, " ~ .")),
          data = train_data,
          ntree = 500,
          mtry = floor(sqrt(ncol(train_data) - 1)),
          importance = TRUE
        )

        pred_t = predict(kmodel, valid_data, type = 'prob')[, 2]
        y_valid = as.numeric(valid_data[[drug]]) - 1

        cauc = roc(y_valid, pred_t, quiet = TRUE)$auc
        kvalid = data.frame(AUC = cauc)

      } else if (learner == "XGB") {
        # ========== XGBoost ==========
        library(xgboost)

        X_train = as.matrix(train_data[, -which(colnames(train_data) == drug)])
        y_train = as.numeric(train_data[[drug]]) - 1

        X_valid = as.matrix(valid_data[, -which(colnames(valid_data) == drug)])
        y_valid = as.numeric(valid_data[[drug]]) - 1

        dtrain = xgb.DMatrix(X_train, label = y_train)
        dvalid = xgb.DMatrix(X_valid, label = y_valid)

        params = list(
          booster = "gbtree",
          objective = "binary:logistic",
          max_depth = 6,
          eta = 0.1,
          subsample = 0.8,
          colsample_bytree = 0.8
        )

        kmodel = xgb.train(
          params = params,
          data = dtrain,
          nrounds = 100,
          watchlist = list(eval = dvalid),
          verbose = 0
        )

        pred_t = predict(kmodel, X_valid)
        cauc = roc(y_valid, pred_t, quiet = TRUE)$auc
        kvalid = data.frame(AUC = cauc)
      }

      # Collect fold results
      kvalid$fold_id = i
      cvSummary[[i]] = kvalid

      # Collect predictions
      pred_t2 = data.frame(
        subj = rownames(valid_data),
        foldid = i,
        bootid = tt,
        prob = pred_t
      )
      cvSummary2 = rbind.data.frame(cvSummary2, pred_t2)
    }

    # Aggregate bootstrap results
    cvtotal[[tt]] = do.call(rbind, cvSummary)
    cvtotal[[tt]]$bootid = tt
    cvtotal2[[tt]] = cvSummary2
  }

  # ============================================================
  # OUTPUT FOR THIS LEARNER
  # ============================================================
  cvtotal = do.call(rbind, cvtotal)
  cvtotal2 = do.call(rbind, cvtotal2)

  prefix_learner = paste0(prefix, '_', learner)

  # Write both cvSummary and predictions to raw_pred subfolder
  raw_pred_dir = file.path(dirname(prefix), 'raw_pred')
  if (!dir.exists(raw_pred_dir)) dir.create(raw_pred_dir, recursive = TRUE)
  write.csv(cvtotal, file.path(raw_pred_dir, paste0(basename(prefix_learner), "_cvSummary.csv")), row.names = FALSE, quote = FALSE)
  write.csv(cvtotal2, file.path(raw_pred_dir, paste0(basename(prefix_learner), "_pred_prob.csv")), row.names = FALSE, quote = FALSE)

  cat('Completed:', prefix_learner, '\n')
  cat('Mean AUC:', round(mean(cvtotal$AUC, na.rm = TRUE), 3), '\n')
}
