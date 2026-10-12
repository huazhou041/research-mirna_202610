#!/usr/bin/env Rscript
library(argparser)

p = arg_parser("Train one learner on combined miRNA data with bootstrap CV and final saved model")
p = add_argument(p, "matrix", help="Input RDS matrix")
p = add_argument(p, "prefix", help="Output prefix")
p = add_argument(p, "--learner", default="GLM", help="Learner: GLM, RF, or XGB")
p = add_argument(p, "--varcut", default=100, help="Number of MAD-ranked features")
p = add_argument(p, "--nboot", default=100, help="Number of bootstrap repeats")
p = add_argument(p, "--kfold", default=10, help="Number of CV folds")
p = add_argument(p, "--cores", default="Slurm", help="future::availableCores method")
argv = parse_args(p)

suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(caret)
  library(pROC)
})

mat = readRDS(argv$matrix)
drug = "recurrence"
mat[[drug]] = as.factor(mat[[drug]])
mat = mat[!is.na(mat[[drug]]), , drop=FALSE]
batch_cols = grep("^batch_", colnames(mat), value=TRUE)   # cohort indicators (always kept; unpenalized in GLM)
genes = setdiff(colnames(mat), c(drug, batch_cols))
origmat = mat

learner = toupper(argv$learner)
varcut = as.integer(argv$varcut)
kboot = as.integer(argv$nboot)
kfold = as.integer(argv$kfold)
Ncores = future::availableCores(method=argv$cores)
options(mc.cores=Ncores)
set.seed(123)

# GLM alpha grid: use the same hyperparams.json grid as the per-source/transfer models
glm_alphas = tryCatch(jsonlite::fromJSON("codes/hyperparams.json", simplifyVector=TRUE),
                      error=function(e) c(0, 0.25, 0.5, 0.75, 1))

FS_top = function(dat, top_n) {
  mad_vals = apply(dat[, genes, drop=FALSE], 2, mad, na.rm=TRUE)
  top = genes[order(mad_vals, decreasing=TRUE)[seq_len(min(top_n, length(mad_vals)))]]
  c(batch_cols, top)   # always include cohort batch covariates alongside top-MAD genes
}

train_predict = function(learner, train_data, valid_data, features) {
  y_train = as.numeric(train_data[[drug]]) - 1
  y_valid = as.numeric(valid_data[[drug]]) - 1

  if (learner == "GLM") {
    library(glmnet)
    alphas = glm_alphas
    x_train = as.matrix(train_data[, features, drop=FALSE])
    x_valid = as.matrix(valid_data[, features, drop=FALSE])
    pf = ifelse(features %in% batch_cols, 0, 1)   # batch covariates unpenalized (always retained)
    models = lapply(alphas, function(a) {
      tryCatch(cv.glmnet(x_train, y_train, family="binomial", alpha=a, nfolds=5, penalty.factor=pf),
               error=function(e) cv.glmnet(x_train, y_train, family="binomial", alpha=a, penalty.factor=pf,
                                           lambda=exp(seq(log(0.001), log(5), length.out=100))))
    })
    cvms = sapply(models, function(m) m$cvm[m$lambda == m$lambda.min])
    model = models[[which.min(cvms)]]
    pred = as.vector(predict(model, x_valid, s="lambda.min", type="response"))
  } else if (learner == "RF") {
    library(randomForest)
    x_train = train_data[, features, drop=FALSE]
    x_valid = valid_data[, features, drop=FALSE]
    model = randomForest(x=x_train,
                         y=train_data[[drug]],
                         ntree=500,
                         mtry=max(1, floor(sqrt(length(features)))),
                         importance=TRUE)
    pred = predict(model, x_valid, type="prob")[, 2]
  } else if (learner == "XGB") {
    library(xgboost)
    x_train = as.matrix(train_data[, features, drop=FALSE])
    x_valid = as.matrix(valid_data[, features, drop=FALSE])
    dtrain = xgb.DMatrix(x_train, label=y_train)
    params = list(booster="gbtree", objective="binary:logistic", eval_metric="auc",
                  max_depth=6, eta=0.1, subsample=0.8, colsample_bytree=0.8)
    model = xgb.train(params=params, data=dtrain, nrounds=100, verbose=0)
    pred = predict(model, x_valid)
  } else {
    stop("Unknown learner: ", learner)
  }

  list(model=model, pred=pred, auc=as.numeric(roc(y_valid, pred, quiet=TRUE, direction = "<")$auc))  # direction="<": avoid autoflip inflating AUC on anti-correlated folds
}

cvtotal = vector("list", kboot)
cvpred = vector("list", kboot)

message(sprintf("[%s] START %s | learner=%s varcut=%d n=%d genes=%d batch=%d | nboot=%d kfold=%d",
                format(Sys.time(), "%H:%M:%S"), basename(argv$matrix), learner, varcut,
                nrow(origmat), length(genes), length(batch_cols), kboot, kfold))
.t0 = Sys.time()
for (tt in seq_len(kboot)) {
  .tt0 = Sys.time()
  mat = origmat
  flds = createFolds(mat[[drug]], k=kfold, list=TRUE, returnTrain=FALSE)
  fold_id = rep(NA_integer_, nrow(mat))
  for (i in seq_len(kfold)) fold_id[flds[[i]]] = i

  features = mclapply(seq_len(kfold), function(i) FS_top(mat[fold_id != i, , drop=FALSE], varcut))
  features[[kfold + 1L]] = FS_top(mat, varcut)

  fold_rows = vector("list", kfold)
  pred_rows = vector("list", kfold)
  for (i in seq_len(kfold)) {
    train_data = mat[fold_id != i, , drop=FALSE]
    valid_data = mat[fold_id == i, , drop=FALSE]
    fit = train_predict(learner, train_data, valid_data, features[[i]])
    fold_rows[[i]] = data.frame(AUC=fit$auc, fold_id=i, bootid=tt)
    pred_rows[[i]] = data.frame(subj=rownames(valid_data), foldid=i, bootid=tt, prob=fit$pred)
  }
  cvtotal[[tt]] = bind_rows(fold_rows)
  cvpred[[tt]] = bind_rows(pred_rows)
  message(sprintf("[%s] boot %3d/%d | %.1fs/boot | mean fold AUC=%.3f | total %.1f min",
                  format(Sys.time(), "%H:%M:%S"), tt, kboot,
                  as.numeric(difftime(Sys.time(), .tt0, units = "secs")),
                  mean(cvtotal[[tt]]$AUC, na.rm = TRUE),
                  as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
}

cvtotal = bind_rows(cvtotal)
cvpred = bind_rows(cvpred)
dir.create(dirname(argv$prefix), recursive=TRUE, showWarnings=FALSE)
dir.create("raw_pred", recursive=TRUE, showWarnings=FALSE)
base = basename(argv$prefix)
write.csv(cvtotal, file.path("raw_pred", sprintf("%s_%s_v%d_cvSummary.csv", base, learner, varcut)), row.names=FALSE, quote=FALSE)
write.csv(cvpred, file.path("raw_pred", sprintf("%s_%s_v%d_pred_prob.csv", base, learner, varcut)), row.names=FALSE, quote=FALSE)

final_features = FS_top(origmat, varcut)
final_fit = train_predict(learner, origmat, origmat, final_features)
final_obj = list(model=final_fit$model,
                 learner=learner,
                 varcut=varcut,
                 features=final_features,
                 training_n=nrow(origmat),
                 apparent_auc=final_fit$auc,
                 created=Sys.time())
saveRDS(final_obj, file.path("res", sprintf("%s_%s_v%d.rds", base, learner, varcut)))

cat(sprintf("Completed %s v%d: mean CV AUC %.3f, apparent AUC %.3f\n",
            learner, varcut, mean(cvtotal$AUC, na.rm=TRUE), final_fit$auc))
