#!/usr/bin/env Rscript
# inputs    -------------------------------------------------------------------
library(argparser)

p = arg_parser("Run training pipeline")
p = add_argument(p, "drug", help = "Drug to learn")
p = add_argument(p, "matrix", help = "Matrix to train on (RDS file)")
p = add_argument(p, "hyperparams", help = "json file with hyper parameters")
p = add_argument(p, "prefix", help = "prefix to write files")
p = add_argument(p, "--max_models", default = 15L,
                 help = "Max number of models to compute")
p = add_argument(p, "--varcut", default = 100,
                 help = "Percent of genes to keep based on variance.")
p = add_argument(p, "--nboot", default = 50,
                 help = "Number of bootstrap.")
p = add_argument(p, "--learner", default = "DL",
                 help = "Algorithm to use for learning (DL, RF, GLM)")
p = add_argument(p, "--test", help = "Path to clinical data")
p = add_argument(p, "--kfold", default = 5L, help = "Number of CV-splits")
p = add_argument(p, "--max_mem", default = "64G", help = "H2O parameter")
p = add_argument(p, "--cores", default = "Slurm",
                 help = "Method to determine cores with future::availableCores")

argv = parse_args(p)
# argv = parse_args(p, c("Bortezomib", "bortezomib_train.rds", "hyperparam_grid_dl.json",
#                        "debug", "--max_models", "2", "--varcut", "50"))

message("Running with the following arguments")
# str(argv)

drug = argv$drug
mat_file = argv$matrix
lift_null_file = argv$lift
hyperparam_file = argv$hyperparams
prefix = argv$prefix
Nmodels = argv$max_models
varcut = argv$varcut
learner = argv$learner
kfold = argv$kfold
kboot=argv$nboot
max_mem = argv$max_mem
test_file = argv$test
cores_method = argv$cores

stopifnot(learner %in% c("DL", "RF", "GLM"))

Ncores = future::availableCores(method = cores_method)
options(mc.cores = Ncores)
cat('*****************\n')
cat(learner,' : ' ,mat_file,'\n')
cat('*****************\n')

# Library    -------------------------------------------------------------------

suppressPackageStartupMessages({
    library(tidyverse)
    library(parallel)
    library(glmnet)
    library(glmnetUtils)
    library(h2o)
	library(pROC)
  library(caret)
})

all_equal = function(x, y) {
    if (length(y) == 1)
        y = rep(y, length(x))
    !is.character(all.equal(x, y, check.attributes = FALSE))
}

is_scaled = function(mat) {
    # result of `stat::scale`
    if (!is.null(attr(mat, "scaled:center")) )
        return(TRUE)
    # zero-mean
    cmeans = colMeans(mat, na.rm = TRUE)
    if (!all_equal(cmeans, 0))
        return(FALSE)
    # unit variance
    csd = apply(mat, 2, sd)
    if (!all_equal(csd, 1))
        return(FALSE)

    return(TRUE)
}

rbind_data_frame = function(...)
    # simple wrapper
    rbind.data.frame(..., make.row.names = FALSE, stringsAsFactors = FALSE)

delist_col = function(x, collapse = '|') {
    if (!is.list(x))
        return(x)
    vapply(x, paste, '', sep = '', collapse = collapse)
}

delist_df = function(x, collapse = '|')
    as.data.frame(lapply(x, delist_col, collapse = collapse))

FS_top = function(mat, top) {
    ix = apply(mat[, genes], 2, mad, na.rm = TRUE)
    if (length(ix)<top)
	top=length(ix)
    ix = order(ix, decreasing = TRUE)[1:top]
    genes = genes[ix]
    #if (!is.na(genes))
if (length(genes)>1)
    c(tissues,genes)
	else
	tissues
}

expand_hyperparams_h2o = function(hyper_params) {
    hyper_params = split(hyper_params, seq(nrow(hyper_params)))
    hyper_params = lapply(hyper_params, as.list)
    names(hyper_params) = NULL
    for (i in seq(length(hyper_params))) {
        for (j in seq(length(hyper_params[[i]]))) {
            if (is.list(hyper_params[[i]][[j]]))
                hyper_params[[i]][[j]] = unlist(hyper_params[[i]][[j]],
                                                recursive = FALSE,
                                                use.names = FALSE)
        }
    }
    hyper_params
}

validate_h2o = function(model, newdata) {
    y_hat = h2o.performance(model, newdata)
    data.frame(MSE = h2o.mse(y_hat), AUC = h2o.auc(y_hat))
}

validate_glmnet = function(model, newdata, s = 'lambda.min') {
  idrug = which(colnames(newdata) == drug)
  #//y_hat = predict(model, as(newdata[, -idrug],'dgCMatrix'), s = s)
  y_hat = predict(model, newdata[, -idrug], s = s)
  print(idrug)
  dy = y_hat - newdata[, idrug]
  d1=roc(as.numeric(newdata[,idrug]),as.numeric(y_hat))
  data.frame(lambda = model[[s]],  # hacky
             MSE = mean(dy * dy), AUC =d1$auc)
}

# Preprocessing -----------------------------------------------------------------

message(paste0('Sourcing ', hyperparam_file, "..."))
hyper_params <- jsonlite::fromJSON(hyperparam_file,
                                   simplifyVector = TRUE,
                                   simplifyDataFrame = FALSE,
                                   simplifyMatrix = FALSE)

message(paste0('Reading ', mat_file, '...'))
#mat = as.matrix(readRDS(mat_file))
mat = readRDS(mat_file)
mat$recurrence=as.factor(mat$recurrence) ##hua zhou 2020.03.13
tissues = grep("^Tissue_", colnames(mat), value = TRUE)
if (any(grepl("_GE", colnames(mat)))) {
    genes = grep("_GE$", colnames(mat), value = TRUE)
} else {
    genes = setdiff(colnames(mat), c(drug, tissues))
}

message('Subseting matrix...')
mat = mat[!is.na(mat[, drug]),  ]
origmat=mat

set.seed(123)

cvtotal=vector("list", kboot)

cvtotal2 = vector("list", kboot)

###doing kboot boostraps
for (tt in 1:kboot) {
# cat('doing ',tt,'\n') 
  # k-fold cross-validation
  mat=origmat
#mat=mat[sample(1:nrow(mat),nrow(mat),replace=F),]
fold_id = cut(1:nrow(mat), breaks = kfold, labels = FALSE)

flds <- createFolds(mat[,1], k = kfold, list = TRUE, returnTrain = FALSE)
for (i in 1:kfold) {
  fold_id[flds[[i]]]=i
}
fold_size = sum(fold_id != 1L)

message(sprintf('Selecting genes based on top %f genes cutoff', varcut))
#top = floor(varcut * length(genes) / 100)
top = varcut
features = mclapply(1:kfold, function(i) FS_top(mat[fold_id != i, ], top))
features[[kfold + 1L]] = FS_top(mat, top)
if (length(features[[1]])<100)
cat(features[[1]])

message('Setting up learning algorithm...')
if (learner == "DL") {
    h2o.init(nthreads = Ncores, max_mem_size = max_mem)
    convert = as.h2o
    learn = h2o.deeplearning
    validate = validate_h2o
    predict = h2o.predict
    hyper_params_df = lapply(hyper_params, expand.grid,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    hyper_params_df = Reduce(rbind_data_frame, hyper_params_df)
    Nmodels = min(Nmodels, nrow(hyper_params_df))
    hyper_params_df = hyper_params_df[sample(nrow(hyper_params_df), Nmodels), ]
    hyper_params = expand_hyperparams_h2o(hyper_params_df)
    hyper_params_df = cbind(data.frame(model_id = 1:Nmodels), hyper_params_df)
    saveModel = function(x, nam) h2o.saveModel(x, paste0(nam, ".h2o"))
} else if (learner == "RF") {
    h2o.init(nthreads = Ncores, max_mem_size = max_mem)
    convert = as.h2o
    learn = h2o.randomForest
    validate = validate_h2o
    predict = h2o.predict
    hyper_params_df = lapply(hyper_params, expand.grid,
                             KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    hyper_params_df = Reduce(rbind_data_frame, hyper_params_df)
    Nmodels = min(Nmodels, nrow(hyper_params_df))
    hyper_params_df = hyper_params_df[sample(nrow(hyper_params_df), Nmodels), ]
    hyper_params = expand_hyperparams_h2o(hyper_params_df)
    hyper_params_df = cbind(data.frame(model_id = 1:Nmodels), hyper_params_df)
    saveModel = function(x, nam) h2o.saveModel(x, paste0(nam, ".h2o"))
} else {
    # GLM
    validate = validate_glmnet
    convert = function(x, ...) {
		x[,1]=(x[,1]==1)+0
			as.matrix(x)}
# learn = function(training_frame, y, x, alpha) {
#        cv.glmnet(training_frame[, x], training_frame[, y],
#                  parallel = TRUE, alpha = alpha, family = "gaussian")
#    }
  #  learn = function(training_frame, y, x, alpha) {
  #      cv.glmnet(training_frame[, x], as.vector(training_frame[, y]),
  #                parallel = TRUE, alpha = alpha, family = "binomial") }
###hua added, prevent failing in lambda default
learn=function(training_frame, y, x, alpha) {
  result=tryCatch( {
    cv.glmnet(training_frame[, x], as.vector(training_frame[, y]),
              parallel = TRUE, alpha = alpha, family = "binomial") },
    error=function(err) {
      cat('something wrong')
      cv.glmnet(training_frame[, x], as.vector(training_frame[, y]),
                parallel = TRUE, alpha = alpha, lambda=exp(seq(log(0.001), log(5), length.out=100)),family = "binomial")
    })
  return(result)
}


    hyper_params = hyper_params
    hyper_params_df = data.frame(model_id = 1:length(hyper_params),
                                 alpha = hyper_params)
    saveModel = function(x, nam) saveRDS(x, paste0(nam, ".rds"))
    Nmodels = length(hyper_params)
}

message('Running cross-validation')

cvSummary = vector("list", kfold)

cvSummary2=data.frame(matrix(NaN,0,4))
names(cvSummary2)=c('subj','foldid','bootid','prob')
for (i in seq_len(kfold)) {
    kmat = mat[, c(drug, features[[i]]), drop = FALSE]
    ix = which(fold_id == i)
    train_hex = convert(kmat[-ix, , drop = FALSE])
    valid_hex = convert(kmat[ix, , drop = FALSE])
 if (length(features[[i]])>1 || learner!='GLM') {
    kmodel = lapply(hyper_params, function(params) {
                    x = features[[i]]
                    if ('mtries' %in% names(params)) {# special case
                        params[['mtries']] = round(params[['mtries']] * length(x))
			  if (params[['mtries']]==0)
                          params[['mtries']]=-1  ##hua added to block mtries=0 error
                    }
                    do.call(learn, c(list(y = drug, x = x, training_frame = train_hex), params))
                  })
    cvms = sapply(kmodel, function(mod) mod$cvm[mod$lambda == mod$lambda.min])
    imod = which.min(cvms)
    kmodel = kmodel[imod]
    
    
    ####get pred prob
    
    
    pred_t=predict(kmodel[[1]],valid_hex[,-which(colnames(valid_hex) == drug)],type='response')
    
    pred_t2=data.frame(subj=rownames(pred_t),foldid=i,bootid=tt,prob=pred_t[,1])
    cvSummary2=rbind.data.frame(cvSummary2,pred_t2)
    
    
    
    kvalid = lapply(kmodel, validate, newdata = valid_hex)
    kvalid = do.call(rbind, kvalid)
    kvalid$fold_id = i
    # kvalid$model_id = imod
    tmp_df = kvalid
    # tmp_df = cbind(data.frame(fold_id = i),
    #                hyper_params_df[imod],
    #                kvalid)
    # rownames(tmp_df) = NULL
    cvSummary[[i]] = tmp_df } else {
    kmodel=glm(as.formula(paste0(drug,'~',features[[i]])),family='binomial',data=data.frame(train_hex))
    cauc=auc(valid_hex[,1],predict(kmodel,newdata=data.frame(valid_hex)))

    kvalid = data.frame(model_id=1,MSE=0,AUC=cauc)
    tmp_df = cbind(data.frame(fold_id = rep(i, 1)),
                   hyper_params_df[1,],
                   kvalid)
       pred_t=predict(kmodel,newdata=data.frame(valid_hex),type='response')
    
    pred_t2=data.frame(subj=names(pred_t),foldid=i,bootid=tt,prob=pred_t)
    cvSummary2=rbind.data.frame(cvSummary2,pred_t2)
    rownames(tmp_df) = NULL
    cvSummary[[i]] = kvalid }

}



# Write results -----------------------------------------------------------
cvSummary = do.call(rbind, cvSummary)
cvSummary$bootid=tt
# cvSummary = delist_df(cvSummary)
# cvSummary = cvSummary[order(cvSummary[["AUC"]],decreasing = T), ]
cvtotal[[tt]]=cvSummary
cvtotal2[[tt]]=cvSummary2
#write_csv(cvSummary, paste0(prefix, "_cv.csv"))
#stop()
# 
}

cvtotal = do.call(rbind, cvtotal)
cvtotal2 = do.call(rbind, cvtotal2)
write_csv(cvtotal, paste0(prefix, "_cvSummary.csv"))

write_csv(cvtotal2, paste0(prefix, "_pred_prob.csv"))



# Train best model on all of the data
#best_model = cvSummary$model_id[nrow(cvSummary)]
train_hex = convert(mat)
params = hyper_params[[1]]
x = features[[kfold + 1L]]
if ('mtries' %in% names(params))
  params[['mtries']] = floor(params[['mtries']] * length(x))
if (length(x)>1 || learner != 'GLM') {
  kmodel = lapply(hyper_params, function(params) {
    x = features[[i]]
    if ('mtries' %in% names(params)) {# special case
      params[['mtries']] = round(params[['mtries']] * length(x))
      if (params[['mtries']]==0)
        params[['mtries']]=-1  ##hua added to block mtries=0 error
    }
    do.call(learn, c(list(y = drug, x = x, training_frame = train_hex), params))
  })
  cvms = sapply(kmodel, function(mod) mod$cvm[mod$lambda == mod$lambda.min])
  imod = which.min(cvms)
  model=kmodel[[imod]] } else
    model=glm(as.formula(paste0(drug,'~',x)),family='binomial',data=data.frame(train_hex))
saveModel(model, prefix)







if (!is.na(test_file)) {
    test = readRDS(test_file)
    patient_id = rownames(test)
    test = convert(test[, x, drop = FALSE])
    y = as.vector(predict(model, test))
    df = tibble(patient_id = patient_id, ic50 = y)
    write_csv(df, paste0(prefix, "_clinicalPrediction.csv"))
}

if (learner %in% c("DL", "RF"))
    h2o.shutdown(FALSE)


cat('done\n',date())


