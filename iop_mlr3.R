####################################
## iop: inequality of opportunity ##
####################################
# base on ml3
library("mlr3verse")
library("DescTools")
library("R6")



# (1) inequality of opportunity evaluation
# parameters:  samples (list; data.frame list with targets as list's names), learner (mlr3 LearnerRegr)
# iop results: iop (list; IOp, Total, IOR results of each data.frame)
iop <- R6Class(
  "inequality of opportunity", 
  public = list(
    
    samples = NULL, 
    learner = NULL, 
    
    targets = c(), 
    task    = list(), 
    pred    = list(), 
    iop     = list(), 
    formula = list(), 
    
    initialize = function(samples, learner) {
      stopifnot(is.list(samples) & !is.null(names(samples)))
      stopifnot(class(learner)[[2]] == "LearnerRegr")
      
      self$samples = samples
      self$learner = learner
      self$targets = names(samples)
      
      # set each task
      for (i in 1:length(samples)) {
        self$task[[self$targets[[i]]]] = as_task_regr(samples[[i]], target = self$targets[[i]])
      }
      
      # calculate IOp by each task
      for (i in 1:length(self$task)) {   # train learner and predict by each task
        self$learner$train(self$task[[i]])
        self$pred[[self$targets[i]]] = self$learner$predict(self$task[[i]])
      }
      self$iop_measure()
    }, 
    
    # (1) calculate IOp, total inequality & IOR
    iop_measure = function() {
      for (i in self$targets) {
        self$iop[[i]] = c(0, 0, 0, 0, 0)
        names(self$iop[[i]]) = c("IOp", "Total", "IOR", "MSE", "N")
        
        self$iop[[i]][["IOp"]]   = log(mean(self$pred[[i]]$response)) - mean(log(self$pred[[i]]$response))
        self$iop[[i]][["Total"]] = log(mean(self$pred[[i]]$truth)) - mean(log(self$pred[[i]]$truth))
        self$iop[[i]][["IOR"]]   = self$iop[[i]][["IOp"]] / self$iop[[i]][["Total"]]
        self$iop[[i]][["MSE"]]   = self$pred[[i]]$score(msr("regr.mse"))
        self$iop[[i]][["N"]]     = length(self$pred[[i]]$truth)
      }
    }
  )
)



# (2) Shapley decomposition (R2)
# parameters: sample (data.frame or tbl), target (characters), learner (mlr3 LearnerRegr), 
#             subgroups (list; subgroup variable name vectors)
Shapley <- R6Class(
  "Shapley", 
  public = list(
    
    sample    = NULL, 
    target    = NA, 
    learner   = NULL, 
    subgroups = NULL, 
    measure   = NULL, 
    
    Shapley_imp = list(), 
    
    initialize = function(sample, target, learner, measure = NULL) {
      stopifnot(is.character(target))
      
      self$sample  = sample
      self$target  = target
      self$learner = learner
      self$subgroups = subgroups
      self$measure = measure
      self$Shapley(self$sample, self$subgroups, self$target)
    }, 
    
    # measure
    cal = function(x, y = NULL, measure) {
      
      if (measure == "R2")   r = cor(x, y)^2
      if (measure == "MLD")  r = log(mean(x)) - mean(log(x))
      if (measure == "MSE")  r = var(x)
      if (measure == "pdb")  r = sum(abs(x - mean(x))) / (length(x)*mean(x))
      if (measure == "ws")   r = sum(abs(x - mean(x))) / (length(x)*2)
      if (measure == "Gini") r = Gini(x, unbiased = FALSE)
      
      return(r)
    }, 
    
    Shapley = function(sample = NULL, subgroups = NULL, target = NA) {
      
      # create all the combinations and formulas
      combinations <- lapply(1:length(subgroups), function(y) combn(subgroups, y))
      combination_names <- lapply(1:length(subgroups), function(y) combn(names(subgroups), y))
      formulas <- list()
      n <- 1
      for (i in 1:length(combinations)) {
        for (j in 1:ncol(combinations[[i]])) {
          formulas[[n]] <- unlist(combinations[[i]][, j])
          n <- n+1
        }
      }
      
      # calculate Shapley variable importance
      for (i in 1:length(subgroups)) {
        self$Shapley_imp[[names(subgroups)[i]]] = 0
      }
      for (i in 1:length(subgroups)) {
        for (j in (length(subgroups)+1):length(formulas)) {
          if (all(subgroups[[i]] %in% formulas[[j]])) {
            formula = list()
            R2 = c(0, 0)
            formula[[1]] = append(target, formulas[[j]])
            formula[[2]] = append(target, setdiff(formulas[[j]], subgroups[[i]]))
            
            task_1 = as_task_regr(sample[, formula[[1]]], target = target)
            task_2 = as_task_regr(sample[, formula[[2]]], target = target)
            
            # with & without subgroups[[i]]
            self$learner$train(task_1)
            pred_1 = self$learner$predict(task_1)
            
            self$learner$train(task_2)
            pred_2 = self$learner$predict(task_2)
            
            # change
            r = c(0, 0)
            if (var(pred_1$response) == 0 | var(pred_2$response) == 0) {
              change = 0
            } else {
              r[1] = self$cal(pred_1$response, pred_1$truth, self$measure)
              r[2] = self$cal(pred_2$response, pred_2$truth, self$measure)
              change = r[1] - r[2]
            }
            if (self$measure == "R2" && change < 0) change = 0
            
            self$Shapley_imp[[names(subgroups)[i]]] = self$Shapley_imp[[names(subgroups)[i]]] + change
          }
        }
        print(sprintf("%s complete %s / %s", class(self$learner)[1], i, length(subgroups)))
      }
      
      imp = unlist(self$Shapley_imp)
      self$Shapley_imp = c()
      for (i in 1:length(imp)) {
        self$Shapley_imp[i] = imp[[i]] / sum(imp)
      }
      names(self$Shapley_imp) = names(subgroups)
    }
  )
)


