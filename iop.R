library("tidyverse")
library("tidymodels")
library("DescTools")
library("R6")



# notes: 
# target: , predictors: vectors
# measure: 
iop <- R6Class(
  "inequality of opportunity", 
  public = list(
    sample     = NULL, 
    target     = NULL, 
    predictors = NULL, 
    bootstrap  = NULL, 
    measure    = NULL, 
    workflow   = NULL, 
    Shapley_imp = NULL, 
    
    pred = NULL, 
    iop  = NULL, 
    
    initialize = function(sample = NULL, target = NULL, predictors = NULL, measure = NULL, 
                          bootstrap = NULL, subgroups = NULL, model = NULL) {
      stopifnot(!is.null(sample))
      
      self$sample      = sample
      self$target      = target
      self$predictors  = predictors
      self$measure     = measure
      
      # create a workflow with linear regression model as default
      self$workflow = workflow() %>% 
        when(!is.null(model) ~ (.) %>% add_model(model), 
             is.null(model)  ~ (.) %>% add_model(
               linear_reg() %>% 
                 set_mode("regression") %>% 
                 set_engine("lm")
             )) %>% 
        add_recipe(
          recipe(cformula(x = self$predictors, y = self$target), data = self$sample)
        ) %>% 
        fit(self$sample)
      
      prediction  = predict(self$workflow, new_data = self$sample %>% select(all_of(self$predictors)))
      self$pred   = self$sample %>% select(self$target) %>% bind_cols(prediction)
      
      self$iop_measure()
    }, 
    
    # measure
    cal = function(x = NULL, y = NULL, measure = self$measure) {
      
      if (measure == "R2")   r = cor(x, y)^2
      if (measure == "MLD")  r = log(mean(x)) - mean(log(x))
      if (measure == "MSE")  r = var(x)
      if (measure == "pdb")  r = sum(abs(x - mean(x))) / (length(x)*mean(x))
      if (measure == "ws")   r = sum(abs(x - mean(x))) / (length(x)*2)
      if (measure == "Gini") r = Gini(x, unbiased = FALSE)
      
      return(r)
    }, 
    
    # calculate IOp, total inequality and IOR
    iop_measure = function(...) {
      
      self$iop = rep(0, times = 3)
      names(self$iop) = c("IOp", "Total", "IOR")
      
      # IOp results
      self$iop[["IOp"]]             = self$cal(self$pred$.pred)
      self$iop[["Total"]]           = self$cal(self$pred[[self$target]])
      self$iop[["IOR"]]             = self$iop[["IOp"]] / self$iop[["Total"]]
    }, 
    
    # create formula from a vector
    cformula = function(x = NULL, y = NULL) {
      y = str_c(y, " ~ ", collapse = "")
      x = str_c(x, collapse = " + ")
      formula = str_c(y, x, collapse = "")
      return(as.formula(formula))
    }, 
    
    # Shapley decomposition
    Shapley = function(sample = NULL, subgroups = NULL) {
      
      # create all the combinations and formulas
      target = self$target
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
            formula[[1]] = formulas[[j]]
            formula[[2]] = setdiff(formulas[[j]], subgroups[[i]])
            
            task_1 = self$sample %>% select(self$target, formula[[1]])
            task_2 = self$sample %>% select(self$target, formula[[2]])
            
            # with & without subgroups[[i]]
            prediction = self$workflow %>% 
              remove_recipe() %>% 
              add_recipe(
                recipe(cformula(x = formula[[1]], y = self$target), data = task_1)
              ) %>% 
              fit(task_1) %>% 
              predict(new_data = task_1)
            pred_1     = self$sample %>% select(self$target) %>% bind_cols(prediction) %>% rename("truth" = self$target)
            
            prediction = self$workflow %>% 
              remove_recipe() %>% 
              add_recipe(
                recipe(cformula(x = formula[[2]], y = self$target), data = task_2)
              ) %>% 
              fit(task_2) %>% 
              predict(new_data = task_2)
            pred_2     = self$sample %>% select(self$target) %>% bind_cols(prediction) %>% rename("truth" = self$target)
            
            # change
            r = c(0, 0)
            if (var(pred_1$.pred) == 0 | var(pred_2$.pred) == 0) {
              change = 0
            } else {
              r[1] = self$cal(pred_1$.pred, pred_1$truth, self$measure)
              r[2] = self$cal(pred_2$.pred, pred_2$truth, self$measure)
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


