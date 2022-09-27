# IOp
Inequality of opportunity (IOp) measurement and Shapley decomposition based on each or groups of variables, with multiple inequality measure indexs.  

It's able to use all the models in `tidymodels` for IOp measurement.

## Usage
```
IOp <- iop$new(sample = data, 
        target = "fincomeavg2", 
        predictors = circumstances, 
        measure = "MLD")
IOp$Shapley(subgroups = subgroups)

IOp$iop
IOp$Shapley_imp
```

## Dependency: 
`tidyverse`, `tidymodels`, `DescTools`


# IOp_mlr3

mlr3 version

## Dependency: 
`mlr3verse`, `DescTools`
