fix_rf_workflow <- function(workflow) {
  specs <- pull_workflow_spec(workflow)
  
  predictors <-
    workflow %>%
    pull_workflow_preprocessor() %>%
    summary() %>%
    filter(role == "predictor") %>%
    pull(variable)
  
  # ensure less mtry than predictors
  specs %<>% update(mtry = !!min(length(predictors) - 2 ))
  
  workflow %>%
    update_model(specs)
}