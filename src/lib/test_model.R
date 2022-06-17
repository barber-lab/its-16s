test_model <- function(model, split, outcome) {
  levels <- model %>% pull_workflow_fit() %>% pluck("lvl")
  
  truth_tbl <-
    split %>%
    testing() %>%
    transmute(
      .truth_class := !!sym(outcome) %>% factor(levels = levels),
      .truth_pred = .truth_class %>% as.integer() - 1
    )
  
  list(
    model %>% predict(new_data = split %>% testing(), type = "prob"),
    model %>% predict(new_data = split %>% testing(), type = "class")
  ) %>%
    bind_cols() %>%
    bind_cols(truth_tbl) %>% {
      # auto invert classifier e.g. to constrain roc_auc to [0.5, 1]
      x <- .
      a <- yardstick::metrics(x, .truth_class, .pred_class, paste0(".pred_", levels[[1]]))
      b <- yardstick::metrics(x, .truth_class, .pred_class, paste0(".pred_", levels[[2]]))
      auc <- function(x) x %>% filter(.metric == "roc_auc") %>% pull(.estimate) %>% first()
      
      if(auc(a) > auc(b)) {
        a
      } else {
        b
      }
    }
}