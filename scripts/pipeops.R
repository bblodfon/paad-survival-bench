library('mlr3pipelines')
library('mlr3proba')
library('paradox')

#' shuffles train set's Surv targets (time, status)
PipeOpSurvShuffle = R6::R6Class('PipeOpSurvShuffle', inherit = PipeOpTaskPreproc,
  public = list(
    initialize = function(id = 'survshuffle', param_vals = list()) {
      p = ps(replace = p_lgl(tags = 'required'))
      p$values = list(replace = FALSE)
      super$initialize(id = id, param_set = p, param_vals = param_vals,
        can_subset_cols = FALSE, task_type = 'TaskSurv'
      )
    }
  ),
  private = list(
    .train_task = function(task) {
      pvals = self$param_set$get_values()
      surv  = task$data(cols = c('time', 'status'))
      if (nrow(surv) > 1) {  # `sample` 'misbehaves' when 1st argument has length 1!
        surv$time   = sample(surv$time,   replace = pvals$replace)
        surv$status = sample(surv$status, replace = pvals$replace)
      }
      # $cbind() overwrites old task columns
      # check if this breaks inside resample()
      task$cbind(surv)
    },
    .predict_task = function(task) task # don't shuffle task during prediction!
  )
)

poss = PipeOpSurvShuffle$new()