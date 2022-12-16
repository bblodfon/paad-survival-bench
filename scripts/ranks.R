library(combinat)
library(tidyverse)

# Construct example data ----
tasks = c('A', 'B', 'C')
tasks = LETTERS[1:8]
n = length(tasks)
n

task_subsets = lapply(1:length(tasks), combinat::combn, x = tasks,
  simplify = FALSE) %>% unlist(recursive = FALSE)

task_subsets = lapply(task_subsets, function(el) {
  paste0(el, collapse = '-')
}) %>% unlist()
task_subsets # ok, all there

ranks = c(5,4,7,2,6,3,1)
ranks = sample(x = 1:(2^n-1), size = 2^n-1, replace = FALSE)
names(ranks) = task_subsets
rank_tbl = ranks %>% enframe(name = 'task_combo', value = 'rank')
rank_tbl %>%
  arrange(rank) # ok, this is our input ranking

# get per-task ranking
rank_sum = sapply(tasks, function(task) {
  rank_tbl %>%
    filter(grepl(pattern = task, x = task_combo)) %>%
    summarise(s = sum(rank)) %>%
    pull(s)
})
rank_sum
sort(rank_sum)
sort(rank(rank_sum))

best_rank_sum  = sum(1:2^(n-1))
best_rank_sum # smaller
worst_rank_sum = sum(seq(to = 2^n-1, length.out = 2^(n-1)))
worst_rank_sum # larger

# normalized score
norm1 = (rank_sum - best_rank_sum)/(worst_rank_sum - best_rank_sum)
sort(norm1) # smaller is best
# norm2 = 1 - norm1
norm2 = (worst_rank_sum - rank_sum)/(worst_rank_sum - best_rank_sum)
sort(norm2, decreasing = TRUE) # higher is best


