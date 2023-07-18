# distribute(workloads, across = ncores)
# Load balancing based on some amount of work per task (workloads)
# to be distributed as evenly as possible across some number of 
# cores (across). Messages to tell what each iteration did. 
# Basically, this follows the "Big Rocks" principle (steven covey):
# - sort the tasks by amount of work required (descending) 
# - iteratively assign each task to the core with the lowest workload
# 
# returns: list with:
#           - core.loads: assigned cumulative loads for each core
#           - core.ids: core assignments for each task
# example: distribute(c(123, 2, 324, 1234, 3212, 32, 3, 500, 23, 234, 3, 3214, 1234), 2)
# 
# MTJ 2023

distribute <- function(workloads, across, messages = F){
  
  sort.order <- order(workloads, decreasing = T)
  workloads <- workloads[sort.order]
  id <- rep(0,length(workloads))
  core.loads <- rep(0, across)

  for (i in 1:length(workloads)){
    # Add tasks by descending workload to the core that has the least
    lazy.core <- which.min(core.loads)
    if (messages){message('Adding ', workloads[i], ' to core ', lazy.core, '.')}
    core.loads[lazy.core] <- core.loads[lazy.core] + workloads[i]
    id[i] <- lazy.core
  }
  
  return(list(core.loads = core.loads,
              core.ids = id[order(sort.order)]))
}
