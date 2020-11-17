library("amap")
source("nsga2.r")


evaluator.multiobjective.clustering <- function(solutions) {
  return( solutions )
}

evaluator.biological.significance <- function() {
}

evaluator.metaheuristics <- function(metaheuristic, meta_params, run_evaluator = evaluator.multiobjective.clustering, runs = 10) {
  results <- sapply( 1:runs, function(n) run_evaluator(do.call(metaheuristic, meta_params)) )
  print( results )
  return( colMeans( results ) )
  #return( rowMeans( results ) )
}

evaluator.friedman <- function() {
}
