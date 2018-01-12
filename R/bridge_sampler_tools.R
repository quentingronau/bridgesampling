
#--------------------------------------------------------------------------
# functions for nimble support
#--------------------------------------------------------------------------


.log_posterior_nimble <- nimble::nimbleFunction(

  # based on code by Perry de Valpine

  ## setup code is executed in R and specializes an instance
  ## of the nimbleFunction to a particular model or nodes
  setup = function(model, nodes) {
    calcNodes <- model$getDependencies(nodes)
  },
  ## run code is called repeatedly and can be converted into C++
  run = function(sample = double(1)) {
    values(model, nodes) <<- sample
    out <- model$calculate(calcNodes)
    return(out)
    returnType(double(0))
  }
)

.nimble_bounds <- function(samples, model, which) {

  if ( ! (which %in% c("lower", "upper")) ) {
    stop('"which" needs to be either "lower" or "upper"\n',  call. = FALSE)
  }

  cn <- colnames(samples)
  bounds <- numeric(length(cn))
  names(bounds) <- cn

  for (i in seq_along(cn)) {
    bounds[[cn[i]]] <- model$getBound(cn[i], which)
  }

  return(bounds)

}
