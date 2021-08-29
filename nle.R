### set-up of system of nonlinear equations and solver ----


nle_setup_solver <- function(x,               # points through which the ordered terms pass
                             d,               # points through which the interaction terms pass
                             fixed_points,    # values for pi(x) or pi(x, d)
                             j,               # values for the last equation in the system
                             iteration = 1,   # number of additional iterations for optim loop
                             model = 'logit', # choice of model
                             order,           # order of the polynomial
                             interaction = 0  # number of interaction terms
){
  
  # dummy list for solutions
  nle_solution <- list()
  
  # loop over all j values
  ## progress bar
  pb = txtProgressBar(min = 1, max = length(j), style = 3) 
  for(p in 1:(length(j))){
    ## set up the polynomial function
    set_up <- function(beta){
      ### logit model
      if(model %in% 'logit'){
        #### equation vector dummy
        eq <- c()
        
        #### set up loop for polynomials
        for(n in 1:(order + 1)){
          
          ##### case of no interaction terms
          if(interaction == 0){
            v <- c()
            for(i in 0:order){
              v1 <- x[n]^i
              v <- c(v, v1)
            }
            
            #### polynomial
            poly <- t(beta) %*% v
            
            #### creating the nth equation
            eqn <- exp(poly)/(1+exp(poly))
            
            #### add to equation vector
            eq <- c(eq, eqn)
            
          }else{
            
            #### case of interaction terms
            for(q in 1:length(d)){
              
              v <- c()
              
              #### different order terms
              for(i in 0:order){
                v1 <- x[n]^i
                v <- c(v, v1)
              }
              
              #### interaction terms
              k <- length(v)
              for(l in 1:k){
                v2 <- v[l]*d[q]
                v <- c(v, v2)
              }
              
              #### polynomial
              poly <- t(beta) %*% v
              
              #### creating the nth equation
              eqn <- exp(poly)/(1+exp(poly))
              
              #### add to equation vector
              eq <- c(eq, eqn)
              
            }
            
          }
          
        }
        
      }
      
      ### combine fixed points with j
      points <- c(fixed_points, j[p])
      
      ### creating sqaured "residuals"
      sol <- crossprod(eq - points)
      
      return(sol)
    }
    
    ## number of betas
    if(interaction != 0){
      nbeta <- order + 1 + (order + 1)*interaction
    }else{
      nbeta <- order + 1
    }
    
    ## solve for the pth j value
    nle_solution_p <- optim(rep(0, nbeta), set_up)
    
    ## additional iterations
    for(r in 1:iteration){
      nle_solution_p <- optim(nle_solution_p$par, set_up)
    }
    
    ## add to list and rename
    nle_solution[[p]] <- nle_solution_p
    names(nle_solution)[p] <- paste0('j =', j[p])
    
    # update progress bar
    setTxtProgressBar(pb,p)
  }
  
  # close progress bar
  close(pb)
  
  return(nle_solution)
}
