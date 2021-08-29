oracle_eHL <- function(y,P,cep, boot = 10, s=.5){
  # y is the realization
  # cep: true qunditional event probability
  
  df <- data.frame("y"=y, "CEP"=cep, "P"= P)
  
  
  #both ifs in the loop is super ineeficient. we dont need the loop for the full sample. 
  
  boot <- foreach(t = 1:boot)%do% {
    
    # split the data set
    
    if(s!=0){
    
    split <- sample(nrow(df), (nrow(df)*s) , replace = F)
    splitted.1 <- df[split,]
    splitted.2 <- df[-split,]
    } 
    
    if(s==0){
      
      splitted.2 <- df
    }

    
    suppressMessages(
      eval.orcale <- splitted.2 %>%  tibble() %>% 
        mutate(q = CEP) %>% 
        mutate(E = (q^y*(1-q)^(1-y))/(P^y*(1-P)^(1-y)) ) %>% 
        summarise(evalue = prod(E))
    )
    
  }
  
  eval.orcale <- mean(do.call(rbind, boot)$evalue)

  return(eval.orcale)
  }

  



s=0
cep=pr
P=dt$prediction

nrow(splitted.2)


