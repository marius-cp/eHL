rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreach);library(tidyverse);library(doParallel)
source("nle.R")
source("eHL.R")
source("oracle_eHL.R")
set.seed(2021)

# Obtain the parameter values 
J<-c(.001,.002,.007,seq(0.015,0.2,0.01), .2)#c(.001,.002,.007,.01,.025,.05,.1,.2)
no_interact <- nle_setup_solver(x = c(-1.5, 3, -3), fixed_points = c(0.05, 0.95),
                                j = J, iteration = 10,
                                order = 2)
list <- lapply(no_interact, '[[', 1)
pmeter <- data.frame(J , matrix(unlist(list), nrow=length(list), byrow=TRUE))
colnames(pmeter) <- c("J", "beta0", "beta1","beta2")          
pmeter

N <-  2^c(9,10,11,12) 

S <- seq(.2,.8,.1)

core.max <- 50
cl <- makeCluster(min(parallel::detectCores()-1, core.max) )
registerDoParallel(cl)
start.time <- Sys.time()
MCsim <- foreach(i = 1:5000, .combine=rbind, .packages = c("tidyverse", "foreach"))%dopar%{

iter.obs <- foreach(i = 1:length(N))%do% {
n <-  N[i]
HLs <- foreach(j = 1:nrow(pmeter))%do% {
  

x1 <- runif(n, min = -3, max = 3)
x2 <- x1^2
xb <- pmeter$beta0[j] + pmeter$beta1[j] * x1 + pmeter$beta2[j] * x2
pr <- exp(xb)/(1 + exp(xb)) #= exp(xb)/(1+exp(xb))
y <- purrr::rbernoulli(n, p = pr)*1
dt <- data.frame(y,x1,x2)


mod.hl <- glm(y~x1, family = binomial(link = "logit"), data = dt)
dt$prediction <- predict(mod.hl, type = "response") # same as fitted(mod.hl)

classic.HL <- ResourceSelection::hoslem.test(dt$y,fitted(mod.hl), g=10)$p.value
cHL <- list("test"= "cHL", "value" = classic.HL, "obs" = n, "J" = J[j], "S" = NA)
p.val.cHL <- do.call(cbind, cHL) # das kostet bestimmt zeit...

oracle.eHL <- oracle_eHL(y=y, cep=pr, P=dt$prediction)
e.val.oracle.eHL <- do.call(cbind,list("test"= "oracle_eHL", "value" = as.numeric(oracle.eHL), "obs" = n, "J" = J[j], "S" = 0.5))


oracle.eHL.s0 <- oracle_eHL(y=y, cep=pr, P=dt$prediction, s=0)
e.val.oracle.eHL.s0 <- do.call(cbind,list("test"= "oracle_eHL", "value" = as.numeric(oracle.eHL.s0), "obs" = n, "J" = J[j], "S" = 0))



iter.split <- foreach(s = 1:length(S))%do% {
eHL.s. <- eHL(y=dt$y, P=dt$prediction, s=S[s], boot=10)
return(list("test"= "eHL", "value" = eHL.s., "obs" = n, "J" = J[j], "S" = S[s]))  
}
e.val.HL <- do.call(rbind, iter.split)# das kostet bestimmt zeit...

list(p.val.cHL,e.val.HL, e.val.oracle.eHL, e.val.oracle.eHL.s0)

#return(list("eHL" = e.val.HL, "cHL"=classic.HL, "obs" = n, "J" = J[j], "S" = S[s])) 
}

do.call(rbind,do.call(Map, c(f = rbind, HLs))) # wahrscb´heinlich ist es am besten nur ein mal do.call aufzurufen

#data.frame(do.call(rbind, HLs)) %>% gather(key="test", value = "val", c(eHL,cHL))
}
}

stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)

# Time difference of 8.669647 hours, 50 cores

saveRDS(do.call(rbind, MCsim), "MCsim_data.rds")


dd <- do.call(rbind, MCsim) %>% data.frame() %>% 
  mutate(J = as.numeric(J),
         value=as.numeric(value),
         obs=as.numeric(obs),
         S=as.numeric(S), 
         test = as.character(test))
dd
