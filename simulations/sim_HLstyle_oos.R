rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreach);library(tidyverse);library(doParallel)
source("../nle.R")
source("../eHL.R")
source("../oracle_eHL.R")

set.seed(123)

# Obtain the parameter values

J<-seq(0,0.10,.002)+.007


no_interact <- nle_setup_solver(
  x = c(-1.5, 3, -3),
  fixed_points = c(0.05, 0.95),
  j = J,
  iteration = 10,
  order = 2
)
list <- lapply(no_interact, '[[', 1)
pmeter_quad <- data.frame(
  J ,
  matrix(unlist(list),
         nrow=length(list),
         byrow=TRUE)
)
colnames(pmeter_quad) <- c("J", "beta0", "beta1","beta2")
pmeter_quad %>% arrange(abs(beta2))


N <- 2^c(10:14)
S <- sort(c(1/4,2/3,1/3,1/2,3/4))   #seq(.2,.8,.1)
setup <- c("quadratic")#,"interaction")

core.max <- 100
cl <- makeCluster(min(parallel::detectCores()-1, core.max))
registerDoParallel(cl)

start.time <- Sys.time()

MC <- 1000
MCsim <- foreach(
  MCi = 1:MC,
  .errorhandling = "pass",
  #.combine=rbind,
  .packages = c("tidyverse", "foreach")
  )%dopar%{
    seed <- MCi^2

    iter.setup  <- foreach(st = 1:length(setup)) %do%{
    set <- setup[st]

    iter.obs <- foreach(i = 1:length(N))%do% {
      n <-  N[i]
      HLs <- foreach(j = 1:nrow(pmeter_quad))%do% {



        if(set=="quadratic"){
          set.seed(seed)
          x1 <- runif(n*2, min = -3, max = 3)
          x2 <- x1^2
          xb <- pmeter_quad$beta0[j] + pmeter_quad$beta1[j] * x1 + pmeter_quad$beta2[j] * x2
          pr <- exp(xb)/(1 + exp(xb)) #= exp(xb)/(1+exp(xb))
          y <- purrr::rbernoulli(n*2, p = pr)*1
          dt <- data.frame(y,x1,x2)
          dt.is <- dt[1:n,]
          dt.oos <- dt[(n+1):(n*2),1:2]



          jj <- pmeter_quad$J[j]

          mod.hl <- glm(y~x1, family = binomial(link = "logit"), data = dt.is)
        }

        if(isTRUE(jj==.007)){
          dt.oos$prediction <- pr[(n+1):(n*2)]
        } else {
          dt.oos$prediction <- predict(mod.hl, type = "response", newdata = dt.oos)
        }

        classic.HL <- ResourceSelection::hoslem.test(dt.oos$y,dt.oos$prediction, g=10)

        # OOS HL pvalue:
        #PVAL = 1 - pchisq(as.numeric(classic.HL$statistic),10)

        cHL <- list("test"= "cHL", "value" = classic.HL$statistic, "obs" = n, "J" = jj, "S" = NA, "setup" = set, "seed" = seed)
        p.val.cHL <- do.call(cbind, cHL)

        oracle.eHL <- oracle_eHL(y=y[(n+1):(n*2)], cep=pr[(n+1):(n*2)], P=dt.oos$prediction)
        e.val.oracle.eHL <- do.call(
          cbind,
          list("test"= "oracle_eHL",
               "value" = as.numeric(oracle.eHL),
               "obs" = n,
               "J" = jj,
               "S" = 0.5,
               "setup" = set,
               "seed" = seed
          )
        )

        oracle.eHL06 <- oracle_eHL(y=y[(n+1):(n*2)], cep=pr[(n+1):(n*2)], P=dt.oos$prediction,s=0.6)
        e.val.oracle.eHL06 <- do.call(
          cbind,
          list("test"= "oracle_eHL",
               "value" = as.numeric(oracle.eHL06),
               "obs" = n,
               "J" = jj,
               "S" = 0.6,
               "setup" = set,
               "seed" = seed
          )
        )


        oracle.eHL.s0 <- oracle_eHL(
          y=y[(n+1):(n*2)], cep=pr[(n+1):(n*2)], P=dt.oos$prediction,
          s=0)
        e.val.oracle.eHL.s0 <- do.call(
          cbind,
          list(
            "test"= "oracle_eHL",
            "value" = as.numeric(oracle.eHL.s0),
            "obs" = n,
            "J" = jj,
            "S" = 0,
            "setup" = set,
            "seed" = seed
          )
        )

        iter.split <- foreach(s = 1:length(S))%do% {
          eHL.s. <- eHL(y=dt.oos$y, P=dt.oos$prediction, s=S[s], boot=10)
          return(list("test"= "eHL", "value" = eHL.s.$HLe, "obs" = n, "J" = jj, "S" = S[s],  "setup" = set,"seed" = seed))
        }
        e.val.HL <- do.call(rbind, iter.split)


        list(p.val.cHL,e.val.HL, e.val.oracle.eHL, e.val.oracle.eHL.s0,e.val.oracle.eHL06)
      }

      do.call(rbind,do.call(Map, c(f = rbind, HLs)))


    }

    do.call(rbind,iter.obs)


    #do.call(rbind,do.call(Map, c(f = rbind, HLs)))
    }


  do.call(rbind,iter.setup) %>% data.frame() %>%
    tibble() %>%
    mutate_at(c('J','S','obs','value','seed'), as.numeric) %>%
    mutate_at(c('test','setup'),as.character)
}
stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)
#Time difference of 8.258534 hours 2000 reps

# length(N)*length(J)*length(S)+(length(N)*length(J))*3
# thats the number of observations in each run

saveRDS(MCsim, "sim_HLsytle_oos.rds")

MCsim <-
  readRDS("sim_HLsytle_oos.rds")

check <- sapply(MCsim, function(x) inherits(x, 'error'))
out <- MCsim[!check]
MCsim[check]
dat <-do.call("rbind", out)
dat
### REMEMBER THAT VALUE OF cHL NOW IS THE TEST STATISTIC--> TRANSFORM TO p-VALUE!!!!


# rejection rates ----
dt <-
  dat %>%
  mutate(
    S=as.factor(S),
    dec = if_else(
      test=="cHL",
      as.numeric(value<=0.05),
      as.numeric(value>=20))
    ) %>%
  group_by(obs,J,test,S,setup) %>%
  summarize(rej = mean(dec)) %>%
  dplyr::filter(
    setup == "quadratic",
    #test == "eHL" && S%in%c("0.2","0.4","0.6","0.8",NA,0,1)
    !(test=="eHL" && S%in%c("0.3","0.5","0.7")),
    J<.105,obs<9000,obs>1000, J>=.007
    )


dt$J %>% unique()
dt$test %>% unique()
dt$S %>% unique()

rr <-
ggplot(
  dt %>%
    filter(test=="eHL"),
  aes(x=J,y=rej, colour = factor(S))
  )+
  geom_line(size=.55) +
  geom_point(
    dt%>% filter(test=="eHL"),
    mapping = aes(x=J,y=rej, shape=factor(S)),
    size=2
  )+
  theme_bw() +
  theme(
    legend.position=c(.7,.1),
    legend.box="horizontal"
  )+
  ylab("Rejection rate")+
  xlab("lack of linearity")+
  scale_color_brewer(palette = "Dark2")+
  geom_hline(
    aes(yintercept = 0.05),
    alpha=.25,
    size = .25,
    linetype="dashed"
  )+
  geom_vline(
    aes(xintercept = 0.007),
    alpha=.25,
    size = .25,
    linetype="dashed")+
  facet_wrap(obs ~ ., scales = "free", ncol = 1)+
  labs(colour = "legend title")+
  geom_line(
    dt %>% filter(test=="oracle_eHL", S==0),
    mapping=aes(J,rej, linetype = "oracle eHL,s=0"), size =.5, colour = "black"
  )+
  geom_line(
    dt %>% filter(test=="oracle_eHL", S==.5),
    mapping= aes(J,rej, linetype = "oracle eHL,s=1"),
    size =.5,
    colour = "mediumblue"
  )+
  geom_line(
    dt %>% filter(test=="cHL"),
    mapping=aes(J,rej, linetype = "classic HL"), size = .5, colour = "darkgray")+
  labs(
    color  = "Proportion of data used to estimate q:",
    linetype = "Test:",
    shape = "Proportion of data used to estimate q:"
  )+
  scale_linetype_manual(
    values=c("dotdash", "solid", "dotted", "dashed")
  )+
  scale_shape_manual(
    values = c(3,8,17,19,4,5,9)
  )


  scale_y_continuous(
    trans="log10"
    #breaks = log(c(0.2,.4)),
    #labels = c(".2",".4")
  )+
  scale_x_continuous(
    trans="log10"
    #breaks = log(c(0.2,.4)),
    #labels = c(".2",".4")
  )

breaks = breaks, labels = comma(breaks, digits = 1)
scale_y_continuous(
  trans="sqrt"
  #breaks = log(c(0.2,.4)),
  #labels = c(".2",".4")
)
scale_y_continuous(
  trans="pseudo_log"
  #breaks = log(c(0.2,.4)),
  #labels = c(".2",".4")
)





# growth rates -----

grow.evals <-
  data.frame(dat) %>%
  mutate(
    J = as.numeric(J),
    value=as.numeric(value),
    obs=as.numeric(obs),
    S=as.factor(S),
    test=as.character(test)) %>%
  filter(
    test %in% c("eHL","oracle_eHL")
  ) %>%
  group_by(obs,J,test,S) %>%
  summarize(
    grow.eval = mean(log(value))
  )%>%
  dplyr::filter(
    !(test=="eHL" & S%in%c("0.3","0.5","0.7")),
    J<.105,obs<9000,obs>1000, J>=.007
  )


  filter(!S%in%c(.8))

gr <-
ggplot(
  grow.evals %>%
    filter(test=="eHL"),
  aes(J, grow.eval, colour = factor(S), linetype=test)
) +
  geom_line(size=.75) +
  geom_point(
    grow.evals%>% filter(test=="eHL"),
    mapping = aes(x=J,y=grow.eval, shape=factor(S)),
    size=3
  )+
  geom_vline(
    aes(xintercept = 0.007),
    alpha=.25,
    size = .7,
    linetype="dashed") +
  geom_line(
    grow.evals %>% filter(test=="oracle_eHL", S == 0),
    mapping=aes(
      J,
      grow.eval,
      linetype = "oracle eHL,s=0"
    ),
    size =0.5,
    colour = "black"
  )+
  geom_line(
    grow.evals %>% filter(test=="oracle_eHL", S == 0.5),
    mapping=aes(
      J,
      grow.eval,
      linetype = "oracle eHL,s=0.5"
    ),
    size =.5,
    colour = "black"
  )+
  theme_bw() +
  theme(legend.position="bottom")+#, legend.box="vertical")+
  labs(
    color  = "Proportion of data for q:",
    linetype = "Test:",
    shape = "Proportion of data for q:")+
  ylab("growth rate e-value")+
  xlab("lack of linearity")+
  scale_color_brewer(palette = "Dark2")+
  scale_shape_manual(
    values = c(3,8,17,19,4,5,9)
  )+
  scale_linetype_manual(values=c("solid","dotted", "dashed"))+
  facet_grid(obs ~ ., scales = "free", )




gr

rr

comdat <-
bind_rows(
  dt %>% mutate(value = rej, id = "a) rejection rate") %>% select(!c(setup,rej)),
  grow.evals %>% mutate(value = grow.eval, id = "b) growth rate") %>% select(!c(grow.eval))
) %>%
  mutate(
    test = ifelse(test=="oracle_eHL", "oracle eHL", test),
    labs = ifelse(test=="cHL", "cHL",paste(test," (s=",S,")",sep = ""))
  ) %>% arrange(labs)


c(comdat$J %>% unique())


c(.007,seq(0.02,0.1,0.02))

ggplot(
  comdat,
  aes(x=J,y=value, colour = labs, linetype=labs,
      )
  )+
  geom_line(size=.55) +
  geom_point(
    comdat %>% filter(J%in%c(.007,seq(0.02,0.1,0.02))),
    mapping = aes(
      x=J,y=value,
      shape=labs,  colour = labs,
      #linetype=factor(labs)
      ),
    size=2
  )+
  theme_bw() +
  theme(
    legend.position="bottom",
    legend.box="horizontal"
  )+
  ylab("rate")+
  xlab("lack of linearity")+
  #scale_color_brewer(palette = "Dark2")+
  ggh4x::facet_nested(obs ~ id, scales = "free", independent ="y")+
  guides(
    colour = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1)
    )+
  labs(
    color  = "Test:",
    linetype = "Test:",
    shape = "Test:"
  )+
  scale_shape_manual(
    values = c(NA,1,4,8,2,NA,NA)
  )+
  scale_linetype_manual(
    values = c("dotted", "dashed", "dashed", "dashed", "dashed","solid","solid")
  )+
  scale_colour_manual(
    values = c("black", "darkgreen", "blue", "red", "coral","#CC00FFFF","darkgray")
  )+
  scale_x_continuous(breaks=c(.007,seq(0.02,0.1,0.02)))

ggsave("Fig_Simulation.pdf", width = 12, height = 8 )
getwd()
