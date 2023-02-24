rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreach);library(tidyverse);library(doParallel)
source("../nle.R")
source("../eHL.R")
source("../oracle_eHL.R")
library(ResourceSelection)
library(calibrationband)
source("../eHL.R")
source("../application/instabilities.R") # flexibel cHL calcualtions
library(kableExtra)
s_=.5

# bagged recal
bag <- readRDS("../application/c3recaliso.rds") %>% tibble() %>% mutate(bagged = 1); bag
# oneshotrecal
ones <- readRDS("../application/c3recaliso_oneshot.rds")%>% mutate(bagged = 0); ones

dat <- bind_rows(bag,ones);dat


core.max <-100
cl <- makeCluster(min(parallel::detectCores()-1, core.max))
registerDoParallel(cl)
start.time <- Sys.time()
MC <- 500
S <- c(1/3,1/2,2/3)
B <- c(100, 500, 1000, 10000)
MCsim <- foreach(
  MCi = 1:MC,
  .errorhandling = "pass",
  #.combine=rbind,
  .packages = c("tidyverse", "foreach"))%dopar%{
  seed <- MCi^2
  set.seed(seed)
  iter.fcast <- foreach(i = unique(dat$bagged))%do% {
    bagged_ <- i
    dattemp1 <- dat %>% dplyr::filter(bagged==i)

    iter.split <- foreach(s = 1:length(S))%do% {
      iter.boot <- foreach(b = 1:length(B))%do% {
        eHL.s. <- eHL(y=dattemp1$Y, P=dattemp1$P, s=S[s], boot=B[b])
        return(list("test"= "eHL", "value" = eHL.s.$HLe, "bagged" = i, "S" = S[s],"B" = B[b], "seed" = seed))
      }
      do.call(rbind, iter.boot)
    }
    eval.HL <- do.call(rbind, iter.split)
    eval.HL
  }
  do.call(rbind, iter.fcast)
}

stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time) # 7 min on mac per run

do.call(rbind, MCsim)


#saveRDS(MCsim, "sim_instability.rds")

res <- do.call(rbind,readRDS("sim_instability.rds")) %>%
  data.frame() %>% tibble() %>%
  mutate_at(c('B','S','value','bagged','seed'), as.numeric) %>%
  mutate_at(c('test'),as.character) %>%
  mutate(
  B = as.factor(B),
  bagged = as.factor(bagged),
  bagged=recode_factor(
      bagged,
      "0" =  "Isotonic recalibration",
      '1' = "Bagged isotonic recalibration",
      .ordered = T)
  ) %>%
  filter(S==.5)

res


bsample_sizes = c(
  '100' =  "italic(B)==100",
  '500' = "italic(B)==500",
  '1000' = "italic(B)==1000",
  '10000' = "italic(B)==10000"
)


helpers <-
  tibble(
    #S=.5,
    bagged = rep(c("Bagged isotonic recalibration", "Isotonic recalibration"),4),
    value = rep(c(6.14,20.04),4),
    #S = rep(c(1/3,1/2,2/3),2),
    B = rep(c(100,500,1000,10000),2)
  )

ggplot(res)+
  # geom_vline(
  #   helpers,
  #   mapping=aes(xintercept=value),
  #   linetype = "dotted"
  #   #mapping = aes(xintrecept=value)
  # )+
  geom_density(aes(x=value, fill=as.factor(B)), alpha=.8)+
  #geom_histogram(aes(x=value, y=after_stat(..density..)), bins = 50, alpha=.75)+
  ggh4x::facet_nested(
    B~bagged, scales = "free",#,independent = "all"
    labeller = labeller(B  = as_labeller(bsample_sizes, label_parsed)))+
  #ggh4x::facet_nested(bagged~.,scales = "free",independent ="all")+#, scales = "free",independent ="all")+
  theme_bw()+
  theme(
    #aspect.ratio = 1,
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    legend.position="none",
    legend.box="horizontal",
    legend.text = element_text(size = 12,
                               colour = "black",
                               margin = margin(r = 12, unit = "pt")
    ),
    legend.title = element_text(size = 12),
    legend.key.width = unit(1.75, "line"),
    legend.spacing.x = unit(0.05, "cm"),
    axis.text.x = element_text(
      colour = "black",
      size = 12
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 12
    ),
    strip.text = element_text(#face = "bold",
                              size = 10
    ),
    #strip.text.y = element_blank(),
    axis.title=element_text(
      #face = "bold",
      size=12
    )
  )+
  xlab("e-value")+
  xlim(c(0,40))+
  labs(fill=" ")+
  geom_vline(xintercept = 20, color="black")




ggsave("Fig_eHLstability.pdf", width = 9.5, height = 6)


res %>% group_by(B, bagged) %>% summarise(rejections =sum(as.numeric(value>20)))
