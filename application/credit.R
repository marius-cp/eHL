rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(ResourceSelection)
library(calibrationband)
library(doParallel)# for loop in eHL function
source("../eHL.R")
source("instabilities.R") # flexibel cHL calcualtions
library(betacal)
library(kableExtra)
bbb=10000 # overall use of bootstrap reps
s_=.5

# main calculations ----

credit <- readxl::read_xls("default of credit card clients.xls", col_names =T)
names(credit) <- as.matrix(credit[1, ])
credit <- credit[-1, ]
credit[] <- lapply(credit, function(x) type.convert(as.character(x)))
credit$default <- credit$`default payment next month`
credit <- credit %>%  select(!c(`default payment next month`)) %>%  select(!ID)

set.seed(1234)
samples = sample(c(1, 1, 2, 2, 3),dim(credit)[1], rep=TRUE)

c1 <- credit[(samples==1),]
c2 <- credit[(samples==2),]
c3 <- credit[(samples==3),]

# estimate logit in c1 data
logit.tr <- glm(default ~  ., data = c1, family=binomial(link='logit'))
#predict in c2
pairs.c2<- tibble(
  "Y" = c2$default, "P" = predict(logit.tr,new =c2, type='response')
) %>%
  arrange(P)



pairs.c2 %>% filter(P==0&Y==1)
pairs.c2 %>% filter(P==1&Y==0)


bagg <- function(dat, B=100, seed=1234,recal="iso"){
  set.seed(seed)
  list <- list()
  n <- dim(dat)[1]
  for (b in 1:B) {
    obs <- sample(c(1:n), n, replace = TRUE)
    bootdat <- dat[obs,] %>% arrange(P)
    if(recal=="iso"){
      # isoreg in val dat
      isofit <- isoreg(x=bootdat$P,y=bootdat$Y)

      isolin <- approxfun(
        y=c(min(isofit$yf),isofit$yf,max(isofit$yf)),
        x=c(0,isofit$x,1),
        method = "linear",
        ties = "mean"
      )
      list[[paste0("isoBagg_", b)]] = isolin
    }

    if(recal=="beta"){
      # beta calibration
      bc <- beta_calibration(bootdat$P,bootdat$Y, parameters = "abm")
      list[[paste0("isoBagg_", b)]] = bc
    }
  }
  return(list)
}


baggedIso <- bagg(dat=pairs.c2, B=100,recal="iso")
#baggedBeta <- bagg(dat=pairs.c2, B=100, seed=123,recal="beta")

#predict in c3
pairs.c3<- tibble(
  "Y" = c3$default, "P" = predict(logit.tr,new =c3, type='response')
) %>%
  arrange(P)


pairs.c3 %>% filter(P==0&Y==1)
pairs.c3 %>% filter(P==1&Y==0)

recalc3iso <- lapply(baggedIso, function(x){x(pairs.c3$P)})
c3recaliso <-
  data.frame(
    do.call("cbind",recalc3iso),
    P=do.call("cbind",recalc3iso) %>% rowMeans(),
    Y = pairs.c3$Y
  ) %>%
  arrange(P)

saveRDS(c3recaliso %>% select(P,Y), "c3recaliso.rds")


# at least one is non-empty if the e-value is INF
c3recaliso %>% filter(P==0&Y==1)
c3recaliso %>% filter(P==1&Y==0)

## Fig 2 ----

pdat1 <- lapply(baggedIso, function(x){x(seq(1,1000)/1000)})
pdat2 <-
  data.frame(
    X=seq(1,1000)/1000,
    do.call("cbind",pdat1),
    P=do.call("cbind",pdat1) %>% rowMeans()
    #Y = pairs.c3$Y
  ) %>% arrange(P)


pdat2 %>%
  pivot_longer(
    cols = isoBagg_1:P,
    values_to = "value"
  ) %>%
  group_by(X) %>%
  summarise(
    miniso = min(value),
    q05 = quantile(value,.99),
    maxiso = max(value),
    q95 = quantile(value,.01)
  ) %>%
  ggplot()+
  geom_ribbon(
    aes(x=X,ymin=q05,ymax=q95), alpha=.2, color=NA,fill="red"
  )+
  theme_bw()+
  coord_fixed()+
  geom_line(pdat2,mapping = aes(x=X,y=P), color="blue")+
  ylab("Recalibrated Probabilities")+
  xlab("Predictions")+
  theme(
    legend.position="bottom",
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
                              size = 12
    ),
    #strip.text.y = element_blank(),
    axis.title=element_text(
      #face = "bold",
      size=12
    )
  )


ggsave("Fig_calibrationcurve.pdf", height = 4, width = 9.5)


# results for re-calibrated model
ehlre <- eHL(y=c3recaliso$Y, P=c3recaliso$P, s=s_, boot=bbb)
ehlre$HLe
chlre<- ResourceSelection::hoslem.test(y=c3recaliso$P,x=c3recaliso$Y, g = 10)
chlre$p.value <- 1 - pchisq(as.numeric(chlre$statistic), 10) # OOS HL


#calibration_bands(x=c3recaliso$P,y=c3recaliso$Y,digits = 3,method = "round", nc=F)

#eHL(y=c3recaliso$Y, P=c3recaliso$P, s=s_, boot=500)$HLe # see: gets more stable with larger boot

# results for non re-calibrated model
ehlnore <- eHL(y=pairs.c3$Y, P=pairs.c3$P, s=s_, boot=bbb)
ehlnore$HLe
chlnore <- ResourceSelection::hoslem.test(y=pairs.c3$P,x=pairs.c3$Y, g = 10)
chlnore$p.value <- 1 - pchisq(as.numeric(chlnore$statistic), 10)
#calibration_bands(x=pairs.c3$P,y=pairs.c3$Y,digits = 3,method = "round", nc=F)

data.frame(
  "eHL" = rbind("plain"=ehlnore$HLe,"post-processed"=ehlre$HLe) %>% round(3),
  "HL" = rbind("plain"=chlnore$p.value,"post-processed"=chlre$p.value) %>% round(3)
)

ggpubr::ggarrange(
  calibration_bands(x=pairs.c3$P,y=pairs.c3$Y,digits = 3,method = "round", nc=F)%>% autoplot()+ ggtitle("plain"),
  calibration_bands(x=c3recaliso$P,y=c3recaliso$Y,digits = 3,method = "round", nc=F) %>% autoplot()+ ggtitle("post-processed")
)


# robustness checks ----
# results for robustness check 1

c1rob <- credit[(samples!=1),]
c2rob <- credit[(samples==3),]

# estimate logit in c1 data
logit.tr.rob <- glm(default ~  ., data = c1rob, family=binomial(link='logit'))
#predict in c2
pairs.c2.rob<- tibble(
  "Y" = c2rob$default, "P" = predict(logit.tr.rob,new =c2rob, type='response')
) %>%
  arrange(P)


ehlrob1<- eHL(y=pairs.c2.rob$Y, P=pairs.c2.rob$P, s=s_, boot=bbb)
ehlrob1$HLe
chlrob1 <- ResourceSelection::hoslem.test(y=pairs.c2.rob$P,x=pairs.c2.rob$Y, g = 10)
chlrob1
chlrob1$p.value <- 1- pchisq(as.numeric(chlrob1$statistic), 10) # OOS HL
#calibration_bands(x=pairs.c2.rob$P,y=pairs.c2.rob$Y,digits = 3,method = "round", nc=F)


isofit <- isoreg(x=pairs.c2$P,y=pairs.c2$Y)
isolin <- approxfun(
  y=c(min(isofit$yf),isofit$yf,max(isofit$yf)),
  x=c(0,isofit$x,1),
  method = "linear",
  ties = "mean"
)

saveRDS(tibble(P=isolin(pairs.c3$P),Y=pairs.c3$Y), "c3recaliso_oneshot.rds")

# at least one is non-empty if the e-value is INF
tibble(P=isolin(pairs.c3$P),Y=pairs.c3$Y) %>% filter(P==0&Y==1)
tibble(P=isolin(pairs.c3$P),Y=pairs.c3$Y) %>% filter(P==1&Y==0)

ehlrob2<- eHL(y=pairs.c3$Y, P=isolin(pairs.c3$P), s=s_, boot=bbb)
ehlrob2$HLe
chlrob2 <- ResourceSelection::hoslem.test(y=isolin(pairs.c3$P),x=pairs.c3$Y, g = 10)
chlrob2$p.value <- 1- pchisq(as.numeric(chlrob2$statistic), 10)
chlrob2
#calibration_bands(x=isolin(pairs.c3$P),y=pairs.c3$Y,digits = 3,method = "round", nc=F)

## Tab 3 -----

tabdat <-
  data.frame(
    "eHL (e-value)" = rbind(
      "plain"=ehlnore$HLe,
      "bagged recalibration"=ehlre$HLe,
      "plain with 4/5 in estimation set" = ehlrob1$HLe,
      "one-shot recalibration" = ehlrob2$HLe
    ) %>% round(3),
    "HL (p-value)" = rbind(
      "plain"=chlnore$p.value,
      "bagged recalibration" = chlre$p.value,
      "plain with 4/5 in estimation set" = chlrob1$p.value,
      "one-shot recalibration" = chlrob2$p.value
    ) %>% round(3)
  )


kbl(
  tabdat,
  booktabs = T,
  digits = 2,
  'latex'
) %>%
  kable_styling()

# \begin{table}
# \centering
# \begin{tabular}[t]{lrr}
# \toprule
# & eHL..e.value. & HL..p.value.\\
# \midrule
# plain & 6.95070e+28 & 0.00\\
# bagged recalibration & 6.14000e+00 & 0.11\\
# plain with 4/5 in estimation set & 9.58592e+22 & 0.00\\
# one-shot recalibration & 2.03500e+01 & 0.67\\
# \bottomrule
# \end{tabular}
# \end{table}

# instabilities ----

# results for re-calibrated model
ehlre <- eHL(y=c3recaliso$Y, P=c3recaliso$P, s=s_, boot=bbb)
#ehlre$HLe
chlre<- ResourceSelection::hoslem.test(y=c3recaliso$P,x=c3recaliso$Y, g = 10)
chlre
chlre$p.value <- 1 - pchisq(as.numeric(chlre$statistic), 10) # OOS HL


#calibration_bands(x=c3recaliso$P,y=c3recaliso$Y,digits = 3,method = "round", nc=F)


# results for non re-calibrated model
ehlnore <- eHL(y=pairs.c3$Y, P=pairs.c3$P, s=s_, boot=bbb)
#ehlnore$HLe
chlnore <- ResourceSelection::hoslem.test(y=pairs.c3$P,x=pairs.c3$Y, g = 10)
chlnore$p.value <- 1 - pchisq(as.numeric(chlnore$statistic), 10) # OOS HL
chlnore

# results for re-calibrated model
ehlre <- eHL(y=c3recaliso$Y, P=c3recaliso$P, s=s_, boot=bbb)
#ehlre$HLe
chlre<- ResourceSelection::hoslem.test(y=c3recaliso$P,x=c3recaliso$Y, g = 10)
chlre$p.value <- 1 - pchisq(as.numeric(chlre$statistic), 10) # OOS HL


q=10

instab_cHL_recal <-
  data.frame(
    type = "bagged recalibration",
    "Q1"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q1", m.bins=q)$PVAL,
    "RS"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "RS", m.bins=q)$PVAL,
    "Q2"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q2", m.bins=q)$PVAL,
    "Q3"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q3", m.bins=q)$PVAL,
    "Q4"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q4", m.bins=q)$PVAL,
    "equal"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "equal", m.bins=q)$PVAL

  )
instab_cHL_recal

instab_eHL_recal <-
  data.frame(
  "s14" = eHL(y=c3recaliso$Y, P=c3recaliso$P, s=1/4, boot=bbb)$HLe,
   "s13" = eHL(y=c3recaliso$Y, P=c3recaliso$P, s=1/3, boot=bbb)$HLe,
  "s12" =eHL(y=c3recaliso$Y, P=c3recaliso$P, s=1/2, boot=bbb)$HLe,
   "s23" =eHL(y=c3recaliso$Y, P=c3recaliso$P, s=2/3, boot=bbb)$HLe,
   "s34" =eHL(y=c3recaliso$Y, P=c3recaliso$P, s=3/4, boot=bbb)$HLe
  )
instab_eHL_recal


cbind(instab_cHL_recal, instab_eHL_recal)



# results for non re-calibrated model
ehlnore <- eHL(y=pairs.c3$Y, P=pairs.c3$P, s=s_, boot=bbb)
#ehlnore$HLe
chlnore <- ResourceSelection::hoslem.test(y=pairs.c3$P,x=pairs.c3$Y, g = 10)
chlnore$p.value <- 1 - pchisq(as.numeric(chlnore$statistic), 10) # OOS HL




chlnore
instab_cHL_norecal <-
  data.frame(
    type = "logit",
    "Q1"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q1", m.bins=q)$PVAL,
    "RS"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "RS", m.bins=q)$PVAL,
    "Q2"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q2", m.bins=q)$PVAL,
    "Q3"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q3", m.bins=q)$PVAL,
    "Q4"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q4", m.bins=q)$PVAL,
    "equal"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "equal", m.bins=q)$PVAL
  )

instab_eHL_norecal <-
  data.frame(
    "s14" = eHL(y=pairs.c3$Y, P=pairs.c3$P, s=1/4, boot=bbb)$HLe,
    "s13" = eHL(y=pairs.c3$Y, P=pairs.c3$P, s=1/3, boot=bbb)$HLe,
    "s12" =eHL(y=pairs.c3$Y, P=pairs.c3$P, s=1/2, boot=bbb)$HLe,
    "s23" =eHL(y=pairs.c3$Y, P=pairs.c3$P, s=2/3, boot=bbb)$HLe,
    "s34" =eHL(y=pairs.c3$Y, P=pairs.c3$P, s=3/4, boot=bbb)$HLe
  )
instab_eHL_norecal

cbind(instab_cHL_norecal, instab_eHL_norecal)




# results robustness 1

ehlrob1$HLe
chlrob1


instab_cHL_inclogit <-
  data.frame(
    type = "increased logit",
    "Q1"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q1", m.bins=q)$PVAL,
    "RS"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "RS", m.bins=q)$PVAL,
    "Q2"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q2", m.bins=q)$PVAL,
    "Q3"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q3", m.bins=q)$PVAL,
    "Q4"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q4", m.bins=q)$PVAL,
    "equal"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "equal", m.bins=q)$PVAL
  )

instab_eHL_inclogit <-
  data.frame(
    "s14" = eHL(y=pairs.c2.rob$Y, P=pairs.c2.rob$P, s=1/4, boot=bbb)$HLe,
    "s13" = eHL(y=pairs.c2.rob$Y, P=pairs.c2.rob$P, s=1/3, boot=bbb)$HLe,
    "s12" =eHL(y=pairs.c2.rob$Y, P=pairs.c2.rob$P, s=1/2, boot=bbb)$HLe,
    "s23" =eHL(y=pairs.c2.rob$Y, P=pairs.c2.rob$P, s=2/3, boot=bbb)$HLe,
    "s34" =eHL(y=pairs.c2.rob$Y, P=pairs.c2.rob$P, s=3/4, boot=bbb)$HLe
  )
instab_eHL_inclogit


cbind(instab_cHL_inclogit,instab_eHL_inclogit)


# results robustness 2
ehlrob2$HLe
chlrob2

q=10


instab_cHL_oneshot <-
  data.frame(
    type = "oneshot",
    "Q1"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q1", m.bins=q)$PVAL,
    "RS"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "RS", m.bins=q)$PVAL,
    "Q2"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q2", m.bins=q)$PVAL,
    "Q3"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q3", m.bins=q)$PVAL,
    "Q4"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q4", m.bins=q)$PVAL,
    "equal"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "equal", m.bins=q)$PVAL
  )
instab_cHL_oneshot

instab_eHL_oneshot <-
  data.frame(
    "s14" =eHL(y=pairs.c3$Y, P=isolin(pairs.c3$P), s=1/4, boot=bbb)$HLe,
    "s13" = eHL(y=pairs.c3$Y, P=isolin(pairs.c3$P), s=1/3, boot=bbb)$HLe,
    "s12" =eHL(y=pairs.c3$Y, P=isolin(pairs.c3$P), s=1/2, boot=bbb)$HLe,
    "s23" =eHL(y=pairs.c3$Y, P=isolin(pairs.c3$P), s=2/3, boot=bbb)$HLe,
    "s34" =eHL(y=pairs.c3$Y, P=isolin(pairs.c3$P), s=3/4, boot=bbb)$HLe
  )
instab_eHL_oneshot


tabdat2 <-
rbind(
  cbind(instab_cHL_recal, instab_eHL_recal),
  cbind(instab_cHL_norecal, instab_eHL_norecal),
  cbind(instab_cHL_inclogit,instab_eHL_inclogit),
  cbind(instab_cHL_oneshot,instab_eHL_oneshot))

kbl(
  tabdat2,
  booktabs = T,
  digits = 2,
  #'latex'
) %>%
  kable_styling()

## Fig 3----

save <- list()
for(q_ in 5:20 ){

  instab_cHL_recal <-
    data.frame(
      type = "bagged recalibration",
      "Q1"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q1", m.bins=q_)$PVAL,
      "RS"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "RS", m.bins=q_)$PVAL,
      "Q2"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q2", m.bins=q_)$PVAL,
      #"Q3"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q3", m.bins=q_)$PVAL,
      "Q4"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "Q4", m.bins=q_)$PVAL,
      "equal"=HL.binning(FC=c3recaliso$P, rlz=c3recaliso$Y, binning.method = "equal", m.bins=q_)$PVAL
      ) %>% mutate(bins = q_)

  instab_cHL_norecal <-
    data.frame(
      type = "logit",
      "Q1"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q1", m.bins=q_)$PVAL,
      "RS"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "RS", m.bins=q_)$PVAL,
      "Q2"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q2", m.bins=q_)$PVAL,
      #"Q3"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q3", m.bins=q_)$PVAL,
      "Q4"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "Q4", m.bins=q_)$PVAL,
      "equal"=HL.binning(FC=pairs.c3$P, rlz=pairs.c3$Y, binning.method = "equal", m.bins=q_)$PVAL
    )%>% mutate(bins = q_)

  instab_cHL_inclogit <-
    data.frame(
      type = "increased logit",
      "Q1"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q1", m.bins=q_)$PVAL,
      "RS"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "RS", m.bins=q_)$PVAL,
      "Q2"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q2", m.bins=q_)$PVAL,
      #"Q3"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q3", m.bins=q_)$PVAL,
      "Q4"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "Q4", m.bins=q_)$PVAL,
      "equal"=HL.binning(FC=pairs.c2.rob$P, rlz=pairs.c2.rob$Y, binning.method = "equal", m.bins=q_)$PVAL
    )%>% mutate(bins = q_)

  instab_cHL_oneshot <-
    data.frame(
      type = "oneshot",
      "Q1"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q1", m.bins=q_)$PVAL,
      "RS"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "RS", m.bins=q_)$PVAL,
      "Q2"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q2", m.bins=q_)$PVAL,
      #"Q3"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q3", m.bins=q_)$PVAL,
      "Q4"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "Q4", m.bins=q_)$PVAL,
      "equal"=HL.binning(FC=isolin(pairs.c3$P), rlz=pairs.c3$Y, binning.method = "equal", m.bins=q_)$PVAL
    )%>% mutate(bins = q_)

  tbs <- rbind(instab_cHL_recal,instab_cHL_norecal,instab_cHL_inclogit,instab_cHL_oneshot)
  save <- rlist::list.append(save,tbs)

}

do.call("rbind",save) %>%
  filter(type%in%c("oneshot", "bagged recalibration"))  %>%
  mutate(
    type=recode_factor(
      type,
      oneshot =  "Isotonic recalibration",
      'bagged recalibration' = "Bagged isotonic recalibration",
      .ordered = T)
    ) %>%
  pivot_longer(
    cols = "Q1":"equal",
    names_to = "test"
    ) %>%
  arrange(desc(type)) %>%
  ggplot()+
  geom_histogram(
    mapping = aes(x =value),
    bins = 20,
    boundary=0,
    fill=alpha("blue",.5),
    inherit.aes=FALSE,
    color="black"
    )+#, y = ..ndensity..
  facet_wrap(type ~ .)+
  theme_bw()+
  theme(
    plot.margin = unit(c(.5,.5,.5,.5), "cm"),
    legend.position="bottom",
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
    strip.text = element_text(face = "bold",
                              size = 12
    ),
    #strip.text.y = element_blank(),
    axis.title=element_text(
      #face = "bold",
      size=12
    )
  )+
  #ylab("density")+
  xlab("p-value")+
  scale_x_continuous(breaks=seq(0.05,1,0.15))+
  coord_fixed(ratio = 1/40)

ggsave("Fig_HLinstabilty.pdf", width = 9.5, height = 4)



## Tab 4 ----
overallinstb_groups <-
do.call("rbind",save) %>%
  #filter(type%in%c("oneshot", "bagged recalibration"))  %>%
  pivot_longer(
    cols = "Q1":"equal",
    names_to = "test"
  ) %>%
  group_by(type) %>%
    summarise(
      min = min(value),
      max = max(value)
    ) %>%
  mutate(
    interval = paste("[",round(min,3),",", round(max,3) ,"]", sep = "")
  ) %>%
  select(!c(min,max))

overallinstb_groups

overallinstb_groups %>%
  mutate(
    type=recode_factor(
      type,
      logit = "Logistic model",
      'increased logit'="Logistic model with increased estimation set",
      oneshot =  "Isotonic recalibration",
      'bagged recalibration' = "Bagged isotonic recalibration",
      .ordered = T)
  ) %>% arrange(type) %>%
  kbl(
    booktabs = T,
    digits = 3,
    'latex'
  ) %>%
  kable_styling()

# \begin{table}
# \centering
# \begin{tabular}[t]{ll}
# \toprule
# type & interval\\
# \midrule
# Logistic model & {}[0,0]\\
# Logistic model with increased estimation set & {}[0,0]\\
# Isotonic recalibration & {}[0,0.912]\\
# Bagged isotonic recalibration & {}[0,0.529]\\
# \bottomrule
# \end{tabular}
# \end{table}


## Tab 5 ----
# all pvalues
tabdat3 <-
  do.call("rbind",save) %>%
  select(type,bins,Q1,RS,Q2,Q4,equal)%>%
  filter(type%in%c("oneshot", "bagged recalibration"))  %>%
  mutate(
    type=recode_factor(
      type,
      oneshot =  "Isotonic recalibration",
      'bagged recalibration' = "Bagged isotonic recalibration",
      .ordered = T)
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = Q1:equal,
    names_glue = "{type}_{.value}"
  ) %>%
  select(
    bins, starts_with("Iso"),everything()
  ) %>%
  mutate(across(`Isotonic recalibration_Q1`:`Bagged isotonic recalibration_equal`, round, 2))

tabdat3

kbl(
  tabdat3,
  booktabs = T,
  #digits = 3,
  'latex',
  col.names = c("g",rep(
    c("Q1","RS","Q2","Q4","E")
    ,2)
  )
) %>%
  kable_styling()%>%
  add_header_above(c("", "Isotonic recalibration" = 5, "Bagged isotonic recalibration" = 5))

# \begin{table}
# \centering
# \begin{tabular}[t]{rrrrrrrrrrr}
# \toprule
# \multicolumn{1}{c}{} & \multicolumn{5}{c}{Isotonic recalibration} & \multicolumn{5}{c}{Bagged isotonic recalibration} \\
# \cmidrule(l{3pt}r{3pt}){2-6} \cmidrule(l{3pt}r{3pt}){7-11}
# g & Q1 & RS & Q2 & Q4 & E & Q1 & RS & Q2 & Q4 & E\\
# \midrule
# 5 & 0.34 & 0.46 & 0.09 & 0 & 0.24 & 0.08 & 0.05 & 0.09 & 0.00 & 0.11\\
# 6 & 0.16 & 0.60 & 0.22 & 0 & 0.31 & 0.10 & 0.21 & 0.18 & 0.26 & 0.33\\
# 7 & 0.24 & 0.56 & 0.01 & 0 & 0.00 & 0.02 & 0.16 & 0.01 & 0.00 & 0.37\\
# 8 & 0.59 & 0.55 & 0.17 & 0 & 0.38 & 0.11 & 0.20 & 0.17 & 0.02 & 0.15\\
# 9 & 0.20 & 0.53 & 0.06 & 0 & 0.02 & 0.08 & 0.20 & 0.07 & 0.00 & 0.26\\
# \addlinespace
# 10 & 0.26 & 0.67 & 0.19 & 0 & 0.36 & 0.20 & 0.11 & 0.22 & 0.00 & 0.10\\
# 11 & 0.27 & 0.33 & 0.08 & 0 & 0.77 & 0.06 & 0.09 & 0.10 & 0.00 & 0.08\\
# 12 & 0.15 & 0.91 & 0.19 & 0 & 0.02 & 0.10 & 0.18 & 0.21 & 0.01 & 0.11\\
# 13 & 0.57 & 0.58 & 0.27 & 0 & 0.60 & 0.16 & 0.17 & 0.31 & 0.00 & 0.00\\
# 14 & 0.22 & 0.87 & 0.01 & 0 & 0.09 & 0.03 & 0.07 & 0.01 & 0.00 & 0.03\\
# \addlinespace
# 15 & 0.60 & 0.68 & 0.04 & 0 & 0.64 & 0.04 & 0.09 & 0.06 & 0.00 & 0.00\\
# 16 & 0.80 & 0.28 & 0.17 & 0 & 0.11 & 0.29 & 0.37 & 0.20 & 0.00 & 0.01\\
# 17 & 0.86 & 0.45 & 0.11 & 0 & 0.02 & 0.25 & 0.07 & 0.14 & 0.00 & 0.00\\
# 18 & 0.36 & 0.63 & 0.14 & 0 & 0.10 & 0.19 & 0.19 & 0.17 & 0.00 & 0.00\\
# 19 & 0.48 & 0.73 & 0.38 & 0 & 0.01 & 0.40 & 0.53 & 0.42 & 0.00 & 0.10\\
# \addlinespace
# 20 & 0.59 & 0.83 & 0.35 & 0 & 0.61 & 0.42 & 0.30 & 0.39 & 0.00 & 0.02\\
# \bottomrule
# \end{tabular}
# \end{table}


## Tab 6 ----
# table for stability of eHL in application
tabdat2 %>%
  select(type, s14:s34) %>%
  mutate(
  type=recode_factor(
    type,
    logit = "Logistic model",
    'increased logit'="Logistic model with increased estimation set",
    oneshot =  "Isotonic recalibration",
    'bagged recalibration' = "Bagged isotonic recalibration",
    .ordered = T)
  ) %>% arrange(type) %>%
  kbl(
    booktabs = T,
    digits = 1,
    'latex'
  ) %>%
  kable_styling()

# \begin{table}
# \centering
# \begin{tabular}[t]{lrrrrr}
# \toprule
# type & s14 & s13 & s12 & s23 & s34\\
# \midrule
# Logistic model & 9.212073e+35 & 5.227884e+33 & 2.294016e+28 & 3.639274e+21 & 7.363704e+17\\
# Logistic model with increased estimation set & 1.054202e+29 & 6.960192e+26 & 1.071366e+22 & 4.524770e+16 & 3.266091e+13\\
# Isotonic recalibration & 8.700000e+00 & 1.560000e+01 & 2.060000e+01 & 1.420000e+01 & 9.900000e+00\\
# Bagged isotonic recalibration & 2.700000e+00 & 4.200000e+00 & 5.000000e+00 & 4.300000e+00 & 3.700000e+00\\
# \bottomrule
# \end{tabular}
# \end{table}
