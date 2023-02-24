rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreach);library(tidyverse);library(doParallel)
source("../nle.R")
source("../eHL.R")
source("../oracle_eHL.R")
MCsim <-
  readRDS("sim_HLsytle_oos.rds")

check <- sapply(MCsim, function(x) inherits(x, 'error'))
out <- MCsim[!check]
MCsim[check]
dat <-do.call("rbind", out)
dat
# rejection rates ----
dt <-
  dat %>%
  mutate(
    value = ifelse(test =="cHL",  (1 - pchisq(value, 9)),value),
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
    !(test=="eHL" & S%in%c("0.25","0.75")),
    !(test=="oracle_eHL" & S=="0.6"),
    obs %in%c(1024,4096),
    J>=.007
  )


# grwoth rates ----

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
    !(test=="eHL" & S%in%c("0.25","0.75")),
    !(test=="oracle_eHL" & S=="0.6"),
    obs %in%c(1024,4096),
    J>=.007
  )


comdat <-
  bind_rows(
    dt %>% mutate(value = rej, id = "a) rejection rates") %>% select(!c(setup,rej)),
    grow.evals %>% mutate(value = grow.eval, id = "b) growth rates (only for e-values)") %>% select(!c(grow.eval))
  ) %>%
  mutate(
    S_ =recode_factor(
      S,
      "0.333333333333333"="1/3",
      "0.5" = "1/2",
      "0.666666666666667"="2/3",
      .ordered = T
    )
  ) %>%
  mutate(
    test = ifelse(test=="oracle_eHL", "oracle eHL", test),
    labs = ifelse(test=="cHL", "cHL",paste(test," (s=",S_,")",sep = "")) %>% as.factor(),
    J=J-.007
  ) %>%
  mutate(
  labs1 = recode_factor(
      labs,
      "cHL"="cHL      ",
      "eHL (s=1/3)" = "eHL (s=1/3)",
      "eHL (s=1/2)" = "eHL (s=1/2)",
      "eHL (s=2/3)" = "eHL (s=2/3)      ",
      "oracle eHL (s=0)" = "oracle eHL (s=0)",
      "oracle eHL (s=1/2)" = "oracle eHL (s=1/2)",
      .ordered = T
    ))



  #%>% arrange(S)

sample_sizes = c(
  '1024' =  "n==1024",
  '2048' = "n==2048",
  '4096' = "n==4096",
  '8192' = "n==8192"
  )


ggplot(
  comdat,
  aes(
    x=(J),
    y=(value),
    colour = labs1, linetype=labs1,
  )
  )+
  geom_line(size=.55, linewidth=1.005) +
  theme_bw()+
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
    strip.text = element_text(face = "bold",
                              size = 12
    ),
    #strip.text.y = element_blank(),
    axis.title=element_text(
      #face = "bold",
      size=12
    )
  )+
  ylab("rate")+
  xlab("lack of linearity")+
  #scale_color_brewer(palette = "Dark2")+
  ggh4x::facet_nested(obs ~ id, scales = "free",independent ="all",
                      labeller = labeller(obs  = as_labeller(sample_sizes, label_parsed)))+
  guides(
    colour = guide_legend(nrow = 1),
    shape = guide_legend(nrow = 1),
    linetype = guide_legend(nrow = 1)
  )+
  labs(
    color  = "",
    linetype = "",
    shape = ""
  )+
  scale_shape_manual(
    values = c(NA,1,4,8,2,NA,NA,NA)
  )+
  scale_linetype_manual(
    values = c(11,"solid","solid","solid",21,21)
  )+
  scale_colour_manual(
    values = c("black", "#F1BB7B", "#FD6467", "#972D15","#7AA6DC","blue")
  )+
  geom_line(
    comdat %>% filter(labs1=="eHL (s=1/2)"),
    mapping = aes(x=J,y=value), color="#FD6467",show.legend =FALSE,
    size=.55, linewidth=1.005, linetype = "solid"
  )
  #scale_x_continuous(breaks=seq(0,0.1,0.015))


ggsave("Fig_Simulation.pdf", width = 9.5, height = 7)



