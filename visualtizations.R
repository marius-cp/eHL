rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse);library(dplyr);library(kableExtra);library(ggforce);library(scales)

res <- readRDS("MCsim_quad_misspec_boot10_Sseq_05022021.rds")

# power plot -------------------------------------------------------------------

dt <- 
  data.frame(res) %>% 
  mutate(J = as.numeric(J),
         value=as.numeric(value),
         obs=as.numeric(obs),
         S=as.numeric(S),
         test=as.character(test),
         dec = if_else(test=="cHL", 
                       as.numeric(value<=0.05),
                       as.numeric(value>=20))) %>% 
  group_by(obs,J,test, S) %>% 
  summarize(rej = mean(dec)) %>% 
  filter(!S%in%c(.8,.7), obs>512 & J <=.105 | obs == 512) 

dt



ggplot(dt %>% filter(test=="eHL") , aes(x=J,y=rej, colour = factor(S), linetype=test))+
  geom_line(size=.55) +
  geom_point(dt%>% filter(test=="eHL"), mapping = aes(x=J,y=rej, shape=factor(S)), size=2)+
  theme_bw() + 
  theme(legend.position=c(.7,.1), legend.box="horizontal")+
  ylab("Rejection rate")+
  xlab("lack of linearity")+
  scale_color_brewer(palette = "Dark2")+
  geom_hline(aes(yintercept = 0.05), alpha=.25, size = .25, linetype="dashed")+
  geom_vline(aes(xintercept = 0.007), alpha=.25, size = .25, linetype="dashed") +
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x)))+
  facet_wrap(obs ~ ., scales = "free", ncol = 1)+
  labs(colour = "legend title") +
  geom_line(dt %>% filter(test=="oracle_eHL", S==0), 
            mapping=aes(J,rej, linetype = "oracle eHL,s=0"), size =.5, colour = "black")+
  geom_line(dt %>% filter(test=="oracle_eHL", S==.5), 
            mapping=aes(J,rej, linetype = "oracle eHL,s=1"), size =.5, colour = "mediumblue")+
  geom_line(dt %>% filter(test=="cHL"), 
            mapping=aes(J,rej, linetype = "classic HL"), size = .5, colour = "darkgray")+
  #labs(caption = "eHL test under differnt choices for $s$, the sample split. The dasehd line is the classical HL test based on 10 equally populated bins. The dotted line is the orcale eHL.")+
  labs(color  = "Proportion of data used to estimate q:", 
       linetype = "Test:", 
       shape = "Proportion of data used to estimate q:")+
  scale_linetype_manual(values=c("dotdash", "solid", "dotted", "dashed"))+
  scale_shape_manual(values = c(3,8,17,19,4))

ggsave("../plots/MCsim_interaction_boot10_Sseq_17042021.pdf", width = 9, height = 12 )

# growth rate plot -------------------------------------------------------------

grow.evals <- 
  data.frame(res) %>% 
  mutate(J = as.numeric(J),
         value=as.numeric(value),
         obs=as.numeric(obs),
         S=as.numeric(S),
         test=as.character(test)) %>% 
  filter(test %in% c("eHL","oracle_eHL")) %>% 
  group_by(obs,J,test, S) %>% 
  summarize(grow.eval = mean(log(value)))%>% 
  filter(!S%in%c(.8))


ggplot(grow.evals%>% filter(test=="eHL"), aes(J, grow.eval, colour = factor(S), linetype=test)) +
  geom_line(size=.75) +
  geom_point(grow.evals%>% filter(test=="eHL"), mapping = aes(x=J,y=grow.eval, shape=factor(S)), size=3)+
  geom_vline(aes(xintercept = 0.007), alpha=.25, size = .7, linetype="dashed") +
  geom_line(grow.evals %>% filter(test=="oracle_eHL", S == 0), 
           mapping=aes(J,grow.eval, linetype = "oracle eHL,s=0"), size =0.5, colour = "black")+
  geom_line(grow.evals %>% filter(test=="oracle_eHL", S == 0.5), 
            mapping=aes(J,grow.eval, linetype = "oracle eHL,s=0.5"), size =.5, colour = "black")+
theme_bw() +
  theme(legend.position="bottom")+#, legend.box="vertical")+
  labs(color  = "Proportion of data for q:", 
       linetype = "Test:", 
       shape = "Proportion of data for q:")+
  ylab("growth rate e-value")+
  xlab("lack of linearity")+
  scale_color_brewer(palette = "Dark2")+
  scale_shape_manual(values = c(20,4,18,15,8,3))+
  scale_linetype_manual(values=c("solid","dotted", "dashed"))+
  facet_wrap(obs ~ ., scales = "free", ncol = 2)

ggsave("../plots/MCsim_quadmiss_boot10_Sseq_10052021_GROW_eval.pdf", 
       width = 9, height = 8 )


# tab --------------------------------------------------------------------------
    
sizes <- dt %>% ungroup() %>%  filter(J == 0.007) %>% select(obs, test, rej, S) %>% 
  spread(test, rej)

sizes 


 tbl <-
cbind(  
sizes %>% select(-"oracle_eHL", -"eHL") %>% filter(S%in%c(NA)) %>% select(-S),
pivot_wider( sizes %>% select(-"cHL", -"eHL") %>% filter(complete.cases(.))
             , names_from= S, values_from = oracle_eHL) %>% select(-obs),
pivot_wider( sizes %>% select(-"cHL", -"oracle_eHL") %>% filter(complete.cases(.))
             , names_from= S, values_from = eHL) %>% select(-obs)
)

 kbl(data.table::data.table(tbl)*100, booktabs = T, digits = 2,'latex') %>% 
    kable_styling()%>% 
    add_header_above(c(" " = 1, " " = 1, "oracle eHL" = 2, "eHL" = 5))
  
  