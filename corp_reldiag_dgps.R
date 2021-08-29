rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse);library(reliabilitydiag);library(pdp)
source("nle.R")
J<-c(.001,.002,.007,0.06,.1,.2  )
no_interact <- nle_setup_solver(x = c(-1.5, 3, -3), 
                                fixed_points = c(0.05, 0.95),
                                j = J, iteration = 10,
                                order = 2)
list <- lapply(no_interact, '[[', 1)
pmeter <- data.frame(J , matrix(unlist(list), nrow=length(list), byrow=TRUE))
colnames(pmeter) <- c("J", "beta0", "beta1","beta2")          
pmeter
decomp <- list()
for (i in 1:length(J)) {
set.seed(2021)
j=pmeter$J[i]
n <- 50000
x1 <- runif(n, min = -3, max = 3)
x2 <- x1^2
xb <- pmeter$beta0[i] + pmeter$beta1[i] * x1 + pmeter$beta2[i] * x2
pr <- exp(xb)/(1+exp(xb))
u <- runif(n)
y <- 1*(u < pr)
dt <- data.frame(y,x1,x2)
mod.hl <- glm(y~x1, family = binomial(link = "logit"), data = dt)
dt$prediction <- predict(mod.hl, type = "response")
emp.cep <- predict(glm(y~x1, family = binomial(link = "logit"), data = dt), 
                   type = "response")
plot <- reliabilitydiag(X=emp.cep, y=dt$y, region.position = NA)
p <- autoplot(plot,  type = "miscalibration")
p <- p + 
  xlab("Prediction") +
  ggtitle(paste("j=", pmeter$J[i], sep=""))+
  theme(plot.margin=unit(c(1,1,1,1),"mm"))+
  annotate(
    "text",
    x = .15,
    y = .94,
    label = sprintf("MCB = .%03d",
                    round(summary(plot)$miscalibration * 1000)),
    color = "red"
  ) +
  annotate(
    "text",
    x = .15,
    y = .88,
    label = sprintf("DSC = .%03d",
                    round(summary(plot)$discrimination * 1000))
  ) +
  annotate(
    "text",
    x = .15,
    y = .82,
    label = sprintf("UNC = .%03d",
                    round(summary(plot)$uncertainty * 1000))
  )
decomp <- rlist::list.append(decomp, p)
}
pp <- grid.arrange(grobs=decomp, ncol = 3)
pp
ggsave("../plots/CORPreldiagsDGP.pdf", plot = pp,width = 9, height = 6 )









