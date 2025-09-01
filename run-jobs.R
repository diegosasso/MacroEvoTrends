rstudioapi::jobRunScript(
  path = "validate-hmm-multi.R",
  name = "Compute AIC Job",
  workingDir = getwd(),
  importEnv = TRUE
)

library(dplyr)

#---- A1
dir <- 'R-hmm/results'
file <- 'AIC_A1-ch100-reg1.rds'
path <- file.path(dir, file)
AIC_A1 <- readRDS(path)
AIC_A1 <- as_tibble(AIC_A1)
class(AIC_A1)


hist(AIC_A1$Q.r1 - AIC_A1$Q.r2.asym)
sum((AIC_A1$Q.r1 - AIC_A1$Q.r2.asym) > 2)

hist(AIC_A1$Q.r1 - AIC_A1$Q.r2.sym)
sum((AIC_A1$Q.r1 - AIC_A1$Q.r2.sym) > 2)

#---- A2
dir <- 'R-hmm/results'
file <- 'AIC_A2-ch100-reg2-asym.rds'
path <- file.path(dir, file)
AIC_A2 <- readRDS(path)
AIC_A2 <- as_tibble(AIC_A2)
AIC_A2


hist(AIC_A2$Q.r2.asym - AIC_A2$Q.r1)
sum((AIC_A2$Q.r2.asym - AIC_A2$Q.r1) < -2)

hist(AIC_A2$Q.r2.asym - AIC_A2$Q.r2.sym)
sum((AIC_A2$Q.r2.asym - AIC_A2$Q.r2.sym) < -2)

hist(AIC_A2$Q.r1- AIC_A2$Q.r2.sym)


#---- A3
dir <- 'R-hmm/results'
file <- 'AIC_A3-ch100-reg2-sym.rds'
path <- file.path(dir, file)
AIC_A3 <- readRDS(path)
AIC_A3 <- as_tibble(AIC_A3)
AIC_A3


hist(AIC_A3$Q.r2.sym - AIC_A3$Q.r1)
hist(AIC_A3$Q.r2.sym - AIC_A3$Q.r2.asym)


