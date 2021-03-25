library(data.table)
library(dplyr)
library(ggplot2)

source("glm_functions.R")

set.seed(1)

df <- data.frame("Factor1"=runif(10^6,18,100)
                 ,"Factor2"=runif(10^6)
                 ,"Factor3"=round(runif(10^6,1,10),0)
                 ,"Factor4"=runif(10^6))

df$Factor2 <- ifelse(df$Factor2<0.2,"Detached"
                     ,ifelse(df$Factor2<0.7,"Flat","Semi"))
df$Factor4 <- ifelse(df$Factor4<0.1,"Unoccupied"
                     ,ifelse(df$Factor4<0.9,"Occupied","Weekend"))
df$Factor2 <- relevel(factor(df$Factor2),ref="Flat")
df$Factor4 <- relevel(factor(df$Factor4),ref="Occupied")

df$gamma_mean <- exp(0.21 + 
                       0.02*df$Factor1 + 
                       - 0.6*ifelse(df$Factor2=="Detached",1,0) + 0.3*ifelse(df$Factor2=="Semi",1,0) +
                       0.24*df$Factor3 +
                       -0.05*ifelse(df$Factor2=="Unoccupied",1,0) + 0.00001*ifelse(df$Factor2=="Weekend",1,0))
df$gamma_response <- rgamma(nrow(df),shape=15.21837,scale=df$gamma_mean*15.21837)

df$poisson_mean <- exp(0.15 + 
                         0.089*df$Factor1 + 
                         - 0.4*ifelse(df$Factor2=="Detached",1,0) + 0.0008*ifelse(df$Factor2=="Semi",1,0) +
                         0.13*df$Factor3 +
                         -0.5*ifelse(df$Factor2=="Unoccupied",1,0) + 0.01*ifelse(df$Factor2=="Weekend",1,0))
df$poisson_response <- rpois(nrow(df),df$poisson_mean)


start.time <- Sys.time()
test <- glm(gamma_response ~ Factor1+Factor2+Factor3+Factor4,family = Gamma(link="log"), data=df)
round(Sys.time() - start.time,2)

start.time <- Sys.time()
simple_glm(data.table(df),"gamma_response",c("Factor1","Factor2","Factor3","Factor4"),"gamma")
round(Sys.time() - start.time,2)

summary(test)

start.time <- Sys.time()
test <- glm(poisson_response ~ Factor1+Factor2+Factor3+Factor4,family = poisson(link="log"), data=df)
round(Sys.time() - start.time,2)

start.time <- Sys.time()
simple_glm(data.table(df),"poisson_response",c("Factor1","Factor2","Factor3","Factor4"),"poisson")
round(Sys.time() - start.time,2)

summary(test)
