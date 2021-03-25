
#page 40 in http://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf
#Generalized Linear Models Second Edition by P. McCullagh and J.A. Nelder FRS
#This link is also helpful http://h2o-release.s3.amazonaws.com/h2o/master/1241/docs-website/datascience/glm.html
gammaglm <- function(data,response,factors){
  factors_new <- c()
  for(i in 1:length(factors)){
    if(!is.numeric(data[[factors[i]]])){
      fac_new_len <- length(factors_new)
      ux <- unique(data[[factors[i]]])
      cat_fac <- ux[!(ux==ux[which.max(tabulate(match(data[[factors[i]]], ux)))][1])]
      factors_new <- c(factors_new,paste0(factors[i],gsub(" ", "", cat_fac, fixed = TRUE)))
      for(j in 1:length(cat_fac)){
        data[[factors_new[fac_new_len+j]]] <- ifelse(data[[factors[i]]]==cat_fac[j],1,0)
      }
    }else{
      factors_new <- c(factors_new,factors[i])
    }
  }
  X <- as.matrix(cbind("Intercept"=rep(1,nrow(data)),data[,factors_new,with=F]))
  b <- c(log(mean(data[[response]])),rep(0,length(factors_new)))
  y <- matrix(data[[response]],ncol=1,byrow=F)
  b_prev <- b + 2
  while(sum(abs((b/b_prev) - 1)>0.001)>0){
    b_prev <- b
    n <- X%*%b
    mu <- exp(n)
    z <- n + (y[,1] - mu)/mu
    
    X_D <- t(X)
    
    X2 <- X_D%*%X
    X2I <- solve(X2)
    
    y2 <- X_D%*%z
    answer <- X2I%*%y2
    b <- c(answer)
  }
  n <- X%*%b
  mu <- exp(n)
  disp_par <- (1/(nrow(data)-ncol(X))) * sum(((y[,1]-mu)^2)/(mu^2))
  est_se <- sqrt(diag(X2I*disp_par))
  t_value <- b/est_se
  p_value <- 2*pt(-abs(t_value),nrow(data)-ncol(X))
  
  summary_table <- data.frame(factor=c("Intercept",factors_new)
                              ,estimate=b
                              ,standard_error = est_se
                              ,t_value=t_value
                              ,p_value=ifelse(p_value<2*10^-16,2*10^-16,p_value)
                              ,Sig=ifelse(p_value<=0.001,"***",ifelse(p_value<=0.01,"**",ifelse(p_value<=0.05,"*",ifelse(p_value<=0.1,".",""))))
                              ,row.names = NULL)
  print(summary_table)
  cat("\n","Dispersion Parameter is ",disp_par)
}


poissonglm <- function(data,response,factors){
  factors_new <- c()
  for(i in 1:length(factors)){
    if(!is.numeric(data[[factors[i]]])){
      fac_new_len <- length(factors_new)
      ux <- unique(data[[factors[i]]])
      cat_fac <- ux[!(ux==ux[which.max(tabulate(match(data[[factors[i]]], ux)))][1])]
      factors_new <- c(factors_new,paste0(factors[i],gsub(" ", "", cat_fac, fixed = TRUE)))
      for(j in 1:length(cat_fac)){
        data[[factors_new[fac_new_len+j]]] <- ifelse(data[[factors[i]]]==cat_fac[j],1,0)
      }
    }else{
      factors_new <- c(factors_new,factors[i])
    }
  }
  X <- as.matrix(cbind("Intercept"=rep(1,nrow(data)),data[,factors_new,with=F]))
  b <- c(log(mean(data[[response]])),rep(0,length(factors_new)))
  y <- matrix(data[[response]],ncol=1,byrow=F)
  b_prev <- b + 2
  WX <- X
  Wz <- y
  while(sum(abs((b/b_prev) - 1)>0.001)>0){
    b_prev <- b
    n <- X%*%b
    mu <- exp(n)
    z <- n + (y[,1] - mu)/mu
    for(i in 1:nrow(data)){
      Wz[i] <- z[i]*mu[i]
      for(j in 1:ncol(X)){
        WX[i,j] <- X[i,j]*mu[i]
      }
    }
    
    X_D <- t(X)
    
    X2 <- X_D%*%WX
    X2I <- solve(X2)
    
    y2 <- X_D%*%Wz
    answer <- X2I%*%y2
    b <- c(answer)
  }
  n <- X%*%b
  mu <- exp(n)
  disp_par <- 1
  est_se <- sqrt(diag(X2I*disp_par))
  z_value <- b/est_se
  p_value <- 2*pnorm(-abs(z_value))
  
  summary_table <- data.frame(factor=c("Intercept",factors_new)
                              ,estimate=b
                              ,standard_error = est_se
                              ,z_value=z_value
                              ,p_value=ifelse(p_value<2*10^-16,2*10^-16,p_value)
                              ,Sig=ifelse(p_value<=0.001,"***",ifelse(p_value<=0.01,"**",ifelse(p_value<=0.05,"*",ifelse(p_value<=0.1,".",""))))
                              ,row.names = NULL)
  print(summary_table)
  cat("\n","Dispersion Parameter is ",disp_par)
}


