if(!exists("prepare_X_matrix", mode="function")){
  source("prepare_X_matrix.R")
}

#page 40 in http://www.utstat.toronto.edu/~brunner/oldclass/2201s11/readings/glmbook.pdf
#Generalized Linear Models Second Edition by P. McCullagh and J.A. Nelder FRS
#This link is also helpful http://h2o-release.s3.amazonaws.com/h2o/master/1241/docs-website/datascience/glm.html
simple_glm <- function(data,response,factors,distribution){

    X <- prepare_X_matrix(data,factors)
    beta <- c(log(mean(data[[response]])),rep(0,ncol(X)-1))
    y <- matrix(data[[response]],ncol=1,byrow=F)
    beta_prev <- beta + 2
    WX <- X
    X_t <- t(X)
    if(distribution=="gamma"){
        X_tWX <- X_t%*%WX
        X_tWXI <- solve(X_tWX)      
    }
    while(sum(abs((beta/beta_prev) - 1)>0.001)>0){  
        beta_prev <- beta
        eta <- X%*%beta
        mu <- exp(eta)
        z <- eta + (y[,1] - mu)/mu
        if(distribution=="poisson"){
            for(i in 1:nrow(data)){
                z[i] <- z[i]*mu[i]
                for(j in 1:ncol(X)){
                  WX[i,j] <- X[i,j]*mu[i]
                }
            }
            X_tWX <- X_t%*%WX
            X_tWXI <- solve(X_tWX)
        }
        X_tWz <- X_t%*%z
        beta <- c(X_tWXI%*%X_tWz)
    }
    eta <- X%*%beta
    mu <- exp(eta)
    disp_par <- ifelse(distribution=="poisson"
                       ,1
                       ,(1/(nrow(data)-ncol(X))) * sum(((y[,1]-mu)^2)/(mu^2)))
    beta_se <- sqrt(diag(X_tWXI*disp_par))
    stat_value <- beta/beta_se
    if(distribution=="poisson"){
        p_value <- 2*pnorm(-abs(stat_value))
    }else{
        p_value <- 2*pt(-abs(stat_value),nrow(X)-ncol(X))
    }
    summary_table <- data.frame(factor=c("Intercept",colnames(X)[2:ncol(X)])
                                ,estimate=signif(beta,3)
                                ,standard_error = signif(beta_se,3)
                                ,stat_value=round(stat_value,3)
                                ,p_value=signif(ifelse(p_value<2*10^-16,2*10^-16,p_value),3)
                                ,Sig=ifelse(p_value<=0.001,"***",ifelse(p_value<=0.01,"**",ifelse(p_value<=0.05,"*",ifelse(p_value<=0.1,".",""))))
                                ,row.names = NULL)
    colnames(summary_table)[colnames(summary_table)=="stat_value"] <- ifelse(distribution=="poisson","z-value","t-value")
    print(summary_table)
    cat("\n","Dispersion Parameter is ",disp_par)
}
