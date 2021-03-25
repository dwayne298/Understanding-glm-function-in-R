# One-hot encode categorical variables
prepare_X_matrix <- function(data,factors){
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
    return(X)
}
