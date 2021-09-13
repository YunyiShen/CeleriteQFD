majority <- function(x){
    as.numeric(names(which.max(table(x))))
}

error_decoding <- function(decoding, states){
    res <- data.frame(matrix(NA,1,7))
    colnames(res) <- c("n_detected","n_injected","TP","FP","FN","SEN","SPC")

    states <- states > 1
    decoding <- decoding > 1

    location_inj <- which(states)
    res$n_injected <- (length(location_inj)!=0)+sum((location_inj[-1]-location_inj[-length(location_inj)])!=1)
    location_decoding <- which(decoding)
    res$n_detected <- (length(location_decoding)!=0)+sum((location_decoding[-1]-location_decoding[-length(location_decoding)])!=1)

    location_TP <- which( states & decoding)
    res$TP <- (length(location_TP)!=0)+sum((location_TP[-1]-location_TP[-length(location_TP)])!=1)
    location_FP <- which(!states & decoding)
    
    res$FP <- (length(location_FP)!=0)+sum((location_FP[-1]-location_FP[-length(location_FP)])!=1)
    res$FN <- res$n_injected-res$TP
    res$SEN <- res$TP/res$n_injected
    res$SPC <- res$TP/res$n_detected
    return(res)

}

