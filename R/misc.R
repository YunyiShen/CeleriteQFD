majority <- function(x){
    as.numeric(names(which.max(table(x))))
}