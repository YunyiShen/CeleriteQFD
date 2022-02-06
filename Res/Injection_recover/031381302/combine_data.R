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
  res$TP <- count_success(location_TP, location_inj)
  location_FP <- which(!states & decoding)
  
  res$FP <- count_success(location_FP, location_decoding)
  res$FN <- res$n_injected-res$TP
  res$SEN <- res$TP/res$n_injected
  res$SPC <- res$TP/(res$TP+res$FP)
  return(res)
  
}

error_decoding_QFD <- function(decoding, states){
  res <- data.frame(matrix(NA,1,7))
  colnames(res) <- c("n_detected","n_injected","TP","FP","FN","SEN","SPC")
  
  states <- states == 2
  decoding <- decoding > 1
  
  location_inj <- which(states)
  res$n_injected <- (length(location_inj)!=0)+sum((location_inj[-1]-location_inj[-length(location_inj)])!=1)
  location_decoding <- which(decoding)
  res$n_detected <- (length(location_decoding)!=0)+sum((location_decoding[-1]-location_decoding[-length(location_decoding)])!=1)
  
  location_TP <- which( states & decoding)
  res$TP <- count_success(location_TP, location_inj)
  location_FP <- which(!states & decoding)
  
  res$FP <- count_success(location_FP, location_decoding)
  res$FN <- res$n_injected-res$TP
  res$SEN <- res$TP/res$n_injected
  res$SPC <- res$TP/(res$TP+res$FP)
  return(res)
  
}

count_success <- function(location_both, location_base){
  n_true <- (length(location_base)!=0)+sum((location_base[-1]-location_base[-length(location_base)])!=1)
  if(n_true==1 & length(location_both)>0) return(1)
  starting <- location_base[-1]
  ending <- location_base[-length(location_base)]
  tt <- which((location_base[-1]-location_base[-length(location_base)])!=1)
  
  starting <- starting[tt]
  starting <- c(location_base[1], starting)
  ending <- ending[tt]
  ending <- c(ending, location_base[length(location_base)])
  
  success <- 0
  for(i in 1:n_true){
    temp <- (location_both >= starting[i]) & (location_both<=ending[i])
    if(sum(temp)>0){
      success <- success + 1
    }
    
  }
  
  return(success)
  
}

clean_CHTC_data <- function(thedir){
  loss <- read.csv(paste(thedir,"inj_rec_loss.csv",sep = "/"),row.names = 1)
  QFD <- read.csv(paste(thedir,"inj_rec_QFD.csv",sep = "/"),row.names = 1)
  QFD <- as.numeric(QFD)
  gt <- read.csv(paste(thedir,"inj_rec_gtstate.csv",sep = "/"),row.names = 1)
  gt <- as.numeric(gt)
  sigmarule <- read.csv(paste(thedir,"inj_rec_3sigma.csv",sep = "/"),row.names = 1)
  sigmarule <- as.numeric(sigmarule)
  
  loss[1,9:15] <- error_decoding(QFD,gt[-1])
  loss[2,9:15] <- error_decoding(sigmarule, gt)
  loss$res_dir <- thedir
  
  return(loss)
  
}

all_res_dir <- list.dirs()[-c(1,2)]

all_res_df <- lapply(all_res_dir, clean_CHTC_data)

all_res_df1 <- Reduce(rbind, all_res_df)
write.csv(all_res_df1,"tiny_flares_full.csv", row.names = F)


library(reshape2)
library(ggplot2)

colnames(all_res_df1)[14:15] <- c("Sensitivity","Specificity")
plot_data <- melt(all_res_df1[,], measure.vars = colnames(all_res_df1)[9:15])




ggplot(data = plot_data[plot_data$variable %in% c("Sensitivity","Specificity"),], aes(x=method,y = value)) + 
  geom_point( alpha=.1, size=.5, position = "jitter")+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~variable) + 
  #ylab("Small flare detection") +
  ylab("") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        #axis.text.y = element_text(angle=180),
        plot.margin = margin(.15, .15, .15, .15, "cm")) + 
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  coord_flip()

ggsave("../injection-recovery-small-flare.pdf")



