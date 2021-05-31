rm(list(ls()))
library('igraph')
library('data.table')
library('GenomicRanges')
library('chaser')
wd<-setwd("C:/Users/Delouche Alexis/Documents/pro/scolarité/L3 Poitiers 2020-2021/stage")
datapath <- paste0(wd,"/data/PCHiC_peak_matrix_cutoff5_Mono.tsv")
nwlist <- data.frame()

data_PCHiC=read.table(datapath, sep='\t',header=TRUE)
summary(data_PCHiC)

data_PCHiC=data_PCHiC[which(data_PCHiC$Mon>=5),]
data_PCHiC=data_PCHiC[which(data_PCHiC$oeName!="."),]

# data_PCHiCproc1=data.frame(data_PCHiC[,1], data_PCHiC[,2], data_PCHiC[,3])
# data_PCHiCproc2=data.frame(data_PCHiC[,5], data_PCHiC[,6], data_PCHiC[,7])
# G1Name=data.frame(paste0(data_PCHiC[,1],':', data_PCHiC[,2],'-', (data_PCHiC[,3])),data_PCHiC[,4])
# G2Name=data.frame(paste0(data_PCHiC[,5],':', data_PCHiC[,6],'-', (data_PCHiC[,7])),data_PCHiC[,8])
# 
# colnames(data_PCHiCproc1)=colnames(data_PCHiCproc2)=c('chr', 'start', 'end')
# colnames(G1Name) <- colnames(G2Name) <- c('coordinates','Name')


chromnetObjct <- chaser::make_chromnet(data_PCHiC[,c(1,2,3,6,7,8)])

# nwlist <- rbind( G1Name,G2Name)

# annotation <- Reduce(rbind,nwlist)
# network <- function(chromnetObjct){
#   return(graph.data.frame(chromnetObjct$edgesdf[,7:8]))
# }


featfn<- read.table(paste0(wd,"/data/Monocytes_Exp_EV.txt"), header=TRUE)
chromnetObjct <- chaser::load_features(chromnetObjct, featfn, type="features_table", missingv=0)


net_chas<-chaser::chas(chromnetObjct)
net_chas
random_chromnetObjct<-chaser::randomize(chromnetObjct,nrandom=10)
random_net_chas=lapply(random_chromnetObjct,chaser::chas)
random_net_chas


moy <- mean(c(unlist(random_net_chas),net_chas))
SD <- sd(c(unlist(random_net_chas),net_chas))
zscore <- (net_chas-moy)/SD
zscore

#listeDeZscore <- list()
# moy <- mean(unlist(random_net_chas))
# SD <- sd(unlist(random_net_chas))
# zscore <- (chromnetObjct[[j]]-moy)/SD
# for (j in length(random_net_chas)){
#   zscore <- (chromnetObjct[[j]]-moy)/SD
#   listeDeZscore[[j]] <- zscore
# }
# listeDeZscore

# for(i in length(random_net_chas)){
#   print(names(listeDeZscore[[i]]) <- names(random_net_chas[[1]][i]))
# }

# ZscoreFtRnaExpNetBcells <- c(zscore,length(random_net_chas))
# hist(ZscoreFtRnaExpNetBcells)

plotchas<-function(feat, chas, randfeat,randchas, xmin, xmax, ymin, ymax, title){
  
  if(missing(xmin)) {
    xmin=min(unlist(lapply(randfeat, colMeans)))
  } 
  if(missing(xmax)){
    xmax=max(unlist(lapply(randfeat, colMeans)))
  }   
  if(missing(ymin)) {
    ymin=min(unlist(randchas))
  } 
  if(missing(ymax)){
    ymax=max(unlist(randchas))
  }   
  if(missing(ymax)){
    title='ChAs Plot'
  }   
  cols=rainbow(length(chas))
  par(oma=c(1,1,1,1))
  plot( colMeans(feat),chas, pch=20,xlab='global gene expression', ylab='ChAs', cex.axis=1.5, cex.lab=1.5, xlim=c(xmin, xmax), ylim=c(ymin, ymax), col=cols, main=title)
  #text(  colMeans(feat),chas,labels=colnames(feat), pos=3,  cex=0.5, srt=90)
  for (i in 1:length(randchas)){
    points( colMeans(randfeat[[i]]),randchas[[i]],  col=cols)
  }
  #text(  colMeans(randfeat[[i]]),randchas[[i]],chas,labels=colnames(feat), pos=3,  cex=0.5, srt=90)
}

plotchas(chromnetObjct, net_chas, random_chromnetObjct, random_net_chas,xmin=1.45, xmax=2.0, ymin=-0.05, ymax=0.4, title='Genes assortativity for gene age')