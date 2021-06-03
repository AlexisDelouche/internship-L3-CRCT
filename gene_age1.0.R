rm(list(ls()))
library('igraph')
library('data.table')
library('GenomicRanges')
library('chaser')
#defining the working directory
wd<-setwd("C:/Users/Delouche Alexis/Documents/pro/scolarité/L3 Poitiers 2020-2021/stage")
#defining the path to the network
datapath <- paste0(wd,"/data/PCHiC_peak_matrix_cutoff5_Mono.tsv")
#making a table with the network info
data_PCHiC=read.table(datapath, sep='\t',header=TRUE)

#here we filter the table to keep only the lines were the monocyte columns as a score above 5
data_PCHiC=data_PCHiC[which(data_PCHiC$Mon>=5),]
#here we filter data to not keep the lines were the "other End Name" column as a value"."
data_PCHiC=data_PCHiC[which(data_PCHiC$oeName!="."),]

#Now we load the network to make a chromnet object (Chr start end//Chr start en)
chromnetObjct <- chaser::make_chromnet(data_PCHiC[,c(1,2,3,6,7,8)])

#we read the feature table
featfn<-read.table(paste0(wd,"/data/HSgeneAges.tab"), header=TRUE)
#we select the chr start end and feature columns
featfn<-data.frame(featfn[,c(2,3,4,9)])
#we load the translate table (gene age are associated to a number)
corresp_table<-read.table(paste0(wd,"/data/table_age_gene_num.txt"),header=TRUE)


# we create a new data frame chr start en and feature where feature is the same number 
#for a same gene age 
featfn_new<-merge(featfn,corresp_table)
#removing the first column so that the data frame can be used with "feature table " method from
#from load_feture
featfn_new<-featfn_new[,-1]
#converting the features into characters
featfn_new$value<-as.character(featfn_new$value)

# loading the features
chromnetObjct <- chaser::load_features(chromnetObjct, featfn_new, type="features_table", missingv=18,auxfun = "min")

#calculting the chas score for gene age 
net_chas<-chaser::chas(chromnetObjct,method='categorical')
#randomizing chromnet to asses the previous calcul of chas 
random_chromnetObjct<-chaser::randomize(chromnetObjct,nrandom=10)
random_net_chas=lapply(random_chromnetObjct,chas, method='categorical')


#calcul for the zscore and plot
moy <- mean(c(unlist(random_net_chas),net_chas))
SD <- sd(c(unlist(random_net_chas),net_chas))
zscore <- (net_chas-moy)/SD
zscore
barplot(zscore, col="red")


#function to plot assortativity and randomized assortativities
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
  #plotting the mean value of gene age? Why? What does it tell us?  
  plot(colMeans(feat),chas, pch=20,xlab='Gene age', ylab='ChAs', cex.axis=1.5, cex.lab=1.5, xlim=c(xmin, xmax), ylim=c(ymin, ymax), col=cols, main=title)
  text(colMeans(feat),chas,labels=colnames(feat), pos=3,  cex=0.5)
  for (i in 1:length(randchas)){
    points( colMeans(randfeat[[i]]$features),unlist(randchas[[i]]),  col=cols)
  }
}
#making the features as number again to be able to plot
chromnetObjct$features <- apply(chromnetObjct$features, MARGIN = c(1,2), FUN = as.integer)
# same for random
for (i in 1:length(random_chromnetObjct)){
  random_chromnetObjct[[i]]$features <- apply(random_chromnetObjct[[i]]$features,MARGIN = c(1,2), FUN = as.integer)
}
plotchas(chromnetObjct$features, net_chas, random_chromnetObjct, random_net_chas,xmin=4, xmax=7, ymin=-0.05, ymax=0.3, title='Genes assortativity for gene age')