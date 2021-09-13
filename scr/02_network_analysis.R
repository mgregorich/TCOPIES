# -----------------------------------------------------------
# Author: Mariella Gregorich
# Date: 
# Info: Network Analysis
#---------------------------------------------------------------


rm(list=ls())

library(dplyr)
library(reshape2) 
library(ggplot2)          #Plotting package
library(stringdist)       # Calculate Levenshtein distance
library(igraph)
library(openxlsx)
library(tableone)
library(gtools)
library(ggrepel)
library(future.apply)
library(Matrix)
library(parallel)
library(data.table)
library(xtable)
library(stringr)
library(ineq)

source("helper_functions.R")

# Data cleaning parameter
ambiguity.ratio=2
fold.change=5
frequency=0.0000000001

# Network generation parameters
threshold_lv=1 # Levenshtein threshold
pheno ="CD4"  # CD4 or CD8
only.allo=F
large.network=F

# Foldernames
distance=paste0("LV",threshold_lv)
if(only.allo){foldername="allo"
}else{foldername="10ktop"}

# Define paths
path <- "../Daten/All_rerun_corr/"
clean.path <- paste0(path, "clean_tcopies_fold", fold.change,"_freq",frequency,"_ar",ambiguity.ratio,"/")
path <- clean.path

foldername_sim <- paste0("network_tcopies_",pheno,"_",foldername,"_",distance,"_final")
out.path <- "../Output/"
if(dir.exists( paste0(out.path,foldername_sim,"/"))){
  do.call(file.remove, list(list.files(paste0(out.path,foldername_sim,"/"), full.names = TRUE)))
}else{dir.create(paste(out.path,foldername_sim, sep=""))}
out.path <- paste0(out.path,foldername_sim,"/")



# -------------------------------------------------- Functions ---------------------------------------------------

delete.isolates <- function(g, mode = 'all') {
  N.iso=0
  isolates <- names(which(igraph::degree(g) < 1))
  while(length(isolates) !=0){
    N.iso <- N.iso + length(isolates)
    g <- delete.vertices(g, isolates)
    isolates <- names(which(igraph::degree(g) < 1))
  }
  
  return(list("graph"=g, "isolates"=N.iso))
}



compute_row <- function(i, v, threshold){
  res = stringdist(v[i], v[-(1:i)], method = "lv")
  node.pairs <- data.frame("Node.1"=v[i], "Node.2"=v[-(1:i)], "Leven"=0)
  res_threshold = which(res %in% c(1:threshold))
  res_threshold
  
  if(length(res_threshold)>0){
    node.pairs[res_threshold,] <- data.frame(Node.1 = v[i], Node.2 = v[i+res_threshold], Leven = res[res_threshold])
  }else{
   #res <-data.frame(Node.1 = NULL, Node.2 = NULL, Leven =NULL)
      }
  return(node.pairs)
}

compute_levenshtein <- function(v, threshold, workers = 5) {
  n = length(v) 
  
  if (workers > 1) {
    plan(multisession, workers = workers)  
    on.exit(plan(sequential))
  }
  res <- future_lapply(1:(n-1), function(i) compute_row(i, v, threshold))
  res[[length(res)+1]] <- c("Node.1"=v[1], "Node.2"=v[1], "Leven"=0)
  res[[length(res)+1]] <- c("Node.1"=v[n], "Node.2"=v[n], "Leven"=0)
  
  adjlist <- as.data.table(do.call(rbind, res))

  return(adjlist)  
}

get_eff <- function(i,g){return((1/(length(V(g)) - 1))*sum(1/distances(g, V(g)[i])[-i]))}


plot_network <- function(filename, only.allo=F, pheno, threshold, large.network=F){
  # @param patID ... Id of patient (e.g. R106)
  # @param CDtype ... "CD4" or "CD8"
  # @param stage ... "BL", "low" or "late"
  # @param only.converts ... plot only connected vertices
  # TEST: filename=tcopies.samplenames[3]; only.allo=F; threshold=1
  print(filename)
  
  # Load in files
  df <- read.table(paste0(path,filename),header=T, stringsAsFactors = F)
  nonprodclones <- df$cdr3aa[stringr::str_detect(df$cdr3aa,"[_*]")]
  df <- df[!(df$cdr3aa %in% nonprodclones),]
  sample.info <- extract_nameinfo(filename=filename)
  
  # Prune data to most abundant clones 
  df <- plyr::ddply(df[,c("cdr3aa","count","allo")], "cdr3aa", plyr::numcolwise(sum))
  df[df$allo >= 1,]$allo <- 1
  if(sample.info[4]%in%"low1"){df <- df[df$allo==1,]}
  total.clones <- length(unique(df$cdr3aa))
  sample.allo <- length(unique(df[df$allo==1,]$cdr3aa))
  df <- mark_newclones(df, group.name=sample.info["patID"], phenotype=pheno)
  df <- df[order(df$count, decreasing = T),]
  df <- df[str_length(df$cdr3aa) >= 9 ,]
  df <- df[1:min(10000,nrow(df)),]

  # Cut first and last 3 aas
  df$cdr3aa.full <- df$cdr3aa
  df$cdr3aa <- str_sub(df$cdr3aa, 4, -4)

  # 1. Adjancency matrix with Levenshtein condition for the AA sequences
  # Compute Levenshtein distance
  if(large.network){
    adj.list <- compute_levenshtein(unique(df$cdr3aa), threshold = threshold)
    if(any(adj.list$leven>1)){adj.list$leven[adj.list$leven>1]<-1}

    # Create adjacency matrix
    amat <- xtabs(Leven~Node.1+Node.2, adj.list, sparse = T)
    amat[is.na(amat)] <- 0
    diag(amat) <- 0
  }else{
    amat <- Matrix((stringdistmatrix(df$cdr3aa,df$cdr3aa, method="lv")==1)*1, sparse=T)
    colnames(amat)<-rownames(amat) <- df$cdr3aa
  }

  #saveRDS(adj.list, file=paste0(out.path, "adjmat_",sample.info["patID"],"_",sample.info["CDtype"],"_",sample.info["stage"],".Rda"))
  
  # Create graph
  net <- graph_from_adjacency_matrix(amat, mode="undirected", diag=F)
  
  # Convert adjacency matrix to igraph object
  e.density <- edge_density(net)*100
  total.iso <- sum(igraph::degree(net)==0)
  con.clones <- (nrow(df[df$cdr3aa %in% names(which(igraph::degree(net)>0)),])/nrow(df))*100
  allo.iso <- sum(df[!(df$cdr3aa %in% names(V(net))),]$allo)
  allo.con <- (nrow(df[df$cdr3aa %in% names(which(igraph::degree(net)>0)) & df$allo==1,])/nrow(df[df$allo==1,]))*100
  allo.clones <- nrow(df[df$allo==1,])

  
  # # Plot the network
  # if(sum(amat==1)!=0){
  #   ## PLOTTING THE TCR REPERTOIRE NETWORK
  #   # Set color of nodes and edges
  #   net.red <- delete.vertices(net, names(which(igraph::degree(net)==0)))
  #   E(net.red)$color <- 'black'
  #   vertex.color <- data.frame(clone=names(V(net.red)),group=df[which(names(V(net.red)) %in% df$cdr3aa),]$allo)
  #   
  #   # Layout algorithm
  #   l <- layout_with_fr(net.red, dim=3)
  #   
  #   # Plot
  #   file_pattern <- paste0(sample.info["patID"],"_",sample.info["CDtype"],"_",sample.info["stage"])
  #   plotname <- paste0("plot_",file_pattern)
  #   png(file = paste0(out.path, plotname,".png", sep = ""), width = 7*600, height = 5*600, units = "px", res = 600)
  #   plot(net.red, vertex.label=NA, remove.multiple = T, edge.arrow.size=0.1, vertex.size=2, vertex.color = as.factor(vertex.color$group),
  #        layout=l, main=paste0("PatID: " ,sample.info["patID"],' - TCR Repertoire'))
  #   dev.off()
  # }

  # Clustering
  cl <- clusters(net)
  cl.max <- max(cl$csize)
  cl.mean <- mean(cl$csize)
  cl.median <- median(cl$csize)
  cl.min <- min(cl$csize)
  
  cl.pairs <- sum(cl$csize==2)
  cl.mult <- sum(cl$csize>2)
  ratio.cl.mult.pairs <- cl.mult/cl.pairs
  
  me <- which(cl$csize == max(cl$csize))
  res <- unlist(split(names(cl$membership), cl$membership)[me])
  largest.cluster <- df[df$cdr3aa %in% res,]
  write.xlsx(largest.cluster, paste0(out.path, "data_LargestCluster_",sample.info["patID"],"_",sample.info["CDtype"],"_",sample.info["stage"],".xlsx"))
  LC.allo <- sum(largest.cluster$allo==1)/cl.max 
  LC.newBX <- sum(largest.cluster$newBX==1)/cl.max 
  LC.count <- mean(df$count)
  
  cluster_analysis <- function(df, me_cl, cl, exclude.cl.1=F){
    res <- unlist(split(names(cl$membership), cl$membership)[me_cl])
    cluster.data <- df[df$cdr3aa %in% res,]
    
    perc.allo <- ((sum(cluster.data$allo==1))/nrow(cluster.data))*100
    
    # Exclude clusters of size <=1
    if(exclude.cl.1){
          if(nrow(cluster.data) > 1) {return(perc.allo)
      }else{return(NA)}
    }else{
      return(perc.allo)
    }
  }
  
  plan(multisession, workers=10)
  clusters.perc.allo <- future_sapply(unique(cl$membership), function(x) cluster_analysis(df, x, cl))
  clusters.perc.allo.exclude.small.cl <- future_sapply(unique(cl$membership), function(x) cluster_analysis(df, x, cl, exclude.cl.1=T))
  plan(sequential)
  
  cl.perc.allo <- summary(clusters.perc.allo)
  cl.perc.allo.ex.1 <- summary(clusters.perc.allo.exclude.small.cl)
  
  cluster.info <- c( sample.info[2:5], cl$no, "csize"=summary(cl$csize), "p.allo"=cl.perc.allo, "p.allo.noiso"=cl.perc.allo.ex.1, allo.con, "#allo.subset"=allo.clones, "#allo.total"=sample.allo)
  write.xlsx(cluster.info, paste0(out.path, "table_clusterinfo_",sample.info["patID"],"_",sample.info["CDtype"],"_",sample.info["stage"],".xlsx"))
  
  # General metrics
  mean.cdr3length <- mean(str_length(df$cdr3aa))
  nconnected_nodes <- length(V(net))
  median.degree <- median(igraph::degree(net))
  mean.degree <- mean(igraph::degree(net))
  
  ratio.iso.con <- total.iso/length(V(net))
  ratio.iso.con.allo <- allo.iso/sum(df[match(names(V(net)), df$cdr3aa),]$allo==1)
  ratio.sing.doub <- sum(igraph::degree(net)==1)/sum(igraph::degree(net)==2)
  ratio.mult.sing <- sum(igraph::degree(net)>1)/sum(igraph::degree(net)==1)
  
  clustcoeff.glob <- transitivity(net)
  char.pathlength.unconT <- mean_distance(net, directed = F, unconnected = TRUE)
  char.pathlength.unconF <- mean_distance(net, directed = F, unconnected = FALSE)
  mod <- modularity(net, membership = cl$membership)
  ass <- assortativity.degree(net)
  dia <- diameter(net)
  rad <- radius(net)
  
  # Local network-specific metrics
  ev <- eigen_centrality(net)
  eigen.score <- ev$value
  
  plan(multisession, workers=detectCores()/2)
  local.eff <- future_sapply(1:length(V(net)), function(i) get_eff(i,g=net))
  plan(sequential)
  global.eff <- mean(local.eff)
  
  df.local <- data.frame("cdr3aa"=names(V(net)), "cluster"=cl$membership,"allo"=df[match(names(V(net)),df$cdr3aa),]$allo, 
                          "newBX"=df[match(names(V(net)), df$cdr3aa),]$newBX, "degree"=igraph::degree(net), 
                          "bc"=igraph::betweenness(net), "eigenc"=ev$vector, "localeff"=local.eff)
  df.local <- left_join(df.local, df[,c("cdr3aa", "cdr3aa.full")], by="cdr3aa")
  df.local <- df.local %>%
    relocate(cdr3aa.full, .before=cluster)
  
  median.degree.allo <- median(df.local[df.local$allo==1,]$degree)
  mean.degree.allo <- mean(df.local[df.local$allo==1,]$degree)
  
  # Local graph-theoretical features
  write.xlsx(df.local, paste0(out.path, "table_localfeatures_",sample.info["patID"],"_",sample.info["CDtype"],"_",sample.info["stage"],".xlsx"))
  file_pattern <- paste0(sample.info["patID"],"_",sample.info["CDtype"],"_",sample.info["stage"])
  plotname <- paste0("hist_distr_degree",file_pattern)
  png(file = paste0(out.path, plotname,".png", sep = ""), width = 7*600, height = 5*600, units = "px", res = 600)
  barplot(table(igraph::degree(net)), main="Frequency barplot of nodal network degree")
  dev.off()
  
  
  # Clone with the largest count
  largest.clone <- df[which(df$count==max(df$count)),]$cdr3aa
  lc.degree <- igraph::degree(net)[which(names(V(net))==largest.clone)]
  
  # Gini index measure
  Gini.vdegree <- Gini(df.local$degree)
  Gini.csize <- Gini(cl$csize)
  Gini.bc <- Gini(df.local$bc)
  Gini.localeff <- Gini(df.local$localeff)
  
  
  
  return(data.frame(group=sample.info["patID"], phenotype=sample.info["CDtype"], stage=sample.info["stage"], rejection=sample.info["rejection"],
                    sample.clones=total.clones, reduced.clones=length(df$cdr3aa), clones.network = vcount(net), edges.network=gsize(net),
                    sample.allo=sample.allo, pruned.allo = sum(df$allo==1), allo.net=nrow(df.local[df.local$allo==1,]), 
                    total.iso=total.iso, allo.iso=allo.iso,con.clones=con.clones, allo.con=allo.con,
                    perc.clones.nw=(vcount(net)/total.clones)*100, perc.allo.nw=(nrow(df.local[df.local$allo==1,])/sample.allo)*100,
                    median.degree=median.degree, mean.degree= mean.degree, median.degree.allo=median.degree.allo, mean.degree.allo,
                    mean.cdr3length=mean.cdr3length, 
                    ratio.iso.con=ratio.iso.con, ratio.iso.con.allo=ratio.iso.con.allo, ratio.sing.doub=ratio.sing.doub, ratio.mult.sing,
                    clustering.coefficient=clustcoeff.glob, 
                    cluster.no =cl$no, cl.max=cl.max, cl.mean=cl.mean, cl.median=cl.median, cl.pairs=cl.pairs, cl.mult=cl.mult, 
                    ratio.cl.mult.pairs=ratio.cl.mult.pairs, LC.allo=LC.allo, LC.newBX=LC.newBX, LC.count=LC.count,
                    char.pathlength.unconT=char.pathlength.unconT, char.pathlength.unconF=char.pathlength.unconF, 
                    global.efficiency=global.eff, edge.density=e.density, modularity=mod, assortativity=ass, diameter=dia, radius=rad, 
                    Gini.degree=Gini.vdegree, Gini.cl=Gini.csize, Gini.bc=Gini.bc, Gini.localeff=Gini.localeff))
  
}


# Generate networks for all groups and samples
tcopies.samplenames <- get_RelevantFiles(path, phenotype=pheno)
#tcopies.samplenames <- tcopies.samplenames[!str_detect(tcopies.samplenames, "low")]

gtmetrics <- lapply(tcopies.samplenames, function(x) plot_network(x, only.allo=only.allo, pheno=pheno, threshold = threshold_lv))
gtmetrics <- do.call("rbind",gtmetrics)
gtmetrics$stage  <- factor(gtmetrics$stage , levels = c("BL1","low1", "late1"))
gtmetrics <- gtmetrics[order(gtmetrics$stage),]
write.xlsx(gtmetrics, file=paste0(out.path, "table_tcopies_networks_",pheno,"_full.xlsx"))
#gtmetrics <- read.xlsx(paste0("C:/Users/Mariella/Documents/PhD/Projects/TCOPIES/Output/test_n10000_kk/table_tcopies_networks_n10000_full.xlsx"))

net.bl <- gtmetrics[gtmetrics$stage %in% "BL1",]
net.mlr <- gtmetrics[gtmetrics$stage %in% "low1",]
net.late <- gtmetrics[gtmetrics$stage %in% "late1",]


tbl_bl <- CreateTableOne(colnames(net.bl)[5:length(net.bl)], strata="rejection", data=net.bl, test=T)
tbl_mlr <- CreateTableOne(colnames(net.mlr)[5:length(net.mlr)], strata="rejection", data=net.mlr, test=T)
tbl_late <- CreateTableOne(colnames(net.late)[5:length(net.late)], strata="rejection", data=net.late, test=T)

tbl_netmetrics <- CreateTableOne(colnames(gtmetrics)[5:ncol(gtmetrics)], strata=c("rejection","stage"), data=gtmetrics,test=F, smd=F)
tbl_xlsx <- print(tbl_netmetrics, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.xlsx(tbl_xlsx, file=paste0(out.path, "table_tcopies_networks_",pheno,"_strat.xlsx"),col.names = TRUE, row.names = TRUE)
#xtable(print(tbl_netmetrics))

tbl_netmetrics_stage <- CreateTableOne(colnames(gtmetrics)[5:ncol(gtmetrics)], strata=c("stage"), data=gtmetrics,test=F, smd=F)
tbl_xlsx <- print(tbl_netmetrics_stage, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.xlsx(tbl_xlsx, file=paste0(out.path, "table_tcopies_networks_",pheno,"_strat_stage.xlsx"),col.names = TRUE, row.names = TRUE)
#tbl_netmetrics
#xtable(print(tbl_netmetrics))



