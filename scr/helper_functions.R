# ---------------------------------------
# Author: Mariella Gregorich
# Date: 14/05/2020
# Info: Helper functions
# ----------------------------------------

library(openxlsx)

extract_nameinfo <- function(filename){
  # Extract patient ID, CDtype and the stage at which sample was taken (e.g baseline, late, MLR as low and high)
  # from the text filename
  filename_new <- gsub("_1","1", filename)
  filename_new <- gsub("_2","2", filename_new)
  
  patID <- gsub(".*vdj.(.+)_CD.*", "\\1", filename_new)
  stage <- gsub(".*_(.+).*", "\\1", stringr::str_match(filename_new, "_(.*?)_clones")[2])
  CDtype <- gsub("_", "", stringr::str_match(filename_new, "_CD(.*?)_")[1])
  ifelse(grepl("1",stage,fixed=T), stage<-stage, stage<-paste0(stage,"1"))
  
  tcopies.TCMGRset <- c("R42","R43","R106","R132","R172","R34")
  Rejection <- ifelse(patID %in% tcopies.TCMGRset, 1,  0)
  
  
  info <- c(filename=filename,patID=patID,CDtype=CDtype, stage=stage, rejection=Rejection)
  rownames(info)<-NULL
  return(info)
}

read_sample <- function(filename){
  # Read in one sample
  sample <- read.table(paste0(path,filename),header=T, stringsAsFactors = F)
  return(sample)
}
readCloneTables <- function(directory){
  # Read in of all files in directory
  files <- grep("_clones.txt", dir(directory, full.names = T), fixed = T, val = T)
  names(files) <- sub(".*vdj\\.(.*)_clones.txt", "\\1", files)
  ct <- lapply(files, function(fn){
    read.table(fn, header = T, stringsAsFactors = F)
  })
  return(ct)
}

read_AllSamples <- function(directory){
  # Read in of all files in directory
  files <- get_RelevantFiles(directory)
  rel.files <- paste0(directory, files)
  
  names(rel.files) <- files
  all.ct <- lapply(rel.files, function(fn){
    read.table(fn, header = T, stringsAsFactors = F)
  })
  return(all.ct)
}

get_RelevantFiles <- function(directory, phenotype=NA){
  # Exclude files not relevant for the analysis but included in the directory (MLR-high, MLR replicates and R40 file)
  files <- list.files(directory, pattern="_clones.txt")
  files <- files[!(grepl("high", files)|grepl("low2", files)|grepl("low_2", files)|grepl("low_2", files)|grepl("R40", files))]
  if(phenotype%in%c("CD4","CD8")){files <- files[(grepl(phenotype, files))]}
  return(files)
}


DNAcodontable <- function(ntseq){
  ntseq <- gsub(paste0('.{',(str_length(ntseq)%%3),'}$'), '', ntseq)
  sst <- str_split(ntseq, "")[[1]]
  threeAA <- paste0(sst[c(TRUE, FALSE,FALSE)], sst[c(FALSE, TRUE,FALSE)], sst[c(FALSE,FALSE,TRUE)])
  
  aaseq <- paste(sapply(threeAA, function(x) Biostrings::GENETIC_CODE[[x]]), collapse = "")
  
  return(aaseq)
}

fixSampleNames <- function(ct){
  # Fix file names (ie some contain _1 some do not)
  names(ct) <- ifelse(
    grepl("_\\d$", names(ct)), 
    names(ct), 
    paste(names(ct), "1", sep = "_")
  )
  
  return(ct)
}

# Plots of the results
grouped_dotplot <- function(varname, yname="", title="", data, plot=T){
  df <- data.frame(data[,c("group", "stage", "rejection")])
  df$y <- data[,varname]
  if(length(unique(df$stage))>2){
    df$x <- 'PreTX'
    df[df$stage %in% "low1",]$x <- 'alloreactive'
    df[df$stage %in% "late1",]$x <- 'PostTX'
    labels.names <- c('PreTX','PostTX', "alloreactive")
  }else{
    df$x <- 1
    df[df$stage %in% "late1",]$x <- 2
    labels.names <- c('PreTX','PostTX')}

  plotname <- paste0("plot_comp_",varname)
  if(plot){png(file = paste0(out.path, plotname,".png", sep = ""), width = 7*600, height = 5*600, units = "px", res = 600)}
  g<-ggplot(df, aes(x=x, y=y, label=group))+
    geom_point(aes(x=stage,group=group, color= rejection, shape=rejection)) +
    # geom_label_repel(data=subset(df, stage %in% c("Pre-Transplant")),
    #                  direction    = "y",
    #                  hjust        = 1,
    #                  nudge_x = -0.5,
    #                  segment.size = 0.25,
    #                  box.padding=0.15,
    #                  size = 3,
    #                  segment.color = 'grey50') +    
    geom_line(aes(x=stage, group=group, color= rejection), alpha=1) +
    ggtitle(title) +
    scale_x_discrete(labels = labels.names, guide = guide_axis(angle = 45)) +
    scale_fill_manual(values=c("firebrick2", "deepskyblue4")) +
    scale_color_manual(values=c("firebrick2", "deepskyblue4")) +
    scale_shape_manual(values = c(1,6)) +
    labs(x="",y=yname) +
    theme_bw()
  print(g)
  if(plot){dev.off()}
  return(g)
}

grouped_boxplot <- function(varname, yname="", title="",data, plot=T){
  
  df <- data.frame(data[,c("group", "stage", "rejection")])
  df$y <- data[,varname]
  if(length(unique(df$stage))>2){
    df[df$stage %in% "BL1",]$x <- 'Pre-Transplant'
    df[df$stage %in% "low1",]$x <- 'alloreactive'
    df[df$stage %in% "late1",]$x <- 'Post-Transplant'
    labels.names <- c('Pre-Transplant','Post-Transplant', "alloreactive")
  }else{
    df$x <- 1
    df[df$stage %in% "late1",]$x <- 2
    labels.names <- c('Pre-Transplant','Post-Transplant')}
  
  ggplot(df, aes(x=stage, y=y, fill=rejection)) +
    geom_boxplot(position=position_dodge(0.8), notch=F, width = .4, alpha=0.8) +
    geom_dotplot(binaxis = "y", stackdir = "center",  position=position_dodge(0.8), dotsize=0.7, alpha=1) +
    stat_boxplot(geom ='errorbar', width = 0.2, position=position_dodge(0.8), alpha=0.8) +
    scale_y_continuous(expand=c(0,0), yname) +
    scale_x_discrete(expand=c(0,0),"", labels=labels.names) +
    scale_fill_brewer(palette="Dark2",name = "", labels = c("Controls", "Rejectors")) +
    ggtitle(title) +
    theme_bw() +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5), text=element_text(size=14),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  
}

mark_newclones <- function(df, phenotype="CD4", group.name){
  df$newBX <- 0
  load("TCopies_AA_merged.RData")
  merged <- merged[[group.name]]
  
  
  matches <- grep(phenotype, colnames(merged))
  matches <- unique(grep(paste(c("low", "BL"), collapse="|"), colnames(merged)[matches], value=TRUE))
  merged <- merged[,c("cdr3aa",matches)]
  
  ind <- apply(merged[,2:ncol(merged)], 1, function(x) all(is.na(x)))
  merged <- merged[ !ind, ]
  
  df$newBX[which(is.na(match(df$cdr3aa, merged$cdr3aa)))] <- 1
  
  return(df)
}


########### JSD_analysis.R #################
#calculate sample entropy using shannon.entropy()
#calculate entropy-weighted divergence between any two samples using jensen_shannon()
#make a table of pairwise jensen_shannon divergences using jsdReport(), looking at the top N clones
#print tables of pairwise JSD for all clones, top5000, top1000, top500, and top100 using jsdThresholded()

# data = read.table("fileName.tsv", header=T, sep="\t")
# colnames(data)
# jsdThresholded(data[,c(column1index, column2index, column3index)])

jsdThresholded<-function(data)
{
  write.table(jsdReport(data),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=5000),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=1000),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=500),sep="\t")
  cat("\n")
  write.table(jsdReport(data,topN=100), sep="\t")
}

jsdReport<-function(data, topN=-1)
{
  out<-matrix(nrow=ncol(data),ncol=ncol(data), dimnames=list(colnames(data), colnames(data)))
  for(i in 1:ncol(data)){
    for(j in 1:ncol(data)){
      if(topN==-1){
        out[i,j]<-jensen_shannon(data[,i], data[,j])
      }
      else{
        a<-order(data[,i],decreasing=TRUE)[1:topN]
        b<-order(data[,j],decreasing=TRUE)[1:topN]
        z<-data[union(a,b),]
        out[i,j]<-jensen_shannon(z[,i],z[,j])
      }
    }
  }
  return(out)
}

jensen_shannon <- function(x){
  ## JSD = H(0.5 *(p +q)) - 0.5H(p) - 0.5H(q)
  # H(X) = \sum x_i * log2(x_i)
  #  p = p[p >0 & q >0]
  #  q = q[p>0 & q>0]
  dfP <- read.table(paste0(data.path,x[1]), header=T, stringsAsFactors = F)
  dfQ <- read.table(paste0(data.path,x[2]), header=T, stringsAsFactors = F)
  
  p=dfP$count
  q=dfQ$count
  
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  
  jsd = Hj - 0.5*(Hp+Hq)
  #	cat(Hj, Hp, Hq, jsd, "\n")
  return(jsd)
}

shannon.entropy <- function(p)
{
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
#  p.norm <- p[p>0]/sum(p)
  return(-sum(log2(p.norm)*p.norm))
}
