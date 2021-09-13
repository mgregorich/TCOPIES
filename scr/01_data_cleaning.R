# -----------------------------------------------------------
# Author: Mariella Gregorich, Constantin Aschauer
# Date: 
# Info: Data Cleaning of the TCOPIES; 
#       Scipt performs phenotype cleanup & fold change and frequency for alloreactive clone definition
#---------------------------------------------------------------


rm(list=ls())
pacman::p_load(stringr, openxlsx, reshape2)

# ----------- Initialisation -----------------------

# path = working directory holding the vdj..._clones.txt files, 
# new.path = path were the cleaned data vdj..._clones.txt files should go
path <- "../Daten/All_rerun/"
cor.path <- "../Daten/All_rerun_corr/clean_tcopies_fold5_freq1e-10_ar2/"
new.path <- "../Daten/All_rerun_corr/"
source("helper_functions.R")


# Define ambiguity ratio (ar), fold change (fc) and frequency (freq)
# Fc and freq can be set NULL to solely perform phenotype cleanup
# ar .... ratio in which CD4/CD8 cells have to occurr in to be assigned one of the two phenogroups
# fc .... threshold defined by frequency in unstimulated samples/frequency in unstimulated samples to be considered alloreactive
ambiguity.ratio=2
fold.change=5
frequency=0.0000000001



# ------ FUNCTIONS ----------------------------

retainCloneTableColumns <- function(ct, col.retain){
  # Retain only the columns freq and count
  return(
    lapply(
      ct, 
      function(x){
        return(x[, col.retain])
      }
    )
  )
}

prefixColnames <- function(df, prefix, col.exclude = NULL){
  # Adds a prefix to colnames e.g. typ of sample to freq : CD4_BL_1.count or CD4_BL_1.freq
  prefix.name <- make.names(prefix)
  include <- ! (colnames(df) %in% col.exclude)
  colnames(df)[include] <- paste(prefix, colnames(df)[include], sep = ".")
  return(df)
}

getIndivName <- function(sample.name){
  # Get the paIDs from a sample name
  return(sub("_.*", "", sample.name))
}

removeIndivName <- function(sample.name){
  # Removes the patID from string name of a sample
  return(sub("[^_]*_", "", sample.name))
}

useValAsName <- function(x){
  names(x) <- x
  return(x)
}

reduce <- function(fct, values, initial){
  result <- initial
  for(value in values){
    result <- fct(result, value)
  }
  
  return(result)
}

mergeIndivCloneTables <- function(ct, merge.col = "cdr3nt", col.retain = c("count", "freq")){
  merge.order <- c("CD4_BL_1", "CD8_BL_1", "CD4_low_1", "CD4_low_2", "CD8_low_1", "CD8_low_2", "CD4_late_1", "CD8_late_1")
  
  # Reduce dataframes in ct list to onyl contain the columns: cdr3nt, freq & count
  if(!is.null(col.retain)){
    ct <- retainCloneTableColumns(ct, c(merge.col, col.retain))
  }
  
  # Extract which samples are contained, the patIDs and split into list o patID and associated samples
  sample.names <- names(ct)
  indiv.names <- getIndivName(sample.names)
  indiv.groups <- split(sample.names, indiv.names)
  
  # Groups ... patID, samples ... phenotyp and BL/MLR/late mixture
  # Splits into groups, then samples of each group
  result <- lapply(
    indiv.groups, 
    function(member.names){
      samples.in.group <- ct[member.names]
      names(samples.in.group) <- removeIndivName(names(samples.in.group))
      # Add prefix to the columns count&freq of each sample in the group
      samples.in.group <- lapply(
        useValAsName(names(samples.in.group)), 
        function(sample.name){
          sample <- samples.in.group[[sample.name]]
          sample <- prefixColnames(sample, sample.name, merge.col)
          return(sample)
        }
      )
      # Get all unique CDR3 nucleotide sequences of one individual assessed for all sample types and phenotypes
      all.cdr3nt <- unique(
        unlist(
          lapply(
            samples.in.group, 
            function(sample) sample$cdr3nt
          ),
          use.names = F
        )
      )
      
      # Put them in dataframe
      all.cdr3nt <- data.frame(cdr3nt = all.cdr3nt, stringsAsFactors = F)
      
      # Stop function if one sample of a patient is missing
      stopifnot(merge.order %in% names(samples.in.group))
      
      # Returns dataframe indicating the freq and count of each cdr3nt sequence in all samples
      merged <- reduce(
        function(a, b){
          return(merge(a, b, by = merge.col, all.x = T))
        },
        samples.in.group[merge.order],
        all.cdr3nt
      )
      return(merged)
    }
  )
  # Result is a list for each patient associated with one df specifyig the freq and counts of each clonotype in each sample 
  return(result)
  
}
replaceNA <- function(df, ignore.cols = NULL){
  # Replaces NAs with 0
  replaceNAVect <- function(val){
    return(ifelse(is.na(val), 0, val))
  }
  
  cols <- seq(1, ncol(df))
  cols <- cols[!(colnames(df) %in% ignore.cols)]
  for(col in cols){
    df[, col] <- replaceNAVect(df[, col])
  }
  
  return(df)
}


# - Functions : ## TCR analysis- Aleksandar Obradovic

cleanup <- function(cd4, cd8, ratio = 2) {
  ## remove contaminated clones ie cells that have erroneously been assigned to CD4 or CD8
  ambi = (cd4 > 0 & cd8 > 0 & cd4 / cd8 > 1/ratio & cd4 / cd8 < ratio) 
  cd4exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 <= 1/ratio ) 
  cd8exclude = (cd4 > 0 & cd8 > 0 & cd4 / cd8 >= ratio )  
  return(cbind(ambi | cd4exclude, ambi | cd8exclude))
}

#remove the clones that need to be removed, output a cleaned-up dataframe
exclude <- function(data, excluderows) {
  data = data[excluderows == F, ]
  return(data)
}

#divide every column by its sum to convert counts to frequencies
normalize <- function(data) {
  nc = ncol(data)
  for (i in 1:nc) {
    data[,i] = data[,i ] / sum(data[,i])
  }
  return(data)
}

#list the actual sequences of all alloreactive clones
#before running, must do rownames(data)=data[,1]
#outputs a list, where allo[[1]] is the cd4 alloreactives and allo[[2]] is the cd8 alloreactives

listAlloreactive<-function(cd4,cd8, fold=5, freq1=.00001, ambiguityRatio=5) # rows must be indexed by clone ID
{
  # need only column 1 unstim and column 2 stim
  #cd4 = normalize(cd4)
  #cd8 = normalize(cd8)
  # cd4=df[,c("CD4_BL_1.count","CD4_low.count")]; cd8=df[,c("CD8_BL_1.count","CD8_low.count")]; fold=fc; freq1=freq;ambiguityRatio=ar
  rows1 = cleanup(cd4[,1], cd8[,1], ratio = ambiguityRatio)
  rows2 = cleanup(cd4[,2], cd8[,2], ratio = ambiguityRatio)
  
  cd4 = normalize(exclude(cd4, rows1[,1] | rows2[,1]))
  cd8 = normalize(exclude(cd8, rows1[,2] | rows2[,2]))
  
  allHvGcd4 = cd4[cd4[,2] > freq1 & cd4[,2] > cd4[,1] * fold, ]
  allHvGcd8 = cd8[cd8[,2] > freq1 & cd8[,2] > cd8[,1] * fold, ]
  test.freq = sum(sum(cd4[,2] <= freq1 & cd4[,2] != 0),sum(cd8[,2] <= freq1 & cd8[,2] != 0))

  return(list(rownames(allHvGcd4), rownames(allHvGcd8), test.freq))
}

#return a data frame containing only the alloreactive clones
reactiveClones <- function(data, fold=5, freq=1e-5) {
  # MLR = data[,2], BL = data[,1]
  if(is.null(freq)){freq<-0}
  rclones = data[data[,2] > freq & data[,2] > data[,1] * fold, ]
  return(rclones)
}

# ------------ FUNCTIONS END -----------------------------------------------

# Run the next 5 lines just the first time the script is executed, the load 'merged' for speed
# Read in data, and fix names and create aggregated clone tables
# ct <- readCloneTables(path)
# ct <- fixSampleNames(ct)
# merged <- mergeIndivCloneTables(ct)
# save(merged, file = paste0("TCopies_merged.RData"))


load(file = paste0("TCopies_merged.RData"))

#--------------------------------------

# Extract patient ids
pats <- gsub(".*vdj.(.+).*", "\\1",unique(sapply(list.files(path),getIndivName)))

# Clean up each group
clean_group <- function(group.name, fc, freq, ar=2, cpath=new.path){
  
  print(group.name)
  
  # group.name=pats[1]; fc=fold.change; freq=frequency; ar=ambiguity.ratio; cpath=clean.path
  df <- replaceNA(merged[[group.name]])
  
  df$CD4_low.count=df$CD4_low_1.count+df$CD4_low_2.count
  df$CD8_low.count=df$CD8_low_1.count+df$CD8_low_2.count
  
  # Calculate frequencies of added counts
  df$CD4_low.freq=df$CD4_low.count/sum(df$CD4_low.count)
  df$CD8_low.freq=df$CD8_low.count/sum(df$CD8_low.count)
  
  # Mark clonotype that fulfill certain conditions stated in function cleanup()
  excl.bl = cleanup(df$CD4_BL_1.count, df$CD8_BL_1.count, ratio = ar)
  excl.low = cleanup(df$CD4_low.count , df$CD8_low.count, ratio = ar)

  excl.cd4 <- excl.bl[,1]|excl.low[,1]
  excl.cd8 <- excl.bl[,2]|excl.low[,2]
  
  # CDR3 sequences to exclude in group from one phenotyp
  cd4.cdr3.excl <- df$cdr3nt[excl.cd4]
  cd8.cdr3.excl <- df$cdr3nt[excl.cd8]
  
  # Set clonotypes that are marked by TRUE to 0 since they fufilled the condition either for the BL or the MLR samples
  # Function listAlloreactive() already include cleanup() and exclude()
  listAllo <- listAlloreactive(df[,c("CD4_BL_1.count","CD4_low.count")] ,df[,c("CD8_BL_1.count","CD8_low.count")],fc,freq,ar)
  CD4.allo.rownr <- listAllo[[1]]
  CD8.allo.rownr <- listAllo[[2]]
  
  CD4.allo <- df$cdr3nt[as.numeric(as.character(CD4.allo.rownr))]
  CD8.allo <- df$cdr3nt[as.numeric(as.character(CD8.allo.rownr))]
  
  # Get alloreactive clones in BL sample
  CD4.allo.BL=intersect(df[(df$CD4_BL_1.count>0),]$cdr3nt,CD4.allo)
  CD8.allo.BL=intersect(df[(df$CD8_BL_1.count>0),]$cdr3nt,CD8.allo)
  CD4.allo.Bx=intersect(df[(df$CD4_late_1.count>0),]$cdr3nt,CD4.allo)
  CD8.allo.Bx=intersect(df[(df$CD8_late_1.count>0),]$cdr3nt,CD8.allo)
  
  # There should be no overlap between phenotypes
  if(length(intersect(CD4.allo.BL, CD8.allo.BL))>0){print("Same CDR3 sequences found in CD4 and CD8!")}

  # Extract newly detected clones in BX
  newBX_CD4 <- df$cdr3nt[which(df$CD4_BL_1.count==0 & df$CD4_late_1.count>0)]
  newBX_CD8 <- df$cdr3nt[which(df$CD8_BL_1.count==0 & df$CD8_late_1.count>0)]
  newBX_CD4.ratio <- length(newBX_CD4)/length(df$CD4_late_1.count>0)
  newBX_CD8.ratio <- length(newBX_CD8)/length(df$CD8_late_1.count>0)
  
  
  # Get filenames from directory separated by cd4 and cd8
  files.cd4 <- grep(list.files(path), pattern=paste0(group.name,"_CD4"),value=T)
  files.cd8 <- grep(list.files(path), pattern=paste0(group.name,"_CD8"),value=T)

  # Remove clones that were picked by cleanup(), mark alloreactive clones and save table
  correct_sample <- function(filename, vec_excl, vec_allo, vec_newBX){
    # filename=files.cd4[2];vec_excl=cd4.cdr3.excl; vec_allo=CD4.allo; vec_newBX = newBX_CD4
    print(filename)
    dfs <- read.table(paste0(path,filename),header=T, stringsAsFactors = F)
    
    # Exclude non productive clones and clones erroneously classified CD4/CD8
    nonprodclones <- dfs$cdr3nt[stringr::str_detect(dfs$cdr3nt,"[_*]")]
    dfs <- dfs[!(dfs$cdr3nt %in% vec_excl | dfs$cdr3nt %in% nonprodclones),]
    
    # Determine allo
    dfs$allo <- 0
    dfs[dfs$cdr3nt %in% vec_allo,]$allo <- 1
    if(grepl("late",filename, fixed=T)){
      dfs$newBX <- 0
      dfs[dfs$cdr3nt %in% vec_newBX,]$newBX <- 1}
    dfs$freq <- dfs$count/sum(dfs$count)
    dfs$cdr3nt.length <- str_length(dfs$cdr3nt)
    dfs$cdr3aa.length <- str_length(dfs$cdr3aa)
    
    write.table(dfs, file=paste0(cpath,filename),row.names = F)
    return(length(unique(dfs$cdr3aa)))
  } 
  
  # Apply function correct_sample() to every sample of group
  total_CD4 <- lapply(files.cd4, function(filename) correct_sample(filename,vec_excl=cd4.cdr3.excl, vec_allo=CD4.allo, vec_newBX = newBX_CD4))
  total_CD8 <- lapply(files.cd8, function(filename) correct_sample(filename,vec_excl=cd8.cdr3.excl, vec_allo=CD8.allo, vec_newBX = newBX_CD8))
  
  return(c(excl.cd4=length(cd4.cdr3.excl),allo.cd4=length(CD4.allo), allo.cd4.BL=length(CD4.allo.BL), 
           allo.cd4.Bx=length(CD4.allo.Bx), newBXclones.cd4=length(newBX_CD4), newBXclones.cd4.ratio=newBX_CD4.ratio,
           excl.cd8=length(cd8.cdr3.excl),allo.cd8=length(CD8.allo), allo.cd8.BL=length(CD8.allo.BL), 
           allo.cd8.Bx=length(CD8.allo.Bx), newBXclones.cd8=length(newBX_CD4), newBXclones.cd8.ratio=newBX_CD8.ratio))
}

# Create new directory
dir.create(paste0(new.path, "clean_tcopies_fold", fold.change,"_freq",frequency,"_ar",ambiguity.ratio), showWarnings = F)
clean.path <- paste0(new.path, "clean_tcopies_fold", fold.change,"_freq",frequency,"_ar",ambiguity.ratio,"/")

# Execute clean up
inf.groups <-sapply(pats, function(x) clean_group(x, fc=fold.change, freq=frequency, ar=ambiguity.ratio, cpath=clean.path))
write.xlsx(round(inf.groups,3), file=paste0(clean.path,"info_cleanup.xlsx"),row.names = T, col.names = T)



# tcopies.samplenames <- get_RelevantFiles(path)
# 
# # Clean up info
# extract_samplestats <- function(filename,path=clean.path){
#   print(filename)
#   df <- read.table(paste0(path,filename),header=T, stringsAsFactors = F)
#   df <- plyr::ddply(df[,c("cdr3aa","count","allo")], "cdr3aa", plyr::numcolwise(sum))
#   sample.info <- extract_nameinfo(filename=filename)
#   df <- mark_newclones(df, group.name=sample.info["patID"], phenotype=sample.info["CDtype"])
#   
#   res <- data.frame(group=sample.info["patID"], phenotype=sample.info["CDtype"], 
#                     stage=sample.info["stage"], rejection=sample.info["rejection"],
#                     total.clones.AA=length(unique(df$cdr3aa)), total.allo.AA=sum(df$allo==1),
#                     total.newBX.AA=sum(df$newBX==1))
#   return(res)
#   
# }
# 
# info_samples <- sapply(tcopies.samplenames, extract_samplestats)
# info_samples <- data.frame(matrix(unlist(info_samples),ncol=7,byrow=T))
# colnames(info_samples) <- c("group", "phenotype", "stage", "rejection", "total.clones","allo.clones", "newBX.clones")
# write.xlsx(info_samples, file=paste0(clean.path,"info_sample.xlsx"),row.names = T, col.names = T)
# 

# Corrected merged dataset
# Read in data, and fix names and create aggregated clone tables
# ct <- readCloneTables(cor.path)
# ct <- fixSampleNames(ct)
# merged <- mergeIndivCloneTables(ct)
# save(merged, file = paste0("TCopies_merged_corr1.RData"))

load(file = paste0("TCopies_merged_corr1.RData"))


# Clean up remaining CD4/CD8 mix ups

correct_remaining_mixups <- function(group.name, cor.path){
  # group.name=pats[1]; cor.path=cor.path
  print(group.name)
  df_group <- merged[[group.name]]
  ind <- apply(df_group[,2:ncol(df_group)], 1, function(x) all(is.na(x)))
  df_group <- df_group[ !ind, ]
  
  # Check if CD4 clones in CD8 clones
  df_group_CD4 <- df_group[,grepl("CD4", colnames(df_group))]
  df_group_CD4$CD4 <- apply(df_group_CD4, 1, function(x) any(!is.na(x)))
  df_group_CD4$cdr3nt <- df_group$cdr3nt
  df_group_CD8 <- df_group[,grepl("CD8", colnames(df_group))]
  df_group_CD8$CD8 <- apply(df_group_CD8, 1, function(x) any(!is.na(x)))
  df_group_CD8$cdr3nt <- df_group$cdr3nt
  
  large <- full_join(df_group_CD4[,c("cdr3nt","CD4")], df_group_CD8[,c("cdr3nt","CD8")], by="cdr3nt")
  
  # Get CD4/CD8 overlap cdr3 sequences
  df_freq_dbl <- df_group[which(apply(large[,c("CD4","CD8")],1,sum)==2),]
  df_freq_dbl$cdr3nt <- df_group[which(df_group_CD4$CD4 & df_group_CD8$CD8),]$cdr3nt

  # Get max frequency of overlapping cdr3 sequences in CD4 and CD()
  df_freq_dbl$cd4_max <- apply(df_freq_dbl[,grepl("CD4",colnames(df_freq_dbl))],1, function(x) max(x, na.rm = T))
  df_freq_dbl$cd8_max <- apply(df_freq_dbl[,grepl("CD8",colnames(df_freq_dbl))],1, function(x) max(x, na.rm = T))
  
  # Assign clonotype to CD4 if max(freq) in CD4 samples > 2 * max(freq) in CD8 and vice versa; delete if ambiguous
  df_freq_dbl$CD8 <-(df_freq_dbl$cd8_max > df_freq_dbl$cd4_max *2) 
  df_freq_dbl$CD4 <-(df_freq_dbl$cd4_max > df_freq_dbl$cd8_max *2) 

  print(paste0("Number of ambigious clonotypes = ",sum(apply(df_freq_dbl[,c("CD4","CD8")],1,sum)>2)))
  list_CD4 <- df_freq_dbl[df_freq_dbl$CD4==1,c("cdr3nt")]
  list_CD8 <- df_freq_dbl[df_freq_dbl$CD8==1,c("cdr3nt")]
    
  # Remove clones that were picked by cleanup(), mark alloreactive clones and save table
  correct_mixup <- function(filename, cd_list){
      
      print(filename)
      dfs <- read.table(paste0(cor.path,filename), header=T, stringsAsFactors = F)
      dfs_new <- dfs[!(dfs$cdr3nt %in% cd_list),]
      dfs_new$cdr3nt
      
      write.table(dfs_new, file=paste0(cor.path,filename), row.names = F)
      return(length(unique(dfs$cdr3aa)))
    } 
    
    # Get filenames from directory separated by cd4 and cd8
    files.cd4 <- grep(list.files(cor.path), pattern=paste0(group.name,"_CD4"),value=T)
    files.cd8 <- grep(list.files(cor.path), pattern=paste0(group.name,"_CD8"),value=T)
    
    # Apply function correct_sample() to every sample of group
    total_CD4 <- lapply(files.cd4, function(filename) correct_mixup(filename, cd_list=list_CD4))
    total_CD8 <- lapply(files.cd8, function(filename) correct_mixup(filename, cd_list=list_CD8))
}


inf.groups <- sapply(pats, function(x) correct_remaining_mixups(x, cor.path=cor.path))

# Corrected merged dataset
# Read in data, and fix names and create aggregated clone tables
ct <- readCloneTables(cor.path)
ct <- fixSampleNames(ct)
merged <- mergeIndivCloneTables(ct)

save(merged, file = paste0("TCopies_merged_corr.RData"))



# -------------------- Find and mark alloclones in the tissue samples r34 and r43 ------------------------
tissue.path <- "../Daten/Tissue/"
clean.path <- paste0(new.path, "clean_tcopies_fold", fold.change,"_freq",frequency,"_ar",ambiguity.ratio,"/")

df.34.CD4 <- read.table(paste0(clean.path,"vdj.R34_CD4_low_1_clones.txt"), header=T, stringsAsFactors = F)
df.34.CD8 <- read.table(paste0(clean.path,"vdj.R34_CD8_low_1_clones.txt"), header=T, stringsAsFactors = F)

df.43.CD4 <- read.table(paste0(clean.path,"vdj.R43_CD4_low_1_clones.txt"), header=T, stringsAsFactors = F)
df.43.CD8 <- read.table(paste0(clean.path,"vdj.R43_CD8_low_1_clones.txt"), header=T, stringsAsFactors = F)

filename <- "vdj.R34_bulk_tissue_1_clones.txt"
tissue_34 <- read.table(paste0(tissue.path,filename), header=T, stringsAsFactors = F)

filename <- "vdj.R43_bulk_tissue_1_clones.txt"
tissue_43 <- read.table(paste0(tissue.path,filename), header=T, stringsAsFactors = F)

# allo R34
allo.R34 <- unique(c(df.34.CD4[df.34.CD4$allo==1,]$cdr3nt, df.34.CD8[df.34.CD8$allo==1,]$cdr3nt))

# allo R43
allo.R43 <- unique(c(df.43.CD4[df.43.CD4$allo==1,]$cdr3nt, df.43.CD8[df.43.CD8$allo==1,]$cdr3nt))

tissue_34$allo <- tissue_43$allo <- 0
tissue_34[tissue_34$cdr3nt %in% allo.R34,]$allo <- 1
tissue_43[tissue_43$cdr3nt %in% allo.R43,]$allo <- 1

write.table(tissue_34, file=paste0(tissue.path,"vdj.R34_bulk_tissue_1_clones_allo.txt"), row.names = F)
write.table(tissue_43, file=paste0(tissue.path,"vdj.R43_bulk_tissue_1_clones_allo.txt"), row.names = F)



