
library(phangorn)
library(NELSI)
library(gtools)
library(stringr)

setwd("Documents/Astral")
load("updated.all.trees.Rdata")
load("references.Rdata")

#Write gene trees to file, then pull from directory for Astral
for (i in 1:length(maximalTrsOrd2)) {
  #class(maximalTrsOrd2[[i]]) <- "multiPhylo" #use if not using updated.all.trees.Rdata
  write.tree(maximalTrsOrd2[[i]], paste0(file = "s",i,".trs"))
}
tree_dir <- getwd() 
gene_files <- dir(tree_dir, pattern= ".trs") 
gene_files <- mixedsort(sort(gene_files))

#Estimate species trees via Astral
for (i in 1:length(gene_files)) {
  system(paste0("java -jar astral.5.7.3.jar -i s",i,".trs -o s",i,"_species.tre"))
}

#Pull species tree files 
tree_dir <- getwd() 
species_files <- dir(tree_dir, pattern= "species.tre")
species_files <- mixedsort(sort(species_files)) #Sort numerically, not alphabetically
species_files <- lapply(species_files, read.tree) #Read in as tree for RF.dist

#Calculate branch length variables using gene- and species trees
get_metrics <- function(tre, species_tr, study_num, study_name) {
  snum <- study_num
  sname <- study_name
  tbl <- sum(tre$edge.length)
  mbl <- mean(tre$edge.length)
  coeffVar_bl <- sd(tre$edge.length)/mean(tre$edge.length)
  stemminess <- sum(tre$edge.length[which(tre$edge[,2] > Ntip(tre))], na.rm = T) / sum(tre$edge.length, na.rm = T)
  species_distance <- RF.dist(species_tr, tre, check.labels = TRUE, normalize = TRUE)
  rootTip_lengths <- pathnode(midpoint(tre))
  coeffVar_rtl <- sd(rootTip_lengths$roottotippath)/mean(rootTip_lengths$roottotippath)
  class(tre$node.label) <- "numeric"
  mean_branchSup <- mean(tre$node.label, na.rm = TRUE)
  numberBranches <- Nedge(tre)
  numberTaxa <- Ntip(tre)
  return(c(snum, sname, tbl, mbl, coeffVar_bl, stemminess, coeffVar_rtl, rootTip_lengths, mean_branchSup, numberBranches, numberTaxa, species_distance))
}

topologicalMetrics <- list()
for (i in 1:34) {
  metrics <- lapply(maximalTrsOrd2[[i]], get_metrics, species_tr = species_files[[i]], study = study_numberOrd[[i]], study_name = references[[i]])
  topologicalMetrics <- c(topologicalMetrics, metrics)
}

#Create data frame 
data_df <- as.data.frame(do.call(rbind, lapply(topologicalMetrics, `length<-`, max(lengths(topologicalMetrics)))))
colnames(data_df) <- c("Study.Number", "Reference", "Total.Branch.Lengths", "Mean.Branch.Lengths", "CoV.Branch.Lengths", "Stemminess", "CoV.RootTotip.Lengths", "Mean.Branch.Support", "Number.of.Branches", "Number.of.Taxa", "Species.Distance")

data_df$RF.Distance <- as.numeric(as.character(data_df$RF.Distance))
data_df$Mean.Branch.Lengths <- as.numeric(as.character(data_df$Mean.Branch.Lengths))
data_df$CoV.Branch.Lengths <- as.numeric(as.character(data_df$CoV.Branch.Lengths))
data_df$CoV.RootTotip.Lengths <- as.numeric(as.character(data_df$CoV.RootTotip.Lengths))
data_df$Stemminess <- as.numeric(as.character(data_df$Stemminess))
data_df$Mean.Branch.Support <- as.numeric(as.character(data_df$Mean.Branch.Support))


#Calculate gene distances and add to data_df
get_GDmatrix <- function(studyTrees) {
  pwd = path.dist(studyTrees)
  pwdMatrix <- as.matrix(pwd)
  return(pwdMatrix)
}
GDmatrix <- lapply(maximalTrsOrd2, get_GDmatrix) #returns matrix of gene distances, use later if need for treespace

get_GDmeans <- function(pwdMatrix) {
  diag(pwdMatrix) <- NA
  pwdMeans <- colMeans(pwdMatrix, na.rm = TRUE)
  return(pwdMeans)
}
GDMeans <- lapply(GDmatrix, get_GDmeans) #returns the mean of each gene to another, use for glm analysis

data_df$Gene.Distance=unlist(GDMeans)
write.csv(data_df, file = "data_df.csv")



