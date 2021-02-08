
#Function to generate a list with the names of the taxa that are maximised in each dataset

maximaltaxa_list <- function(genetree.list) {
  setsoftaxa <- list()
  for (i in 1:length(maximalTrsOrd2)) {
    setoftaxa <- maximalTrsOrd2[[i]][[1]]$tip.label
    setsoftaxa[length(setsoftaxa)+1] <- list(setoftaxa)
  }
  retrn(setsoftaxa)
}

#Function to trim the alignments so they only contain the maximal number of taxa for all datasets
trimAlignments <- function(datasets.dir, trimming.variable) {
  
  require(ape)
  
  for (i in 1:length(datasets.dir)) {
    
    alignment.dir <- datasets.dir[i]
    setwd(alignment.dir)
    alignment.files <- dir(alignment.dir)
    trimming.variablei <- trimming.variable[[i]]
    
    for (i in 1:length(alignment.files)) {
      if (isTRUE(grepl("fasta|fa", basename(alignment.files[i])))) {
        alignment <- read.dna(alignment.files[i], format = "fasta")
        if (all(trimming.variablei %in% rownames(alignment)) == TRUE){
          newalignment <- alignment[trimming.variablei,]
          write.dna(newalignment, file = paste0("filepath",basename(alignment.dir)), basename(alignment.files[i]), format = "fasta")
        }
      } else if (isTRUE(grepl("nex|nexus", basename(alignment.files[i])))) {
        alignment <- read.nexus.data(alignment.files[i])
        alignment <- nexus2DNAbin(alignment)
        alignment <- as.matrix(alignment)
        if (all(trimming.variablei %in% rownames(alignment)) == TRUE){
          newalignment <- alignment[trimming.variablei,]
          write.nexus.data(newalignment, file = paste0("filepath",basename(alignment.dir)), basename(alignment.files[i]))
        }
      } else {
        alignment <- read.dna(alignment.files[i]) 
        if (all(trimming.variablei %in% rownames(alignment)) == TRUE){
          newalignment <- alignment[trimming.variablei,]
          write.dna(newalignment, file = paste0("filepath",basename(alignment.dir)), basename(alignment.files[i]))
        }
      }
    }
    print(basename(alignment.dir))
  }
}

#IQ-TREE command
iqtreeCommand <- function(fileName) {
  iqtreepath <- "filepath"
  system(paste0(iqtreepath, " -st DNA -s ",fileName," -m GTR+G -alrt 1000 -quiet"))
}

#Pulls files from a dataset directory
pullfiles <- function(data.directory, filepattern) {
  setwd(data.directory)
  subset.files <- dir(data.directory, pattern = filepattern)
  return(c(directory = basename(data.directory), files = list(subset.files)))
}

#Saves IQ-TREE gene-trees from each dataset as multiphylo object
saveMultiPhylo <- function(directory) {
  
  require(ape)
  
  genetree.info <- pullfiles(directory, "treefile")
  genetrees <- lapply(genetree.info$files, read.tree)
  class(genetrees) <- "multiPhylo"
  ape::write.tree(genetrees, file = here("filepath", paste0(genetree.info$directory,"_allgene.trs")))
  return(genetrees)
}

#Astral command for species tree estimation
runAstral <- function(genetree.file, speciestreeDir) {
  filename <- gsub("_allgene.trs", "", basename(genetree.file))
  system(paste0("java -jar ~filepath ",genetree.file," -o ",speciestreeDir, filename,"_species.tre 2> ",speciestreeDir, filename,".log"))
}

#Function to calculate the branch based metrics from each gene tree
extractBLvariables <- function(tre, species_tr, study_num, study_name) {
  
  require(phangorn)
  require(NELSI)
  
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
  return(c(snum, sname, tbl, mbl, coeffVar_bl, stemminess, coeffVar_rtl, mean_branchSup, numberBranches, numberTaxa, species_distance))
}

#Calculates pairwise distance matrix for gene trees
getGDmatrix <- function(studyTrees) {
  pwd = path.dist(studyTrees)
  pwdMatrix <- as.matrix(pwd)
  return(pwdMatrix)
}

#Calculates mean of pairwise distances for each gene tree
getGDmeans <- function(pwdMatrix) {
  diag(pwdMatrix) <- NA
  pwdMeans <- colMeans(pwdMatrix, na.rm = TRUE)
  return(pwdMeans)
}

#Function to run regression models for individual datasets
runIndvGLMS <- function(df, study) {
  
  indvGDglms <- list()
  indvSDglms <- list()
  
  for (i in 1:30) {
    #gene distance
    glmi <- glm(log(Gene.Distance) ~ (Total.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support), data = subset(df, df$Study == i))
    indvGDglms[[(length(indvGDglms)+1)]] <- glmi
    
    #species distance
    sglmi <- glm(Species.Distance ~ (Total.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support), data= subset(df, df$Study == i))
    indvSDglms[[(length(indvSDglms)+1)]] <- sglmi
  }
  return(c(indvGDglms, indvSDglms))
}

#Function to run regression models across all datasets
runbigGLM <- function(df) {
  
  require(lme4)
  
  #all studies: gene distance (GD)
  GDglmLerB <- lmer(Gene.Distance ~ (Mean.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support) + (1 | Study.Number), weights = Number.of.Branches, data= df)
  GDglmLerT <- lmer(Gene.Distance ~ (Mean.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support) + (1 | Study.Number), weights = Number.of.Taxa, data= df)
  GDglmLer <- lmer(Gene.Distance ~ (Mean.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support) + (1 | Study.Number), data= df)
  
  #all studies: species distance (SD)
  SDglmLerB <- lmer(Species.Distance ~ (Mean.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support) + (1 | Study.Number), weights = Number.of.Branches, data= df)
  SDglmLerT <- lmer(Species.Distance ~ (Mean.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support) + (1 | Study.Number), weights = Number.of.Taxa, data= df)
  SDglmLer <- lmer(Species.Distance ~ (Mean.Branch.Lengths + CoV.RootTotip.Lengths + Stemminess + Mean.Branch.Support) + (1 | Study.Number), data= df)
  
  return(c(GDglm = GDglmLer, GDglmT = GDglmLerT, GDglmB = GDglmLerB, SDglm = SDglmLer, SDglmT = SDglmLerT, SDglmB = SDglmLerB))
}

#Extracts t-statistic from regression model summary
getTstat <- function(glm) {
  coef(summary(glm))[,"t value"]
}

#Creates dataframe with t-statistics and variables
tstatDF<- function(glmslist, glm) {
  
  require(reshape2)
  
  #extract t-stats from individual and all studies glms
  indvglm_tstats <- lapply(glmslist, getTstat) #get t-stats for individual models
  bigglm_tstats <- getTstat(glm) #get t-stats for big model
  
  #extract P value from individual study glms
  pvals <- lapply(glmslist, getPvalue)
  pvals <- lapply(pvals, round, digits = 2)
  
  #concatenate indidvidual and all studies t-stats
  allTstats <- indvglm_tstats
  allTstats[length(allTstats)+1] <- list(bigglm_tstats)
  
  #tstats create dataframe and format
  df <- as.data.frame(do.call(cbind, lapply(allTstats, `length<-`, max(lengths(allTstats)))))
  rownames(df) <- c("Intercept", "Total Branch Lengths", "CV Root to Tip Lengths", "Stemminess", "Mean Branch Support")
  df <- df[-1,] #Remove Intercept variable
  for (i in 1:31) {
    colnames(df)[i] <- paste0(i)
  } #Rename rows to study.number
  df$group <- row.names(df) #Put variables in group for melt transformation
  df.m <- melt(df, id.vars = "group")
  
  #pvalues create df and format
  df <- as.data.frame(do.call(cbind, lapply(pvals, `length<-`, max(lengths(pvals)))))
  rownames(df) <- c("Intercept", "Total Branch Lengths", "CV Root to Tip Lengths", "Stemminess", "Mean Branch Support")
  df <- df[-1,] #Remove Intercept variable
  v31 <- c("NA", "NA", "NA", "NA")
  df$V31=v31
  for (i in 1:31) {
    colnames(df)[i] <- paste0(i)
  } #Rename rows to study.number
  df$group <- row.names(df) #Put variables in group for melt transformation
  df.mPVAL <- melt(df, id.vars = "group")
  
  #merge tstats and pvalues
  df.m$pvals=df.mPVAL$value
  
  #make pvalues T or F based on significance level
  for (i in 1:length(df.m$pvals)) {
    df.m$pvals[i] <- ifelse(df.m$pvals[i] < 0.01, "T", "F")
  }
  return(df.m)
}

#Plots t-statistics
plotViolinplotsCOLOUR <- function(sd.df, gd.df, filename1) {
  
  require(tidyverse)
  require(gridExtra)
  require(scales)
  
  #species distance plot
  SDplot <- ggplot(sd.df, 
                   aes(x=group, y=value, 
                       colour = pvals,
                       shape = factor(ifelse(variable != "31", "Individual Study", "All Studies")))) +
    
    #data 
    geom_violin(data = subset(sd.df, sd.df$variable != '31'), colour = "black")  +
    geom_jitter(data = subset(sd.df, sd.df$variable == '31'), size = 4, colour = "#35B779FF") +
    geom_jitter(data = subset(sd.df, sd.df$variable != '31'), alpha = 0.5, size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed", show.legend = TRUE) +
    
    #axis co-ordinates
    coord_cartesian(ylim = c(-100, 25)) +
    
    #theme and plot labels
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold")) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + #Wrap variable labels
    labs(title = "A. Distance to Species Trees", x="Topological Measurement", y="T Statistic") +
    
    #legend
    theme(legend.position = "right", 
          legend.text = element_text(size = 6)) +
    
    scale_color_discrete(name = "P-values", labels = c("Non-significant", "Significant")) +
    scale_shape_discrete(name = "GLM Analysis", labels = c("Individual study", "All studies")) +
    
    guides(colour = guide_legend(override.aes = list(shape = 16)), 
           shape = guide_legend(override.aes = list(shape = c(16, 8)))) +
    
    #shape and colour points
    scale_shape_manual(name = "GLM Analysis", 
                       labels = c("Individual study", "All studies"),
                       values = c(8, 16)) +
    scale_colour_manual(labels = c("Non-significant", "Significant"), 
                        values = c("#440154FF", "#35B779FF")) 
  
  
  #gene distance plot
  GDplot <- ggplot(gd.df,
                   aes(x=group, y=value, 
                       colour = pvals,
                       shape = factor(ifelse(variable != "31", "Individual Study", "All Studies")))) +
    
    #data
    geom_violin(data = subset(gd.df, gd.df$variable != '31'), colour = "black")  +
    geom_jitter(data = subset(gd.df, gd.df$variable == '31'), size = 4, colour = "#35B779FF") +
    geom_jitter(data = subset(gd.df, gd.df$variable != '31'), alpha = 0.5, size = 4) +
    geom_hline(yintercept = 0, linetype = "dashed", show.legend = TRUE) +
    
    #axis co-ordinates
    coord_cartesian(ylim = c(-100, 25)) +
    
    #shape and colour points
    scale_shape_manual(name = "GLM Analysis", 
                       labels = c("Individual study", "All studies"),
                       values = c(8, 16)) +
    scale_colour_manual(labels = c("Non-significant", "Significant"), 
                        values = c("#440154FF", "#35B779FF")) +
    
    #theme and plot labels
    theme_bw() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + #Wrap variable labels
    theme(plot.title = element_text(size = 10, face = "bold")) +
    labs(title = "B. Distance to Gene Trees", x="Topological Measurement", y="T Statistic") +
    
    #legend
    theme(legend.position = "none") 
  
  #Extract legend for PDF
  mylegend<-generateLegend(SDplot)
  
  #save to PDF
  aspectR <- 2.5
  pdf(file = filename1, onefile = TRUE, height = 7, width = 7*aspectR, useDingbats = FALSE)
  lineplot <- grid.arrange(arrangeGrob(SDplot + theme(legend.position="none"),
                                       GDplot + theme(legend.position="none"),
                                       nrow=1),
                           mylegend, widths=c(70, 22))
  dev.off()
} 

#Calculates psuedo R-squared value for regression models
RsquaredGLM <- function(indvglm, bigglm) {
  
  require(lme4)
  require(MuMIn)
  
  rbig <- lapply(bigglm , MuMIn::r.squaredGLMM)
  rindv <- lapply(indvglm, MuMIn::r.squaredGLMM)
  
  return(c(rindv, rbig))
}

#Plot R-squared values
plotR2 <- function(rv, filename) {
  
  require(viridis)
  
  #extract R conditional from R values
  r2val <-lapply(rv, getRval)
  
  #Extract species distance to data frame
  r2 <- data.frame(r2val[31:60])
  for (i in 1:30) {
    colnames(r2)[i] <- paste0(i)
  }
  
  #Extract gene distance to data frame
  r2gd <- data.frame(r2val[1:30])
  for (i in 1:30) {
    colnames(r2gd)[i] <- paste0(i)
  }
  
  #merge data frames
  r2all <- rbind.data.frame(r2, r2gd, make.row.names = FALSE)
  rownames(r2all) <- c("1Species.Distance", "2Gene.Distance")
  
  #melt data frame
  r2all$group <- row.names(r2all) #Put variables in group for melt transformation
  r2.m <- melt(r2all, id.vars = "group")
  
  #plot R2
  rplot <- ggplot(r2.m, aes(y = value, x = group, colour = group)) +
    geom_violin(colour = "black") +
    geom_text(aes(label=variable), position=position_jitter(0.4)) +
    theme_bw() +
    labs(title = "R-squared Conditional Value for GLMs", x="Models", y="R2 Conditional") + 
    scale_colour_manual(values = c("#440154FF", "#35B779FF")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  #save plot to file
  pdf(file = filename, onefile = TRUE, height = 7, width = 7, useDingbats = FALSE)
  plot(rplot)
  dev.off()
}

#Function to sort gene-trees according to branch metric, create loci subsets and estimate new species trees
getBLVtrees <- function(studydf, blvariable) {
  
  require(dplyr)
  
  #Isolate study from full dataframe
  study.number <- studydf$Study.Number[1]
  
  #Arrange by COV RTTL and branch support
  ascending.df <- studydf %>% arrange(blvariable)
  descending.df <- studydf %>% arrange(desc(blvariable))
  
  #Create 20%, 40%, 60%, 80% loci lists
  totaloci <- nrow(studydf) 
  loci <- c(totaloci*0.2, totaloci*0.4, totaloci*0.6, totaloci*0.8)
  
  #Get a random assortment of loci for control based on number of loci, , gives index position
  randomloci <- c(list(sample(studydf$X.1, loci[1])), list(sample(studydf$X.1, loci[2])), 
                  list(sample(studydf$X.1, loci[3])), list(sample(studydf$X.1, loci[4])))
  
  #Extract loci from dataframe for ascending subset
  asc20 <- ascending.df[1:loci[1],]
  asc40 <- ascending.df[1:loci[2],]
  asc60 <- ascending.df[1:loci[3],]
  asc80 <- ascending.df[1:loci[4],]
  ascloci <- c(list(asc20$X.1), list(asc40$X.1), list(asc60$X.1), list(asc80$X.1))
  
  #Extract loci from dataframe for descending subset
  desc20 <- descending.df[1:loci[1],]
  desc40 <- descending.df[1:loci[2],]
  desc60 <- descending.df[1:loci[3],]
  desc80 <- descending.df[1:loci[4],]
  descloci <- c(list(desc20$X.1), list(desc40$X.1), list(desc60$X.1), list(desc80$X.1))
  
  #Concatenate loci lists 
  allloci <- c(ascloci, descloci, randomloci)
  
  #Make new list of trees with loci list
  studytrs <- list()
  for (i in 1:length(allloci)) {
    trsi <- list(maxtaxgenetrees[[study.number]][allloci[[i]]])
    studytrs[(length(studytrs)+1)] <- trsi
  }
  for (i in 1:length(studytrs)) {class(studytrs[[i]]) <- "multiPhylo"}
  
  #write trees to file as multiPhylo
  for (i in 1:12) {
    ape::write.tree(studytrs[[i]], file = paste0(treDir, ref, "_", fileprefix[i], ".trs"))
  }
}

#Extracts measurable outcomes of subset species trees
extractFilteringVariables <- function(species_tr, ref_tr, 
                                      bl_variable, data_percentage, loci_rank, 
                                      study_name, study_num) {
  
  require(ape)
  require(phangorn)
  
  
  v <- bl_variable
  dpercent <- data_percentage
  lrank <- loci_rank
  snum <- study_num
  sname <- study_name
  
  RFdistance <- RF.dist(species_tr, ref_tr, check.labels = TRUE, normalize = TRUE)
  class(species_tr$node.label) <- "numeric"
  mean_branchSup <- mean(species_tr$node.label, na.rm = TRUE)
  
  return(c(snum, v, dpercent, lrank, sname, mean_branchSup, RFdistance))
}

#Creates data frame of measurable outcomes
createAccuracyMetricsDF <- function(rttlSTrs, mbsSTrs, sdSTrs, referenceTrees) {
  
  #extract co-efficient of variation in root-to-tip distances (RTTL) and create df
  r1 <- extractFilteringAccuracyMetrics(species_trees = rttlSTrs$top, reference_species = referenceTrees, 
                                        bl_vari = "rttl", tbr = "top")
  
  r2 <- extractFilteringAccuracyMetrics(species_trees = rttlSTrs$bottom, reference_species = referenceTrees, 
                                        bl_vari = "rttl", tbr = "bottom")
  
  r3 <- extractFilteringAccuracyMetrics(species_trees = rttlSTrs$random, reference_species = referenceTrees, 
                                        bl_vari = "rttl", tbr = "random")
  Rdf <- c(r1, r2, r3)
  Rdf <- as.data.frame(do.call(rbind, lapply(Rdf, `length<-`, max(lengths(Rdf)))))
  colnames(Rdf) <- c("Study.Number", "Variable", "Percentage.of.Loci", "Loci.Rank", "Reference", "Mean.Branch.Support", "Distance.to.Reference")

  #extract mean branch support (MBS) and create df
  m1 <- extractFilteringAccuracyMetrics(species_trees = mbsSTrs$top, reference_species = referenceTrees, 
                                        bl_vari = "mbs", tbr = "top")
  
  m2 <- extractFilteringAccuracyMetrics(species_trees = mbsSTrs$bottom, reference_species = referenceTrees, 
                                        bl_vari = "mbs", tbr = "bottom")
  
  m3 <- extractFilteringAccuracyMetrics(species_trees = mbsSTrs$random, reference_species = referenceTrees, 
                                        bl_vari = "mbs", tbr = "random")
  Mdf <- c(m1, m2, m3)
  Mdf <- as.data.frame(do.call(rbind, lapply(Mdf, `length<-`, max(lengths(Mdf)))))
  colnames(Mdf) <- c("Study.Number", "Variable", "Percentage.of.Loci", "Loci.Rank", "Reference", "Mean.Branch.Support", "Distance.to.Reference")
  
  #extract control species distance (SD) metrics (Robinson-Folds distance) and create df
  s1 <- extractFilteringAccuracyMetrics(species_trees = sdSTrs$top, reference_species = referenceTrees, 
                                        bl_vari = "sd", tbr = "top")
  
  s2 <- extractFilteringAccuracyMetrics(species_trees = sdSTrs$bottom, reference_species = referenceTrees, 
                                        bl_vari = "sd", tbr = "bottom")
  
  s3 <- extractFilteringAccuracyMetrics(species_trees = sdSTrs$random, reference_species = referenceTrees, 
                                        bl_vari = "sd", tbr = "random")
  Sdf <- c(s1, s2, s3)
  Sdf <- as.data.frame(do.call(rbind, lapply(Sdf, `length<-`, max(lengths(Sdf)))))
  colnames(Sdf) <- c("Study.Number", "Variable", "Percentage.of.Loci", "Loci.Rank", "Reference", "Mean.Branch.Support", "Distance.to.Reference")
  
  return(c(Rdf, Mdf, Sdf))
}

#Extracts measurable outcomes for reference tree
getReferenceTreesMetrics <- function(reference.Tree) {
  
  require(ape)
  require(phangorn)
  
  #Extract metrics
  metrics <- list()
  snumber <- 1:30
  for (i in 1:30) {
    si <- getRefAccuracy(species_tr = reference.Tree[[i]], ref_tr = reference.Tree[[i]], 
                         study_name = references[[i]], study_num = snumber[i])
    metrics[length(metrics)+1] <- list(si)
  }
  
  #Create df and save 
  df <- as.data.frame(do.call(rbind, lapply(metrics, `length<-`, max(lengths(metrics)))))
  colnames(df) <- c("Study.Number", "Variable", "Percentage.of.Loci", "Loci.Rank", "Reference", "Mean.Branch.Support", "Distance.to.Reference")
  write.csv(df, file = paste0("filepath", "filename.csv"))
}

#Generates common legend for multiple plots
generateLegend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#Makes plot with filtering outcomes
filteringViolinplots <- function (rttlDF, mbsDF, sdDF, legend) {
  
  require(ggplot2)
  require(ggpubr)
  require(gridExtra)
  
  #Plot clocklikeness
  r1 <- ggplot(rttlDF, aes(x = Percentage.of.Loci, y = Mean.Branch.Support, colour = mixedsort(sort(Loci.Rank)))) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw() +
    coord_cartesian(ylim = c(0.6, 1)) +
    scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#F9B641FF")) +
    labs(y = "Mean Branch Support (Astral)") +
    theme(axis.title.x=element_blank()) 
  
  r2 <- ggplot(rttlDF, aes(x = Percentage.of.Loci, y = Distance.to.Reference, colour = Loci.Rank)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 0.8)) +
    scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#F9B641FF")) +
    labs(title = "Clocklikeness", y = "Distance to Reference Tree (Robinson-Foulds)") +
    theme(axis.title.x=element_blank()) 
  
  #plot mbs
  m1 <- ggplot(mbsDF, aes(x = Percentage.of.Loci, y = Mean.Branch.Support, colour = Loci.Rank)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw() +
    coord_cartesian(ylim = c(0.6, 1)) +
    scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#F9B641FF")) +
    theme(axis.title.x=element_blank(), 
          axis.title.y = element_blank()) 
  
  m2 <- ggplot(mbsDF, aes(x = Percentage.of.Loci, y = Distance.to.Reference, colour = Loci.Rank)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 0.8)) +
    labs(title = "SH-aLRT Mean Branch Support") +
    scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#F9B641FF")) +
    theme(axis.title.x=element_blank(), 
          axis.title.y = element_blank()) 
  
  #Plot species distance
  s1 <- ggplot(sdDF, aes(x = Percentage.of.Loci, y = Mean.Branch.Support, colour = Loci.Rank)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw() +
    coord_cartesian(ylim = c(0.6, 1)) +
    scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#F9B641FF")) +
    theme(axis.title.x=element_blank(), 
          axis.title.y = element_blank()) 
  
  s2 <- ggplot(sdDF, aes(x = Percentage.of.Loci, y = Distance.to.Reference, colour = Loci.Rank)) +
    geom_boxplot(outlier.shape = NA, show.legend = F) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 0.8)) +
    labs(title = "RF Distance to Species Tree") +
    scale_color_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#F9B641FF")) +
    theme(axis.title.x=element_blank(), 
          axis.title.y = element_blank()) 
  
  #Create group plot
  plot1 <- grid.arrange(arrangeGrob(r2, m2, s2,
                                    r1, m1, s1,
                                    ncol = 3, nrow = 2, 
                                    bottom = "Percentage of Total Loci", 
                                    top = "Loci Rank", 
                                    left = "Phylogenetic Inferece Accuracy Metrics"), 
                        "bottom" = box_legend)
  
  pdf(file = paste0("filepath", "filename.pdf"), onefile = TRUE)
  plot(plot1)
  dev.off()
  return(plot1)
}

#Function to test each subset using ANOVA and TukeyHSD
anovaCommand <- function(df) {
  
  require(car)
  require(dplyr)
  
  #Oneway ANOVA: Distance to Reference
  a20DTR <- aov(Distance.to.Reference ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "20" | Percentage.of.Loci == "100"))
  A20DTR <- Anova(a20DTR, type = "III")
  
  a40DTR <- aov(Distance.to.Reference ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "40" | Percentage.of.Loci == "100"))
  A40DTR <- Anova(a40DTR, type = "III")
  
  a60DTR <- aov(Distance.to.Reference ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "60" | Percentage.of.Loci == "100"))
  A60DTR <- Anova(a60DTR, type = "III")
  
  a80DTR <- aov(Distance.to.Reference ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "80" | Percentage.of.Loci == "100"))
  A80DTR <- Anova(a80DTR, type = "III")
  
  #oneway ANOVA: Species Tree Mean Branch Support
  a20MBS <- aov(Mean.Branch.Support ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "20" | Percentage.of.Loci == "100"))
  A20MBS <- Anova(a20MBS, type = "III")
  
  a40MBS <- aov(Mean.Branch.Support ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "40" | Percentage.of.Loci == "100"))
  A40MBS <- Anova(a40MBS, type = "III")
  
  a60MBS <- aov(Mean.Branch.Support ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "60" | Percentage.of.Loci == "100"))
  A60MBS <- Anova(a60MBS, type = "III")
  
  a80MBS <- aov(Mean.Branch.Support ~ Loci.Rank, data = dplyr::filter(df, Percentage.of.Loci == "80" | Percentage.of.Loci == "100"))
  A80MBS <- Anova(a80MBS, type = "III")
  
  #TukeyHSD: Distance to Reference
  t20DTR <- TukeyHSD(a20DTR)
  t40DTR <- TukeyHSD(a40DTR)
  t60DTR <- TukeyHSD(a60DTR)
  t80DTR <- TukeyHSD(a80DTR)
  
  #TukeyHSD: Species Tree Mean Branch Support
  t20MBS <- TukeyHSD(a20MBS)
  t40MBS <- TukeyHSD(a40MBS)
  t60MBS <- TukeyHSD(a60MBS)
  t80MBS <- TukeyHSD(a80MBS)
  
  return(c(list(A20DTR), list(A40DTR), list(A60DTR), list(A80DTR),
           list(A20MBS), list(A40MBS), list(A60MBS), list(A80MBS),
           list(t20DTR), list(t40DTR), list(t60DTR), list(t80DTR),
           list(t20MBS), list(t40MBS), list(t60MBS), list(t80MBS)))
}
