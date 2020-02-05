This repository is associated with the article 'Phylogenetic signal is associated with the degree of variation in root-to-tip distances.'

It contains the R code and tree files used in the analysis:

  1. MaximalTrs: code used to estimate gene trees with maximal number of taxa.
  
  2. updated.all.trees.Rdata: contains 36,075 gene trees from the 34 studies, estimated using above code, ordered by study from smallest to largest number of loci.
  
  3. speciesTr_estimationRTT.R: code that estimated the species tree for each study using Astral (version 5.7.3) and updated.all.trees.Rdata.
  
  3. data_extractionRTT: code to extract branch length metrics, gene distance and species distance (RF.dist).
  
  4. data_analysisRTT: code used to run glm analysis and plot t-statistics.



