# From a set of trees extract a set that maximises the number of loci and taxa:

maximalTaxonTrs <- function(trs){
	
	allntax <- table(unlist(lapply(trs, function(z) z$tip.label)))
	taxnames <- names(allntax)
      	allntax <- as.numeric(allntax)
      	names(allntax) <- taxnames
      	allntax <- sort(allntax, decreasing = T)
	
	# CHECK THIS FOLLOWING SECTION!!
	taxrankproduct <- 1 * allntax[1]
	for(i in 2:length(allntax)){
	      NtrsWithTargets <- sum(sapply(trs, function(x) all(names(allntax)[1:i] %in% x$tip.label)))
	      taxrankproduct <- c(taxrankproduct, i * NtrsWithTargets)
	}
	names(taxrankproduct) <- names(allntax)
	
	maximalTaxa <- names(allntax)[1:tail(which(taxrankproduct == max(taxrankproduct)), 1)]
	maximalTrs <- list()
	for(i in 1:length(trs)){
	      if(all(maximalTaxa %in% trs[[i]]$tip.label) && length(maximalTaxa) == length(trs[[i]]$tip.label)){
	      	     maximalTrs[[length(maximalTrs) + 1]] <- trs[[i]]
		     names(maximalTrs)[length(maximalTrs)] <- names(trs)[i]
	      } else if(all(maximalTaxa %in% trs[[i]]$tip.label)){
	      	     maximalTrs[[length(maximalTrs) + 1]] <- drop.tip(trs[[i]], trs[[i]]$tip.label[which(!trs[[i]]$tip.label %in% maximalTaxa)])
		     names(maximalTrs)[length(maximalTrs)] <- names(trs)[i]
	      }
	}
	return(maximalTrs)
}