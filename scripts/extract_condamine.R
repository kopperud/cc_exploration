library(ape)

paths <- Sys.glob("condamine/*.Rdata")

print("extracting condamine .Rdata, calculating sampling fraction and saving to file")

for (path in paths){
  load(path)
}

trees <- c(FamilyAmphibiaTrees, FamilyBirdTrees, FamilyCrocoTurtleTrees, 
           FamilyMammalTrees, FamilySquamateTrees)

treenames <- names(trees)

for (treename in treenames){
  tree <- trees[[treename]]
  n_extant <- tree$totalsp
  phy <- tree$tree
  
  write.tree(phy, file = paste0("trees/condamine_etal_2019_", treename, ".tre"))
  rho <- length(phy$tip.label) / n_extant
  cat("rho <- ", rho, file = paste0("trees/condamine_etal_2019_", treename, "_sampling.Rev"))
  cat(".")
}
cat("\n")
print("done")

#hist(sapply(trees, function(e) e$totalsp))
#hist(sapply(trees, function(e) Ntip(e$tree)))


