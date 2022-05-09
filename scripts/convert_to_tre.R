#convert to *.tre
library(ape)


paths <- Sys.glob("trees/*.tree")
newpaths <- gsub(".tree", ".tre", paths)

for (i in seq_along(paths)){
    print(paste0("Renaming: ", paths[i]))
    phy <- read.tree(paths[i])
    write.tree(phy, newpaths[i])
}

paths <- Sys.glob("trees/*.nex")
newpaths <- gsub(".nex", ".tre", paths)

for (i in seq_along(paths)){
    print(paste0("Renaming: ", paths[i]))
    phy <- read.nexus(paths[i])
    write.tree(phy, newpaths[i])
}


print("Testing if all trees can be loaded")

trees <- Sys.glob("trees/*.tre")
for (tree in trees){
    out <- tryCatch(
        {
        phy <- read.tree(tree)
        cat(".")
        },
        error = function(cond) {
        message(paste0("Tree can't load:", tree))
        message("Here is original error:")
        message(cond)
        }
    )
}
cat("\n")
print("done")
