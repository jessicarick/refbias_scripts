library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")
library("phytools")

output <- "053019-output"
results.raxml <- read.csv(paste("output/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

results.raxml$simulation <- as.factor(results.raxml$simulation)

subset <- raxml.trees[results.raxml$simulation == "2" &
                        results.raxml$noref == "REF" &
                        results.raxml$height == "2000000"]
j <- which(ml.tree.info$simulation == 2 & ml.tree.info$height == 2000000)
subset2 <- as.multiPhylo(c(subset,ml.tree[[j]]))
rooted.subset <- root(subset2, "sim_0_0_0", resolve.root=TRUE)

res <- treespace(rooted.subset, nf=10, method="RF")
names(res)
table.image(res$D, nclass=30)
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(6))

plotGroves(res$pco, lab.show=FALSE, lab.cex=1.5)
plotGrovesD3(res$pco, treeNames=paste(c(results.raxml$maf[results.raxml$simulation == "2" &
                                                    results.raxml$noref == "REF" &
                                                    results.raxml$height == "2000000"],"true"),
                                      c(results.raxml$missing[results.raxml$simulation == "2" &
                                                    results.raxml$noref == "REF" &
                                                    results.raxml$height == "2000000"],"true"),
                                      sep = ","),
             symbol_var = c(results.raxml$int[results.raxml$simulation == "2" &
                                            results.raxml$noref == "REF" &
                                            results.raxml$height == "2000000"],"true"))

clusters <- findGroves(res, nclust=6)
names(clusters)
vircols <- c("#440154", "#414487", "#2A788E", "#22A884", "#7AD151", "#FDE725")
plotGrovesD3(res$pco, treeNames=paste(c(results.raxml$maf[results.raxml$simulation == "2" &
                                                                          results.raxml$noref == "REF" &
                                                                          results.raxml$height == "2000000"],"true"),
                                                    c(results.raxml$missing[results.raxml$simulation == "2" &
                                                                              results.raxml$noref == "REF" &
                                                                              results.raxml$height == "2000000"],"true"),
                                                    sep = ","),
             col_lab="cluster", legend_width=150, groups=clusters$groups, colors=vircols,
             point_opacity=0.4,
             point_size=100,
             ellipses=TRUE, fixed=TRUE, menu=TRUE)

colours <- fac2col(clusters$groups, col.pal=virid)
plot3d(clusters$treespace$pco$li[,1],
      clusters$treespace$pco$li[,2],
      clusters$treespace$pco$li[,3],
      col=colours, type="s", size=1.5,
      xlab="", ylab="", zlab="")

medTrees <- medTree(rooted.subset, groups=clusters$groups)
medTrees.clust <- lapply(medTrees, function(e) ladderize(e$trees[[1]]))
par(mfrow=c(2,2))
for (i in 1:length(medTrees.clust)){
  plot(medTrees.clust[[i]],
       main=paste("cluster",i),
       cex=1.5)
}

# Compare median trees from clusters to true tree, as a visual:
for (i in 1:length(medTrees.clust)){
  plotTreeDiff(medTrees.clust[[i]],ml.tree[[j]], use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)
}
