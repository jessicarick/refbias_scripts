library(TreeDist)
library(TreeTools)

spectrum <- hcl.colors(6,'plasma')
results.Q <- results.raxml %>%
  mutate(distQ = numeric(1898),
         distCI = numeric(1898))

for (height in c("SHORT","MED","LONG")) {
  for (i in 1:10){
    trees <- c(ml.tree[ml.tree.info$height == height & ml.tree.info$simulation == i],raxml.trees[results.raxml$height == height & results.raxml$simulation == i])
    distancesQ <- Quartet::QuartetDivergence(Quartet::QuartetStatus(trees), similarity=FALSE)[-1]
    distancesCI <- as.matrix(ClusteringInfoDistance(trees,normalize = TRUE))[-1,1]
    results.Q$distQ[results.Q$height == height & results.Q$simulation == i] <- distancesQ/(2/3)
    results.Q$distCI[results.Q$height == height & results.Q$simulation == i] <- distancesCI
    
    # info <- as.factor(c("ML",results.raxml$maf[results.raxml$height == "SHORT" & results.raxml$simulation == i]))
    # mapping <- cmdscale(distances3, k=5)
    # plot(mapping, asp=1,ann=FALSE,col=scales::alpha(spectrum[info],0.6),pch=16,cex=4)
    # legend('bottomleft',levels(info),
    #        col=spectrum,pch=16)
  }
}


lm(distQ ~ maf + missing + int + maf:int + missing:int, data=results.Q, pch=16, col=as.factor(results.Q$int))
results.Q %>%
  filter() %>%
  ggplot(aes(x=maf,y=distCI,col=as.factor(int))) +
  geom_point() +
  theme_custom() +
  geom_smooth(method="lm") +
  facet_wrap(~height) 



#----------------#
distances1 <- as.matrix(ClusteringInfoDistance(trees))[,1]
distances2 <- RobinsonFoulds(trees)

plot(distances3,distances1)

distancesQ <- Quartet::QuartetDivergence(Quartet::QuartetStatus(trees), similarity=FALSE)
plot(results.raxml$maf[results.raxml$height == "SHORT" & results.raxml$simulation == i],distancesQ,pch=16,col=as.factor(results.raxml$int[results.raxml$height == "SHORT" & results.raxml$simulation == i]))
info <- results.raxml[results.raxml$height == "SHORT" & results.raxml$simulation == i,]
info$Qdist <- distancesQ
summary(lm(Qdist ~ maf + int + maf:int, data=info))
abline(lm(Qdist ~ maf, data=info[info$int == "EXT",]),col="black",lty=2)
abline(lm(Qdist ~ maf, data=info[info$int == "INT",]),col="red",lty=2)

possibleClusters <- 2:6
pamClusters <- lapply(possibleClusters, function(k) cluster::pam(distances1,k=k))
pamSils <- vapply(pamClusters, function(pamcluster) {
  mean(cluster::silhouette(pamcluster)[,3])
}, double(1))
bestPam <- which.max(pamSils)
pamSil <- pamSils[bestPam]
pamCluster <- pamClusters[[bestPam]]$cluster 

htree <- protoclust::protoclust(distances1)
hClusters <- lapply(possibleClusters, function(k) cutree(htree, k=k))
hSils <- vapply(hClusters, function(hcluster) {
  mean(cluster::silhouette(hcluster, distances1)[,3])
},double(1))

bestH <- which.max(hSils)
hSil <- hSils[bestH]
hCluster <- hClusters[[bestH]]

plot(pamSils ~ possibleClusters,
     xlab = 'Number of clusters', ylab = 'Silhouette coefficient',
     ylim = range(c(pamSils, hSils)))
points(hSils ~ possibleClusters, pch = 2)
legend('topright', c('PAM', 'Hierarchical'), pch = 1:2)


cluster <- hClusters[[3 - 1]]
class(htree) <- 'hclust'
plot(htree, labels = FALSE, main = '')
points(seq_along(trees), rep(21.5, length(trees)), pch = 16,
       col = pamCluster)
