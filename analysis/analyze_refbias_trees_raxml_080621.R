############################################
## Script for analyzing refbias output trees from RAxML
## Written by J. Rick, 11 March 2019
## Updated 6 August 2021
## Made to be run on the command line
############################################

suppressMessages(
  c(library(ape),
    library(phytools),
    library(apTreeshape),
    library(phangorn),
    library(tidyverse),
    # library(readtext),
    library(ggpubr),
    # library(CCA)
    # library(GGally)
    # library(vegan)
    # library(LMERConvenienceFunctions)
    # library(lme4)
    # library(MuMIn)
    # library(car)
    library(ggrepel),
    library(ggsci),
    library(TreeDist)))
source("analysis/theme_custom.R")

############################################
## Specifying file names ########
############################################

raxml.trees <- "092321-all-raxml.trees"
raxml.tree.names <- "092321-all-raxml.names"
# astral.trees <- "101819-all-astral.trees"
# astral.tree.names <- "101819-all-astral.names"
ml.trees <- "092321-s_tree.tree"
ml.tree.names <- "092321-s_tree.names"
output <- "092321-output"
refdist <- "092321-refdist.txt"
args <- mget(c("raxml.trees","raxml.tree.names","ml.trees","ml.tree.names","refdist","output"))

############################################
## Reading in command line arguments #######
############################################

#check for the argparse package installed
# if (!('argparse' %in% installed.packages()[, 'Package'])) {
#   install.packages('argparse', repos='http://cran.rstudio.com/')
# }
# 
# suppressPackageStartupMessages(library('argparse'))
# 
# parse.args <- function() {
#   parser <- ArgumentParser()
#   parser$add_argument('--ml.trees', help='file with ML (true) trees')
#   parser$add_argument('--ml.tree.names', help='file with names of ML trees')
#   parser$add_argument('--raxml.trees', help='file with raxml trees')
#   parser$add_argument('--raxml.tree.names', help='file with raxml tree names')
#   parser$add_argument('--astral.trees', help='file with astral trees')
#   parser$add_argument('--astral.tree.names', help='file with astral tree names')
#   parser$add_argument('--refdist', help='file with reference distances')
#   parser$add_argument('-o', '--output', help='output file prefix')
#   return(parser$parse_args())
# }
# args <- parse.args()

####Read in concatenated tree file and ML tree file (make sure the latter is ultrametric)
raxml.trees<-read.tree(paste("output/new/",args$raxml.trees,sep=""))
#astral.trees<-read.tree(paste("output/",args$astral.trees,sep=""))
ml.tree<-read.tree(paste("output/new/",args$ml.trees,sep=""))

###Read in file of tree names
raxml.tree.names<-read.table(paste("output/new/",args$raxml.tree.names,sep=""),stringsAsFactors=FALSE)[,1]
#astral.tree.names<-read.table(paste("output/",args$astral.tree.names,sep=""),stringsAsFactors = FALSE)
ml.tree.names <- read.table(paste("output/new/",args$ml.tree.names,sep=""),stringsAsFactors = FALSE)[,1]

###################
## changing tip labels on ml trees to match the sim trees
## ONLY DO ONCE!!
####################
for (i in 1:length(ml.tree)){
  tips <- ml.tree[[i]]$tip.label
  tips.new <- paste("sim_",tips,"_0_0",sep="")
  ml.tree[[i]]$tip.label <- tips.new
}

###Pull info about ml species trees
## and calculate tree statistics
ml.tree.info <- tibble(simulation = numeric(length(ml.tree.names)),
                           height = character(length(ml.tree.names))) %>%
  mutate(simulation = sapply(ml.tree.names, function(x) as.integer(regmatches(x, regexec('sim([0-9]+)', x))[[1]][2])),
         height = sapply(ml.tree.names, function(x) regmatches(x, regexec('height([A-Z]+)', x))[[1]][2]),
         tree.height = sapply(ml.tree,function(x) max(branching.times(x),na.rm=T)), 
         ingroup.tree.height = sapply(ml.tree,function(x) max(branching.times(drop.tip(x,"sim_0_0_0")))),
         Avg.BLs = sapply(ml.tree,function(x) mean(x$edge.length,na.rm=T)),
         gamma = sapply(ml.tree,function(x) gammaStat(x)),
         ingroup.gamma = sapply(ml.tree,function(x) gammaStat(drop.tip(x,"sim_0_0_0"))),
         colless = sapply(ml.tree,function(x) colless(as.treeshape(x,model="pda"),norm="pda")),
         ingroup.colless = sapply(ml.tree,function(x) colless(as.treeshape(drop.tip(x,"sim_0_0_0"),model="pda"),norm="pda"))
  )

refdist <- read_table2(paste("output/new/",args$refdist,sep=""), col_names=c("simulation","height","int","taxa_ref","avg_dxy"), skip=1) %>%
  group_by(simulation,height,int,taxa_ref) %>%
  mutate(simulation = as.integer(simulation)) %>%
  slice_tail(n=1) %>%
  ungroup()

gt.rf <- read_table2("output/new/092321-LONG-EXT.TTR.gt_rf",col_names=c("sim","gt_rf")) %>%
  mutate(height = "LONG",
         simulation = as.integer(gsub("sim","",sim))) %>%
  add_row(read_table2("output/new/092321-MED-EXT.TTR.gt_rf",col_names=c("sim","gt_rf")) %>%
            mutate(height = "MED"),
          simulation = as.integer(gsub("sim","",sim))) %>%
  add_row(read_table2("output/new/092321-SHORT-EXT.TTR.gt_rf",col_names=c("sim","gt_rf")) %>%
            mutate(height = "SHORT"),
          simulation = as.integer(gsub("sim","",sim))) %>%
  drop_na() %>%
  group_by(height,sim) %>%
  slice_head(n=1)

ml.tree.info2 <- ml.tree.info %>%
  left_join(gt.rf,by=c("simulation","height")) %>%
  filter(simulation > 15 & simulation < 26)


#################
## plots of ml tree characteristics
#################
plot1 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=tree.height,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Tree Height")

plot2 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=Avg.BLs,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Average Branch Lengths")

plot3 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=gt_rf,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Gene Tree Discordance")

plot4 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=ingroup.tree.height,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Ingroup Tree Height")

plot5 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=ingroup.colless,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Ingroup Colless Imbalance")

plot6 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=ingroup.gamma,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Ingroup gamma")

ggarrange(plot1,plot2,plot4,plot3,plot5,plot6,nrow=2,ncol=3)

#####################
## Start of analysis!
#####################

###Create empty data frame and with named, empty columns 
num.trees <- length(raxml.tree.names)
results.raxml <- tibble(tree.num=seq(1:length(raxml.trees)))

## pull info from tree names, and calculate tree statistics
results.raxml <- results.raxml %>%
  # mutate(
  #   simulation = sapply(raxml.tree.names,function(x) as.integer(regmatches(x, regexec('-s([0-9]+)\\_q', x))[[1]][2])),
  #   height = sapply(raxml.tree.names,function(x) as.character(regmatches(x, regexec('([A-Z]+)-s[0-9]', x))[[1]][2])),
  #   method = "raxml",
  #   quality = sapply(raxml.tree.names,function(x) as.integer(regmatches(x, regexec('q(.*?)\\_m', x))[[1]][2])),
  #   missing = sapply(raxml.tree.names,function(x) as.numeric(regmatches(x, regexec('miss([0-9]?\\.?[0-9]+)\\_mac', x))[[1]][2])),
  #   maf = sapply(raxml.tree.names, function(x) as.numeric(regmatches(x, regexec('_mac(0?.?[0-9]+)_sites', x))[[1]][2])),
  #   sites = sapply(raxml.tree.names, function(x) as.integer(regmatches(x, regexec('sites([0-9]*?)\\.', x))[[1]][2])),
  #   taxa_ref = sapply(raxml.tree.names, function(x) regmatches(x,regexec('[A-Z]-([0-9]+_0_0)\\.phylip', x))[[1]][[2]]),
  #   int = sapply(raxml.tree.names, function(x) regmatches(x, regexec('REF.([A-Z]+)\\.', x))[[1]][2]),
  #   noref = "REF"
  # ) %>%
  left_join(refdist,by=c("simulation","height","int","taxa_ref"),copy=TRUE) %>%
  left_join(select(gt.rf,gt_rf:simulation),by=c("simulation","height")) %>%
  mutate(
    tree.height = sapply(raxml.trees, function(x) max(branching.times(root(x,"sim_0_0_0",resolve.root=TRUE)),na.rm=T)),
    ingroup.tree.height = sapply(raxml.trees, function(x) max(branching.times(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0")),na.rm=T)),
    Avg.BLs = sapply(raxml.trees, function(x) mean(root(x,"sim_0_0_0",resolve.root=TRUE)$edge.length,na.rm=T)),
    SD.BLs = sapply(raxml.trees, function(x) sd(root(x,"sim_0_0_0",resolve.root=TRUE)$edge.length,na.rm=T)),
    mean.support = sapply(raxml.trees, function(x) mean(as.numeric(root(x,"sim_0_0_0",resolve.root=TRUE)$node.label[-1]),na.rm=T)),
    sd.support = sapply(raxml.trees, function(x) sd(as.numeric(root(x,"sim_0_0_0",resolve.root=TRUE)$node.label[-1]),na.rm=T))
  ) %>%
  mutate(
    gamma = sapply(raxml.trees, function(x) gammaStat(chronopl(root(x,"sim_0_0_0",resolve.root=TRUE),lambda=1))),
    ingroup.gamma = sapply(raxml.trees, function(x) gammaStat(drop.tip(chronopl(root(x,"sim_0_0_0",resolve.root=TRUE),lambda=1),"sim_0_0_0"))),
    colless = sapply(raxml.trees, function(x) colless(as.treeshape(root(x,"sim_0_0_0",resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.colless = sapply(raxml.trees, function(x) colless(as.treeshape(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0"),model="pda"),norm="pda")),
    sackin = sapply(raxml.trees, function(x) sackin(as.treeshape(root(x,"sim_0_0_0",resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.sackin = sapply(raxml.trees, function(x) sackin(as.treeshape(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0"),model="pda"),norm="pda"))
  ) 

###################
## statistics requiring separation by simulation number
## to compare to ML tree
###################
results.raxml <- results.raxml %>%
  mutate(RF.Dist.ML = numeric(nrow(results.raxml)),
         Q.Dist.ML = numeric(nrow(results.raxml)),
         CI.Dist.ML = numeric(nrow(results.raxml)),
         std.gamma = numeric(nrow(results.raxml)),
         std.colless = numeric(nrow(results.raxml)),
         std.sackin = numeric(nrow(results.raxml)),
         std.sites = numeric(nrow(results.raxml)),
         std.ingroup.gamma = numeric(nrow(results.raxml)),
         std.ingroup.colless = numeric(nrow(results.raxml)),
         std.ingroup.sackin = numeric(nrow(results.raxml)))
for (i in 1:nrow(results.raxml)){
  s <- as.integer(results.raxml$simulation[i])
  h <- as.character(results.raxml$height[i])
  j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
  
  ml.tree.pruned <- drop.tip(ml.tree[[j]],results.raxml$taxa_ref[i])
  
  if(results.raxml$taxa_ref[i] %in% raxml.trees[[i]]$tip.label){
    trees.both <- c(raxml.trees[[i]],ml.tree[[j]])
  } else {
    trees.both <- c(raxml.trees[[i]],ml.tree.pruned)
  }
  
  if(results.raxml$taxa_ref[i] %in% raxml.trees[[i]]$tip.label){
    trees.both.root <- c(root(raxml.trees[[i]],"sim_0_0_0",resolve.root=TRUE),
                         root(ml.tree[[j]],"sim_0_0_0",resolve.root=TRUE))
  } else {
    trees.both.root <- c(root(raxml.trees[[i]],"sim_0_0_0",resolve.root=TRUE),
                         root(ml.tree.pruned,"sim_0_0_0",resolve.root=TRUE))
  }
  
  ##Calculating distances to ML tree
  results.raxml$RF.Dist.ML[i]<- TreeDist::InfoRobinsonFoulds(trees.both,normalize=TRUE)[1]
  #results.raxml$weighted.rf[i] <- wRF.dist(raxml.trees[[i]],ml.tree[[j]]) # this always seems to come out wonky
  results.raxml$Q.Dist.ML[i] <- Quartet::QuartetDivergence(Quartet::QuartetStatus(trees.both), similarity=FALSE)[-1]/(2/3) #2/3 is max Q
  results.raxml$CI.Dist.ML[i] <- as.matrix(ClusteringInfoDistance(trees.both,normalize = TRUE))[-1,1] #normalized relative to max CI for the tree

  ## normalized gamma and imbalance
  results.raxml$std.gamma[i] <- gammaStat(chronos(trees.both.root[[1]],lambda=1,quiet=TRUE))[1] - gammaStat(trees.both.root[[2]])[1]
  results.raxml$std.colless[i] <- colless(as.treeshape(trees.both.root[[1]], model="pda"),norm="pda")[1] - colless(as.treeshape(trees.both.root[[2]],model="pda"),norm="pda")
  results.raxml$std.sackin[i] <- sackin(as.treeshape(trees.both.root[[1]], model="pda"),norm="pda") - sackin(as.treeshape(trees.both.root[[2]],model="pda"),norm="pda")

  results.raxml$std.sites[i] <- results.raxml$sites[i] / max(results.raxml$sites[results.raxml$simulation == s])

   if("sim_0_0_0" %in% trees.both.root[[1]]$tip.label){
     ingroup <- drop.tip(trees.both.root[[1]],"sim_0_0_0")
   } else {
     ingroup <- trees.both.root[[1]]
   }

   if("sim_0_0_0" %in% trees.both.root[[2]]$tip.label){
     ml.ingroup <- drop.tip(trees.both.root[[2]],"sim_0_0_0")
   } else {
     ml.ingroup <- trees.both.root[[2]]
   }

  results.raxml$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda"),norm="pda") - colless(as.treeshape(ml.ingroup,model="pda"),norm="pda")
  results.raxml$std.ingroup.gamma[i] <- gammaStat(chronos(ingroup,lambda=1))[1] - gammaStat(chronos(ml.ingroup,lambda=1))[1]
  results.raxml$std.ingroup.sackin[i] <- sackin(as.treeshape(ingroup, model="pda"),norm="pda") - sackin(as.treeshape(ml.ingroup,model="pda"),norm="pda")
}

summary(results.raxml)  
write.csv(results.raxml,
          file=paste("output/new/",args$output,"-raxml.csv",sep=""),
          quote=FALSE,row.names=TRUE,na="NA")

######################
## univariate plots
#####################
results.raxml$maf <- as.factor(results.raxml$maf)
plot1 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot2 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=Avg.BLs,
                    fill=NULL)) +
  geom_boxplot(aes(alpha=0.9),outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot3 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=ingroup.tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot4 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=ingroup.colless,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot5 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=ingroup.gamma,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot6 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=RF.Dist.ML,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.1,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot7 <- ggplot(data = results.raxml, 
                aes(x=height,
                    y=Q.Dist.ML,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.1,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,nrow=1) 


# subsamp.plot1 <- ggplot(data = results.raxml.subsamp, 
#                 aes(x=height,
#                     y=tree.height,
#                     fill=NULL)) +
#   geom_boxplot(alpha=0.9,outlier.shape=NA) +
#   geom_jitter(aes(col=height),width=0.1,height=0.1,alpha=0.4) +
#   theme_custom() +
#   ylim(0,0.5) +
#   theme(legend.position="none",
#         axis.title.x = element_blank())
# 
# subsamp.plot2 <- ggplot(data = results.raxml.subsamp, 
#                 aes(x=height,
#                     y=Avg.BLs,
#                     fill=NULL)) +
#   geom_boxplot(aes(alpha=0.9),outlier.shape=NA) +
#   geom_jitter(aes(col=height),width=0.1,height=0.1,alpha=0.4) +  theme_custom()  +
#   theme(legend.position="none",
#         axis.title.x = element_blank())
# 
# subsamp.plot3 <- ggplot(data = results.raxml.subsamp, 
#                 aes(x=height,
#                     y=ingroup.gamma,
#                     fill=NULL)) +
#   geom_boxplot(alpha=0.9,outlier.shape=NA) +
#   geom_jitter(aes(col=height),width=0.1,height=0.1,alpha=0.4) +  theme_custom()  +
#   theme(legend.position="none",
#         axis.title.x = element_blank()) 
# 
# ggarrange(subsamp.plot1,subsamp.plot2,subsamp.plot3,nrow=1) 
# 

####################BEGIN ANALYSIS ################################
####Calculate Matrix of RF distances for all trees with the same species tree in multi tree object
#results.raxml <- read.csv(file=paste("output/new/",output,"-raxml.csv",sep=""),row.names=1,na="NA",header=TRUE)
pdf(paste("output/new/",args$output,".pdf",sep=""))
# for (h in unique(as.character(results.raxml$height))){
#   for (i in unique(as.numeric(as.character(results.raxml$simulation)))){
#   #for (i in seq(16,18,by=1)){
#     j <- which(ml.tree.info$simulation == i & ml.tree.info$height == h)
#     subset <- which(results.raxml$simulation == i & results.raxml$height == h & results.raxml$noref == "REF")
#     if (length(subset) != 0) {
#       trees.subset <- c(raxml.trees[subset],ml.tree[[j]])
#     } else {
#       print(paste("no trees for sim",i," for height ",h)) 
#       next
#     }
#     rf_matrix <- multiRF(trees.subset)
#     
#     ####PcOA of RF distance matrix for plotting trees in tree space
#     rf_pcoa <- pcoa(rf_matrix)
#     
#     rf_df <- data.frame(axis1=rf_pcoa$vectors[,1],axis2=rf_pcoa$vectors[,2])
#     biplot <- ggplot(rf_df,aes(axis1,axis2))
#     plot <- biplot +
#       geom_jitter(aes(color=as.factor(c(results.raxml$int[subset],"ML_tree")),shape=as.factor(c(results.raxml$method[subset],"ML_tree"))),alpha=0.5,size=5,height=2,width=2)+
#       theme_bw()+
#       theme(legend.text = element_text(size=rel(1.5)),
#             legend.title = element_blank(),
#             plot.margin = unit(c(6,5.5,20,10),"points"),
#             line = element_line(size=1),
#             axis.title = element_text(size=rel(1.5)),
#             axis.text = element_text(size=rel(1)),
#             panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank())+
#       ggtitle(paste("PCoA on RF Distances, Sim",i,", height ",h, sep=""))+
#       xlab(paste("PCoA Axis 1 (",round(rf_pcoa$values$Relative_eig[1]*100,1),"%)",sep=""))+
#       ylab(paste("PCoA Axis 2 (",round(rf_pcoa$values$Relative_eig[2]*100,1),"%)",sep=""))+
#       #geom_label_repel(label=c(paste("MAF",results.raxml$maf[subset],"MISS",results.raxml$missing[subset],sep=" "),"truth"),size=rel(1))+
#       scale_color_manual(labels = c("EXT","INT","truth"),values=c("#009980", "#006699", "black"))+
#       scale_shape_manual(labels = c("raxml","ml"),values=c(1,16,17))
#     print(plot)
#     
#     plot2 <- biplot +
#       geom_jitter(aes(color=as.factor(c(results.raxml$maf[subset],6)),shape=as.factor(c(results.raxml$missing[subset],"true"))),alpha=0.5,size=5,height=2,width=2)+
#       theme_bw()+
#       theme(legend.text = element_text(size=rel(1.5)),
#             legend.title = element_blank(),
#             plot.margin = unit(c(6,5.5,20,10),"points"),
#             line = element_line(size=1),
#             axis.title = element_text(size=rel(1.5)),
#             axis.text = element_text(size=rel(1)),
#             panel.grid.major = element_blank(), 
#             panel.grid.minor = element_blank())+
#       ggtitle(paste("PCoA on RF Distances, Sim",i,", height ",h, sep=""))+
#       xlab(paste("PCoA Axis 1 (",round(rf_pcoa$values$Relative_eig[1]*100,1),"%)",sep=""))+
#       ylab(paste("PCoA Axis 2 (",round(rf_pcoa$values$Relative_eig[2]*100,1),"%)",sep=""))+
#       #geom_label_repel(label=c(paste("MISS",results.raxml[subset,3],"Q",results.raxml[subset,2],sep=","),"truth"),size=rel(1))+
#       scale_color_manual(values=c("#009980", "#006699","magenta","lightblue","orange","green","black","pink","turquoise"))+
#       scale_shape_manual(values=c(0,1,2,3,16,17))
#     print(plot2)
#     
#     #heatmap(rf_matrix,labCol=FALSE,labRow=paste(results.raxml[subset,]$maf,results.raxml[subset,]$missing,results.raxml[subset,]$int,sep=","),cex.lab=0.5)
#   }
# }
# dev.off()

## Plot trees in PCoA space by simulation
results.raxml <- results.raxml[results.raxml$simulation > 15 & results.raxml$simulation < 26,]
for (h in unique(as.character(results.raxml$height))){
  for (i in unique(as.numeric(as.character(results.raxml$simulation)))){
    #for (i in seq(16,30,by=1)){
    j <- which(ml.tree.info$simulation == i & ml.tree.info$height == h)
    subset <- which(results.raxml$simulation == i & results.raxml$height == h & results.raxml$noref == "REF")
    if (length(subset) != 0) {
      trees.subset <- c(raxml.trees[subset],ml.tree[[j]])
    } else {
      print(paste("no trees for sim",i," for height ",h))
      next
    }
    
    #trees.subset <- sapply(trees.subset,function(x) root(x,outgroup,resolve.root=TRUE))
    #trees.subset.root <- root(trees.subset,outgroup,resolve.root=TRUE)
    class(trees.subset) <- "multiPhylo"
    rf.dist <- InfoRobinsonFoulds(trees.subset,normalize=TRUE)
    rf.pcoa <- pcoa(rf.dist)
    rf.pcoa.plot <- rf.pcoa$vectors %>%
      as_tibble() %>%
      ggplot(aes(x=Axis.1,y=Axis.2)) +
      geom_jitter(aes(col=factor(c(results.raxml$maf[subset],"true_tree")),
                     #text=results.raxml$missing[results.raxml$simulation == s],
                     shape=factor(c(results.raxml$missing[subset],"true_tree")),
                     size=factor(c(results.raxml$int[subset],"true_tree"))),
                 alpha=0.6) +
      theme_bw(base_size=12, base_family="Arial") +
      xlab("PCoA Axis 1") +
      ylab("PCoA Axis 2") +
      theme(
        panel.background  = element_blank(),
        plot.background = element_rect(fill="white", color=NA), 
        legend.background = element_rect(fill="transparent", colour=NA),
        legend.key = element_rect(fill="transparent", colour=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=18, family="Arial"),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank()
      ) +
      ggtitle(paste0("Full Sim PCoA, Simulation ",i,", ",h)) +
      geom_hline(yintercept=0,lty=2,col="gray50") +
      geom_vline(xintercept=0,lty=2,col="gray50") +
      # biplot arrows for MAF
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.character(results.raxml$maf[subset])), rf.pcoa$vectors[,1][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
        yend = cor(as.numeric(as.character(results.raxml$maf[subset])), rf.pcoa$vectors[,2][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1),
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1, 
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) + 
      annotate(geom = "text",
               x = cor(as.numeric(as.character(results.raxml$maf[subset])), rf.pcoa$vectors[,1][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
               y = cor(as.numeric(as.character(results.raxml$maf[subset])), rf.pcoa$vectors[,2][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
               label = "MAF",
               hjust=0,vjust=0
      ) +
      # biplot arrows for Missing
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.character(results.raxml$missing[subset])), rf.pcoa$vectors[,1][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
        yend = cor(as.numeric(as.character(results.raxml$missing[subset])), rf.pcoa$vectors[,2][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1),
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1, 
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) + 
      annotate(geom = "text",
               x = cor(as.numeric(as.character(results.raxml$missing[subset])), rf.pcoa$vectors[,1][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
               y = cor(as.numeric(as.character(results.raxml$missing[subset])), rf.pcoa$vectors[,2][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
               label = "Missing",
               hjust=1,vjust=1
      ) +
      # biplot arrows for INT
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.factor(results.raxml$int[subset])), rf.pcoa$vectors[,1][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
        yend = cor(as.numeric(as.factor(results.raxml$int[subset])), rf.pcoa$vectors[,2][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1),
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1, 
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) + 
      annotate(geom = "text",
               x = cor(as.numeric(as.factor(results.raxml$int[subset])), rf.pcoa$vectors[,1][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
               y = cor(as.numeric(as.factor(results.raxml$int[subset])), rf.pcoa$vectors[,2][1:length(subset)]) * 0.5 * sqrt(length(subset) - 1), 
               label = "INT",
               hjust=0,vjust=0
      )
    
    
    print(rf.pcoa.plot)
    #plotly::ggplotly(rf.pcoa.plot)
  }
}

dev.off()
