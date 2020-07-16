############################################
## Script for analyzing refbias output trees from RAxML
## Written by J. Rick, 11 March 2019
## Made to be run on the command line
############################################

suppressMessages(
  c(library(ape),
  library(phytools),
  library(apTreeshape),
  library(phangorn),
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
  library(ggsci)))

############################################
## For debugging-- specifying files ########
############################################

#raxml.trees <- "103019-all-raxml.trees"
#raxml.tree.names <- "103019-all-raxml.names"
#astral.trees <- "101819-all-astral.trees"
#astral.tree.names <- "101819-all-astral.names"
#ml.trees <- "103019-s_tree.tree"
#ml.tree.names <- "103019-s_tree.names"
#output <- "103019-output"
#refdist <- "103019-refdist.txt"
#args <- mget(c("raxml.trees","raxml.tree.names","astral.trees","astral.tree.names","ml.trees","ml.tree.names","refdist","output"))

############################################
## Reading in command line arguments #######
############################################

#check for the argparse package installed
if (!('argparse' %in% installed.packages()[, 'Package'])) {
  install.packages('argparse', repos='http://cran.rstudio.com/')
}

suppressPackageStartupMessages(library('argparse'))

parse.args <- function() {
  parser <- ArgumentParser()
  parser$add_argument('--ml.trees', help='file with ML (true) trees')
  parser$add_argument('--ml.tree.names', help='file with names of ML trees')
  parser$add_argument('--raxml.trees', help='file with raxml trees')
  parser$add_argument('--raxml.tree.names', help='file with raxml tree names')
  parser$add_argument('--astral.trees', help='file with astral trees')
  parser$add_argument('--astral.tree.names', help='file with astral tree names')
  parser$add_argument('--refdist', help='file with reference distances')
  parser$add_argument('-o', '--output', help='output file prefix')
  return(parser$parse_args())
}
args <- parse.args()


####Read in concatenated tree file and ML tree file (make sure the latter is ultrametric)
raxml.trees<-read.tree(paste("output/new/",args$raxml.trees,sep=""))
#astral.trees<-read.tree(paste("output/",args$astral.trees,sep=""))
ml.tree<-read.tree(paste("output/new/",args$ml.trees,sep=""))

###Read in file of tree names
raxml.tree.names<-read.table(paste("output/new/",args$raxml.tree.names,sep=""),stringsAsFactors=FALSE)
#astral.tree.names<-read.table(paste("output/",args$astral.tree.names,sep=""),stringsAsFactors = FALSE)
ml.tree.names <- read.table(paste("output/new/",args$ml.tree.names,sep=""),stringsAsFactors = FALSE)

###Pull info about ml species trees
ml.tree.info <- data.frame(simulation = numeric(nrow(ml.tree.names)),
		height=character(nrow(ml.tree.names)))
for (i in 1:nrow(ml.tree.names)){
  ml.tree.info$simulation[i] <- as.integer(regmatches(ml.tree.names[i,1], regexec('sim([0-9]+)', ml.tree.names[i,1]))[[1]][2])
  ml.tree.info$height[i] <- as.factor(regmatches(ml.tree.names[i,1], regexec('height([A-Z]+)', ml.tree.names[i,1]))[[1]][2])
}

refdist <- read.table(paste("output/new/",args$refdist,sep=""),header=T,stringsAsFactors = FALSE)

###################
## changing tip labels on ml trees to match the sim trees
## ONLY DO ONCE!!
####################
for (i in 1:length(ml.tree)){
  tips <- ml.tree[[i]]$tip.label
  tips.new <- paste("sim_",tips,"_0_0",sep="")
  ml.tree[[i]]$tip.label <- tips.new
}

for (i in 1:length(raxml.trees)){
  if ("reference" %in% raxml.trees[[i]]$tip.label){
    print(i)
    new <- drop.tip(raxml.trees[[i]],"reference")
    raxml.trees[[i]] <- new
  }
}

###################
## changing tip labels on astral trees to match the sim trees
## ONLY DO ONCE!!
####################
# for (i in 1:length(astral.trees)){
#   tips <- astral.trees[[i]]$tip.label
#   tips.new <- unlist(strsplit(unlist(strsplit(tips,"aln_")),".sorted.bam"))
#   astral.trees[[i]]$tip.label <- tips.new
# }

#####################
## Start of analysis!
#####################

###Create empty data frame and with named, empty columns 
num.trees <- nrow(raxml.tree.names)
results.raxml<-data.frame(simulation=integer(num.trees),
                    height=character(num.trees),
                    method=character(num.trees),
                    quality=integer(num.trees), 
                    missing=numeric(num.trees), 
                    maf=numeric(num.trees), 
                    int=character(num.trees), 
                    noref=character(num.trees),
                    taxa_ref=character(num.trees),
                    refdist=numeric(num.trees),
                    sites=integer(num.trees), 
                    std.sites=numeric(num.trees),
                    gamma=integer(num.trees), 
                    colless=integer(num.trees), 
                    sackin=integer(num.trees), 
                    tree.height=integer(num.trees), 
                    Avg.BLs=integer(num.trees), 
                    SD.BLs=integer(num.trees), 
                    RF.Dist.ML=integer(num.trees), stringsAsFactors=FALSE)

for (i in 1:length(raxml.trees)) {
  
  ####Use regex and sub to extract number for filtering and sim parameters and add them to data frame
  results.raxml[i,1] <- as.integer(regmatches(raxml.tree.names[i,1], regexec('-s([0-9]+)\\_q', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$height[i] <- as.character(regmatches(raxml.tree.names[i,1], regexec('([0-9]+)-s[0-9]', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$method[i] <- as.character("raxml")
  results.raxml$quality[i] <- as.integer(regmatches(raxml.tree.names[i,1], regexec('q(.*?)\\_m', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$missing[i] <- as.numeric(regmatches(raxml.tree.names[i,1], regexec('miss(.*?)\\_m', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$maf[i] <- as.numeric(regmatches(raxml.tree.names[i,1], regexec('maf(0.*?)\\.[A-Z]', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$sites[i] <- as.integer(regmatches(raxml.tree.names[i,1], regexec('sites([0-9]*?)\\.', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$int[i] <- (regmatches(raxml.tree.names[i,1], regexec('.([A-Z]+)-', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$noref[i] <- (regmatches(raxml.tree.names[i,1], regexec('[0-9].([A-Z]+)\\.[A-Z]', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$taxa_ref[i] <- regmatches(raxml.tree.names[i,1],regexec('[A-Z]-([0-9]+_0_0)\\.phylip', raxml.tree.names[i,1]))[[1]][[2]]
  results.raxml$refdist[i] <- refdist$avg_dist[refdist$sim == results.raxml$simulation[i] &
                                               refdist$tree_height == results.raxml$height[i] &
                                               refdist$int == results.raxml$int[i]][1]
}

results.raxml$int <- case_when(results.raxml$taxa_ref == "0_0_0" ~ "EXT",
                               TRUE ~ "INT")
results.raxml$noref <- "REF"

for (i in 1:length(raxml.trees)){
  ###Root on outgroup and make ultrametric
  if(!("sim_0_0_0" %in% raxml.trees[[i]]$tip.label)){
    rooted <- raxml.trees[[i]]
  } else {
    rooted <- root(raxml.trees[[i]],"sim_0_0_0",resolve.root=TRUE)
  }
  
  ultra <- chronopl(rooted, lambda=1, iter.max=100)
  
  ###Replace NA branch lengths
  #raxml.trees[[i]]$edge.length[is.na(raxml.trees[[i]]$edge.length)] <- 0
  
  ###Gamma stat
  results.raxml$gamma[i] <- gammaStat(ultra)[1]

  ###Gamma on ingroup only
  if(!("sim_0_0_0" %in% raxml.trees[[i]]$tip.label)){
    ingroup <- ultra
  } else {
    ingroup <- drop.tip(ultra,"sim_0_0_0")
  }

  results.raxml$ingroup.gamma[i] <- gammaStat(ingroup)[1]

  ###Colless or Ic stat
  results.raxml$colless[i] <- colless(as.treeshape(ultra, model="pda"),norm="pda")
  results.raxml$ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda"),norm="pda")

  ###Sackin stat
  results.raxml$sackin[i] <- sackin(as.treeshape(ultra, model="pda"))
  results.raxml$ingroup.sackin[i] <- sackin(as.treeshape(ingroup, model="pda"))

  ###Total tree height
  results.raxml$tree.height[i] <- max(branching.times(rooted),na.rm=T)

  ###Mean BLs
  results.raxml$Avg.BLs[i] <- mean(rooted$edge.length,na.rm=T)

  ###SD BLs
  results.raxml$SD.BLs[i] <- sd(rooted$edge.length,na.rm=T)

  ####Add mean and SD of BS support values
  results.raxml$mean.support[i] <- mean(as.numeric(rooted$node.label[-1]),na.rm=T)
  results.raxml$sd.support[i] <- sd(as.numeric(rooted$node.label[-1]),na.rm=T)
  
}

###################
## statistics requiring separation by simulation number
## to compare to ML tree
###################

for (i in 1:length(raxml.trees)){
  s <- as.integer(results.raxml$simulation[i])
  h <- as.character(results.raxml$height[i])
  j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
  
  ml.tree.pruned <- drop.tip(ml.tree[[j]],results.raxml$taxa_ref[i])
  
  if(results.raxml$taxa_ref[i] %in% raxml.trees[[i]]$tip.label){
    trees.both <- as.multiPhylo(c(raxml.trees[[i]],ml.tree[[j]]))
  } else {
    trees.both <- as.multiPhylo(c(raxml.trees[[i]],ml.tree.pruned))
  }

  if(results.raxml$taxa_ref[i] %in% raxml.trees[[i]]$tip.label){
    trees.both.root <- as.multiPhylo(c(root(raxml.trees[[i]],"sim_0_0_0",resolve.root=TRUE),
                                       root(ml.tree[[j]],"sim_0_0_0",resolve.root=TRUE)))
  } else {
    trees.both.root <- as.multiPhylo(c(root(raxml.trees[[i]],"sim_0_0_0",resolve.root=TRUE),
                                       root(ml.tree.pruned,"sim_0_0_0",resolve.root=TRUE)))
  }
  
  ###RF distance to ML tree
  results.raxml$RF.Dist.ML[i]<- multiRF(trees.both,multi2di=TRUE)[[2]]
  results.raxml$weighted.rf[i] <- wRF.dist(raxml.trees[[i]],ml.tree[[j]])

  ## normalized gamma and imbalance
  results.raxml$std.gamma[i] <- gammaStat(trees.both.root[[1]])[1] - gammaStat(trees.both.root[[2]])[1]
  results.raxml$std.colless[i] <- colless(as.treeshape(trees.both.root[[1]], model="pda")) - colless(as.treeshape(trees.both.root[[2]],model="pda"))

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
        
  results.raxml$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda")) - colless(as.treeshape(ml.ingroup,model="pda"))
  results.raxml$std.ingroup.gamma[i] <- gammaStat(ingroup)[1] - gammaStat(ml.ingroup)[1]
}

summary(results.raxml)
write.csv(results.raxml,file=paste("output/",args$output,"-raxml.csv",sep=""),quote=FALSE,row.names=TRUE,na="NA")

#########################
####Clean up dataframe
# #########################
# results.raxml$int <- as.factor(results.raxml$int)
# results.raxml$noref <- as.factor(results.raxml$noref)
# results.raxml$simulation <- as.factor(results.raxml$simulation)
# results.raxml$missing <- as.factor(results.raxml$missing)
# results.raxml$method <- as.factor(results.raxml$method)
# 

####As a start, visualize correlations between all variables
#pdf(paste("output/pairs-",args$output,".pdf",sep=""))
#psych::pairs.panels(results.raxml[,c("missing","maf","int","gamma","colless","sackin","tree.height","Avg.BLs","SD.BLs","ingroup.gamma","ingroup.colless","ingroup.sackin","mean.support","sd.support")],cex.cor=2)
#dev.off()

####################BEGIN ANALYSIS ################################
####Calculate Matrix of RF distances for all trees with the same species tree in multi tree object
results.raxml <- read.csv(file=paste("output/",args$output,"-raxml.csv",sep=""),row.names=1,na="NA",header=TRUE)
pdf(paste("output/",args$output,".pdf",sep=""))
for (h in unique(as.character(results.raxml$height))){
  for (i in unique(as.numeric(as.character(results.raxml$simulation)))){
    
    j <- which(ml.tree.info$simulation == i & ml.tree.info$height == h)
    subset <- which(results.raxml$simulation == i & results.raxml$height == h & results.raxml$noref == "REF")
    if (length(subset) != 0) {
      trees.subset <- c(raxml.trees[subset],ml.tree[[j]])
    } else {
      print(paste("no trees for sim",i," for height ",h)) 
      next
    }
    rf_matrix<-multiRF(trees.subset)
    
    ####PcOA of RF distance matrix for plotting trees in tree space
    rf_pcoa <- pcoa(rf_matrix)
    
    rf_df <- data.frame(axis1=rf_pcoa$vectors[,1],axis2=rf_pcoa$vectors[,2])
    biplot <- ggplot(rf_df,aes(axis1,axis2))
    plot <- biplot +
      geom_jitter(aes(color=as.factor(c(results.raxml$int[subset],3)),shape=as.factor(c(results.raxml$method[subset],3))),alpha=0.5,size=5,height=2,width=2)+
      theme_bw()+
      theme(legend.text = element_text(size=rel(1.5)),
            legend.title = element_blank(),
            plot.margin = unit(c(6,5.5,20,10),"points"),
            line = element_line(size=1),
            axis.title = element_text(size=rel(1.5)),
            axis.text = element_text(size=rel(1)),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
      ggtitle(paste("PCoA on RF Distances, Sim",i,", height ",h, sep=""))+
      xlab(paste("PCoA Axis 1 (",round(rf_pcoa$values$Relative_eig[1]*100,1),"%)",sep=""))+
      ylab(paste("PCoA Axis 2 (",round(rf_pcoa$values$Relative_eig[2]*100,1),"%)",sep=""))+
      #geom_label_repel(label=c(paste("MAF",results.raxml$maf[subset],"MISS",results.raxml$missing[subset],sep=" "),"truth"),size=rel(1))+
      scale_color_manual(labels = c("EXT","INT","truth"),values=c("#009980", "#006699", "black"))+
      scale_shape_manual(labels = c("raxml","ml"),values=c(1,16,17))
    print(plot)
   
    plot2 <- biplot +
      geom_jitter(aes(color=as.factor(c(results.raxml$maf[subset],6)),shape=as.factor(c(results.raxml$miss[subset],"true"))),alpha=0.5,size=5,height=2,width=2)+
      theme_bw()+
      theme(legend.text = element_text(size=rel(1.5)),
            legend.title = element_blank(),
            plot.margin = unit(c(6,5.5,20,10),"points"),
            line = element_line(size=1),
            axis.title = element_text(size=rel(1.5)),
            axis.text = element_text(size=rel(1)),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
      ggtitle(paste("PCoA on RF Distances, Sim",i,", height ",h, sep=""))+
      xlab(paste("PCoA Axis 1 (",round(rf_pcoa$values$Relative_eig[1]*100,1),"%)",sep=""))+
      ylab(paste("PCoA Axis 2 (",round(rf_pcoa$values$Relative_eig[2]*100,1),"%)",sep=""))+
      #geom_label_repel(label=c(paste("MISS",results.raxml[subset,3],"Q",results.raxml[subset,2],sep=","),"truth"),size=rel(1))+
      scale_color_manual(values=c("#009980", "#006699","magenta","lightblue","orange","green","black","pink","turquoise"))+
      scale_shape_manual(values=c(0,1,2,3,16,17))
    print(plot2)
   
    #heatmap(rf_matrix,labCol=FALSE,labRow=paste(results.raxml[subset,]$maf,results.raxml[subset,]$miss,results.raxml[subset,]$int,sep=","),cex.lab=0.5)
  }
}
dev.off()

#write.csv(results.raxml,file=paste("output/",args$output,"-raxml.csv",sep=""),quote=FALSE,row.names=TRUE,col.names=TRUE,na="NA")
