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
  library(readtext),
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

raxml.trees <- "031419-all-raxml.trees"
raxml.tree.names <- "031419-all-raxml.names"
astral.trees <- "031419-all-astral.trees"
astral.tree.names <- "031419-all-astral.names"
ml.trees <- "031419-s_tree.trees"
ml.tree.names <- "031419-s_tree.names"
output <- "031419-output"
args <- mget(c("raxml.trees","raxml.tree.names","astral.trees","astral.tree.names","ml.trees","ml.tree.names","output"))

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
  parser$add_argument('-o', '--output', help='output file prefix')
  return(parser$parse_args())
}
args <- parse.args()


####Read in concatenated tree file and ML tree file (make sure the latter is ultrametric)
raxml.trees<-read.tree(paste("output/",args$raxml.trees,sep=""))
#astral.trees<-read.tree(paste("output/",args$astral.trees,sep=""))
ml.tree<-read.tree(paste("output/",args$ml.trees,sep=""))

###Read in file of tree names
raxml.tree.names<-read.table(paste("output/",args$raxml.tree.names,sep=""),stringsAsFactors=FALSE)
#astral.tree.names<-read.table(paste("output/",args$astral.tree.names,sep=""),stringsAsFactors = FALSE)
ml.tree.names <- read.table(paste("output/",args$ml.tree.names,sep=""),stringsAsFactors = FALSE)

###Pull info about ml species trees
ml.tree.info <- data.frame(simulation = numeric(nrow(ml.tree.names)),height=numeric(nrow(ml.tree.names)))
for (i in 1:nrow(ml.tree.names)){
  ml.tree.info$simulation[i] <- as.integer(regmatches(ml.tree.names[i,1], regexec('sim([0-9]+)', ml.tree.names[i,1]))[[1]][2])
  ml.tree.info$height[i] <- as.integer(regmatches(ml.tree.names[i,1], regexec('height([0-9]+)', ml.tree.names[i,1]))[[1]][2])
}

###################
## changing tip labels on ml trees to match the sim trees
## ONLY DO ONCE!!
####################
for (i in 1:length(ml.tree)){
  tips <- ml.tree[[i]]$tip.label
  tips.new <- paste("sim_",tips,"_0_0",sep="")
  ml.tree[[i]]$tip.label <- tips.new
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

results.raxml<-data.frame(simulation=integer(),
                    height=integer(),
                    method=character(),
                    quality=integer(), 
                    missing=numeric(), 
                    maf=numeric(), 
                    int=character(), 
                    noref=character(),
                    taxa_ref=character(),
                    sites=integer(), 
                    std.sites=numeric(),
                    gamma=integer(), 
                    colless=integer(), 
                    sackin=integer(), 
                    tree.height=integer(), 
                    Avg.BLs=integer(), 
                    SD.BLs=integer(), 
                    RF.Dist.ML=integer())
results.raxml$int <- as.character(results.raxml$int)
results.raxml$noref <- as.character(results.raxml$noref)
results.raxml$taxa_ref <- as.character(results.raxml$taxa_ref)
results.raxml$method <- as.character(results.raxml$method)

for (i in 1:length(raxml.trees)) {
  
  ####Use regex and sub to extract number for filtering and sim parameters and add them to data frame
  results.raxml[i,1]<-as.integer(regmatches(raxml.tree.names[i,1], regexec('-s([0-9]+)\\_q', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$height[i] <-as.integer(regmatches(raxml.tree.names[i,1], regexec('([0-9]+)-s[0-9]', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$method[i]<- as.character("raxml")
  results.raxml$quality[i]<-as.integer(regmatches(raxml.tree.names[i,1], regexec('q(.*?)\\_m', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$missing[i]<-as.numeric(regmatches(raxml.tree.names[i,1], regexec('miss(.*?)\\_m', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$maf[i]<-as.numeric(regmatches(raxml.tree.names[i,1], regexec('maf(0.*?)\\.[A-Z]', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$sites[i]<-as.integer(regmatches(raxml.tree.names[i,1], regexec('sites([0-9]*?)\\.', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$int[i]<-(regmatches(raxml.tree.names[i,1], regexec('.([A-Z]+).filtered', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$noref[i]<-(regmatches(raxml.tree.names[i,1], regexec('[0-9].([A-Z]+)\\.[A-Z]', raxml.tree.names[i,1]))[[1]][2])
  results.raxml$taxa_ref[i]<-regmatches(raxml.tree.names[i,1],regexec('ref-([sim[0-9]+_0_0)\\.phylip', raxml.tree.names[i,1]))[[1]][[2]]
}

# for (i in 1:length(astral.trees)) {
# ####Use regex and sub to extract number for filtering and sim parameters and add them to data frame
#   j <- i + length(raxml.trees)
#   results[j,1]<-as.integer(regmatches(astral.tree.names[i,1], regexec('s([0-9]+)\\_q', astral.tree.names[i,1]))[[1]][2])
#   results$height[j] <-as.integer(regmatches(astral.tree.names[i,1], regexec('([0-9]+)\\_s[0-9]', astral.tree.names[i,1]))[[1]][2])
#   results$method[j]<- as.character("astral")
#   results$quality[j]<-as.integer(regmatches(astral.tree.names[i,1], regexec('q(.*?)\\_m', astral.tree.names[i,1]))[[1]][2])
#   results$missing[j]<-as.numeric(regmatches(astral.tree.names[i,1], regexec('miss(.*?)\\_m', astral.tree.names[i,1]))[[1]][2])
#   results$maf[j]<-as.numeric(regmatches(astral.tree.names[i,1], regexec('maf(0.*?)\\.[A-Z]', astral.tree.names[i,1]))[[1]][2])
#   results$sites[j]<-as.integer(regmatches(astral.tree.names[i,1], regexec('sites([0-9]*?)\\.', astral.tree.names[i,1]))[[1]][2])
#   results$int[j]<-(regmatches(astral.tree.names[i,1], regexec('.([A-Z]+).ref', astral.tree.names[i,1]))[[1]][2])
#   results$noref[j]<-(regmatches(astral.tree.names[i,1], regexec('[0-9].([A-Z]+)', astral.tree.names[i,1]))[[1]][2])
#   results$taxa_ref[j]<-regmatches(astral.tree.names[i,1], regexec('ref-([0-9]+_0_0)\\.astral', astral.tree.names[i,1]))[[1]][2]
# }

# combined.trees <- as.multiPhylo(c(raxml.trees,astral.trees))
# combined.names <- rbind(raxml.tree.names,astral.tree.names)
#row.names(results) <- as.character(combined.names)

# rm(raxml.trees)
# rm(astral.trees)

for (i in 1:length(raxml.trees)){
  ###Root on outgroup?
  #rooted <- root(raxml.trees[[i]],results.raxml$taxa_ref[i])
  
  ###Replace NA branch lengths
  raxml.trees[[i]]$edge.length[is.na(raxml.trees[[i]]$edge.length)] <- 0
  
  ###Gamma stat
  results.raxml$gamma[i]<-gammaStat(raxml.trees[[i]])[1]
  
  ###Gamma on ingroup only
  ingroup <- drop.tip(raxml.trees[[i]],"sim_0_0_0")
  results.raxml$ingroup.gamma[i] <- gammaStat(ingroup)[1]
  
  ###Colless or Ic stat
  results.raxml$colless[i]<-colless(as.treeshape(raxml.trees[[i]], model="pda"))
  results.raxml$ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda"))
  
  ###Sackin stat
  results.raxml$sackin[i]<-sackin(as.treeshape(raxml.trees[[i]], model="pda"))
  results.raxml$ingroup.sackin[i]<-sackin(as.treeshape(ingroup, model="pda"))
  
  ###Total tree height
  results.raxml$tree.height[i]<-max(branching.times(raxml.trees[[i]]),na.rm=T)
  
  ###Mean BLs
  results.raxml$Avg.BLs[i]<-mean(raxml.trees[[i]]$edge.length,na.rm=T)
  
  ###SD BLs
  results.raxml$SD.BLs[i]<-sd(raxml.trees[[i]]$edge.length,na.rm=T)
  
  ####Add mean and SD of BS support values?
  results.raxml$mean.support<-mean(raxml.trees[[i]]$node.label,na.rm=T)
  results.raxml$sd.support<-sd(raxml.trees[[i]]$node.label,na.rm=T)
  
}

###################
## statistics requiring separation by simulation number
## to compare to ML tree
###################

for (i in 1:length(raxml.trees)){
  if (results.raxml$taxa_ref[i] == "0_0_0"){
    if (results.raxml$noref[i] == "REF"){
      s <- as.integer(results.raxml$simulation[i])
      h <- as.integer(results.raxml$height[i])
      j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
      trees.both <- as.multiPhylo(c(raxml.trees[[i]],ml.tree[[j]]))

      ###RF distance to ML tree
      results.raxml$RF.Dist.ML[i]<- multiRF(trees.both,multi2di=TRUE)[[2]]
      results.raxml$weighted.rf[i] <- wRF.dist(raxml.trees[[i]],ml.tree[[j]],check.labels=FALSE)
      
      ## normalized gamma
      results.raxml$std.gamma[i] <- gammaStat(raxml.trees[[i]])[1] - gammaStat(ml.tree[[j]])[1]

      results.raxml$std.colless[i] <- colless(as.treeshape(raxml.trees[[i]], model="pda")) - colless(as.treeshape(ml.tree[[j]],model="pda"))

      results.raxml$std.sites[i] <- results.raxml$sites[i] / max(results.raxml$sites[results.raxml$simulation == s])
      
      ingroup <- drop.tip(raxml.trees[[i]],"sim_0_0_0")
      ml.ingroup <- drop.tip(ml.tree[[j]],"sim_0_0_0")
      results.raxml$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda")) - colless(as.treeshape(ml.ingroup,model="pda"))
    } else {
      s <- as.integer(results.raxml$simulation[i])
      h <- as.integer(results.raxml$height[i])
      j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
      ml.tree.pruned <- drop.tip(ml.tree[[j]],"sim_0_0_0")
      trees.both <- as.multiPhylo(c(raxml.trees[[i]],ml.tree.pruned))
      
      ###RF distance to ML tree
      results.raxml$RF.Dist.ML[i]<- multiRF(trees.both,multi2di=TRUE)[[2]]
      results.raxml$weighted.rf[i] <- wRF.dist(raxml.trees[[i]],ml.tree.pruned,check.labels=FALSE)
       
      ## normalized gamma
      results.raxml$std.gamma[i] <- gammaStat(raxml.trees[[i]])[1] - gammaStat(ml.tree.pruned)[1]
      
      results.raxml$std.colless[i] <- colless(as.treeshape(raxml.trees[[i]], model="pda")) - colless(as.treeshape(ml.tree.pruned,model="pda"))
      
      results.raxml$std.sites[i] <- results.raxml$sites[i] / max(results.raxml$sites[results.raxml$simulation == s])
      
      ingroup <- drop.tip(raxml.trees[[i]],"sim_0_0_0")
      ml.ingroup <- drop.tip(ml.tree.pruned,"sim_0_0_0")
      results.raxml$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda")) - colless(as.treeshape(ml.ingroup,model="pda"))
    }
  } else if (results.raxml$taxa_ref[i] != "0_0_0"){
    if (results.raxml$noref[i] == "REF"){
      s <- as.integer(results.raxml$simulation[i])
      h <- as.integer(results.raxml$height[i])
      j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
      trees.both <- as.multiPhylo(c(raxml.trees[[i]],ml.tree[[j]]))
      
      ###RF distance to ML tree
      results.raxml$RF.Dist.ML[i] <- multiRF(trees.both,multi2di=TRUE)[[2]]
      results.raxml$weighted.rf[i] <- wRF.dist(raxml.trees[[i]],ml.tree[[j]],check.labels=FALSE)
       
      ## normalized gamma
      results.raxml$std.gamma[i] <- gammaStat(raxml.trees[[i]])[1] - gammaStat(ml.tree[[j]])[1]
      
      results.raxml$std.colless[i] <- colless(as.treeshape(raxml.trees[[i]], model="pda")) - colless(as.treeshape(ml.tree[[j]],model="pda"))
      
      results.raxml$std.sites[i] <- results.raxml$sites[i] / max(results.raxml$sites[results.raxml$simulation == s])
      
      ingroup <- drop.tip(raxml.trees[[i]],"sim_0_0_0")
      ml.ingroup <- drop.tip(ml.tree[[j]],"sim_0_0_0")
      results.raxml$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda")) - colless(as.treeshape(ml.ingroup,model="pda"))
    } else {
      s <- as.integer(results.raxml$simulation[i])
      h <- as.integer(results.raxml$height[i])
      j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
      ml.tree.pruned <- drop.tip(ml.tree[[j]],paste("sim_",results.raxml$taxa_ref[i]))
      trees.both <- as.multiPhylo(c(raxml.trees[[i]],ml.tree.pruned))
      
      ###RF distance to ML tree
      results.raxml$RF.Dist.ML[i]<- multiRF(trees.both,multi2di=TRUE)[[2]]
      results.raxml$weighted.rf[i] <- wRF.dist(raxml.trees[[i]],ml.tree.pruned,check.labels=FALSE)
                           
      ## normalized gamma
      results.raxml$std.gamma[i] <- gammaStat(raxml.trees[[i]])[1] - gammaStat(ml.tree.pruned)[1]
      
      results.raxml$std.colless[i] <- colless(as.treeshape(raxml.trees[[i]], model="pda")) - colless(as.treeshape(ml.tree.pruned,model="pda"))
      
      results.raxml$std.sites[i] <- results.raxml$sites[i] / max(results.raxml$sites[results.raxml$simulation == s])
      
      ingroup <- drop.tip(raxml.trees[[i]],"sim_0_0_0")
      ml.ingroup <- drop.tip(ml.tree.pruned,"sim_0_0_0")
      results.raxml$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda")) - colless(as.treeshape(ml.ingroup,model="pda"))
    }
  }
}

summary(results.raxml)
write.csv(results.raxml,file=paste("output/",args$output,"-raxml.csv",sep=""),quote=FALSE,row.names=TRUE,col.names=TRUE,na="NA")

#########################
####Clean up dataframe
# #########################
results.raxml$int <- as.factor(results.raxml$int)
results.raxml$noref <- as.factor(results.raxml$noref)
results.raxml$simulation <- as.factor(results.raxml$simulation)
results.raxml$missing <- as.factor(results.raxml$missing)
results.raxml$method <- as.factor(results.raxml$method)
# 

####As a start, visualize correlations between all variables
#psych::pairs.panels(results.raxml[results.raxml$maf %in% c(0,0.01,0.05),c("quality","missing","maf","int","noref","gamma","colless","sackin","tree.height","Avg.BLs","SD.BLs")],cex.cor=2)


####################BEGIN ANALYSIS ################################
####Calculate Matrix of RF distances for all trees with the same species tree in multi tree object
pdf(paste("output/",args$output,".pdf",sep=""))
for (h in unique(as.numeric(results.raxml$height))){
  for (i in unique(as.numeric(results.raxml$simulation))){
    
    j <- which(ml.tree.info$simulation == i & ml.tree.info$height == h)
    subset <- which(results.raxml$simulation == i & results.raxml$noref == "REF" & results.raxml$height == h)
    trees.subset <- c(raxml.trees[subset],ml.tree[[j]])
    
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
      geom_label_repel(label=c(paste("MAF",results.raxml$maf[subset],"MISS",results.raxml$missing[subset],sep=" "),"truth"),size=rel(1))+
      scale_color_manual(labels = c("EXT","INT","truth"),values=c("#009980", "#006699", "black"))+
      scale_shape_manual(labels = c("astral","raxml","ml"),values=c(1,16,17))
    print(plot)
    
    plot2 <- biplot +
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
      #geom_label_repel(label=c(paste("MISS",results.raxml[subset,3],"Q",results.raxml[subset,2],sep=","),"truth"),size=rel(1))+
      scale_color_manual(labels = c("EXT","INT","truth"),values=c("#009980", "#006699", "black"))+
      scale_shape_manual(labels = c("astral","raxml","ml"),values=c(1,16,17))
    print(plot2)
    
    heatmap(rf_matrix,labCol=FALSE,labRow=paste(results.raxml[subset,]$maf,results.raxml[subset,]$int,sep=","))
  }
}
dev.off()

#write.csv(results.raxml,file=paste("output/",args$output,"-raxml.csv",sep=""),quote=FALSE,row.names=TRUE,col.names=TRUE,na="NA")
