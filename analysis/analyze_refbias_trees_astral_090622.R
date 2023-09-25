############################################
## Script for analyzing refbias output trees from RAxML
## Written by J. Rick, 11 March 2019
## Updated 6 August 2021
## Updated 18 Feb 2022 to use here package
## updated 27 sept 2022 to work with astral trees
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
    library(TreeDist),
    library(ggdist),
    library(here)))
here::i_am("analysis/analyze_refbias_trees_astral_090622.R")
source(here("analysis","theme_custom.R"))

############################################
## Specifying file names ########
############################################

raxml.trees <- "092321-all-raxml.trees"
raxml.tree.names <- "092321-all-raxml.names"
astral.trees <- "090622-all-astral.trees"
astral.tree.names <- "090622-all-astral.names"
ml.trees <- "092321-s_tree.tree"
ml.tree.names <- "092321-s_tree.names"
output <- "092321-output"
refdist <- "092321-refdist.txt"
args <- mget(c("raxml.trees","raxml.tree.names","astral.trees","astral.tree.names",
               "ml.trees","ml.tree.names","refdist","output"))

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
raxml.trees<-read.tree(here("output","new",args$raxml.trees))
astral.trees<-read.tree(here("output","new",args$astral.trees))
ml.tree<-read.tree(here("output","new",args$ml.trees))

###Read in file of tree names
raxml.tree.names<-read.table(here("output","new",args$raxml.tree.names),stringsAsFactors=FALSE)[,1]
astral.tree.names<-read.table(here("output","new",args$astral.tree.names),stringsAsFactors = FALSE)[,1]
ml.tree.names <- read.table(here("output","new",args$ml.tree.names),stringsAsFactors = FALSE)[,1]

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

refdist <- read_table(here("output","new",args$refdist), col_names=c("simulation","height","int","taxa_ref","avg_dxy"), skip=1) %>%
  group_by(simulation,height,int) %>%
  mutate(simulation = as.integer(simulation)) %>%
  slice_tail(n=1) %>%
  ungroup()

gt.rf <- read_table(here("output","new","092321-LONG-EXT.TTR.gt_rf"),col_names=c("sim","gt_rf")) %>%
  mutate(height = "LONG",
         simulation = as.integer(gsub("sim","",sim))) %>%
  filter(!(simulation %in% c(23,24,25))) %>%
  add_row(read_table2(here("output","new","092321-MED-EXT.TTR.gt_rf"),col_names=c("sim","gt_rf")) %>%
            mutate(height = "MED",
                   simulation = as.integer(gsub("sim","",sim)))) %>% 
  add_row(read_table2(here("output","new","092321-SHORT-EXT.TTR.gt_rf"),col_names=c("sim","gt_rf")) %>%
            mutate(height = "SHORT",
                   simulation = as.integer(gsub("sim","",sim))) %>%
            filter(!(simulation %in% c(23,24,25)))) %>%
  add_row(read_table2(here("output","new","100922-SHORT-EXT.TTR.gt_rf"),col_names=c("sim","gt_rf")) %>%
            mutate(height = "SHORT",
                   simulation = as.integer(gsub("sim","",sim)))) %>%
  add_row(read_table2(here("output","new","100922-LONG-EXT.TTR.gt_rf"),col_names=c("sim","gt_rf")) %>%
            mutate(height = "LONG",
                   simulation = as.integer(gsub("sim","",sim)))) %>%
  drop_na() %>%
  group_by(height,sim) %>%
  slice_head(n=1) %>%
  ungroup()

ml.tree.info2 <- ml.tree.info %>%
  left_join(gt.rf,by=c("simulation","height")) %>%
  filter(simulation > 15 & simulation < 26)


#################
## plots of ml tree characteristics
#################
ml.tree.info2 %>%
  pivot_longer(cols = c(ingroup.tree.height,Avg.BLs,gt_rf,ingroup.gamma,ingroup.colless),names_to="params",values_to="value") %>%
  mutate(params = factor(params, levels=c("ingroup.tree.height","Avg.BLs","gt_rf","ingroup.colless","ingroup.gamma"))) %>%
  ggplot(aes(x=height,y=value,col=height)) +
  geom_boxplot(outlier.shape=NA, aes(fill=height),col="black",alpha=0.5) +
  geom_jitter(aes(fill=height),col="black",width=0.05,height=0,pch=21,size=2)+
  theme_custom() +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_blank()) +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_y_continuous(position="right") +
  scale_fill_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int,aesthetics = c("fill","col")) +
  facet_grid(vars(params),scales="free_y",
             labeller = as_labeller(c(ingroup.tree.height = "Tree Height",
                                      ingroup.gamma = "Ingroup\ngamma", 
                                      ingroup.colless = "Ingroup\nimbalance",
                                      gt_rf = "Gene tree\ndiscordance",
                                      Avg.BLs = "Average\nbranch lengths")))

plot1 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=tree.height,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Tree Height") +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  facet_wrap(~simulation)

plot2 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=Avg.BLs,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Average Branch Lengths") +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int)

plot3 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=gt_rf,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Gene Tree Discordance") +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int)

plot4 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=ingroup.tree.height,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Ingroup Tree Height") +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int)

plot5 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=ingroup.colless,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Ingroup Colless Imbalance") +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int)

plot6 <- ml.tree.info2 %>%
  ggplot(aes(x=height,y=ingroup.gamma,fill=NULL)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(col=height),width=0.1,height=0.1) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  ylab("Ingroup gamma") +
  scale_x_discrete(labels=c("Low ILS","Med ILS","High ILS")) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int)

ggarrange(plot1,plot2,plot4,plot3,plot5,plot6,nrow=2,ncol=3)

#####################
## Start of analysis!
#####################

## RAXML TREES FIRST
###Create empty data frame and with named, empty columns 

####################
## need to remove short & long runs from sim23-25, and replace them with re-done simulations
new.trees<-read.tree(here("output","new","100922-all-raxml.trees"))
new.tree.names<-read.table(here("output","new","100922-all-raxml.names"))

t23.short <- results.raxml$tree.num[results.raxml$simulation == 23 & results.raxml$height == "SHORT"]
t23.long <- results.raxml$tree.num[results.raxml$simulation == 23 & results.raxml$height == "LONG"]
t24.long <- results.raxml$tree.num[results.raxml$simulation == 24 & results.raxml$height == "LONG"]
t24.short <- results.raxml$tree.num[results.raxml$simulation == 24 & results.raxml$height == "SHORT"]
t25.short <- results.raxml$tree.num[results.raxml$simulation == 25 & results.raxml$height == "SHORT"]
t25.long <- results.raxml$tree.num[results.raxml$simulation == 25 & results.raxml$height == "LONG"]
remove <- c(t23.short,t23.long,t24.short,t24.long,t25.short,t25.long)

raxml.trees.new <- c(raxml.trees[-remove],new.trees)
raxml.names.new <- c(raxml.tree.names[-remove],new.tree.names$V1)

raxml.trees <- raxml.trees.new
raxml.tree.names <- raxml.names.new
######################

num.trees <- length(raxml.tree.names)
results.raxml <- tibble(tree.num=seq(1:length(raxml.trees)))

## pull info from tree names, and calculate tree statistics
results.raxml <- results.raxml %>%
  mutate(
    simulation = sapply(raxml.tree.names,function(x) as.integer(regmatches(x, regexec('-s([0-9]+)\\_q', x))[[1]][2])),
    height = sapply(raxml.tree.names,function(x) as.character(regmatches(x, regexec('([A-Z]+)-s[0-9]', x))[[1]][2])),
    method = "raxml",
    quality = sapply(raxml.tree.names,function(x) as.integer(regmatches(x, regexec('q(.*?)\\_m', x))[[1]][2])),
    missing = sapply(raxml.tree.names,function(x) as.numeric(regmatches(x, regexec('miss([0-9]?\\.?[0-9]+)\\_mac', x))[[1]][2])),
    maf = sapply(raxml.tree.names, function(x) as.numeric(regmatches(x, regexec('_mac(0?.?[0-9]+)_sites', x))[[1]][2])),
    sites = sapply(raxml.tree.names, function(x) as.integer(regmatches(x, regexec('sites([0-9]*?)\\.', x))[[1]][2])),
    taxa_ref = sapply(raxml.tree.names, function(x) regmatches(x,regexec('[A-Z]-([0-9]+_0_0)\\.phylip', x))[[1]][[2]]),
    int = sapply(raxml.tree.names, function(x) regmatches(x, regexec('REF.([A-Z]+)\\.', x))[[1]][2]),
    noref = "REF"
  ) %>%
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

raxml.node.desc <- results.raxml %>%
  mutate(avg.node.descendants = sapply(raxml.trees, function(x) mean(sapply(seq(length(x$tip.label),length(x$tip.label)+x$Nnode), function(y) length(Descendants(x,y,"tips")[[1]])))))

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
          file=here("output","new",paste(args$output,"-raxml.csv",sep="")),
          quote=FALSE,row.names=TRUE,na="NA")

###################################
## ASTRAL TREES
###################################

####################
## need to remove short & long runs from sim23-25, and replace them with re-done simulations
new.trees.astr<-read.tree(here("output","new","100922-all-astral.trees"))
new.tree.names.astr<-read.table(here("output","new","100922-all-astral.names"))

t23.short.astr <- which(sapply(astral.tree.names,function(x) grepl("sim23",x)) & sapply(astral.tree.names,function(x) grepl("SHORT",x)))
t23.long.astr <- which(sapply(astral.tree.names,function(x) grepl("sim23",x)) & sapply(astral.tree.names,function(x) grepl("LONG",x)))
t24.long.astr <- which(sapply(astral.tree.names,function(x) grepl("sim24",x)) & sapply(astral.tree.names,function(x) grepl("LONG",x)))
t24.short.astr <- which(sapply(astral.tree.names,function(x) grepl("sim24",x)) & sapply(astral.tree.names,function(x) grepl("SHORT",x)))
t25.short.astr <- which(sapply(astral.tree.names,function(x) grepl("sim25",x)) & sapply(astral.tree.names,function(x) grepl("SHORT",x)))
t25.long.astr <- which(sapply(astral.tree.names,function(x) grepl("sim25",x)) & sapply(astral.tree.names,function(x) grepl("LONG",x)))
remove.astr <- c(t23.short.astr,t23.long.astr,t24.short.astr,t24.long.astr,t25.short.astr,t25.long.astr)

astral.trees.new <- c(astral.trees[-remove.astr],new.trees.astr)
astral.names.new <- c(astral.tree.names[-remove.astr],new.tree.names.astr$V1)

astral.trees <- astral.trees.new
astral.tree.names <- astral.names.new
######################


###Create empty data frame and with named, empty columns 
num.trees <- length(astral.tree.names)
results.astral <- tibble(tree.num=seq(1:length(astral.trees)))

## pull info from tree names, and calculate tree statistics
results.astral <- results.astral %>%
  mutate(
    simulation = sapply(astral.tree.names,function(x) as.integer(regmatches(x, regexec('-sim([0-9]+)\\_m', x))[[1]][2])),
    height = sapply(astral.tree.names,function(x) as.character(regmatches(x, regexec('([A-Z]+)-sim[0-9]', x))[[1]][2])),
    method = "astral",
    quality = 40,
    missing = sapply(astral.tree.names,function(x) as.numeric(regmatches(x, regexec('_miss([0-9]?\\.?[0-9]+)\\_mac', x))[[1]][2])),
    maf = sapply(astral.tree.names, function(x) as.numeric(regmatches(x, regexec('_mac(0?.?[0-9]+)_[A-Z]', x))[[1]][2])),
    #sites = sapply(astral.tree.names, function(x) as.integer(regmatches(x, regexec('sites([0-9]*?)\\.', x))[[1]][2])),
    #taxa_ref = sapply(astral.tree.names, function(x) regmatches(x,regexec('[A-Z]-([0-9]+_0_0)\\.phylip', x))[[1]][[2]]),
    int = sapply(astral.tree.names, function(x) regmatches(x, regexec('_([A-Z]+)', x))[[1]][2]),
    noref = "REF"
  ) %>%
  left_join(select(results.raxml,c(simulation:height,missing:maf,taxa_ref,int)),by=c("simulation","height","missing","maf","int")) %>%
  left_join(refdist,by=c("simulation","height","int","taxa_ref"),copy=TRUE) %>%
  #left_join(select(gt.rf,gt_rf:simulation),by=c("simulation","height")) %>%
  mutate(
    #tree.height = sapply(astral.trees, function(x) max(branching.times(root(x,"sim_0_0_0",resolve.root=TRUE)),na.rm=T)),
    #ingroup.tree.height = sapply(astral.trees, function(x) max(branching.times(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0")),na.rm=T)),
    #Avg.BLs = sapply(astral.trees, function(x) mean(root(x,"sim_0_0_0",resolve.root=TRUE)$edge.length,na.rm=T)),
    #SD.BLs = sapply(astral.trees, function(x) sd(root(x,"sim_0_0_0",resolve.root=TRUE)$edge.length,na.rm=T)),
    mean.support = sapply(astral.trees, function(x) mean(as.numeric(root(x,"sim_0_0_0",resolve.root=TRUE)$node.label[-1]),na.rm=T)),
    sd.support = sapply(astral.trees, function(x) sd(as.numeric(root(x,"sim_0_0_0",resolve.root=TRUE)$node.label[-1]),na.rm=T))
  ) %>%
   mutate(
  #   gamma = sapply(astral.trees, function(x) gammaStat(chronopl(root(x,"sim_0_0_0",resolve.root=TRUE),lambda=1))),
  #   ingroup.gamma = sapply(astral.trees, function(x) gammaStat(drop.tip(chronopl(root(x,"sim_0_0_0",resolve.root=TRUE),lambda=1),"sim_0_0_0"))),
     colless = sapply(astral.trees, function(x) colless(as.treeshape(root(x,"sim_0_0_0",resolve.root=TRUE),model="pda"),norm="pda")),
     ingroup.colless = sapply(astral.trees, function(x) colless(as.treeshape(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0"),model="pda"),norm="pda")),
     sackin = sapply(astral.trees, function(x) sackin(as.treeshape(root(x,"sim_0_0_0",resolve.root=TRUE),model="pda"),norm="pda")),
     ingroup.sackin = sapply(astral.trees, function(x) sackin(as.treeshape(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0"),model="pda"),norm="pda"))
   ) 

astral.node.desc <- results.astral %>%
  mutate(avg.node.descendants = sapply(astral.trees, function(x) mean(sapply(seq(length(x$tip.label),length(x$tip.label)+x$Nnode), function(y) length(Descendants(x,y,"tips")[[1]])))))

###################
## statistics requiring separation by simulation number
## to compare to ML tree
###################
results.astral <- results.astral %>%
  mutate(RF.Dist.ML = numeric(nrow(results.astral)),
         Q.Dist.ML = numeric(nrow(results.astral)),
         CI.Dist.ML = numeric(nrow(results.astral)),
         std.gamma = numeric(nrow(results.astral)),
         std.colless = numeric(nrow(results.astral)),
         std.sackin = numeric(nrow(results.astral)),
         std.sites = numeric(nrow(results.astral)),
         std.ingroup.gamma = numeric(nrow(results.astral)),
         std.ingroup.colless = numeric(nrow(results.astral)),
         std.ingroup.sackin = numeric(nrow(results.astral)),
         ml.tree = integer(nrow(results.astral)))
for (i in 1:nrow(results.astral)){
  s <- as.integer(results.astral$simulation[i])
  h <- as.character(results.astral$height[i])
  j <- which(ml.tree.info$simulation == s & ml.tree.info$height == h)
  
  results.astral$ml.tree[i] <- j
  
  ml.tree.pruned <- drop.tip(ml.tree[[j]],results.astral$taxa_ref[i])
  
  #if(results.astral$taxa_ref[i] %in% astral.trees[[i]]$tip.label){
    trees.both <- c(astral.trees[[i]],ml.tree[[j]])
  #} else {
  #  trees.both <- c(astral.trees[[i]],ml.tree.pruned)
  #}
  
  if(results.astral$taxa_ref[i] %in% astral.trees[[i]]$tip.label){
    trees.both.root <- c(root(astral.trees[[i]],"sim_0_0_0",resolve.root=TRUE),
                         root(ml.tree[[j]],"sim_0_0_0",resolve.root=TRUE))
  } else {
    trees.both.root <- c(root(astral.trees[[i]],"sim_0_0_0",resolve.root=TRUE),
                         root(ml.tree.pruned,"sim_0_0_0",resolve.root=TRUE))
  }
  
  ##Calculating distances to ML tree
  results.astral$RF.Dist.ML[i]<- TreeDist::InfoRobinsonFoulds(trees.both,normalize=TRUE)[1]
  #results.astral$weighted.rf[i] <- wRF.dist(astral.trees[[i]],ml.tree[[j]]) # this always seems to come out wonky
  results.astral$Q.Dist.ML[i] <- Quartet::QuartetDivergence(Quartet::QuartetStatus(trees.both), similarity=FALSE)[-1]/(2/3) #2/3 is max Q
  results.astral$CI.Dist.ML[i] <- as.matrix(ClusteringInfoDistance(trees.both,normalize = TRUE))[-1,1] #normalized relative to max CI for the tree
  
  ## normalized gamma and imbalance
  #results.astral$std.gamma[i] <- gammaStat(chronos(trees.both.root[[1]],lambda=1,quiet=TRUE))[1] - gammaStat(trees.both.root[[2]])[1]
  results.astral$std.colless[i] <- colless(as.treeshape(trees.both.root[[1]], model="pda"),norm="pda")[1] - colless(as.treeshape(trees.both.root[[2]],model="pda"),norm="pda")
  results.astral$std.sackin[i] <- sackin(as.treeshape(trees.both.root[[1]], model="pda"),norm="pda") - sackin(as.treeshape(trees.both.root[[2]],model="pda"),norm="pda")
  
  #results.astral$std.sites[i] <- results.astral$sites[i] / max(results.astral$sites[results.astral$simulation == s])
  
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
  
  results.astral$std.ingroup.colless[i] <- colless(as.treeshape(ingroup, model="pda"),norm="pda") - colless(as.treeshape(ml.ingroup,model="pda"),norm="pda")
  #results.astral$std.ingroup.gamma[i] <- gammaStat(chronos(ingroup,lambda=1))[1] - gammaStat(chronos(ml.ingroup,lambda=1))[1]
  results.astral$std.ingroup.sackin[i] <- sackin(as.treeshape(ingroup, model="pda"),norm="pda") - sackin(as.treeshape(ml.ingroup,model="pda"),norm="pda")
}

summary(results.astral)  
write.csv(results.astral %>% filter(RF.Dist.ML < 0.9),
          file=here("output","new",paste(args$output,"-astral.csv",sep="")),
          quote=FALSE,row.names=TRUE,na="NA")

######################
## univariate plots
#####################
results.all <- results.raxml %>%
  add_row(results.astral %>% select(-ml.tree)) %>%
  filter(RF.Dist.ML < 0.9)
results.all$maf <- as.factor(results.all$maf)
plot1 <- ggplot(data = results.all, 
                aes(x=height,
                    y=tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot2 <- ggplot(data = results.all, 
                aes(x=height,
                    y=Avg.BLs,
                    fill=NULL)) +
  geom_boxplot(aes(alpha=0.9),outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot3 <- ggplot(data = results.all, 
                aes(x=height,
                    y=ingroup.tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot4 <- ggplot(data = results.all, 
                aes(x=height,
                    y=ingroup.colless,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~method)

plot5 <- ggplot(data = results.all, 
                aes(x=height,
                    y=std.ingroup.gamma,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank())

plot6 <- ggplot(data = results.all, 
                aes(x=height,
                    y=RF.Dist.ML,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~method)

plot7 <- ggplot(data = results.all, 
                aes(x=height,
                    y=Q.Dist.ML,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.1,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~method)

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
#results.raxml <- read.csv(file=here("output","new",paste0(output,"-raxml.csv")),row.names=1,na="NA",header=TRUE)
pdf(here("output","new",paste(args$output,".pdf",sep="")))
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
results.raxml <- results.raxml[results.raxml$simulation > 15 & results.raxml$simulation < 26 & results.raxml$RF.Dist.ML < 0.9,]
raxml.trees2 <- raxml.trees[results.raxml$tree.num]
param_loadings <- tibble(sim=integer(),height=character(),method=character(),param=character(),
                         pct_var_exp=numeric(),
                         x_corr=numeric(),y_corr=numeric(),
                         x_corr_sig=numeric(),y_corr_sig=numeric())
for (m in c("raxml","astral")){
  for (h in unique(as.character(results.raxml$height))){
  for (i in unique(as.numeric(as.character(results.astral$simulation)))){
    #for (i in seq(16,30,by=1)){
    j <- which(ml.tree.info$simulation == i & ml.tree.info$height == h)
    
    if (m == "raxml") {
      subset_r <- results.raxml$tree.num[results.raxml$simulation == i & results.raxml$height == h & results.raxml$noref == "REF"]
      subset_a <- NULL
      
      if (length(subset_r) !=0) {
        trees.subset <- c(raxml.trees[subset_r],ml.tree[[j]])
      } else {
        print(paste("no raxml trees for sim",i," for height ",h))
        next
      }
    }
    if (m == "astral") {
      subset_a <- results.astral$tree.num[results.astral$simulation == i & results.astral$height == h]
      subset_r <- NULL
      if (length(subset_a) !=0) {
        trees.subset <- c(astral.trees[subset_a],ml.tree[[j]])
      } else {
        print(paste("no astral trees for sim",i," for height ",h))
        next
      }
    }
    
    #trees.subset <- sapply(trees.subset,function(x) root(x,outgroup,resolve.root=TRUE))
    #trees.subset.root <- root(trees.subset,outgroup,resolve.root=TRUE)
    class(trees.subset) <- "multiPhylo"
    #rf.dist <- as.matrix(Quartet::QuartetDivergence(Quartet::QuartetStatus(trees.subset), similarity=FALSE)[-1]/(2/3))
    #rf.dist <- as.matrix(ClusteringInfoDist(trees.subset,normalize = TRUE))
    rf.dist <- InfoRobinsonFoulds(trees.subset,normalize=TRUE)
    rf.pcoa <- pcoa(rf.dist)
    rf.pcoa.plot <- rf.pcoa$vectors %>%
      as_tibble() %>%
      ggplot(aes(x=Axis.1,y=Axis.2)) +
      geom_jitter(aes(col=factor(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a],"true"),levels=c("0","1","2","3","4","5","10","true")),
                      fill=factor(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a],"true"),levels=c("0","1","2","3","4","5","10","true")),
                     #text=results.raxml$missing[results.raxml$simulation == s],
                     #shape=factor(c(results.raxml$missing[subset])),
                     #size=factor(c(results.raxml$method[subset_r],results.astral$method[results.astral$tree.num %in% subset_a])),
                     shape=factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a],"true")),
                     ),
                  width=0.03,height=0.03,cex=5) +
      theme_bw(base_size=12, base_family="Arial") +
      #scale_size_manual(values=c(3,6,6)) +
      scale_color_viridis_d(direction=-1) +
      scale_fill_viridis_d(direction=-1,alpha=0.8) +
      scale_shape_manual(values=c(21,23,25)) +
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
        xend = cor(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.2,
        yend = cor(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.2,
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1,
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) +
      # annotate(geom = "text",
      #          x = cor(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.15-0.001,
      #          y = cor(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.15-0.002,
      #          label = "MAF",
      #          hjust=1,vjust=1
      # ) +
      # biplot arrows for Missing
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.2,
        yend = cor(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.2,
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1,
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) +
      # annotate(geom = "text",
      #          x = cor(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.15+0.001,
      #          y = cor(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.15-0.002,
      #          label = "Missing",
      #          hjust=0,vjust=0
      # ) +
      # biplot arrows for INT
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.2,
        yend = cor(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.2,
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1,
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) #+
      # annotate(geom = "text",
      #          x = cor(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.15+0.001,
      #          y = cor(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]) * 0.5 * sqrt(length(trees.subset) - 2)*0.15+0.002,
      #          label = "INT",
      #          hjust=0,vjust=1
      # )
    
    #print(rf.pcoa.plot)
    
    #plotly::ggplotly(rf.pcoa.plot)
    vectors <- tibble(sim = i, height = h,method=m,
                      param = c("maf","missing","int"),
                      pct_var_exp = sum(rf.pcoa$values$Relative_eig[1:2]),
                      x_corr = c(cor(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]),
                                 cor(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1]),
                                 cor(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1])
                      ),
                      y_corr = c(cor(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]),
                                cor(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1]),
                                cor(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1])
                      ),
                      x_corr_sig = c(cor.test(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1])$p.value,
                                     cor.test(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1])$p.value,
                                     cor.test(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)-1])$p.value
                      ),
                      y_corr_sig = c(cor.test(as.numeric(as.character(c(results.raxml$maf[results.raxml$tree.num %in% subset_r],results.astral$maf[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1])$p.value,
                                     cor.test(as.numeric(as.character(c(results.raxml$missing[results.raxml$tree.num %in% subset_r],results.astral$missing[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1])$p.value,
                                     cor.test(as.numeric(as.factor(c(results.raxml$int[results.raxml$tree.num %in% subset_r],results.astral$int[results.astral$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)-1])$p.value
                      ))
    param_loadings <- param_loadings %>%
      add_row(vectors)
  }
  }
}

dev.off()

## proportion of significant correlations
param_loadings %>%
  filter(sim < 26) %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE)) %>%
  group_by(param, method) %>% 
  summarize(n_PC1 = sum(sig_x),
            n_PC2 = sum(sig_y),
            mean_PC1 = mean(abs(x_corr)),
            mean_PC2 = mean(abs(y_corr)))


## plot of parameter correlations for all sims
# param_loadings %>% 
#   filter(sim < 26) %>% 
#   mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
#          height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
#          param2 = recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
#   pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
#   ggplot(aes(x=value)) +
#   geom_density(data = . %>% filter(corr == "abs_x_corr"),fill="gray80",adjust=0.75) +
#   geom_density(data = . %>% filter(corr == "abs_y_corr"),aes(y=-..density..),fill="gray30",adjust=0.75) +
#   facet_wrap(~param2) +
#   # ggscatter(x="abs_x_corr",y="abs_y_corr",
#   #           color="param",fill="param",
#   #           facet.by=c("height2","param2"),
#   #           size=3,alpha=0.75) + 
#   #xlab("PCoA 1 Parameter Correlation") + 
#   #ylab("PCoA 2 Parameter Correlation") +
#   theme_custom() +
#   theme(legend.position = "none",
#         strip.text = element_text(size=rel(1.3)),
#         axis.title = element_text(size=rel(1.5)),
#         axis.text = element_text(size=rel(1.5)))
  
param_loadings_raxml <- param_loadings %>% mutate(method="raxml")
param_loadings_raxml2 <- param_loadings %>% mutate(method="raxml_ci_dist")
param_loadings_astral <- param_loadings %>% mutate(method="astral")
param_loadings_astral2 <- param_loadings %>% mutate(method="astral_ci_dist")
param_loadings %>% 
  #add_row(param_loadings_raxml2) %>%
  #add_row(param_loadings_astral) %>%
  #add_row(param_loadings_astral2) %>%
  filter(sim < 26) %>% 
  mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
         height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
         param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
  #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_cat = case_when(sig_x & sig_y ~ "both",
                             sig_x & !sig_y ~ "PC1",
                             sig_y & !sig_x ~ "PC2",
                             TRUE ~ "neither")) %>%
  ggplot(aes(x=abs_x_corr,y=abs_y_corr)) +
  geom_point(aes(col=sig_cat,fill=sig_cat),size=6,pch=21) + 
  facet_grid(cols=vars(param2),rows=vars(method)) +
  # ggscatter(x="abs_x_corr",y="abs_y_corr",
  #           color="param",fill="param",
  #           facet.by=c("height2","param2"),
  #           size=3,alpha=0.75) + 
  #xlab("PCoA 1 Parameter Correlation") + 
  #ylab("PCoA 2 Parameter Correlation") +
  scale_color_manual(values=PNWColors::pnw_palette("Sunset2",4)) +
  scale_fill_manual(values=scales::alpha(PNWColors::pnw_palette("Sunset2",4),0.8)) +
  theme_custom() +
  theme(legend.position = "none",
        strip.text = element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.5)),
        axis.text = element_text(size=rel(1.2))) +
  xlab("PCoA Axis 1 Correlation") +
  ylab("PCoA Axis 2 Correlation")


# param_loadings_sim %>% 
#   filter(sim < 26) %>% 
#   mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
#          height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
#          param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
#   #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
#   mutate(dom = case_when(x_corr > y_corr ~ "PC1",
#                          y_corr > x_corr ~ "PC2"),
#          sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
#                            x_corr_sig >= 0.01 ~ FALSE,
#                            TRUE ~ FALSE),
#          sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
#                            y_corr_sig >= 0.01 ~ FALSE,
#                            TRUE ~ FALSE),
#          sig_cat = case_when(sig_x & sig_y ~ "both",
#                              sig_x & !sig_y ~ "PC1",
#                              sig_y & !sig_x ~ "PC2",
#                              TRUE ~ "neither")) %>%
#   ggplot() +
#   geom_histogram(aes(x=abs_x_corr,y=..density..),binwidth=0.05) +
#   geom_histogram(aes(x=abs_y_corr,y=-..density..),binwidth=0.05) +
#   facet_wrap(~param2) +
#   theme_custom() +
#   geom_hline(yintercept=0)
# 
# param_loadings_sim %>% 
#   filter(sim < 26) %>% 
#   mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
#          height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
#          param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
#   #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
#   mutate(dom = case_when(x_corr > y_corr ~ "PC1",
#                          y_corr > x_corr ~ "PC2"),
#          sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
#                            x_corr_sig >= 0.01 ~ FALSE,
#                            TRUE ~ FALSE),
#          sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
#                            y_corr_sig >= 0.01 ~ FALSE,
#                            TRUE ~ FALSE),
#          sig_cat = case_when(sig_x & sig_y ~ "both",
#                              sig_x & !sig_y ~ "PC1",
#                              sig_y & !sig_x ~ "PC2",
#                              TRUE ~ "neither")) %>%
#   pivot_longer(cols=starts_with("abs"),names_to="corr_type",values_to="corr_value") %>%
#   mutate(corr_sig = case_when(corr_type == "abs_x_corr" ~ sig_x,
#                               corr_type == "abs_y_corr" ~ sig_y,
#                               TRUE ~ FALSE)) %>%
#   ggplot(aes(x=corr_value,y=param)) +
#   ggdist::stat_dots(aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=1.1) +
#   theme_custom() +
#   theme(panel.grid.major=element_line(),
#         legend.position=c(0.8,0.8)) +
#   xlab("PCoA Correlation Coefficient") +
#   ylab("Parameter") +
#   scale_fill_manual(values=c("gray50","black"),aesthetics=c("fill","color")) +
#   coord_cartesian(ylim=c(0.9,4.1),xlim=c(-0.01,1.01),expand=FALSE)
  
## PLOT OF ALL LOADINGS TOGETHER
## (other dataframes come from analyze_refbias_trees_emp scripts)
param_loadings.all <- param_loadings.cichlids %>% 
  add_row(param_loadings.lates) %>%
  add_row(param_loadings) 

write.csv()
  
#param_loadings_raxml %>%
## ORIGINAL FIG. 4B ################
param_loadings.all %>%
  #add_row(param_loadings_raxml) %>%
  filter(method == "raxml") %>%
  mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
         #height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
         param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
  #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_cat = case_when(sig_x & sig_y ~ "both",
                             sig_x & !sig_y ~ "PC1",
                             sig_y & !sig_x ~ "PC2",
                             TRUE ~ "neither")) %>%
  pivot_longer(cols=starts_with("abs"),names_to="corr_type",values_to="corr_value") %>%
  mutate(corr_sig = case_when(corr_type == "abs_x_corr" ~ sig_x,
                              corr_type == "abs_y_corr" ~ sig_y,
                              TRUE ~ FALSE)) %>%
  ggplot(aes(x=corr_value,y=height)) +
  #ggdist::stat_dots(data=. %>% filter(height != "lates" & height != "cichlids"),aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=0.6,binwidth=0.01) +
  #ggdist::stat_dots(data=. %>% filter(height == "lates" | height == "cichlids"),aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=0.6,binwidth=0.01,side="bottom") +
  ggdist::geom_dots(aes(shape=corr_sig,col=corr_type, fill=corr_type, group=param),alpha=0.8,scale=1.1,binwidth=0.04,dotsize=1.2,layout="bin") +
  stat_pointinterval(position=position_nudge(y=-0.1),point_interval="median_qi") +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1.0),labels=c("0","0.25","0.50","0.75","1.0")) +
  theme_custom() +
  theme(panel.grid.major=element_line(),
        legend.position="right",
        strip.text = element_text(size=rel(1.5))) +
  xlab("PCoA Correlation Coefficient") +
  ylab("") +
  facet_wrap(~param2) +
  scale_shape_manual(values=c(4,19)) +
  scale_fill_manual(values=c("gray60","#000000"),aesthetics=c("fill","color")) +
  geom_hline(yintercept=2.8)

## NEW FIG. 4B ################
plot4b <- param_loadings.all %>%
  #add_row(param_loadings_raxml) %>%
  filter(method == "raxml") %>%
  mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
         #height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
         param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
  #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_cat = case_when(sig_x & sig_y ~ "both",
                             sig_x & !sig_y ~ "PC1",
                             sig_y & !sig_x ~ "PC2",
                             TRUE ~ "neither")) %>%
  pivot_longer(cols=starts_with("abs"),names_to="corr_type",values_to="corr_value") %>%
  mutate(corr_sig = case_when(corr_type == "abs_x_corr" ~ sig_x,
                              corr_type == "abs_y_corr" ~ sig_y,
                              TRUE ~ FALSE)) %>%
  ggplot(aes(x=corr_value,y=height)) +
  #ggdist::stat_dots(data=. %>% filter(height != "lates" & height != "cichlids"),aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=0.6,binwidth=0.01) +
  #ggdist::stat_dots(data=. %>% filter(height == "lates" | height == "cichlids"),aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=0.6,binwidth=0.01,side="bottom") +
  ggdist::stat_slab(alpha=0.5,
                    normalize="all",trim=FALSE,
                    adjust=1.2) +
  stat_pointinterval(data = . %>% filter(corr_type=="abs_x_corr"), 
                     aes(col=corr_type, fill=corr_type), position=position_nudge(y=-0.1),
                     point_interval="median_qi") +
  stat_pointinterval(data = . %>% filter(corr_type=="abs_y_corr"), 
                     aes(col=corr_type, fill=corr_type), position=position_nudge(y=-0.2),
                     point_interval="median_qi") +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1.0),labels=c("0","0.25","0.50","0.75","1.0")) +
  coord_cartesian(xlim=c(0,1),expand=FALSE) +
  theme_custom() +
  theme(panel.grid.major=element_line(),
        legend.position="right",
        strip.text = element_text(size=rel(1.5))) +
  xlab("PCoA Correlation Coefficient") +
  ylab("") +
  facet_grid(cols=vars(param2)) +
  scale_shape_manual(values=c(4,19)) +
  scale_fill_manual(values=c("gray60","#000000","gray80"),aesthetics=c("fill","color")) +
  geom_hline(yintercept=2.6)

param_loadings.all %>%
  #add_row(param_loadings_raxml) %>%
  filter(method == "astral") %>%
  mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
         #height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
         param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
  #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_cat = case_when(sig_x & sig_y ~ "both",
                             sig_x & !sig_y ~ "PC1",
                             sig_y & !sig_x ~ "PC2",
                             TRUE ~ "neither")) %>%
  pivot_longer(cols=starts_with("abs"),names_to="corr_type",values_to="corr_value") %>%
  mutate(corr_sig = case_when(corr_type == "abs_x_corr" ~ sig_x,
                              corr_type == "abs_y_corr" ~ sig_y,
                              TRUE ~ FALSE)) %>%
  ggplot(aes(x=corr_value,y=height)) +
  #ggdist::stat_dots(data=. %>% filter(height != "lates" & height != "cichlids"),aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=0.6,binwidth=0.01) +
  #ggdist::stat_dots(data=. %>% filter(height == "lates" | height == "cichlids"),aes(fill=corr_sig,col=corr_sig),alpha=0.8,scale=0.6,binwidth=0.01,side="bottom") +
  ggdist::geom_dots(aes(shape=corr_sig,col=corr_type, fill=corr_type, group=param),alpha=0.8,scale=1.1,binwidth=0.04,dotsize=1.2,layout="bin") +
  stat_pointinterval(position=position_nudge(y=-0.1),point_interval="median_qi") +
  scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1.0),labels=c("0","0.25","0.50","0.75","1.0")) +
  theme_custom() +
  theme(panel.grid.major=element_line(),
        legend.position="right",
        strip.text = element_text(size=rel(1.5))) +
  xlab("PCoA Correlation Coefficient") +
  ylab("") +
  facet_wrap(~param2) +
  scale_shape_manual(values=c(4,19)) +
  scale_fill_manual(values=c("gray60","#000000"),aesthetics=c("fill","color")) +
  geom_hline(yintercept=2.8)

## another try at Fig 4b #######
param_loadings.all %>%
  #add_row(param_loadings_raxml) %>%
  filter(method == "raxml") %>%
  mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
         #height2 = recode_factor(.$height,LONG="Low ILS",MED="Medium ILS",SHORT="High ILS",.ordered=TRUE),
         param2 = dplyr::recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
  #pivot_longer(cols=starts_with("abs"),names_to="corr") %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_cat = factor(case_when(sig_x & sig_y ~ "both",
                             sig_x & !sig_y ~ "PC1",
                             sig_y & !sig_x ~ "PC2",
                             TRUE ~ "neither"),
                          levels=c("neither","PC2","PC1","both"))) %>%
  pivot_longer(cols=starts_with("abs"),names_to="corr_type",values_to="corr_value") %>%
  mutate(corr_sig = case_when(corr_type == "abs_x_corr" ~ sig_x,
                              corr_type == "abs_y_corr" ~ sig_y,
                              TRUE ~ FALSE)) %>%
  #group_by(height,method,param,sig_cat) %>%
  #summarize(n=n()) %>%
  ggplot(aes(x=height,fill=sig_cat)) +
  geom_bar(position="fill") +
  facet_grid(rows=vars(param2)) +
  theme_custom() +
  ylab("Proportion") +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size=rel(1.5))) +
  #coord_flip() +
  scale_fill_manual(values=c("gray80",viridis(5,option="mako")[2:5]))
  