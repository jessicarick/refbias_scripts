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
    library(TreeDist),
    library(Quartet)))
source("analysis/theme_custom.R")

############################################
## Specifying file names ########
############################################

raxml.trees <- "082021-lates-emp-batch.trees"
raxml.tree.names <- "082021-lates-emp-tree.names"
output <- "082021-lates-emp-output"
refdist <- "082021-emp-refdist.txt"
sites <- "082021-SNPs-emp"
args <- mget(c("raxml.trees","raxml.tree.names","refdist","sites","output"))

####Read in concatenated tree file 
raxml.trees<-read.tree(paste("output/new/",args$raxml.trees,sep=""))

###Read in file of tree names
raxml.tree.names<-read.table(paste("output/new/",args$raxml.tree.names,sep=""),stringsAsFactors=FALSE)[,1]

###Import reference distance data
# refdist <- read_table2(paste("output/new/",args$refdist,sep=""), col_names=c("simulation","height","int","taxa_ref","avg_dxy")) %>%
#   group_by(simulation,height,int,taxa_ref) %>%
#   slice(1) %>%
#   ungroup()

sites <- read_csv(paste0("output/new/",args$sites),col_names=c("tree","snps")) %>% 
  mutate(simulation=as.integer(gsub(".*s([0-9]+)_q40.*","\\1",tree)),
         missing=as.numeric(gsub(".*miss([0-9]\\.?[0-9]*)\\_.*","\\1",tree)),
         maf=as.numeric(gsub(".*_maf([0-9]\\.?[0-9]*)_.*","\\1",tree)),
         int=gsub(".*_([A-Z]+)\\.noInv.*","\\1",tree))

#####################
## Start of analysis!
#####################

###Create empty data frame and with named, empty columns 
num.trees <- length(raxml.tree.names)
results.raxml <- tibble(tree.num=seq(1:length(raxml.trees)))

###Fix names in RAxML trees
fix_labels <- function(x) {
  labs <- x$tip.label
  new_labs <- gsub(".*aln_(.*)\\.sorted\\.bam","\\1",labs)
  x$tip.label <- new_labs
  return(x)
}

raxml.trees <- lapply(raxml.trees, function(x) fix_labels(x))

## pull info from tree names, and calculate tree statistics
results.raxml <- results.raxml %>%
  mutate(
    simulation = sapply(raxml.tree.names,function(x) as.integer(regmatches(x, regexec('_s([0-9]+)\\_q', x))[[1]][2])),
    method = "raxml",
    quality = sapply(raxml.tree.names,function(x) as.integer(regmatches(x, regexec('q(.*?)\\_m', x))[[1]][2])),
    missing = sapply(raxml.tree.names,function(x) as.numeric(regmatches(x, regexec('miss(0.*?)\\_maf', x))[[1]][2])),
    maf = sapply(raxml.tree.names, function(x) as.numeric(regmatches(x, regexec('maf(0.*?)\\.REF', x))[[1]][2])),
    int = sapply(raxml.tree.names, function(x) regmatches(x, regexec('REF.([A-Z]+)\\.emp', x))[[1]][2]),
    noref = "REF"
  ) %>%
  mutate(
    tree.height = sapply(raxml.trees, function(x) max(branching.times(root(x,"SRR3140997",resolve.root=TRUE)),na.rm=T)),
    ingroup.tree.height = sapply(raxml.trees, function(x) max(branching.times(drop.tip(root(x,"SRR3140997",resolve.root=TRUE),"SRR3140997")),na.rm=T)),
    Avg.BLs = sapply(raxml.trees, function(x) mean(root(x,"SRR3140997",resolve.root=TRUE)$edge.length,na.rm=T)),
    SD.BLs = sapply(raxml.trees, function(x) sd(root(x,"SRR3140997",resolve.root=TRUE)$edge.length,na.rm=T)),
    mean.support = sapply(raxml.trees, function(x) mean(as.numeric(root(x,"SRR3140997",resolve.root=TRUE)$node.label[-1]),na.rm=T)),
    sd.support = sapply(raxml.trees, function(x) sd(as.numeric(root(x,"SRR3140997",resolve.root=TRUE)$node.label[-1]),na.rm=T))
  ) %>%
  mutate(
    gamma = sapply(raxml.trees, function(x) gammaStat(chronopl(root(x,"SRR3140997",resolve.root=TRUE),lambda=1))),
    ingroup.gamma = sapply(raxml.trees, function(x) gammaStat(drop.tip(chronopl(root(x,"SRR3140997",resolve.root=TRUE),lambda=1),"SRR3140997"))),
    colless = sapply(raxml.trees, function(x) colless(as.treeshape(root(x,"SRR3140997",resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.colless = sapply(raxml.trees, function(x) colless(as.treeshape(drop.tip(root(x,"SRR3140997",resolve.root=TRUE),"SRR3140997"),model="pda"),norm="pda")),
    sackin = sapply(raxml.trees, function(x) sackin(as.treeshape(root(x,"SRR3140997",resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.sackin = sapply(raxml.trees, function(x) sackin(as.treeshape(drop.tip(root(x,"SRR3140997",resolve.root=TRUE),"SRR3140997"),model="pda"),norm="pda"))
  ) 

results.raxml <- results.raxml %>%
  left_join(sites)

write.csv(results.raxml,file=paste("output/new/",args$output,"-raxml.csv",sep=""),quote=FALSE,row.names=TRUE,na="NA")

######################
## univariate plots
#####################
results.raxml$maf <- as.factor(results.raxml$maf)
results.raxml$missing <- as.factor(results.raxml$missing)
  
results.raxml$simulation <- as.factor(results.raxml$simulation)

plot1 <- ggplot(data = results.raxml, 
                aes(x=maf,
                    y=tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +
  theme_custom() +
  ylim(0,0.5) +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot2 <- ggplot(data = results.raxml, 
                aes(x=maf,
                    y=Avg.BLs,
                    fill=NULL)) +
  geom_boxplot(aes(alpha=0.9),outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot3 <- ggplot(data = results.raxml, 
                aes(x=maf,
                    y=ingroup.tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot4 <- ggplot(data = results.raxml, 
                aes(x=maf,
                    y=colless,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot5 <- ggplot(data = results.raxml, 
                aes(x=maf,
                    y=gamma,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

ggarrange(plot1,plot2,plot3,plot4,plot5,ncol=1) 

results.raxml %>% 
  ggplot(aes(x=as.numeric(as.character(maf)),y=tree.height,shape=int)) + 
    geom_point(aes(col=as.factor(missing))) + 
    geom_smooth(aes(group=missing,col=missing),span=2) + 
    theme_custom()

results.raxml %>% 
  ggplot(aes(x=as.factor(as.character(maf)),y=ingroup.colless,shape=int)) + 
  geom_boxplot(aes(fill=as.factor(int))) +
  #geom_jitter(aes(col=as.factor(missing))) + 
  #geom_smooth(aes(group=missing,col=missing),span=2) + 
  theme_custom()

## Plot trees in PCoA space by simulation
pdf("lates_emp_trees_pcoa.pdf")
for (s in 1:10){
  trees.subset <- raxml.trees[results.raxml$simulation == s]
  class(trees.subset) <- "multiPhylo"
  trees.subset <- root(trees.subset,"SRR3140997",resolve.root=TRUE)
  rf.dist <- multiRF(trees.subset)
  rf.pcoa <- pcoa(rf.dist)
  rf.pcoa.plot <- rf.pcoa$vectors %>%
    as_tibble() %>%
    ggplot(aes(x=Axis.1,y=Axis.2)) +
    geom_point(aes(col=factor(results.raxml$maf[results.raxml$simulation == s]),
                              #text=results.raxml$missing[results.raxml$simulation == s],
                   shape=factor(results.raxml$missing[results.raxml$simulation == s]),
                   size=factor(results.raxml$int[results.raxml$simulation == s])),alpha=0.9) +
    theme_bw(base_size=12, base_family="Open Sans Light") +
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", color=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="Open Sans"),
      axis.text = element_text(size=14),
      legend.text = element_text(size=14),
      legend.title = element_blank()
    ) +
    ggtitle(paste0("Lates empirical PCoA, Simulation ",s)) +
    geom_hline(yintercept=0,lty=2,col="gray50") +
    geom_vline(xintercept=0,lty=2,col="gray50") +
    # biplot arrows for MAF
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.character(results.raxml$maf[results.raxml$simulation == s])), rf.pcoa$vectors[,1]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
      yend = cor(as.numeric(as.character(results.raxml$maf[results.raxml$simulation == s])), rf.pcoa$vectors[,2]) * 0.8 * sqrt(nrow(results.raxml) - 1),
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1, 
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) + 
    annotate(geom = "text",
             x = cor(as.numeric(as.character(results.raxml$maf[results.raxml$simulation == s])), rf.pcoa$vectors[,1]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
             y = cor(as.numeric(as.character(results.raxml$maf[results.raxml$simulation == s])), rf.pcoa$vectors[,2]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
             label = "MAF",
             hjust=0,vjust=0
    ) +
    # biplot arrows for Missing
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.character(results.raxml$missing[results.raxml$simulation == s])), rf.pcoa$vectors[,1]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
      yend = cor(as.numeric(as.character(results.raxml$missing[results.raxml$simulation == s])), rf.pcoa$vectors[,2]) * 0.8 * sqrt(nrow(results.raxml) - 1),
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1, 
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) + 
    annotate(geom = "text",
             x = cor(as.numeric(as.character(results.raxml$missing[results.raxml$simulation == s])), rf.pcoa$vectors[,1]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
             y = cor(as.numeric(as.character(results.raxml$missing[results.raxml$simulation == s])), rf.pcoa$vectors[,2]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
             label = "Missing",
             hjust=1,vjust=1
    ) +
    # biplot arrows for INT
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.factor(results.raxml$int[results.raxml$simulation == s])), rf.pcoa$vectors[,1]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
      yend = cor(as.numeric(as.factor(results.raxml$int[results.raxml$simulation == s])), rf.pcoa$vectors[,2]) * 0.8 * sqrt(nrow(results.raxml) - 1),
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1, 
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) + 
    annotate(geom = "text",
             x = cor(as.numeric(as.factor(results.raxml$int[results.raxml$simulation == s])), rf.pcoa$vectors[,1]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
             y = cor(as.numeric(as.factor(results.raxml$int[results.raxml$simulation == s])), rf.pcoa$vectors[,2]) * 0.8 * sqrt(nrow(results.raxml) - 1), 
             label = "INT",
             hjust=0,vjust=0
    )
    
  
  print(rf.pcoa.plot)
  #plotly::ggplotly(rf.pcoa.plot)
}
dev.off()

rf.pcoa$vectors %>%
  as_tibble() %>%
  plot_ly(x = ~Axis.1, y = ~Axis.2, z = ~Axis.3, 
  color = factor(results.raxml$int[results.raxml$simulation == s]),
                 hoverinfo = "text",
                 hovertext = paste("MAF :", results.raxml$maf[results.raxml$simulation == s],
                                   "<br> Missing :", results.raxml$missing[results.raxml$simulation == s],
                                   "<br> INT :", results.raxml$int[results.raxml$simulation == s])
) %>%
  add_markers() %>%
  layout(
    scene = list(xaxis = list(title = 'Axis 1'),
                 yaxis = list(title = 'Axis 2'),
                 zaxis = list(title = 'Axis 3'))
  )
