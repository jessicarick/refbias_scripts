############################################
## Script for analyzing refbias output trees from RAxML
## Written by J. Rick, 11 March 2019
## Updated 6 August 2021
## Updated 18 Feb 2022 to use here package
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
    library(Quartet),
    library(here)))
#here::i_am("analysis/analyze_refbias_trees_emp_080621.R")
source(here("analysis","theme_custom.R"))

############################################
## Specifying file names ########
############################################

emp.trees <- "010822-lates-emp-batch.trees"
emp.tree.names <- "010822-lates-emp-tree.names"
output <- "010822-lates-emp-output"
refdist <- "010822-emp-refdist.txt"
sites <- "010822-SNPs-emp-lates"
outgroup <- "SRR3140997.Lcal"
args <- mget(c("emp.trees","emp.tree.names","refdist","sites","outgroup","output"))

# emp.trees <- "031623-lates-emp-batch.trees"
# emp.tree.names <- "031623-lates-emp-tree.names"
# output <- "031623-lates-emp-output"
# refdist <- "031623-emp-refdist.txt"
# sites <- "031623-SNPs-emp-lates"
# outgroup <- "SRR3140997.Lcal"
# args <- mget(c("emp.trees","emp.tree.names","refdist","sites","outgroup","output"))

####Read in concatenated tree file 
emp.trees<-read.tree(here("output","new",args$emp.trees))

###Read in file of tree names
emp.tree.names<-read.table(here("output","new",args$emp.tree.names),stringsAsFactors=FALSE)[,1][1:length(emp.trees)]

###Import reference distance data
# refdist <- read_table2(paste("output/new/",args$refdist,sep=""), col_names=c("simulation","height","int","taxa_ref","avg_dxy")) %>%
#   group_by(simulation,height,int,taxa_ref) %>%
#   slice(1) %>%
#   ungroup()

sites <- read_csv(here("output","new",args$sites),col_names=c("tree","snps")) %>% 
  mutate(simulation=as.integer(gsub(".*s([0-9]+)_q40.*","\\1",tree)),
         missing=as.numeric(gsub(".*miss([0-9]\\.?[0-9]*)\\_.*","\\1",tree)),
         maf=as.numeric(gsub(".*_mac([0-9]+).*","\\1",tree)),
         int=gsub(".*\\_([A-Z]+)\\.noInv.*","\\1",tree))
dxy <- read_csv(here("output","new","010822-emp-dxy.csv")) %>%
  filter(height == "Lates")

#####################
## Start of analysis!
#####################

###Create empty data frame and with named, empty columns 
num.trees <- length(emp.tree.names)
results.emp <- tibble(tree.num=seq(1:length(emp.trees)))

###Fix names in RAxML trees
fix_labels <- function(x) {
  labs <- x$tip.label
  new_labs <- gsub(".*aln_(.*)\\.sorted\\.bam","\\1",labs)
  x$tip.label <- new_labs
  return(x)
}

emp.trees <- lapply(emp.trees, function(x) fix_labels(x))

## pull info from tree names, and calculate tree statistics
results.emp <- results.emp %>%
  mutate(
    simulation = sapply(emp.tree.names,function(x) as.integer(regmatches(x, regexec('_s([0-9]+)\\_q', x))[[1]][2])),
    method = "raxml",
    height = "Lates",
    quality = sapply(emp.tree.names,function(x) as.integer(regmatches(x, regexec('_q(40)\\_miss', x))[[1]][2])),
    missing = sapply(emp.tree.names,function(x) as.numeric(regmatches(x, regexec('miss(0.*?)\\_mac', x))[[1]][2])),
    maf = sapply(emp.tree.names, function(x) as.numeric(regmatches(x, regexec('mac([0-9]+)\\.REF', x))[[1]][2])),
    int = sapply(emp.tree.names, function(x) regmatches(x, regexec('REF.([A-Z]+)\\.emp', x))[[1]][2]),
    noref = "REF"
  ) %>%
  mutate(
    tree.height = sapply(emp.trees, function(x) max(branching.times(root(x,outgroup,resolve.root=TRUE)),na.rm=T)),
    ingroup.tree.height = sapply(emp.trees, function(x) max(branching.times(drop.tip(root(x,outgroup,resolve.root=TRUE),"SRR3140997")),na.rm=T)),
    Avg.BLs = sapply(emp.trees, function(x) mean(root(x,outgroup,resolve.root=TRUE)$edge.length,na.rm=T)),
    SD.BLs = sapply(emp.trees, function(x) sd(root(x,outgroup,resolve.root=TRUE)$edge.length,na.rm=T)),
    mean.support = sapply(emp.trees, function(x) mean(as.numeric(root(x,outgroup,resolve.root=TRUE)$node.label[-1]),na.rm=T)),
    sd.support = sapply(emp.trees, function(x) sd(as.numeric(root(x,outgroup,resolve.root=TRUE)$node.label[-1]),na.rm=T))
  ) %>%
  mutate(
    gamma = sapply(emp.trees, function(x) gammaStat(chronopl(root(x,outgroup,resolve.root=TRUE),lambda=1))),
    ingroup.gamma = sapply(emp.trees, function(x) gammaStat(drop.tip(chronopl(root(x,outgroup,resolve.root=TRUE),lambda=1),"SRR3140997"))),
    colless = sapply(emp.trees, function(x) colless(as.treeshape(root(x,outgroup,resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.colless = sapply(emp.trees, function(x) colless(as.treeshape(drop.tip(root(x,outgroup,resolve.root=TRUE),"SRR3140997"),model="pda"),norm="pda")),
    sackin = sapply(emp.trees, function(x) sackin(as.treeshape(root(x,outgroup,resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.sackin = sapply(emp.trees, function(x) sackin(as.treeshape(drop.tip(root(x,outgroup,resolve.root=TRUE),"SRR3140997"),model="pda"),norm="pda"))
  ) 

results.emp <- results.emp %>%
  left_join(sites %>% group_by(tree) %>% slice_head(n=1),by=c("simulation","missing","maf","int")) %>%
  left_join(dxy,by=c("simulation"="sim","int"="int","height"="height"))

results.raxml.lates <- results.emp %>%
  group_by(simulation,method,missing,maf,int) %>%
  slice_head(n=1) %>%
  ungroup()
#write.csv(results.raxml.lates,file=here("output","new",paste(args$output,"-raxml-031623.csv",sep="")),quote=FALSE,row.names=TRUE,na="NA")
results.emp <- read_csv(here("output","new",paste(args$output,"-raxml.csv",sep="")))[,-1]
#results.emp <- results.emp.lates

######################
## univariate plots
#####################
results.emp$maf <- as.factor(results.emp$maf)
results.emp$missing <- as.factor(results.emp$missing)
results.emp$dxy <- as.numeric(results.emp$dxy)
results.emp$simulation <- as.factor(results.emp$simulation)

plot1 <- ggplot(data = results.emp, 
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

plot2 <- ggplot(data = results.emp, 
                aes(x=maf,
                    y=Avg.BLs,
                    fill=NULL)) +
  geom_boxplot(aes(alpha=0.9),outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot3 <- ggplot(data = results.emp, 
                aes(x=maf,
                    y=ingroup.tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot4 <- ggplot(data = results.emp, 
                aes(x=maf,
                    y=colless,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

plot5 <- ggplot(data = results.emp, 
                aes(x=maf,
                    y=gamma,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~int)

ggarrange(plot1,plot2,plot3,plot4,plot5,ncol=1) 

results.emp %>% 
  ggplot(aes(x=as.numeric(as.character(maf)),y=ingroup.gamma,shape=int)) + 
    geom_point(aes(col=as.factor(missing))) + 
    geom_smooth(aes(group=missing,col=missing),span=2) + 
    theme_custom()

results.emp %>% 
  ggplot(aes(x=as.factor(as.character(maf)),y=ingroup.colless,shape=int)) + 
  geom_boxplot(aes(fill=as.factor(int))) +
  #geom_jitter(aes(col=as.factor(missing))) + 
  #geom_smooth(aes(group=missing,col=missing),span=2) + 
  theme_custom()

###################################
## ASTRAL TREES
###################################

####################
# astral.orig.trees.emp<-read.tree(here("output","new","030123-lates-emp-ASTRAL-batch.trees"))
# astral.orig.tree.names.emp<-read.table(here("output","new","030123-lates-emp-ASTRAL-tree.names"))[,1][1:length(astral.trees.emp)]
# astral.extra.trees.emp <- read.tree(here("output","new","030423-lates-emp-ASTRAL-batch.trees"))
# astral.extra.tree.names.emp <- read.table(here("output","new","030423-lates-emp-ASTRAL-tree.names"))[,1]
# 
# astral.all.tree.names <- tibble(name = astral.orig.tree.names.emp) %>%
#   add_column(tree_num = seq(1:length(astral.orig.tree.names.emp))) %>%
#   add_column(set = "orig") %>%
#   add_row(tibble(name=astral.extra.tree.names.emp) %>% 
#             add_column(tree_num = seq(1:length(astral.extra.tree.names.emp))) %>%
#             add_column(set = "extra")) %>%
#   group_by(name) %>%
#   slice_tail(n = 1) %>%
#   ungroup() 
# orig_trees <- astral.all.tree.names %>%
#   filter(set == "orig")
# astral.trees.emp <- c(astral.orig.trees.emp[orig_trees$tree_num],astral.extra.trees.emp)
# astral.tree.names.emp <- c(astral.orig.tree.names.emp[orig_trees$tree_num],astral.extra.tree.names.emp)

astral.trees.emp <- read.tree(here("output","new","031623-lates-emp-ASTRAL-batch.trees"))
astral.tree.names.emp <- read.table(here("output","new","031623-lates-emp-ASTRAL-tree.names"))[,1]

sites.astral <- read_csv(here("output","new","031623-SNPs-emp-lates"),col_names=c("tree","snps")) %>% 
  mutate(simulation=as.integer(gsub(".*s([0-9]+)_q40.*","\\1",tree)),
         missing=as.numeric(gsub(".*miss([0-9]\\.?[0-9]*)\\_.*","\\1",tree)),
         maf=as.numeric(gsub(".*_mac([0-9]+).*","\\1",tree)),
         int=gsub(".*\\.([A-Z]+)\\.noInv.*","\\1",tree)) %>%
  group_by(simulation,missing,maf,int) %>%
  slice(1) # keep only one entry per combination, in case of duplicates

# fix taxon names
astral.trees.emp <- lapply(astral.trees.emp, function(x) fix_labels(x))

###Create empty data frame and with named, empty columns 
num.trees.emp <- length(astral.tree.names.emp)
results.astral.emp <- tibble(tree.num=seq(1:length(astral.trees.emp)))

## pull info from tree names, and calculate tree statistics
results.astral.emp <- results.astral.emp %>%
  mutate(
    simulation = sapply(astral.tree.names.emp,function(x) as.integer(regmatches(x, regexec('sim([0-9]+)\\_m', x))[[1]][2])),
    height = "Lates",
    method = "astral",
    quality = 40,
    missing = sapply(astral.tree.names.emp,function(x) as.numeric(regmatches(x, regexec('_miss([0-9]?\\.?[0-9]+)\\_mac', x))[[1]][2])),
    maf = sapply(astral.tree.names.emp, function(x) as.numeric(regmatches(x, regexec('_mac(0?.?[0-9]+)_[A-Z]', x))[[1]][2])),
    #sites = sapply(astral.tree.names.emp, function(x) as.integer(regmatches(x, regexec('sites([0-9]*?)\\.', x))[[1]][2])),
    #taxa_ref = sapply(astral.tree.names.emp, function(x) regmatches(x,regexec('[A-Z]-([0-9]+_0_0)\\.phylip', x))[[1]][[2]]),
    int = sapply(astral.tree.names.emp, function(x) as.character(regmatches(x, regexec('[0-9]_([A-Z]+)', x))[[1]][2])),
    noref = "REF"
  ) %>% 
  left_join(select(results.emp,c(simulation,quality:noref)),by=c("simulation","missing","maf","int","quality","noref")) %>%
  left_join(dxy,by=c("simulation"="sim","height"="height","int"="int"),copy=TRUE) %>%
  left_join(sites.astral %>% select(-tree),by=c("simulation","missing","maf","int")) %>%
  mutate(
    #tree.height = sapply(astral.trees, function(x) max(branching.times(root(x,"sim_0_0_0",resolve.root=TRUE)),na.rm=T)),
    #ingroup.tree.height = sapply(astral.trees, function(x) max(branching.times(drop.tip(root(x,"sim_0_0_0",resolve.root=TRUE),"sim_0_0_0")),na.rm=T)),
    #Avg.BLs = sapply(astral.trees, function(x) mean(root(x,"sim_0_0_0",resolve.root=TRUE)$edge.length,na.rm=T)),
    #SD.BLs = sapply(astral.trees, function(x) sd(root(x,"sim_0_0_0",resolve.root=TRUE)$edge.length,na.rm=T)),
    mean.support = sapply(astral.trees.emp, function(x) mean(as.numeric(root(x,outgroup,resolve.root=TRUE)$node.label[-1]),na.rm=T)),
    sd.support = sapply(astral.trees.emp, function(x) sd(as.numeric(root(x,outgroup,resolve.root=TRUE)$node.label[-1]),na.rm=T))
  ) %>%
 mutate(
#   gamma = sapply(astral.trees, function(x) gammaStat(chronopl(root(x,"sim_0_0_0",resolve.root=TRUE),lambda=1))),
#   ingroup.gamma = sapply(astral.trees, function(x) gammaStat(drop.tip(chronopl(root(x,"sim_0_0_0",resolve.root=TRUE),lambda=1),"sim_0_0_0"))),
   colless = sapply(astral.trees.emp, function(x) colless(as.treeshape(root(x,outgroup,resolve.root=TRUE),model="pda"),norm="pda")),
   ingroup.colless = sapply(astral.trees.emp, function(x) colless(as.treeshape(drop.tip(root(x,outgroup,resolve.root=TRUE),outgroup),model="pda"),norm="pda")),
   sackin = sapply(astral.trees.emp, function(x) sackin(as.treeshape(root(x,outgroup,resolve.root=TRUE),model="pda"),norm="pda")),
   ingroup.sackin = sapply(astral.trees.emp, function(x) sackin(as.treeshape(drop.tip(root(x,outgroup,resolve.root=TRUE),outgroup),model="pda"),norm="pda"))
 ) 

astral.node.desc.emp <- results.astral.emp %>%
  mutate(avg.node.descendants = sapply(astral.trees.emp, function(x) mean(sapply(seq(length(x$tip.label),length(x$tip.label)+x$Nnode), function(y) length(Descendants(x,y,"tips")[[1]])))))


summary(results.astral.emp)  
results.astral.emp.uniq <- results.astral.emp %>%
  group_by(simulation,height,method,quality,missing,maf,int,noref) %>%
  slice_tail(n=1)
write.csv(results.astral.emp.uniq,
           file=here("output","new",paste(args$output,"-astral.csv",sep="")),
           quote=FALSE,row.names=TRUE,na="NA")

######################
## univariate plots
#####################
results.all <- results.emp %>%
  add_row(results.astral.emp)  
results.all$maf <- as.factor(results.all$maf)
plot1 <- ggplot(data = results.all, 
                aes(x=int,
                    y=tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +
  theme_custom() +
  theme(legend.position="none",
        axis.title.x = element_blank())+
  facet_wrap(~method)

plot2 <- ggplot(data = results.all, 
                aes(x=int,
                    y=ingroup.colless,
                    fill=NULL)) +
  geom_boxplot(aes(alpha=0.9),outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~method)

plot3 <- ggplot(data = results.all, 
                aes(x=int,
                    y=ingroup.sackin,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=maf),width=0.1,height=0.01,alpha=0.4) +  
  theme_custom()  +
  theme(legend.position="none",
        axis.title.x = element_blank()) +
  facet_wrap(~method)

ggarrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,nrow=1) 




###################
## Plot trees in PCoA space by simulation
# pdf(here("output","new","lates_emp_trees_pcoa.pdf"))
# param_loadings <- tibble(sim=integer(),height=character(),param=character(),x_corr=numeric(),y_corr=numeric(),
#                          x_corr_sig=numeric(),y_corr_sig=numeric())
# h <- "lates"
# for (s in unique(results.emp.lates$simulation)){
#   trees.subset <- c(emp.trees[results.emp.lates$simulation == s],astral.trees.emp[results.astral.emp$simulation == s])
#   class(trees.subset) <- "multiPhylo"
#   trees.subset <- root(trees.subset,outgroup,resolve.root=TRUE)
#   rf.dist <- InfoRobinsonFoulds(trees.subset,normalize=TRUE)
#   rf.pcoa <- pcoa(rf.dist)
#   rf.pcoa.plot <- rf.pcoa$vectors %>%
#     as_tibble() %>%
#     ggplot(aes(x=Axis.1,y=Axis.2)) +
#     geom_jitter(aes(col=as.factor(c(results.emp.lates$method[results.emp.lates$simulation == s],
#                                     results.astral.emp$method[results.astral.emp$simulation == s])),
#                    fill=as.factor(c(results.emp.lates$maf[results.emp.lates$simulation == s],
#                                   results.astral.emp$maf[results.astral.emp$simulation == s])),
#                               #text=results.emp.lates$missing[results.emp.lates$simulation == s],
#                    shape=factor(c(results.emp.lates$int[results.emp.lates$simulation == s],
#                                   results.astral.emp$int[results.astral.emp$simulation == s]))),
#                    size=5,width=0.01,height=0.01) +
#     #theme_bw(base_size=12, base_family="Open Sans Light") +
#     theme_custom() +
#     theme(
#       panel.background  = element_blank(),
#       plot.background = element_rect(fill="white", color=NA), 
#       legend.background = element_rect(fill="transparent", colour=NA),
#       legend.key = element_rect(fill="transparent", colour=NA),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="Open Sans"),
#       axis.text = element_text(size=14),
#       legend.text = element_text(size=14),
#       legend.title = element_blank()
#     ) +
#     #scale_color_manual(values=PNWColors::pnw_palette("Sunset2",7,type="continuous")) +
#     scale_color_viridis_d(direction=-1) +
#     scale_fill_viridis_d(alpha=0.8,direction=-1) +
#     scale_shape_manual(values=c(21,23)) +
#     ggtitle(paste0("Lates empirical PCoA, Simulation ",s)) +
#     geom_hline(yintercept=0,lty=2,col="gray50") +
#     geom_vline(xintercept=0,lty=2,col="gray50") +
#     # biplot arrows for MAF
#     geom_segment(
#       x = 0, y = 0,
#       xend = cor(as.numeric(as.character(c(results.emp.lates$maf[results.emp.lates$simulation == s],
#                                            results.astral.emp$maf[results.astral.emp$simulation == s]))), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#       yend = cor(as.numeric(as.character(c(results.emp.lates$maf[results.emp.lates$simulation == s],
#                                            results.astral.emp$maf[results.astral.emp$simulation == s]))), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015,
#       lineend = "round", # See available arrow types in example above
#       linejoin = "round",
#       size = 1, 
#       arrow = arrow(length = unit(0.2, "inches")),
#       colour = "black" # Also accepts "red", "blue' etc
#     ) + 
#     # annotate(geom = "text",
#     #          x = cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#     #          y = cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#     #          label = "MAC",
#     #          hjust=0,vjust=0
#     # ) +
#     # biplot arrows for Missing
#     geom_segment(
#       x = 0, y = 0,
#       xend = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#       yend = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015,
#       lineend = "round", # See available arrow types in example above
#       linejoin = "round",
#       size = 1, 
#       arrow = arrow(length = unit(0.2, "inches")),
#       colour = "black" # Also accepts "red", "blue' etc
#     ) + 
#     # annotate(geom = "text",
#     #          x = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#     #          y = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#     #          label = "Missing",
#     #          hjust=1,vjust=1
#     # ) +
#     # biplot arrows for INT
#     geom_segment(
#       x = 0, y = 0,
#       xend = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#       yend = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015,
#       lineend = "round", # See available arrow types in example above
#       linejoin = "round",
#       size = 1, 
#       arrow = arrow(length = unit(0.2, "inches")),
#       colour = "black" # Also accepts "red", "blue' etc
#     # ) +
#     # annotate(geom = "text",
#     #          x = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#     #          y = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
#     #          label = "INT",
#     #          hjust=0,vjust=0
#     )
#     
#   
#   print(rf.pcoa.plot)
#   #plotly::ggplotly(rf.pcoa.plot)
#   vectors <- tibble(sim = s, height = h,
#                     param = c("maf","missing","int"),
#                     x_corr = c(cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]),
#                                cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]),
#                                cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1])
#                     ),
#                     y_corr = c(cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]),
#                                cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]),
#                                cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2])
#                     ),
#                     x_corr_sig = c(cor.test(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value
#                     ),
#                     y_corr_sig = c(cor.test(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value
#                     ))
#   param_loadings <- param_loadings %>%
#     add_row(vectors)
# }
# dev.off()
# 
# ## plots of parameter correlations
# param_loadings %>% 
#   mutate(abs_x_corr = abs(x_corr),
#          abs_y_corr = abs(y_corr)) %>%
#   ggscatter(x="abs_x_corr",y="abs_y_corr",
#             color="param",fill="param",
#             facet.by=c("param","height"),
#             size=3,alpha=0.75) + 
#   xlab("PCoA 1 Parameter Correlation") + 
#   ylab("PCoA 2 Parameter Correlation") +
#   theme(legend.position = "none",
#         strip.text = element_text(size=rel(1.3)),
#         axis.title = element_text(size=rel(1.5)),
#         axis.text = element_text(size=rel(1.5)))
# 
# param_loadings %>% 
#   mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
#          param2 = recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
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
#   ggplot(aes(x=abs_x_corr,y=abs_y_corr)) +
#   geom_point(aes(col=sig_cat,fill=sig_cat,shape=as.factor(height)),size=6) + 
#   scale_shape_manual(values=c(21,23)) +
#   facet_wrap(~param2) +
#   # ggscatter(x="abs_x_corr",y="abs_y_corr",
#   #           color="param",fill="param",
#   #           facet.by=c("height2","param2"),
#   #           size=3,alpha=0.75) + 
#   #xlab("PCoA 1 Parameter Correlation") + 
#   #ylab("PCoA 2 Parameter Correlation") +
#   scale_color_manual(values=PNWColors::pnw_palette("Sunset2",4)) +
#   scale_fill_manual(values=scales::alpha(PNWColors::pnw_palette("Sunset2",4),0.8)) +
#   theme_custom() +
#   theme(legend.position = "right",
#         strip.text = element_text(size=rel(1.3)),
#         axis.title = element_text(size=rel(1.5)),
#         axis.text = element_text(size=rel(1.5))) +
#   xlab("PCoA Axis 1 Correlation") +
#   ylab("PCoA Axis 2 Correlation")

results.emp <- results.raxml.lates
# results.astral.emp <- results.astral.lates
# results.astral.emp.uniq <- results.astral.emp %>%
#   group_by(simulation,height,method,quality,missing,maf,int,noref) %>%
#   slice_tail(n=1)
param_loadings.emp <- tibble(sim=integer(),height=character(),method=character(),param=character(),
                             pct_var_exp=numeric(),x_corr=numeric(),y_corr=numeric(),
                         x_corr_sig=numeric(),y_corr_sig=numeric())
h <- "Lates"
for (m in c("astral","raxml")) {
  for (i in unique(as.numeric(as.character(results.emp$simulation)))){
    #for (i in seq(16,30,by=1)){
    #j <- which(ml.tree.info$simulation == i & ml.tree.info$height == h)
    
    # subset_r <- results.emp$tree.num[results.emp$simulation == i & results.emp$height == h & results.emp$noref == "REF"]
    # subset_a <- results.astral.emp$tree.num[results.astral.emp$simulation == i & results.astral.emp$height == h]
    # trees.subset <- c(emp.trees[subset_r],astral.trees.emp[subset_a])
    
    if (m == "raxml"){
      subset_r <- results.emp$tree.num[results.emp$simulation == i & results.emp$height == h & results.emp$noref == "REF"]
      subset_a <- NULL
      
      if (length(subset_r) !=0) {
        trees.subset <- c(emp.trees[results.emp$tree.num %in% subset_r])
      } else {
        print(paste("no trees for sim",i," for height ",h))
        next
      }
    } 
    
    if (m == "astral"){
      subset_a <- results.astral.emp.uniq$tree.num[results.astral.emp.uniq$simulation == i & results.astral.emp.uniq$height == h]
      subset_r <- NULL
      
      if (length(subset_a) !=0) {
        trees.subset <- astral.trees.emp[results.astral.emp$tree.num %in% subset_a]
      } else {
        print(paste("no trees for sim",i," for height ",h))
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
      geom_jitter(aes(col=factor(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]),levels=c("0","1","2","3","4","5","10")),
                      fill=factor(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]),levels=c("0","1","2","3","4","5","10")),
                      #text=results.emp$missing[results.emp$simulation == s],
                      #shape=factor(c(results.emp$missing[subset])),
                      size=factor(c(results.emp$method[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$method[results.astral.emp.uniq$tree.num %in% subset_a])),
                      shape=factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a])),
      ),
      width=0.005,height=0.005) +
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
        xend = cor(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 1)*0.015,
        yend = cor(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 1)*0.015,
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1,
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) +
      annotate(geom = "text",
               x = cor(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015-0.001,
               y = cor(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015-0.002,
               label = "MAF",
               hjust=1,vjust=1
      ) +
      # biplot arrows for Missing
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
        yend = cor(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1,
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) +
      annotate(geom = "text",
               x = cor(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015+0.001,
               y = cor(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015-0.002,
               label = "Missing",
               hjust=0,vjust=0
      ) +
      # biplot arrows for INT
      geom_segment(
        x = 0, y = 0,
        xend = cor(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
        yend = cor(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
        lineend = "round", # See available arrow types in example above
        linejoin = "round",
        size = 1,
        arrow = arrow(length = unit(0.2, "inches")),
        colour = "black" # Also accepts "red", "blue' etc
      ) +
      annotate(geom = "text",
               x = cor(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015+0.001,
               y = cor(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015+0.002,
               label = "INT",
               hjust=0,vjust=1
      )
    
    #print(rf.pcoa.plot)
    #plotly::ggplotly(rf.pcoa.plot)
    vectors <- tibble(sim = i, height = h, method = m,
                      param = c("maf","missing","int"),
                      pct_var_exp = sum(rf.pcoa$values$Relative_eig[1:2]),
                      x_corr = c(cor(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]),
                                 cor(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]),
                                 cor(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])
                      ),
                      y_corr = c(cor(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]),
                                 cor(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]),
                                 cor(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])
                      ),
                      x_corr_sig = c(cor.test(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])$p.value,
                                     cor.test(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])$p.value,
                                     cor.test(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])$p.value
                      ),
                      y_corr_sig = c(cor.test(as.numeric(as.character(c(results.emp$maf[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$maf[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])$p.value,
                                     cor.test(as.numeric(as.character(c(results.emp$missing[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$missing[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])$p.value,
                                     cor.test(as.numeric(as.factor(c(results.emp$int[results.emp$tree.num %in% subset_r],results.astral.emp.uniq$int[results.astral.emp.uniq$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])$p.value
                      ))
    param_loadings.emp <- param_loadings.emp %>%
      add_row(vectors)
  }
}

dev.off()
param_loadings_astral.lates <- param_loadings.emp %>% mutate(method="astral")
param_loadings_raxml.lates <- param_loadings.emp %>% mutate(method="raxml")
param_loadings.lates <- param_loadings_astral.lates %>%
  add_row(param_loadings_raxml.lates)

## proportion of significant correlations
param_loadings.emp %>%
  #add_row(param_loadings_raxml.lates) %>%
  mutate(dom = case_when(x_corr > y_corr ~ "PC1",
                         y_corr > x_corr ~ "PC2"),
         sig_x = case_when(x_corr_sig < 0.01 & x_corr > y_corr ~ TRUE,
                           x_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE),
         sig_y = case_when(y_corr_sig < 0.01 & y_corr > x_corr ~ TRUE,
                           y_corr_sig >= 0.01 ~ FALSE,
                           TRUE ~ FALSE)) %>%
  group_by(param,method) %>%
  summarize(n_PC1 = sum(sig_x),
            n_PC2 = sum(sig_y),
            mean_PC1 = mean(abs(x_corr)),
            mean_PC2 = mean(abs(y_corr)))


p_pcoa_param_lates <- param_loadings.emp %>% 
  #add_row(param_loadings_raxml2.lates) %>%
  #add_row(param_loadings_raxml.lates) %>%
  #add_row(param_loadings_astral2.lates) %>%
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
        axis.text = element_text(size=rel(1.1))) +
  xlab("PCoA Axis 1 Correlation") +
  ylab("PCoA Axis 2 Correlation")


#################
## Plot extreme trees for imbalance
#################

# most imbalanced = tree 204,206,207, sim 5
results.emp.lates[results.emp.lates$ingroup.colless == max(results.emp.lates$ingroup.colless),]

# least imbalanced = tree 132, sim 3
results.emp.lates[results.emp.lates$ingroup.colless == min(results.emp.lates$ingroup.colless),]

# sim 5 trees -- min = 202, max = 204,206,207
sim5_info <- results.emp.lates[results.emp.lates$simulation == 5,]
emp.trees[[202]]$tip.label %in% emp.trees[[204]]$tip.label

t1.lates <- chronos(drop.tip(root(emp.trees[[202]],outgroup,resolve.root = TRUE),outgroup))
t2.lates <- chronos(drop.tip(root(emp.trees[[204]],outgroup,resolve.root = TRUE),outgroup))
lates_metadat <- read_csv("~/Dropbox/i/projects/lates_popgen/scripts/lates_popgen_github_copy/data/lates_all_metadata.csv")
t1.lates.spp <- t1.lates$tip.label %>%
  as_tibble() %>%
  left_join(lates_metadat %>% select(ind,sampling_loc,ent_ID),
            by=c("value" = "ind"))

lates.imb.cophylo <- cophylo(t1.lates,
             t2.lates,
             rotate=TRUE,
             assoc=cbind(t1.lates$tip.label,t1.lates$tip.label))
cols <- PNWColors::pnw_palette("Sunset2",4,type="continuous")
plot(lates.imb.cophylo, type="phylogram", use.edge.length=FALSE,
     length.line=0.5, link.type="curved",link.lwd=3, link.lty="solid",
     link.col=cols[as.factor(t1.lates.spp$ent_ID)],
     assoc=cbind(t1.lates$tip.label,t1.lates$tip.label))             

#################
## Plot extreme trees for gamma
#################
# most imbalanced = tree 99,101,102, sim 2
results.emp.lates[results.emp.lates$ingroup.gamma == max(results.emp.lates$ingroup.gamma),]

# least imbalanced = tree 625, sim 7
results.emp.lates[results.emp.lates$ingroup.gamma == min(results.emp.lates$ingroup.gamma),]

# sim 2 trees -- min = 100, max = 99,101,102
sim2_info <- results.emp.lates[results.emp.lates$simulation == 2,]


t1.lates <- chronos(drop.tip(root(emp.trees[[102]],outgroup,resolve.root = TRUE),outgroup),lambda=1)
t2.lates <- chronos(drop.tip(root(emp.trees[[100]],outgroup,resolve.root = TRUE),outgroup),lambda=1)
lates_metadat <- read_csv("~/Dropbox/i/projects/lates_popgen/scripts/lates_popgen_github_copy/data/lates_all_metadata.csv")
t1.lates.spp <- t1.lates$tip.label %>%
  as_tibble() %>%
  left_join(lates_metadat %>% select(ind,sampling_loc,ent_ID),
            by=c("value" = "ind"))

lates.imb.cophylo <- cophylo(t1.lates,
                             t2.lates,
                             rotate=TRUE,
                             assoc=cbind(t1.lates$tip.label,t1.lates$tip.label))
cols <- PNWColors::pnw_palette("Sunset2",4,type="continuous")
plot(lates.imb.cophylo, type="phylogram", use.edge.length=FALSE,
     length.line=0.5, link.type="curved",link.lwd=3, link.lty="solid",
     link.col=cols[as.factor(t1.lates.spp$ent_ID)],
     assoc=cbind(t1.lates$tip.label,t1.lates$tip.label))   

