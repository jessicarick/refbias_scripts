############################################
## Script for analyzing refbias output trees from RAxML
## Written by J. Rick, 11 March 2019
## Updated 6 August 2021
## UPdated 18 Feb 2022 to use the here package
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
here::i_am("analysis/analyze_refbias_trees_emp_cichlids_080521.R")
source(here("analysis","theme_custom.R"))

############################################
## Specifying file names ########
############################################

emp.trees <- "010822-cichlids-emp-batch.trees"
emp.tree.names <- "010822-cichlids-emp-tree.names"
output <- "010822-cichlids-emp-output"
refdist <- "010822-emp-refdist.txt"
sites <- "010822-cichlids-emp-SNPs"
outgroup <- "SRR14041074.Onil"
args <- mget(c("emp.trees","emp.tree.names","refdist","sites","outgroup","output"))

####Read in concatenated tree file 
emp.trees<-read.tree(here("output","new",args$emp.trees))

###Read in file of tree names
emp.tree.names<-read.table(here("output","new",args$emp.tree.names),stringsAsFactors=FALSE)[,1]

# emp.tree.names.keep <- as_tibble(emp.tree.names) %>%
#   rownames_to_column() %>%
#   group_by(value) %>%
#   slice_tail(n=1)
# 
# emp.trees.keep <- emp.trees[as.integer(emp.tree.names.keep$rowname)]
# 
# emp.tree.names <- emp.tree.names.keep$value
# emp.trees <- emp.trees.keep

###Import reference distance data
# refdist <- read_table2(paste("output/new/",args$refdist,sep=""), col_names=c("simulation","height","int","taxa_ref","avg_dxy")) %>%
#   group_by(simulation,height,int,taxa_ref) %>%
#   slice(1) %>%
#   ungroup()

sites <- read_csv(here("output","new",args$sites),col_names=c("tree","snps")) %>% 
  mutate(simulation=as.integer(gsub(".*s([0-9]+)_q40.*","\\1",tree)),
         missing=as.numeric(gsub(".*miss([0-9]\\.?[0-9]*)\\_.*","\\1",tree)),
         maf=as.numeric(gsub(".*_mac([0-9]+)_.*","\\1",tree)),
         int=gsub(".*_([A-Z]+)\\.noInv.*","\\1",tree))
dxy <- read_csv(here("output","new","010822-emp-dxy.csv")) %>%
  filter(height == "Cichlids")

#####################
## Start of analysis!
#####################

###Create empty data frame and with named, empty columns 
num.trees <- length(emp.tree.names)
results.emp <- tibble(tree.num=seq(1:num.trees))

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
    height = "Cichlids",
    quality = sapply(emp.tree.names,function(x) as.integer(regmatches(x, regexec('_q(40)\\_miss0', x))[[1]][2])),
    missing = sapply(emp.tree.names,function(x) as.numeric(regmatches(x, regexec('miss(0.*?)\\_mac', x))[[1]][2])),
    maf = sapply(emp.tree.names, function(x) as.numeric(regmatches(x, regexec('mac([0-9]+)\\.REF', x))[[1]][2])),
    int = sapply(emp.tree.names, function(x) regmatches(x, regexec('REF.([A-Z]+)\\.emp', x))[[1]][2]),
    noref = "REF"
  ) %>%
  mutate(
    tree.height = sapply(emp.trees, function(x) max(branching.times(root(x,outgroup,resolve.root=TRUE)),na.rm=T)),
    ingroup.tree.height = sapply(emp.trees, function(x) max(branching.times(drop.tip(root(x,outgroup,resolve.root=TRUE),outgroup)),na.rm=T)),
    Avg.BLs = sapply(emp.trees, function(x) mean(root(x,outgroup,resolve.root=TRUE)$edge.length,na.rm=T)),
    SD.BLs = sapply(emp.trees, function(x) sd(root(x,outgroup,resolve.root=TRUE)$edge.length,na.rm=T)),
    mean.support = sapply(emp.trees, function(x) mean(as.numeric(root(x,outgroup,resolve.root=TRUE)$node.label[-1]),na.rm=T)),
    sd.support = sapply(emp.trees, function(x) sd(as.numeric(root(x,outgroup,resolve.root=TRUE)$node.label[-1]),na.rm=T))
  ) %>%
  mutate(
    gamma = sapply(emp.trees, function(x) gammaStat(chronopl(root(x,outgroup,resolve.root=TRUE),lambda=1))),
    ingroup.gamma = sapply(emp.trees, function(x) gammaStat(drop.tip(chronopl(root(x,outgroup,resolve.root=TRUE),lambda=1),outgroup))),
    colless = sapply(emp.trees, function(x) colless(as.treeshape(root(x,outgroup,resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.colless = sapply(emp.trees, function(x) colless(as.treeshape(drop.tip(root(x,outgroup,resolve.root=TRUE),outgroup),model="pda"),norm="pda")),
    sackin = sapply(emp.trees, function(x) sackin(as.treeshape(root(x,outgroup,resolve.root=TRUE),model="pda"),norm="pda")),
    ingroup.sackin = sapply(emp.trees, function(x) sackin(as.treeshape(drop.tip(root(x,outgroup,resolve.root=TRUE),outgroup),model="pda"),norm="pda"))
  ) 

results.emp <- results.emp %>%
  left_join(sites)  %>%
  left_join(dxy,by=c("simulation"="sim","int"="int"))

#write.csv(results.emp.cichlids,file=here("output","new",paste0(args$output,"-raxml.csv")),quote=FALSE,row.names=TRUE,na="NA")
results.emp.cichlids <- read_csv(here("output","new",paste0(args$output,"-raxml.csv")))[,-1]
results.emp <- results.emp.cichlids

######################
## univariate plots
#####################
results.emp$maf <- as.factor(results.emp$maf)
results.emp$missing <- as.factor(results.emp$missing)
  
results.emp$simulation <- as.factor(results.emp$simulation)

plot1 <- ggplot(data = results.emp, 
                aes(x=maf,
                    y=tree.height,
                    fill=NULL)) +
  geom_boxplot(alpha=0.9,outlier.shape=NA) +
  geom_jitter(aes(col=simulation),width=0.1,height=0.01,alpha=0.4) +
  theme_custom() +
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
  #filter(missing != 0) %>%
  ggplot(aes(x=as.numeric(as.character(maf)),y=ingroup.gamma)) + 
    geom_point(aes(col=as.factor(missing))) + 
    geom_smooth(aes(group=missing,col=missing),span=2) + 
    facet_wrap(~int) +
    xlab("MAF") +
    theme_custom()

results.emp %>%
  ggplot(aes(x=maf,y=ingroup.gamma,fill=missing)) +
  geom_jitter(aes(x=maf,col=missing),width=0.1) +
  ggdist::stat_halfeye(adjust = 0.6, # this changes the intervals used for calculating the density plot (bigger=smoother)
                       justification = -0.25,
                       .width = 0,
                       width = 0.8, 
                       point_colour = NA,
                       limits=c(0,NA),
                       alpha=0.5,
                       trim=TRUE, # whether to trim to data range-- looks nicer without trimming, but can be misleading
                       normalize="xy", # this normalizes each line, rather than having them all on the same scale
                       scale=1) +  # change to adjust height (scale=1 means that they fill the whole line height)
  theme_custom() +
  coord_flip() +
  theme(panel.grid.major = element_line(color="gray80", linetype="dashed"))
  

###################################
## ASTRAL TREES
###################################

####################
astral.trees.emp<-read.tree(here("output","new","022023-cichlids-emp-ASTRAL-batch.trees"))
astral.tree.names.emp<-read.table(here("output","new","022023-cichlids-emp-ASTRAL-tree.names"))[,1]

# fix taxon names
astral.trees.emp <- lapply(astral.trees.emp, function(x) fix_labels(x))

###Create empty data frame and with named, empty columns 
num.trees.emp <- length(astral.tree.names.emp)
results.astral.emp <- tibble(tree.num=seq(1:length(astral.trees.emp)))

## pull info from tree names, and calculate tree statistics
results.astral.emp <- results.astral.emp %>%
  mutate(
    simulation = sapply(astral.tree.names.emp,function(x) as.integer(regmatches(x, regexec('sim([0-9]+)\\_m', x))[[1]][2])),
    height = "Cichlids",
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
  #left_join(select(gt.rf,gt_rf:simulation),by=c("simulation","height")) %>%
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
# write.csv(results.astral.emp,
#           file=here("output","new",paste(args$output,"-astral.csv",sep="")),
#           quote=FALSE,row.names=TRUE,na="NA")

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
# pdf(here("output","new","cichlids_emp_trees_pcoa.pdf"))
# param_loadings <- tibble(sim=integer(),height=character(),param=character(),x_corr=numeric(),y_corr=numeric(),
#                          x_corr_sig=numeric(),y_corr_sig=numeric())
# h <- "cichlids"
# for (s in unique(results.emp.cichlids$simulation)){
#   trees.subset <- c(emp.trees[results.emp.cichlids$simulation == s],astral.trees.emp[results.astral.emp$simulation == s])
#   class(trees.subset) <- "multiPhylo"
#   trees.subset <- root(trees.subset,outgroup,resolve.root=TRUE)
#   rf.dist <- InfoRobinsonFoulds(trees.subset,normalize=TRUE)
#   rf.pcoa <- pcoa(rf.dist)
#   rf.pcoa.plot <- rf.pcoa$vectors %>%
#     as_tibble() %>%
#     ggplot(aes(x=Axis.1,y=Axis.2)) +
#     geom_jitter(aes(col=as.factor(c(results.emp.cichlids$method[results.emp.cichlids$simulation == s],
#                                     results.astral.emp$method[results.astral.emp$simulation == s])),
#                    fill=as.factor(c(results.emp.cichlids$maf[results.emp.cichlids$simulation == s],
#                                   results.astral.emp$maf[results.astral.emp$simulation == s])),
#                               #text=results.emp.cichlids$missing[results.emp.cichlids$simulation == s],
#                    shape=factor(c(results.emp.cichlids$int[results.emp.cichlids$simulation == s],
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
#     ggtitle(paste0("cichlids empirical PCoA, Simulation ",s)) +
#     geom_hline(yintercept=0,lty=2,col="gray50") +
#     geom_vline(xintercept=0,lty=2,col="gray50") +
#     # biplot arrows for MAF
#     geom_segment(
#       x = 0, y = 0,
#       xend = cor(as.numeric(as.character(c(results.emp.cichlids$maf[results.emp.cichlids$simulation == s],
#                                            results.astral.emp$maf[results.astral.emp$simulation == s]))), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#       yend = cor(as.numeric(as.character(c(results.emp.cichlids$maf[results.emp.cichlids$simulation == s],
#                                            results.astral.emp$maf[results.astral.emp$simulation == s]))), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015,
#       lineend = "round", # See available arrow types in example above
#       linejoin = "round",
#       size = 1, 
#       arrow = arrow(length = unit(0.2, "inches")),
#       colour = "black" # Also accepts "red", "blue' etc
#     ) + 
#     # annotate(geom = "text",
#     #          x = cor(as.numeric(as.character(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#     #          y = cor(as.numeric(as.character(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#     #          label = "MAC",
#     #          hjust=0,vjust=0
#     # ) +
#     # biplot arrows for Missing
#     geom_segment(
#       x = 0, y = 0,
#       xend = cor(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#       yend = cor(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015,
#       lineend = "round", # See available arrow types in example above
#       linejoin = "round",
#       size = 1, 
#       arrow = arrow(length = unit(0.2, "inches")),
#       colour = "black" # Also accepts "red", "blue' etc
#     ) + 
#     # annotate(geom = "text",
#     #          x = cor(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#     #          y = cor(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#     #          label = "Missing",
#     #          hjust=1,vjust=1
#     # ) +
#     # biplot arrows for INT
#     geom_segment(
#       x = 0, y = 0,
#       xend = cor(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#       yend = cor(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015,
#       lineend = "round", # See available arrow types in example above
#       linejoin = "round",
#       size = 1, 
#       arrow = arrow(length = unit(0.2, "inches")),
#       colour = "black" # Also accepts "red", "blue' etc
#     # ) +
#     # annotate(geom = "text",
#     #          x = cor(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#     #          y = cor(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.cichlids) - 1)*0.015, 
#     #          label = "INT",
#     #          hjust=0,vjust=0
#     )
#     
#   
#   print(rf.pcoa.plot)
#   #plotly::ggplotly(rf.pcoa.plot)
#   vectors <- tibble(sim = s, height = h,
#                     param = c("maf","missing","int"),
#                     x_corr = c(cor(as.numeric(as.character(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]),
#                                cor(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1]),
#                                cor(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1])
#                     ),
#                     y_corr = c(cor(as.numeric(as.character(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]),
#                                cor(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2]),
#                                cor(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2])
#                     ),
#                     x_corr_sig = c(cor.test(as.numeric(as.character(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])])$p.value
#                     ),
#                     y_corr_sig = c(cor.test(as.numeric(as.character(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.character(results.emp.cichlids$missing[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])])$p.value,
#                                    cor.test(as.numeric(as.factor(results.emp.cichlids$int[results.emp.cichlids$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.cichlids$maf[results.emp.cichlids$simulation == s])])$p.value
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


param_loadings.emp <- tibble(sim=integer(),height=character(),param=character(),x_corr=numeric(),y_corr=numeric(),
                             x_corr_sig=numeric(),y_corr_sig=numeric())
h <- "Cichlids"
m <- "astral"
for (i in unique(as.numeric(as.character(results.emp$simulation)))){
  if (m == "raxml"){
    subset_r <- results.emp$tree.num[results.emp$simulation == i & results.emp$height == h & results.emp$noref == "REF"]
    subset_a <- NULL
    
    if (length(subset_r) !=0) {
      trees.subset <- c(emp.trees[subset_r])
    } else {
      print(paste("no trees for sim",i," for height ",h))
      next
    }
  } 
  
  if (m == "astral"){
    subset_a <- results.astral.emp$tree.num[results.astral.emp$simulation == i & results.astral.emp$height == h]
    subset_r <- NULL
    
    if (length(subset_a) !=0) {
      trees.subset <- c(astral.trees.emp[subset_a])
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
    geom_jitter(aes(col=factor(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]),levels=c("0","1","2","3","4","5","10","true_tree")),
                    fill=factor(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]),levels=c("0","1","2","3","4","5","10","true_tree")),
                    #text=results.emp$missing[results.emp$simulation == s],
                    #shape=factor(c(results.emp$missing[subset])),
                    #size=factor(c(results.emp$method[subset_r],results.astral.emp$method[results.astral.emp$tree.num %in% subset_a])),
                    shape=factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a])),
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
      xend = cor(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 1)*0.015,
      yend = cor(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 1)*0.015,
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1,
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) +
    annotate(geom = "text",
             x = cor(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015-0.001,
             y = cor(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015-0.002,
             label = "MAF",
             hjust=1,vjust=1
    ) +
    # biplot arrows for Missing
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
      yend = cor(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1,
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) +
    annotate(geom = "text",
             x = cor(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015+0.001,
             y = cor(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015-0.002,
             label = "Missing",
             hjust=0,vjust=0
    ) +
    # biplot arrows for INT
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
      yend = cor(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015,
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1,
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) +
    annotate(geom = "text",
             x = cor(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015+0.001,
             y = cor(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]) * 0.5 * sqrt(length(trees.subset) - 2)*0.015+0.002,
             label = "INT",
             hjust=0,vjust=1
    )
  
  print(rf.pcoa.plot)
  #plotly::ggplotly(rf.pcoa.plot)
  vectors <- tibble(sim = i, height = h,
                    param = c("maf","missing","int"),
                    x_corr = c(cor(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]),
                               cor(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)]),
                               cor(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])
                    ),
                    y_corr = c(cor(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]),
                               cor(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)]),
                               cor(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])
                    ),
                    x_corr_sig = c(cor.test(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])$p.value,
                                   cor.test(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])$p.value,
                                   cor.test(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,1][1:length(trees.subset)])$p.value
                    ),
                    y_corr_sig = c(cor.test(as.numeric(as.character(c(results.emp$maf[subset_r],results.astral.emp$maf[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])$p.value,
                                   cor.test(as.numeric(as.character(c(results.emp$missing[subset_r],results.astral.emp$missing[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])$p.value,
                                   cor.test(as.numeric(as.factor(c(results.emp$int[subset_r],results.astral.emp$int[results.astral.emp$tree.num %in% subset_a]))), rf.pcoa$vectors[,2][1:length(trees.subset)])$p.value
                    ))
  param_loadings.emp <- param_loadings.emp %>%
    add_row(vectors)
}


dev.off()
param_loadings_astral.cichlids <- param_loadings.emp %>% mutate(method="astral")
param_loadings_raxml.cichlids <- param_loadings.emp %>% mutate(method="raxml")

## proportion of significant correlations
param_loadings_astral.cichlids %>%
  add_row(param_loadings_raxml.cichlids) %>%
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


p_pcoa_param_cich <- param_loadings_astral.cichlids %>% 
  #add_row(param_loadings_raxml2.cichlids) %>%
  add_row(param_loadings_raxml.cichlids) %>%
  #add_row(param_loadings_astral2.cichlids) %>%
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

ggarrange(p_pcoa_param_lates,p_pcoa_param_cich,ncol=1,labels=c("a) Lates","b) Cichlids"),
          label.x=0, label.y=1, font.label = list(size=14,family="Open Sans"))

p1 <- p_pcoa_param_lates +
  ggtitle('b) Lates')
p2 <- p_pcoa_param_cich +
  ggtitle('a) Tropheines')

p2/p1 
