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
here::i_am("analysis/analyze_refbias_trees_emp_080621.R")
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

####Read in concatenated tree file 
emp.trees<-read.tree(here("output","new",args$emp.trees))

###Read in file of tree names
emp.tree.names<-read.table(here("output","new",args$emp.tree.names),stringsAsFactors=FALSE)[,1]

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
  left_join(sites,by=c("simulation","missing","maf","int")) %>%
  left_join(dxy,by=c("simulation"="sim","int"="int"))

results.emp.lates <- results.emp
write.csv(results.emp.lates,file=here("output","new",paste(args$output,"-raxml.csv",sep="")),quote=FALSE,row.names=TRUE,na="NA")

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

## Plot trees in PCoA space by simulation
pdf(here("output","new","lates_emp_trees_pcoa.pdf"))
param_loadings <- tibble(sim=integer(),height=character(),param=character(),x_corr=numeric(),y_corr=numeric(),
                         x_corr_sig=numeric(),y_corr_sig=numeric())
h <- "lates"
for (s in unique(results.emp.lates$simulation)){
  trees.subset <- emp.trees[results.emp.lates$simulation == s]
  class(trees.subset) <- "multiPhylo"
  trees.subset <- root(trees.subset,outgroup,resolve.root=TRUE)
  rf.dist <- InfoRobinsonFoulds(trees.subset,normalize=TRUE)
  rf.pcoa <- pcoa(rf.dist)
  rf.pcoa.plot <- rf.pcoa$vectors %>%
    as_tibble() %>%
    ggplot(aes(x=Axis.1,y=Axis.2)) +
    geom_jitter(aes(col=as.factor(results.emp.lates$maf[results.emp.lates$simulation == s]),
                   fill=as.factor(results.emp.lates$maf[results.emp.lates$simulation == s]),
                              #text=results.emp.lates$missing[results.emp.lates$simulation == s],
                   shape=factor(results.emp.lates$int[results.emp.lates$simulation == s])),size=5,width=0.01,height=0.01) +
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
    #scale_color_manual(values=PNWColors::pnw_palette("Sunset2",7,type="continuous")) +
    scale_color_viridis_d(direction=-1) +
    scale_fill_viridis_d(alpha=0.8,direction=-1) +
    scale_shape_manual(values=c(21,23)) +
    ggtitle(paste0("Lates empirical PCoA, Simulation ",s)) +
    geom_hline(yintercept=0,lty=2,col="gray50") +
    geom_vline(xintercept=0,lty=2,col="gray50") +
    # biplot arrows for MAF
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
      yend = cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015,
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1, 
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) + 
    # annotate(geom = "text",
    #          x = cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
    #          y = cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
    #          label = "MAC",
    #          hjust=0,vjust=0
    # ) +
    # biplot arrows for Missing
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
      yend = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015,
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1, 
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    ) + 
    # annotate(geom = "text",
    #          x = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
    #          y = cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
    #          label = "Missing",
    #          hjust=1,vjust=1
    # ) +
    # biplot arrows for INT
    geom_segment(
      x = 0, y = 0,
      xend = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
      yend = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015,
      lineend = "round", # See available arrow types in example above
      linejoin = "round",
      size = 1, 
      arrow = arrow(length = unit(0.2, "inches")),
      colour = "black" # Also accepts "red", "blue' etc
    # ) +
    # annotate(geom = "text",
    #          x = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
    #          y = cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]) * 0.5 * sqrt(nrow(results.emp.lates) - 1)*0.015, 
    #          label = "INT",
    #          hjust=0,vjust=0
    )
    
  
  print(rf.pcoa.plot)
  #plotly::ggplotly(rf.pcoa.plot)
  vectors <- tibble(sim = s, height = h,
                    param = c("maf","missing","int"),
                    x_corr = c(cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]),
                               cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1]),
                               cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1])
                    ),
                    y_corr = c(cor(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]),
                               cor(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2]),
                               cor(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2])
                    ),
                    x_corr_sig = c(cor.test(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
                                   cor.test(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
                                   cor.test(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,1][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value
                    ),
                    y_corr_sig = c(cor.test(as.numeric(as.character(results.emp.lates$maf[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
                                   cor.test(as.numeric(as.character(results.emp.lates$missing[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value,
                                   cor.test(as.numeric(as.factor(results.emp.lates$int[results.emp.lates$simulation == s])), rf.pcoa$vectors[,2][1:length(results.emp.lates$maf[results.emp.lates$simulation == s])])$p.value
                    ))
  param_loadings <- param_loadings %>%
    add_row(vectors)
}
dev.off()

## plots of parameter correlations
param_loadings %>% 
  mutate(abs_x_corr = abs(x_corr),
         abs_y_corr = abs(y_corr)) %>%
  ggscatter(x="abs_x_corr",y="abs_y_corr",
            color="param",fill="param",
            facet.by=c("param","height"),
            size=3,alpha=0.75) + 
  xlab("PCoA 1 Parameter Correlation") + 
  ylab("PCoA 2 Parameter Correlation") +
  theme(legend.position = "none",
        strip.text = element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.5)),
        axis.text = element_text(size=rel(1.5)))

param_loadings %>% 
  mutate(abs_x_corr = abs(x_corr),abs_y_corr = abs(y_corr),
         param2 = recode(.$param,int="Reference Genome",maf="Minor Allele Count",missing="Missing Data")) %>%
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
  geom_point(aes(col=sig_cat,fill=sig_cat,shape=as.factor(height)),size=6) + 
  scale_shape_manual(values=c(21,23)) +
  facet_wrap(~param2) +
  # ggscatter(x="abs_x_corr",y="abs_y_corr",
  #           color="param",fill="param",
  #           facet.by=c("height2","param2"),
  #           size=3,alpha=0.75) + 
  #xlab("PCoA 1 Parameter Correlation") + 
  #ylab("PCoA 2 Parameter Correlation") +
  scale_color_manual(values=PNWColors::pnw_palette("Sunset2",4)) +
  scale_fill_manual(values=scales::alpha(PNWColors::pnw_palette("Sunset2",4),0.8)) +
  theme_custom() +
  theme(legend.position = "right",
        strip.text = element_text(size=rel(1.3)),
        axis.title = element_text(size=rel(1.5)),
        axis.text = element_text(size=rel(1.5))) +
  xlab("PCoA Axis 1 Correlation") +
  ylab("PCoA Axis 2 Correlation")
