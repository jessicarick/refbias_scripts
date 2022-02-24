###############################
## node descendants analysis
## fall 2021
## j. rick
###############################

#rf.test <- InfoRobinsonFoulds(trees.subset[[1]],trees.subset[[71]],normalize=TRUE,reportMatching=TRUE)
#VisualizeMatching(Func=RobinsonFouldsMatching,root(sp.tree,"sim_0_0_0",resolve.root=TRUE),root(orig.tree,"sim_0_0_0",resolve.root=TRUE))

# import data
here::i_am("analysis/analyze_node_descendants.R")

output <- "092321-output"
results.raxml <- read.csv(here("output","new",paste0(output,"-raxml.csv")),header=TRUE,row.names=1,sep=",")
raxml.trees <- read.tree(here("output","new",gsub("-output","-all-raxml.trees",output)))
raxml.trees2 <- raxml.trees[results.raxml$tree.num]

ml.tree<-read.tree(here("output","new",gsub("output","s_tree.tree",output)))
ml.tree.names <- read.table(here("output","new",gsub("output","s_tree.names",output)),stringsAsFactors = FALSE)[,1]

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
ml.tree.info <- tibble(simulation = numeric(length(ml.tree.names)),
                       height = character(length(ml.tree.names))) %>%
  mutate(simulation = sapply(ml.tree.names, function(x) as.integer(regmatches(x, regexec('sim([0-9]+)', x))[[1]][2])),
         height = sapply(ml.tree.names, function(x) regmatches(x, regexec('height([A-Z]+)', x))[[1]][2]))

##################
### FUNCTIONS ----
##################
# function to match nodes and calculate number of descendants/support for the nodes lost or retained
lostnodes <- function(i) {
  tree1 <- raxml.trees2[[i]]
  tree1.info <- results.raxml[i,]
  orig.info <- results.raxml[results.raxml$simulation == tree1.info$simulation &
                               results.raxml$height == tree1.info$height &
                               results.raxml$int == tree1.info$int &
                               results.raxml$maf == 0 & results.raxml$missing == 0,] %>%
    slice(1) # some are duplicated because of runs that got preempted, so only selecting first
  orig.tree <- raxml.trees2[[which(results.raxml$simulation == tree1.info$simulation &
                                    results.raxml$height == tree1.info$height &
                                    results.raxml$int == tree1.info$int &
                                    results.raxml$maf == 0 & results.raxml$missing == 0)[1]]]
  sp.tree <- ml.tree[[which(ml.tree.info$simulation == tree1.info$simulation &
                       ml.tree.info$height == tree1.info$height)]]
  
  match <- matchNodes(root(sp.tree,"sim_0_0_0",resolve.root=TRUE),root(tree1,"sim_0_0_0",resolve.root=TRUE),method="descendants") %>%
    as_tibble()
  bs <- matchNodes(root(sp.tree,"sim_0_0_0",resolve.root=TRUE),root(orig.tree,"sim_0_0_0",resolve.root=TRUE),method="descendants") %>%
    as_tibble() %>%
    left_join(tibble(node = seq(1:(root(orig.tree,"sim_0_0_0",resolve.root=TRUE,edge.label=TRUE)$Nnode))+Ntip(root(orig.tree,"sim_0_0_0",resolve.root=TRUE,edge.label=TRUE)),
                     bs = root(orig.tree,"sim_0_0_0",resolve.root=TRUE,edge.label=TRUE)$node.label),
                      by =c("tr2" = "node"),copy=TRUE) %>%
    mutate(bs = case_when(bs == "Root" ~ "100",
                          TRUE ~ bs))

  match  %>%
    as_tibble() %>%
    mutate(lost = case_when(is.na(tr2) ~ TRUE,
                            !is.na(tr2) ~ FALSE),
           descendants = sapply(tr1, function(x) length(Descendants(root(sp.tree,"sim_0_0_0",resolve.root=TRUE),x,"tips")[[1]])),
           bs_support = as.numeric(bs$bs),
           desc_category = case_when(descendants < 3 ~ "veryfew",
                                     descendants > 2 & descendants < 10 ~ "few",
                                     descendants > 40 ~ "many",
                                     TRUE ~ "moderate")) %>%
    # group_by(lost,desc_category) %>%
    # summarize(mean_desc = mean(descendants,na.rm=T),
    #           mean_supp = mean(mean_support,na.rm=T),
    #           #max_desc = max(descendants, na.rm=T),
    #           num = n(),
    #           mean_nodedepth = mean(tr1),
    #           median_nodedepth= median(tr1)) %>%
    add_column(tree_num = tree1.info$tree.num,.before = TRUE)
  
}

# same function, but for empirical data (no known tree for comparison)
lostnodes.emp <- function(i) {
  tree1 <- raxml.trees2[[i]]
  tree1.info <- results.raxml[i,]
  orig.info <- results.raxml[results.raxml$simulation == tree1.info$simulation &
                               results.raxml$int == tree1.info$int &
                               results.raxml$maf == 0 & results.raxml$missing == 0,] %>%
    slice(1)
  if (nrow(orig.info) == 1) {
    orig.tree <- raxml.trees2[[which(results.raxml$simulation == tree1.info$simulation &
                                       results.raxml$int == tree1.info$int &
                                       results.raxml$maf == 0 & results.raxml$missing == 0)[1]]]
    
    match <- matchNodes(orig.tree,tree1,method="descendants")
    match  %>%
      as_tibble() %>%
      mutate(lost = case_when(is.na(tr2) ~ TRUE,
                              TRUE ~ FALSE),
             descendants = sapply(tr1, function(x) length(Descendants(orig.tree,x,"tips")[[1]])),
             bs_support = case_when(tr1 == 52 ~ as.integer(100),
                                    TRUE ~ as.integer(orig.tree$node.label)),
             desc_category = case_when(descendants < 3 ~ "veryfew",
                                       descendants > 2 & descendants < 10 ~ "few",
                                       descendants > 40 ~ "many",
                                       TRUE ~ "moderate")) %>%
      # group_by(lost,desc_category) %>%
      # summarize(mean_desc = mean(descendants,na.rm=T),
      #           mean_supp = mean(mean_support,na.rm=T),
      #           #max_desc = max(descendants, na.rm=T),
      #           num = n(),
      #           mean_nodedepth = mean(tr1),
      #           median_nodedepth= median(tr1)) %>%
      add_column(tree_num = tree1.info$tree.num,.before = TRUE)
  }
  
}

###################################
## RUNNING FUNCTIONS ON SIMULATED DATA----
###################################
library(tidyverse)
library(ggdist)
library(viridis)
library(patchwork)

lostnode.info <- tibble(tree_num=integer(), tr1=integer(), tr2=integer(),
                        lost=logical(), descendants=integer(), bs_support=numeric(),
                        desc_category=character())
for (i in 1:nrow(results.raxml)){
  sum <- lostnodes(i)
  lostnode.info <- lostnode.info %>%
    add_row(sum)
}

lostnode.summary <- lostnode.info %>%
  #drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num")) %>%
  mutate(bs_support = as.numeric(bs_support)) %>%
  #filter(maf !=0) %>%
  #mutate(missing=as.factor(missing)) %>%
  group_by(lost,maf,height,simulation) %>%
  summarize(mean_supp = mean(bs_support,na.rm=T),
            med_supp = median(bs_support,na.rm=T),
            mean_desc = mean(descendants,na.rm=T),
            med_desc = median(descendants,na.rm=T),
            mean_nodedepth = mean(tr1,na.rm=T),
            med_nodedepth = median(tr1,na.rm=T),
            n_nodes = n()
            #mean_med_nodedepth = mean(median_nodedepth,na.rm=T),
            #med_med_nodedepth = median(median_nodedepth, na.rm=T)
            ) %>%
  ungroup()
lostnode.info2 <- lostnode.info %>%
  #drop_na() %>%
  mutate(bs_support = as.numeric(bs_support)) %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num")) %>%
  filter(simulation < 26) 

calc_boot_desc <- function(nnodes,mean,quantile=c(0.025,0.975)) {
  bs <- replicate(100,sample(lostnode.info$descendants,nnodes,replace=FALSE))
  if(is.null(nrow(bs))) {
    bs_means <- bs
  } else {
    bs_means <- colMeans(bs)
  }
  sd_bs <- sd(bs_means)
  mean_bs <- mean(bs_means)
  z <- (mean-mean_bs)/sd_bs
  quant <- quantile(bs_means,quantile)
  data.frame(zscore=z,lcl=quant[1],ucl=quant[2],row.names=NULL)
}
calc_boot_bs <- function(x,mean,quantile=c(0.025,0.975)) {
  bs <- replicate(100,sample(lostnode.info$bs_support[!is.na(lostnode.info$bs_support)],x,replace=FALSE))
  if(is.null(nrow(bs))) {
    bs_means <- bs
  } else {
    bs_means <- colMeans(bs,na.rm=T)
  }
  sd_bs <- sd(bs_means,na.rm=T)
  mean_bs <- mean(bs_means,na.rm=T)
  z <- (mean-mean_bs)/sd_bs
  quant <- quantile(bs_means,quantile)  
  data.frame(zscore=z,lcl=quant[1],ucl=quant[2],row.names=NULL)
}

lostnode.info2$bs_support <- as.numeric(lostnode.info2$bs_support)
null_lostnodes <- lostnode.info2 %>%
  #filter(!is.na(bs_support)) %>%
  group_by(lost,height,maf,missing,int) %>%
  summarize(n_nodes = n(),
            mean_desc = mean(descendants,na.rm=T),
            mean_bs = mean(bs_support,na.rm=T)) %>%
  #filter(lost == TRUE) %>%
  mutate(nodestraps = calc_boot_desc(n_nodes,mean_desc,c(0.025,0.975)),
         sig_nodes = case_when(nodestraps$zscore < -1.96 ~ "lower",
                         nodestraps$zscore > 1.96 ~ "higher",
                         TRUE ~ "ns")) %>%
  mutate(bsstraps = calc_boot_bs(n_nodes,mean_bs,c(0.025,0.975)),
         sig_bs = case_when(bsstraps$zscore < -1.96 ~ "lower",
                               bsstraps$zscore > 1.96 ~ "higher",
                               TRUE ~ "ns")) %>%
  unnest(cols=c(nodestraps,bsstraps),names_sep = ".") 
null.node.plot <- null_lostnodes %>% 
  ggplot() + 
  geom_density_ridges(aes(x=nodestraps.zscore,y=as.factor(maf),fill=lost),alpha=0.5) + 
  geom_point(aes(x=nodestraps.zscore,y=as.factor(maf),col=lost),pch="|") +
  theme_custom() + 
  geom_vline(xintercept=c(2,-2),linetype="dashed") +
  xlab("Node Descendants z-score") +
  ylab("Minor Allele Count") +
  facet_grid(~int,scales="free_x") +
  scale_fill_manual(values=c("#0A9396","#9B2226","gray80"),aesthetics = c("color","fill"),
                    labels=c("Retained","Lost")) +
  theme(legend.position=c(0.5,0.9),
        legend.text = element_text(size=rel(1.5)))

null.bs.plot <- null_lostnodes %>% 
  filter(!is.na(mean_bs) & bsstraps.zscore > -50) %>%
  ggplot() + 
  geom_density_ridges(aes(x=bsstraps.zscore,y=as.factor(maf),fill=lost),alpha=0.5) + 
  geom_point(aes(x=bsstraps.zscore,y=as.factor(maf),col=lost),pch="|") + 
  theme_custom() + 
  geom_vline(xintercept=c(2,-2),linetype="dashed") +
  xlab("Nodal Bootstrap Support z-score") +
  ylab("Minor Allele Count") +
  facet_grid(~int,scales="free_x")  +
  scale_fill_manual(values=c("#0A9396","#9B2226","gray80"),aesthetics = c("color","fill"),
                    labels=c("Retained","Lost")) +
  theme(legend.position=c(0.5,0.9),
        legend.text = element_text(size=rel(1.5)))
ggarrange(null.node.plot,null.bs.plot,
          labels="AUTO", ncol=1,
          font.label = list(family="Open Sans",size=25))

null.node.plot.miss <- null_lostnodes %>% 
  ggplot() + 
  geom_density_ridges(aes(x=nodestraps.zscore,y=as.factor(missing),fill=lost),alpha=0.5) + 
  geom_point(aes(x=nodestraps.zscore,y=as.factor(missing),col=lost),pch="|") +
  theme_custom() + 
  geom_vline(xintercept=c(2,-2),linetype="dashed") +
  xlab("Node Descendants z-score") +
  ylab("Missing Data") +
  facet_grid(~int,scales="free_x") +
  scale_fill_manual(values=c("#0A9396","#9B2226","gray80"),aesthetics = c("color","fill"),
                    labels=c("Retained","Lost")) +
  theme(legend.position=c(0.2,0.9),
        legend.text = element_text(size=rel(1.5)))

null.bs.plot.miss <- null_lostnodes %>% 
  filter(!is.na(mean_bs) & bsstraps.zscore > -50) %>%
  ggplot() + 
  geom_density_ridges(aes(x=bsstraps.zscore,y=as.factor(missing),fill=lost),alpha=0.5) + 
  geom_point(aes(x=bsstraps.zscore,y=as.factor(missing),col=lost),pch="|") + 
  theme_custom() + 
  geom_vline(xintercept=c(2,-2),linetype="dashed") +
  xlab("Nodal Bootstrap Support z-score") +
  ylab("Missing Data") +
  facet_grid(~int,scales="free_x")  +
  scale_fill_manual(values=c("#0A9396","#9B2226","gray80"),aesthetics = c("color","fill"),
                    labels=c("Retained","Lost")) +
  theme(legend.position=c(0.2,0.9),
        legend.text = element_text(size=rel(1.5)))

ggarrange(null.node.plot,null.bs.plot,
          labels="AUTO", ncol=1,
          font.label = list(family="Open Sans",size=25))
ggarrange(null.node.plot.miss,null.bs.plot.miss,
          labels="AUTO", ncol=1,
          font.label = list(family="Open Sans",size=25))


# lostnode_ci <- lostnode.info2 %>%
#   filter(!is.na(bs_support)) %>%
#   group_by(lost,simulation,height,missing,maf,int) %>%
#   summarize(n_nodes = n(),
#             mean_desc = mean(descendants,na.rm=T),
#             mean_bs = mean(bs_support,na.rm=T)) %>%
#   #filter(lost == TRUE) %>%
#   group_by(height,maf,simulation,lost) %>%
#   summarize(mean_n = mean(n_nodes),
#             total_n = sum(n_nodes),
#             mean_desc = mean(mean_desc),
#             mean_bs = mean(mean_bs,na.rm=T))
# 
# all_bs <- tibble(height = character(),
#                  maf = integer(),
#                  sim = integer(),
#                  nnode = numeric(),
#                  lost = logical(),
#                  mean_bs_desc = numeric(),
#                  mean_bs_boot = numeric(),
#                  zscore_desc = numeric(),
#                  zscore_boot = numeric())
# for (i in 1:nrow(lostnode_ci)) {
#   nnode <- ceiling(lostnode_ci$mean_n[i])
#   height <- lostnode_ci$height[i]
#   maf <- lostnode_ci$maf[i]
#   sim <- lostnode_ci$simulation[i]
#   lost <- lostnode_ci$lost[i]
#     
#   sp.tree <- ml.tree[[which(ml.tree.info$simulation == sim &
#                               ml.tree.info$height == height)]]
#   sp.desc <- sapply(seq(1:sp.tree$Nnode)+length(sp.tree$tip.label), function(x) length(Descendants(root(sp.tree,"sim_0_0_0",resolve.root=TRUE),x,"tips")[[1]]))
#   
#   bs <- replicate(100,sample(sp.desc,nnode,replace=TRUE))
#   bs_mean <- if (is.null(nrow(bs))) {
#     bs 
#   } else {
#     colMeans(bs)
#   }
#   
#   orig.tree <- raxml.trees2[[which(results.raxml$simulation == sim &
#                               results.raxml$height == height &
#                               results.raxml$int == int &
#                               results.raxml$maf == 0 & results.raxml$missing ==0)]]
#   bs_boot <- replicate(100,sample(as.numeric(orig.tree$node.label),nnode,replace=TRUE))
#   bs_mean_boot <- if (is.null(nrow(bs_boot))) {
#     bs_boot 
#   } else {
#     colMeans(bs_boot,na.rm=T)
#   }
#   
#   dat_bs <- tibble(height = height,
#                    maf = maf,
#                    sim = sim,
#                    lost = lost,
#                    nnode = nnode,
#                    mean_bs_desc = mean(bs_mean,na.rm=T),
#                    mean_bs_boot = mean(bs_mean_boot,na.rm=T),
#                    zscore_desc = (lostnode_ci$mean_desc[i] - mean(bs_mean)) / sd(bs_mean),
#                    zscore_boot = (lostnode_ci$mean_bs[i] - mean(bs_mean_boot,na.rm=T)) / sd(bs_mean_boot,na.rm=T))
#   all_bs <- all_bs %>%
#       add_row(dat_bs)
#   
# }
# 
# all_bs %>% 
#   ggplot() + 
#   geom_density_ridges(aes(x=zscore_desc,y=as.factor(maf),fill=lost),alpha=0.5) + 
#   theme_custom() + 
#   geom_vline(xintercept=c(2,-2),linetype="dashed") +
#   xlab("Node Descendants z-score") +
#   ylab("Minor Allele Count")
# 
# lostnode_ci_miss <- null_lostnodes %>%
#   filter(lost == TRUE) %>%
#   group_by(height,missing,simulation) %>%
#   summarize(mean_desc_lcl = mean(bs_desc_lcl),
#             mean_desc_ucl = mean(bs_desc_ucl),
#             mean_bs_lcl = mean(bs_boot_lcl),
#             mean_bs_ucl = mean(bs_boot_ucl),
#             mean_n = mean(n_nodes),
#             total_n = sum(n_nodes),
#             mean_desc = mean(mean_desc),
#             mean_bs = mean(mean_bs,na.rm=T))
# all_bs_miss <- tibble(height = character(),
#                       missing = numeric(),
#                       sim = integer(),
#                       nnode = numeric(),
#                       mean_bs_desc = numeric(),
#                       mean_bs_boot = numeric())
# for (i in 1:nrow(lostnode_ci_miss)) {
#   nnode <- ceiling(lostnode_ci_miss$mean_n[i])
#   height <- lostnode_ci_miss$height[i]
#   missing <- lostnode_ci_miss$missing[i]
#   sim <- lostnode_ci_miss$simulation[i]
#   
#   sp.tree <- ml.tree[[which(ml.tree.info$simulation == sim &
#                               ml.tree.info$height == height)]]
#   sp.desc <- sapply(seq(1:sp.tree$Nnode)+length(sp.tree$tip.label), function(x) length(Descendants(root(sp.tree,"sim_0_0_0",resolve.root=TRUE),x,"tips")[[1]]))
#   
#   bs <- replicate(100,sample(sp.desc,nnode,replace=FALSE))
#   bs_mean <- if (is.null(nrow(bs))) {
#     bs 
#   } else {
#     colMeans(bs)
#   }
#   
#   bs_boot <- replicate(100,sample(lostnode.info2$bs_support,nnode,replace=FALSE))
#   bs_mean_boot <- if (is.null(nrow(bs_boot))) {
#     bs_boot 
#   } else {
#     colMeans(bs_boot)
#   }
#   
#   dat_bs <- tibble(height = height,
#                    missing = missing,
#                    sim = sim,
#                    nnode = nnode,
#                    mean_bs_desc = bs_mean,
#                    mean_bs_boot = bs_mean_boot)
#   all_bs_miss <- all_bs_miss %>%
#     add_row(dat_bs)
# }

# all_bs %>% 
#   ggplot(aes(x=mean_bs_desc,y=as.factor(maf))) + 
#   geom_jitter(data=lostnode_ci,aes(x=mean_desc,y=as.factor(maf)),width=0.1,height=0.05,alpha=0.2) +
#   stat_halfeye(aes(),position=position_nudge(y=0.1),
#                adjust=0.8,normalize="groups",interval_alpha=0,point_alpha=0,alpha=0.5,
#                scale=0.8) + 
#   stat_halfeye(data=lostnode_ci,aes(x=mean_desc,y=as.factor(maf),fill=as.factor(maf)),position=position_nudge(y=0.1),
#                adjust=0.8,normalize="groups",interval_alpha=0,point_alpha=0,alpha=0.5,
#                scale=0.8) +
#   stat_summary(position=position_nudge(y=0.1), fun.data="mean_cl_boot",
#                geom="errorbar", width=0.1) + 
#   stat_summary(data=lostnode_ci,aes(x=mean_desc,y=as.factor(maf),col=as.factor(maf)), 
#                fun.data="mean_cl_boot", size=1.2) + 
#   theme_custom() +
#   theme(panel.grid.major = element_line(color="gray80",linetype="dashed"),
#         legend.position = "none") +
#   xlab("Mean Node Descendants") +
#   ylab("Minor Alelle Count Threshold") +
#   facet_wrap(~height)
# 
# all_bs_miss %>% 
#   ggplot(aes(x=mean_bs_boot,y=as.factor(missing))) + 
#   geom_jitter(data=lostnode_ci_miss,aes(x=mean_bs,y=as.factor(missing)),width=0.1,height=0.05,alpha=0.2) +
#   stat_halfeye(aes(),position=position_nudge(y=0.1),
#                adjust=0.8,normalize="groups",interval_alpha=0,point_alpha=0,alpha=0.5,
#                scale=0.8) + 
#   stat_halfeye(data=lostnode_ci_miss,aes(x=mean_bs,y=as.factor(missing),fill=as.factor(missing)),position=position_nudge(y=0.1),
#                adjust=0.8,normalize="groups",interval_alpha=0,point_alpha=0,alpha=0.5,
#                scale=0.8) +
#   stat_summary(position=position_nudge(y=0.1), fun.data="mean_cl_boot",
#                geom="errorbar", width=0.1) + 
#   stat_summary(data=lostnode_ci_miss,aes(x=mean_bs,y=as.factor(missing),col=as.factor(missing)), 
#                fun.data="mean_cl_boot", size=1.2) + 
#   theme_custom() +
#   theme(panel.grid.major = element_line(color="gray80",linetype="dashed"),
#         legend.position = "none") +
#   xlab("Mean Node Descendants") +
#   ylab("Missing Data Threshold") +
#   facet_wrap(~height)

# null_lostnodes %>%
#   ggplot(aes(x=bs_desc_lcl,y=as.factor(maf)),alpha=0.5) +
#   geom_density_ridges() +
#   geom_density_ridges(aes(x=bs_desc_ucl,y=as.factor(maf)),fill="turquoise",alpha=0.5) +
#   facet_wrap(~height)

lostnode.info2 %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=maf,y=bs_support)) +
  stat_summary(aes(color=height),alpha=0.8,size=1.5,fun.data = "mean_se") + # basic nonparametric bootstrap for obtaining confidence limits for the population mean without assuming normality. 
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=mean_supp,col=as.factor(maf),group=height),lty=2) +
  #geom_point(alpha=0.3) +
  theme_custom() +
  geom_errorbar(data=null_lostnodes, aes(x=maf,ymin=bsstraps.lcl,ymax=bsstraps.ucl,col=height), inherit.aes=FALSE, width=0.2) +
  #facet_grid(height~int,scales="free") +
  xlab("Minor Allele Count") +
  ylab("Node Descendants") +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position = "none")+
  scale_color_manual(labels=c("Lost"),
                     values=c("#0A9396","#9B2226","gray80")) +
  scale_x_continuous(breaks = c(seq(0, 5, by = 1),10)) +
  facet_wrap(~height)

p1 <- lostnode.info2 %>% 
  #filter(missing == 0) %>%
  #group_by(tree_num,maf,lost) %>%
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=bs_support)) +
  stat_summary(aes(color=lost),alpha=0.8,size=1.5,fun.data = "mean_se") + # basic nonparametric bootstrap for obtaining confidence limits for the population mean without assuming normality. 
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=mean_supp,col=as.factor(maf),group=height),lty=2) +
  theme_custom() +
  geom_hline(yintercept = mean(lostnode.info2$bs_support,na.rm=T),lty=2) +
  # stat_summary(data=all_bs,aes(x=maf,y=mean_bs_boot),fun.data="mean_qi",
  #              geom="errorbar", width=0.1) + 
  facet_grid(~int,scales="free") +
  xlab("Minor Allele Count") +
  ylab("Node Bootstrap Support") +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position = "none")+
  scale_color_manual(values=c("#0A9396","#9B2226","gray80")) +
  scale_x_continuous(breaks = c(seq(0, 5, by = 1),10))
p2 <- lostnode.info2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=descendants)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=1.5,fun.data = "mean_se") +
  # stat_summary(data=all_bs,aes(x=maf,y=mean_bs_desc),fun.data="mean_cl_boot",
  #              geom="errorbar", width=0.1) + 
  #geom_density(alpha=0.5) +
  theme_custom() +
  geom_hline(yintercept = mean(lostnode.info2$descendants,na.rm=T),lty=2) +
  facet_grid(~int) +
  #coord_cartesian(xlim=c(0,100)) +
  ylab("Mean Node Descendants") +
  xlab("Minor Allele Count") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position="none")+
  scale_color_manual(values=c("#0A9396","#9B2226","gray80")) +
  scale_x_continuous(breaks = c(seq(0, 5, by = 1),10))
p1.miss <- lostnode.info2 %>% 
  #filter(missing == 0) %>%
  #group_by(tree_num,maf,lost) %>%
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=missing,y=bs_support)) +
  stat_summary(aes(color=lost),alpha=0.8,size=1.5,fun.data = "mean_se") + # basic nonparametric bootstrap for obtaining confidence limits for the population mean without assuming normality. 
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=mean_supp,col=as.factor(maf),group=height),lty=2) +
  theme_custom() +
  geom_hline(yintercept = mean(lostnode.info2$bs_support,na.rm=T),lty=2) +
  # stat_summary(data=all_bs_miss,aes(x=missing,y=mean_bs_boot),fun.data="mean_cl_boot",
  #              geom="errorbar", width=0.1) + 
  facet_grid(~int,scales="free") +
  xlab("Missing Data") +
  ylab("Node Bootstrap Support") +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position = "none")+
  scale_color_manual(values=c("#0A9396","#9B2226","gray80")) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,0.9))
p2.miss <- lostnode.info2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=missing,y=descendants)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=1.5,fun.data = "mean_se") +
  # stat_summary(data=all_bs_miss,aes(x=missing,y=mean_bs_desc),fun.data="mean_cl_boot",
  #              geom="errorbar", width=0.1) + 
  #geom_density(alpha=0.5) +
  theme_custom() +
  geom_hline(yintercept = mean(lostnode.info2$descendants,na.rm=T),lty=2) +
  facet_grid(~int) +
  #coord_cartesian(xlim=c(0,100)) +
  ylab("Mean Node Descendants") +
  xlab("Missing Data") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position="none")+
  scale_color_manual(values=c("#0A9396","#9B2226","gray80")) +
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,0.9))
p2.2 <- lostnode.info2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=tr1)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=1.5) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  #facet_grid(~int) +
  #coord_cartesian(xlim=c(0,100)) +
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(values=c("#0A9396","#9B2226")) +
  scale_x_continuous(breaks = c(seq(0, 5, by = 1),10))
ggarrange(p1,p2,common.legend=TRUE, nrow=2)
ggarrange(p1.miss,p2.miss,common.legend=TRUE, nrow=2)
ggarrange(p1,p1.miss,p2,p2.miss,common.legend=TRUE)

p3 <- lostnode.summary %>%  
  #mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_supp,fill=maf,col=maf)) +
  geom_density(alpha=0.5) +
  # geom_vline(data=lostnode.summary %>% filter(lost == FALSE),
  #            aes(xintercept=mean_supp,col=maf,group=height),lty=2) +
  theme_custom() +
  facet_grid(~height) +
  #coord_cartesian(xlim=c(50,100)) +
  xlab("Mean Node BS Support") +
  theme(strip.text = element_text(size=rel(1.2)))
p4 <- lostnode.info2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_desc,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(0,100)) +
  xlab("Mean Node Descendants") +
  geom_vline(data=lostnode.summary %>% filter(lost == FALSE),
             aes(xintercept=med_mean_desc,col=missing,group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))
p.density <- ggarrange(p1,p2,nrow=1,
          common.legend=TRUE,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          legend="none")

## saved at 1200x900px
ggarrange(p2,null.node.plot,p1,null.bs.plot,
          labels=c("A","","B",""), ncol=2, nrow=2, widths=c(1,1.6),
          font.label = list(family="Open Sans",size=25))

## saved at 1200x900px
ggarrange(p2.miss,null.node.plot.miss,p1.miss,null.bs.plot.miss,
          labels=c("A","","B",""), ncol=2, nrow=2, widths=c(1,1.6),
          font.label = list(family="Open Sans",size=25))
ggarrange(p2,null.node.plot,p1,null.bs.plot,
          p2.miss,null.node.plot.miss,p1.miss,null.bs.plot.miss,
          labels="AUTO", ncol=2, nrow=4, widths=c(1,1.2,1,1.2),
          font.label = list(family="Open Sans",size=25), common.legend=TRUE,
          legend="right",legend.grob=get_legend(null.node.plot))


lostnode.info2 %>%
  filter(maf !=0 ) %>%
  ggplot(aes(x=maf,y=descendants)) +
  stat_halfeye(data= . %>% filter(lost==TRUE),aes(fill=height,col=height),
               slab_alpha=0.6,interval_alpha=0,adjust=0.6,trim=TRUE,
               position=position_nudge(x=0.03)) +
  #geom_jitter(data= . %>% filter(lost==TRUE),aes(fill=height,col=height)) +
  # stat_halfeye(data= . %>% filter(lost==FALSE),aes(fill=height,col=height),
  #              slab_alpha=0.6,interval_alpha=0,side="left",adjust=0.6,trim=TRUE,
  #              position=position_nudge(x=-0.03)) +
  theme_custom() +
  coord_flip()

p.bars2 <- lostnode.info2 %>%
  #filter(height != "LONG") %>%
  #filter(maf != 0) %>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  group_by(lost,maf,desc_category) %>%
  summarize(sum_num = n()) %>% 
  group_by(desc_category,maf) %>% 
  mutate(prop=sum_num/sum(sum_num)) %>%
  ggplot(aes(x=maf)) +
  #stat_halfeye(aes(fill=desc_category,col=desc_category),slab_alpha=0.6,interval_alpha=0,adjust=1.2,trim=FALSE) +
  geom_bar(data=. %>% filter(lost == FALSE), aes(y=sum_num, fill=desc_category), position="dodge", stat="identity", alpha=0.7) +
  geom_bar(data=. %>% filter(lost == TRUE), aes(y=-sum_num, fill=desc_category), position="dodge", stat="identity", alpha=0.7) +
  geom_hline(yintercept=0) +
  theme_custom() +
  #ylim(-2800,11000) +
  ylab("Number of Nodes") +
  xlab("Minor Alelle Count") +
  annotate(geom="text",x=0.5,y=-4000,label="LOST",hjust=0,size=6,family="Open Sans Light") +
  annotate(geom="text",x=0.5,y=10500,label="RETAINED",hjust=0,size=6, family="Open Sans Light") +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position=c(0.7,0.4),
        legend.direction="horizontal") +
  geom_text(data=. %>% filter(lost == TRUE), 
            aes(x=maf, y=-sum_num, label=paste0(round(prop*100,2),"%")), col="black",
            position=position_nudge(x=c(-0.4,-0.1,0.1,0.3),y=-200),angle=-60,hjust=0,
            family="Open Sans") +
  #facet_grid(~lost) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8) +
  scale_color_viridis(discrete=TRUE,alpha=0.8) +
  scale_x_continuous(breaks = c(seq(0, 5, by = 1),10)) + guides(fill = guide_legend(nrow = 2)) 

p.bars <- lostnode.info2 %>%
  #filter(height != "LONG") %>%
  filter(missing != 0) %>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  group_by(lost,missing,desc_category,height) %>%
  summarize(sum_num = n()) %>% 
  group_by(desc_category,missing,height) %>% 
  mutate(prop=sum_num/sum(sum_num)) %>%
  complete(expand(.,height,desc_category,missing,lost),fill=list(prop=0,sum_num=0)) %>%
  ggplot(aes(x=missing)) +
  #stat_halfeye(aes(fill=desc_category,col=desc_category),slab_alpha=0.6,interval_alpha=0,adjust=1.2,trim=FALSE) +
  #geom_bar(data=. %>% filter(lost == FALSE), aes(y=sum_num, fill=desc_category), position="dodge", stat="identity", alpha=0.7) +
  geom_line(data=. %>% filter(lost == TRUE), aes(y=prop, col=desc_category), lty=3, stat="identity", alpha=0.7) +
  geom_point(data=. %>% filter(lost == TRUE), aes(y=prop, col=desc_category), pch=19, stat="identity", alpha=0.7, size=5) +
  geom_hline(yintercept=0) +
  theme_custom() +
  #ylim(-2800,11000) +
  ylab("Proportion of Nodes Lost") +
  xlab("Missing Data") +
  #annotate(geom="text",x=0.5,y=-4000,label="LOST",hjust=0,size=6,family="Open Sans Light") +
  #annotate(geom="text",x=0.5,y=10500,label="RETAINED",hjust=0,size=6, family="Open Sans Light") +
  theme(strip.text = element_text(size=rel(1.2)),
        legend.position="top",
        legend.direction="horizontal",
        legend.text = element_text(size=rel(1.4))) +
  # geom_text(data=. %>% filter(lost == TRUE), 
  #           aes(x=maf, y=-sum_num, label=paste0(round(prop*100,2),"%")), col="black",
  #           position=position_nudge(x=c(-0.4,-0.1,0.1,0.3),y=-200),angle=-60,hjust=0,
  #           family="Open Sans") +
  facet_grid(~height) +
  #scale_fill_viridis(discrete=TRUE,alpha=0.8,begin=0.1) +
  scale_color_viridis(discrete=TRUE,alpha=0.8,end=0.9) +
  scale_x_continuous(breaks = c(seq(0, 5, by = 1),10)) + guides(fill = guide_legend(nrow = 2)) 

p.bars.prior <- lostnode.info2 %>%
  #filter(height != "LONG")%>%
  mutate(desc_category=factor(desc_category,
                              levels=c("veryfew","few","moderate","many"))) %>%
  filter(maf == 0 & missing == 0) %>%
  ggplot() +
  geom_histogram(aes(x=descendants,y=..density..,fill=desc_category),bins=50) +
  #geom_histogram(aes(x=mean_supp,y=..density..),fill="black",alpha=0.5,bins=50) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8,end=0.9) +
  xlab("Mean Node Descendants") +
  theme_custom() +
  theme(legend.position="none")

p.bars.lab <- ggarrange(p.bars,p.bars.prior,labels=c("D","E"), 
          common.legend=TRUE,
          font.label=list(size=24,font.family="Open Sans"),
          legend="top",widths=c(3,1))

p.bars.lab <- p.bars +
  annotation_custom(
    ggplotGrob(p.bars.prior), 
    xmin = 5.5, xmax = 10, ymin=0.01,ymax=0.13
  )
p.bars.lab.lab <- ggarrange(p.bars.lab,labels=c("C"),
          font.label=list(size=24,font.family="Open Sans"))

# exported as 1800x1000px
combined.fig <- p.density / p.bars.lab +
  plot_layout(heights=c(1,1.5))

ggarrange(p.density,p.bars.lab,ncol=1,
          heights = c(1,1.5))
  
## modeling?
lostnode.summary2 <- lostnode.info2 %>%
  group_by(tree_num,lost,maf,missing,height,simulation,int) %>%
  summarize(mean_bs = mean(bs_support,na.rm=T),
         mean_desc = mean(descendants, na.rm=T),
         count = n(),
         mean_depth = mean(tr1))
model1<-lm(mean_bs ~ maf, data=lostnode.summary2)
model2<-lm(mean_desc ~ maf + int, data=lostnode.summary2)
summary(model1)
summary(model2)
MuMIn::r.squaredGLMM(model1)
MuMIn::r.squaredGLMM(model2)

psum1 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_desc,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +    stat_smooth(method="lm") +
    theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Descendants")
psum2 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_depth,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +
  stat_smooth(method="lm") +
  theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth")
psum3 <- lostnode.summary2 %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=mean_desc,y=mean_bs,col=as.factor(maf))) +
  geom_jitter() +
  theme_custom() 
psum4 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=count,fill=lost)) +
  geom_bar(stat="identity") +
  theme_custom() +
  ylab("Number of Nodes") +
  xlab("Minor Allele Count") +
  scale_fill_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
  
ggarrange(psum1,psum2,psum4,nrow=1,common.legend=TRUE)

lostnode.summary2 %>%
  mutate(height = recode(height,SHORT="High ILS",MED="Med ILS",LONG="Low ILS")) %>%
  #drop_na() %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=as.integer(maf),y=mean_desc,col=height)) +
  geom_jitter(width=0.05,alpha=0.5,height=0.01) +
  stat_summary(aes(x=maf,y=mean_desc,col=height),geom="line",lty=2) +
  stat_summary(aes(x=maf,y=mean_desc,fill=height),col="black",geom="point",size=5,pch=f23) +
  #geom_point(aes(x=maf,y=mean_retained_bs),col="black")+
  #facet_grid(vars(int)) +
  theme_custom() +
  ylab("Mean descendants of lost nodes") +
  xlab("Minor Allele Count") +
  theme(legend.position=c(0.9,0.9),
        legend.text=element_text(size=rel(1.5)))

desc <- sapply(seq(1:101), function(x) length(Descendants(root(sp.tree,"sim_0_0_0",resolve.root=TRUE),x,"tips")[[1]]))
desc <- as_tibble(desc) %>%
  mutate(categ = factor(case_when(desc < 3 ~ "veryfew",
                        desc > 2 & desc < 10 ~ "few",
                        desc > 40 ~ "many",
                        TRUE ~ "moderate"),levels=c("veryfew","few","moderate","many")))
ggtree(sp.tree) +
  geom_nodepoint(aes(col=desc$categ),size=6) +
  scale_color_viridis(discrete=TRUE,alpha=0.8,end=0.9) +
  theme(text = element_text(family="Open Sans Light", size=32),
        legend.title=element_blank(),
        legend.position=c(0.2,0.85))

lostnode.info2 %>%
  filter(maf == 0 & missing == 0 & int ==  "INT") %>%
  ggplot(aes(x=descendants,y=bs_support)) +
  geom_jitter(aes(col=desc_category), alpha=0.5, size=5, pch=19) +
  theme_custom() +
  scale_color_viridis(discrete=TRUE,alpha=0.8,end=0.9) +
  theme(legend.position = "none") +
  xlab("Number of Node Descendants") +
  ylab("Node Bootstrap Support")
  


## hex plots of descendants vs bs support

hexp.int <- lostnode.info2 %>% 
  filter(simulation < 26) %>% 
  filter(int == "INT") %>% 
  ggplot(aes(x=descendants,y=bs_support)) + 
    #geom_hex(bins=30) + 
  geom_jitter(alpha=0.3) +
    theme_custom() + 
    theme(panel.grid.major=element_line(color="gray90",linetype="dashed"),
          strip.background = element_rect(fill="white"),
          strip.text = element_text(size=rel(1.5)),
          legend.title=element_text(size=rel(1.2))) + 
    facet_grid(rows=vars(lost),cols=vars(maf), 
               labeller = labeller(lost = c(`FALSE` = "RETAINED", `TRUE` = "LOST"),
                                   maf = c(`0` = "MAC = 0", `1` = "MAC = 1", 
                                           `2` = "MAC = 2", `3` = "MAC = 3", 
                                           `4` = "MAC = 4", `5` = "MAC = 5", 
                                           `10` = "MAC = 10"))) + 
    scale_fill_continuous(type = "viridis",name="Number\nof Nodes") + 
    ylab("Bootstrap Support") + 
    xlab("Number of Descendant Tips")
hexp.ext <- lostnode.info2 %>% 
  filter(simulation < 26) %>% 
  filter(int == "EXT") %>% 
  ggplot(aes(x=descendants,y=bs_support)) + 
  #geom_hex(bins=30) + 
  geom_jitter(alpha=30) +
  theme_custom() + 
  theme(panel.grid.major=element_line(color="gray90",linetype="dashed"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1.2))) + 
  facet_grid(rows=vars(lost),cols=vars(maf), 
             labeller = labeller(lost = c(`FALSE` = "RETAINED", `TRUE` = "LOST"),
                                 maf = c(`0` = "MAC = 0", `1` = "MAC = 1", 
                                         `2` = "MAC = 2", `3` = "MAC = 3", 
                                         `4` = "MAC = 4", `5` = "MAC = 5", 
                                         `10` = "MAC = 10"))) + 
  scale_fill_continuous(type = "viridis",name="Number\nof Nodes") + 
  ylab("Bootstrap Support") + 
  xlab("Number of Descendant Tips")
hexp.both <- lostnode.info2 %>% 
  #filter(simulation < 26) %>% 
  ggplot(aes(x=descendants,y=bs_support)) + 
  #geom_hex(bins=30) + 
  geom_jitter(alpha=0.3) +
  theme_custom() + 
  theme(panel.grid.major=element_line(color="gray90",linetype="dashed"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=rel(1.5)),
        legend.title=element_text(size=rel(1.2))) + 
  facet_grid(rows=vars(lost),cols=vars(maf), 
             labeller = labeller(lost = c(`FALSE` = "RETAINED", `TRUE` = "LOST"),
                                 maf = c(`0` = "MAC = 0", `1` = "MAC = 1", 
                                         `2` = "MAC = 2", `3` = "MAC = 3", 
                                         `4` = "MAC = 4", `5` = "MAC = 5", 
                                         `10` = "MAC = 10"))) + 
  scale_fill_continuous(type = "viridis",name="Number\nof Nodes") + 
  ylab("Bootstrap Support") + 
  xlab("Number of Descendant Tips")

ggarrange(hexp.int,hexp.ext,ncol=1,labels="AUTO",
          common.legend=TRUE,
          font.label=list(size=24,font.family="Open Sans"),
          legend="right")

##########################
## empirical trees----
##########################
# lostnode.info <- tibble(tree_num=integer(), desc_category=character(),lost = logical(), mean_desc=numeric(), mean_supp=numeric(),
#                         num=integer(), mean_nodedepth=numeric(), median_nodedepth=numeric())
lostnode.info.emp <- tibble(tree_num=integer(), tr1=integer(), tr2=integer(),
                        lost=logical(), descendants=integer(), bs_support=numeric(),
                        desc_category=character())
for (i in 1:nrow(results.raxml)){
  sum <- lostnodes.emp(i)
  lostnode.info.emp <- lostnode.info.emp %>%
    add_row(sum)
}

lostnode.summary.emp <- lostnode.info.emp %>%
  drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num")) %>%
  #filter(maf !=0) %>%
  # mutate(missing=as.factor(missing)) %>%
  group_by(lost,maf) %>%
  summarize(mean_supp = mean(bs_support,na.rm=T),
            med_supp = median(bs_support,na.rm=T),
            mean_desc = mean(descendants,na.rm=T),
            med_desc = median(descendants,na.rm=T),
            #mean_mean_nodedepth = mean(mean_nodedepth,na.rm=T),
            #med_mean_nodedepth = median(mean_nodedepth,na.rm=T),
            #mean_med_nodedepth = mean(median_nodedepth,na.rm=T),
            #med_med_nodedepth = median(median_nodedepth, na.rm=T)
  ) %>%
  ungroup()
lostnode.info.emp2 <- lostnode.info.emp %>%
  #drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num"))

p1.emp <- lostnode.info.emp2 %>% 
  #group_by(tree_num,maf,lost) %>%
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=bs_support)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=mean_supp,col=as.factor(maf),group=height),lty=2) +
  theme_custom() +
  #facet_grid(~height,scales="free") +
  #coord_cartesian(ylim=c(99,100)) +
  xlab("Minor Allele Count") +
  ylab("Node Bootstrap Support") +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
p2.emp <- lostnode.info.emp2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=descendants)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  #facet_grid(~height) +
  #coord_cartesian(xlim=c(0,100)) +
  ylab("Mean Node Descendants") +
  xlab("Minor Allele Count") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
p2.2.emp <- lostnode.info.emp2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=tr1)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  #facet_grid(~height) +
  #coord_cartesian(xlim=c(0,100)) +
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
ggarrange(p1.emp,p2.emp,p2.2.emp, common.legend=TRUE, nrow=1)

p3.emp <- lostnode.info.emp2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_supp,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  geom_vline(data=lostnode.summary.emp %>% filter(lost == FALSE),
             aes(xintercept=mean_mean_supp,col=missing,group=height),lty=2) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(50,100)) +
  xlab("Mean Node BS Support") +
  theme(strip.text = element_text(size=rel(1.2)))
p4.emp <- lostnode.info.emp2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_desc,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(0,100)) +
  xlab("Mean Node Descendants") +
  geom_vline(data=lostnode.summary.emp %>% filter(lost == FALSE),
             aes(xintercept=med_mean_desc,col=missing,group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))
p.density.emp <- ggarrange(p1.emp,p2.emp,p2.2.emp,nrow=1,
                       common.legend=TRUE,
                       labels="AUTO",
                       font.label=list(size=24,font.family="Open Sans"),
                       legend="right")

lostnode.info.emp2 %>%
  filter(maf !=0 ) %>%
  ggplot(aes(x=maf,y=descendants)) +
  stat_halfeye(data= . %>% filter(lost==TRUE),
               slab_alpha=0.6,interval_alpha=0,adjust=0.6,trim=TRUE,
               position=position_nudge(x=0.03)) +
  #geom_jitter(data= . %>% filter(lost==TRUE),aes(fill=height,col=height)) +
  # stat_halfeye(data= . %>% filter(lost==FALSE),aes(fill=height,col=height),
  #              slab_alpha=0.6,interval_alpha=0,side="left",adjust=0.6,trim=TRUE,
  #              position=position_nudge(x=-0.03)) +
  theme_custom() +
  coord_flip()

p.bars.emp <- lostnode.info.emp2 %>%
  #filter(height != "LONG") %>%
  filter(maf != 0) %>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  group_by(lost,maf,desc_category) %>%
  summarize(sum_num = n()) %>% 
  #group_by(desc_category) %>% 
  #mutate(total_categ = sum(sum_num)) %>% ungroup() %>% mutate(prop=sum_num/total_categ) %>%
  ggplot(aes(x=maf)) +
  #stat_halfeye(aes(fill=desc_category,col=desc_category),slab_alpha=0.6,interval_alpha=0,adjust=1.2,trim=FALSE) +
  geom_bar(data=. %>% filter(lost == FALSE), aes(y=sum_num, fill=desc_category), position="dodge", stat="identity") +
  geom_bar(data=. %>% filter(lost == TRUE), aes(y=-sum_num, fill=desc_category), position="dodge", stat="identity") +
  geom_hline(yintercept=0) +
  theme_custom() +
  #ylim(-2800,11000) +
  ylab("Number of Nodes") +
  xlab("Minor Alelle Count") +
  annotate(geom="text",x=0.5,y=-1300,label="LOST",hjust=0,size=6,family="Open Sans Light") +
  annotate(geom="text",x=0.5,y=1200,label="RETAINED",hjust=0,size=6, family="Open Sans Light") +
  theme(strip.text = element_text(size=rel(1.2))) +
  #facet_grid(~lost) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8)

p.bars.prior.emp <- lostnode.info.emp2 %>%
  #filter(height != "LONG")%>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  filter(maf == 0 & missing == 0 & lost == FALSE) %>%
  ggplot() +
  geom_histogram(aes(x=descendants,y=..density..,fill=desc_category),bins=50) +
  #geom_histogram(aes(x=mean_supp,y=..density..),fill="black",alpha=0.5,bins=50) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8) +
  xlab("Mean Node Descendants") +
  theme_custom() +
  theme(legend.position="none")

p.bars.lab.emp <- ggarrange(p.bars.emp,p.bars.prior.emp,labels=c("D","E"), 
                        common.legend=TRUE,
                        font.label=list(size=24,font.family="Open Sans"),
                        legend="right",widths=c(3,1)) 

p.bars.lab.emp <- p.bars.emp +
  annotation_custom(
    ggplotGrob(p.bars.prior.emp), 
    xmin = 5.5, xmax = 9, ymin=3000,ymax=11000
  )


#library(patchwork)
p.density.emp / p.bars.lab.emp +
  plot_layout(heights=c(1,1.5))

## modeling?
lostnode.summary2 <- lostnode.info.emp2 %>%
  group_by(tree_num,lost,maf,missing,simulation,int) %>%
  summarize(mean_bs = mean(bs_support,na.rm=T),
            mean_desc = mean(descendants, na.rm=T),
            count = n(),
            mean_depth = mean(tr1))
model1<-lm(mean_bs ~ maf, data=lostnode.summary2)
model2<-lm(mean_desc ~ maf + int, data=lostnode.summary2)
summary(model1)
summary(model2)
MuMIn::r.squaredGLMM(model1)
MuMIn::r.squaredGLMM(model2)
plot(model2)

psum1 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_desc,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +    stat_smooth(method="lm") +
  theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Descendants")
psum2 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_depth,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +
  stat_smooth(method="lm") +
  theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth")
psum3 <- lostnode.summary2 %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=mean_desc,y=mean_bs,col=as.factor(maf))) +
  geom_jitter() +
  theme_custom() 
psum4 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=count,fill=lost)) +
  geom_bar(stat="identity") +
  theme_custom() +
  ylab("Number of Nodes") +
  xlab("Minor Allele Count") +
  scale_fill_manual(labels=c("Retained","Lost"),
                    values=c("#0A9396","#9B2226"))

ggarrange(psum1,psum2,psum4,nrow=1,common.legend=TRUE)

lostnode.summary2 %>%
  drop_na() %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=as.integer(maf),y=mean_desc)) +
  geom_jitter(width=0.001) +
  stat_summary(aes(x=maf,y=mean_desc),col="black",geom="line",lty=2) +
  stat_summary(aes(x=maf,y=mean_desc),col="black",geom="point",size=5,pch=4) +
  #geom_point(aes(x=maf,y=mean_retained_bs),col="black")+
  facet_grid(vars(int)) +
  theme_custom() +
  ylab("Mean descendants of lost nodes") +
  xlab("Minor Allele Count") +
  theme(legend.position="none")
