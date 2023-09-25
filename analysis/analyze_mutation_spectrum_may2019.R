###############################
## mutation spectrum analysis
## fall 2021
## j. rick
##
## updated feb 2022
###############################

library(tidyverse)
library(ggdist)
library(PNWColors)
library(here)

here::i_am("analysis/analyze_mutation_spectrum_may2019.R")
source(here("analysis","theme_custom.R"))
source(here("analysis","pairwise_ks_test.R"))
pal=rev(pnw_palette("Sunset",100))


## Pulling info on simulated number of sites
mut.files <- list.files(path=here("output","new"),pattern="092321-varSites*")

mut.all <- data.frame(gene=integer(),
                      variants=integer(),
                      sim_num=integer(),
                      height=integer(),
                      int_ext=character(),
                      stringsAsFactors = FALSE)

for (i in 1:length(mut.files)){
  varsites <- scan(paste("output/new/",mut.files[i],sep=""))
  sim <- regmatches(mut.files[i], regexec('sim([0-9]+)', mut.files[i]))[[1]][2]
  height <- regmatches(mut.files[i], regexec('-varSites-([A-Z]+)\\-sim', mut.files[i]))[[1]][2]
  int <- as.character(regmatches(mut.files[i], regexec('-sim[0-9]+-([A-Z]+)', mut.files[i]))[[1]][2])
  mut.df <- data.frame(gene=seq(1,length(varsites[-1])),
                       variants=varsites[-1],
                       sim_num=sim,
                       tree_height=height,
                       int_ext=int)
  mut.all <- rbind(mut.all,mut.df)
}

mut.all$gene_name <- paste0("gene",sprintf("%04d",mut.all$gene))
mut.all$sim_num <- as.integer(as.character(mut.all$sim_num))

#mut.all$tree_height <- as.factor(mut.all$tree_height)
ggpubr::ggdensity(data=mut.all[mut.all$sim_num > 15 & mut.all$sim_num < 26,],
                  x="variants",
                  color="tree_height",
                  fill="tree_height",
                  xlab="Simulated Variant Sites",
                  ylab="Frequency",
                  add="mean",rug=TRUE)

## Pulling information for post-filtering SNPs

var.list <- list.files(path="../../results/results_092321/",pattern="092321-SNPs*")

varsites.all <- data.frame()

for (i in 1:length(var.list)){
  varsites <- read_csv(paste("../../results/results_092321/",var.list[i],sep=""),
                         #row.names=NULL,
                         col_names=c("gene","iteration","nSNP"),
                         #stringsAsFactors = FALSE,
                         na=c("","NA")) %>%
    mutate(height = regmatches(var.list[i], regexec('SNPs-([A-Z]+)\\-', var.list[i]))[[1]][2],
           int = regmatches(var.list[i], regexec('SNPs-[A-Z]+-([A-Z]+)', var.list[i]))[[1]][2],
           sim = gsub("^s([0-9]+)\\_.*","\\1",iteration),
           qual = gsub(".*q([0-9]+)\\_.*","\\1",iteration),
           miss = gsub(".*miss(.*)\\_mac.*","\\1",iteration),
           maf = gsub(".*mac(.*)\\.[A-Z]+\\.REF.*","\\1",iteration),
           noref = gsub(".*\\.([A-Z]+)\\..*","\\1",iteration)) %>%
    mutate(sim = as.integer(as.character(sim))) %>%
    mutate(nSNP = case_when(is.na(nSNP) ~ 0,
                            TRUE ~ nSNP))

  varsites.all <- varsites.all %>%
    bind_rows(varsites)
}

summary(varsites.all)

varsites.all %>%
  mutate(maf = factor(maf,levels=c("0","1","2","3","4","5","10"))) %>% 
  ggplot() +
  stat_halfeye(aes(x=nSNP,y=maf,fill=height,color=height), alpha=0.5,
               normalize="groups", adjust=1.2, p_limits=c(0.025,0.975),
               trim=TRUE,expand=FALSE,point_interval="mean_qi",n=500) +
  ylab("Minor Allele Count Threshold") +
  xlab("SNPs per Locus") +
  theme_custom() 

varsites.all %>%
  mutate(maf = factor(maf,levels=c("0","1","2","3","4","5","10"))) %>%
  group_by(height,maf,int) %>%
  summarize(mean_snps = mean(nSNP,na.rm=TRUE)) %>% view()

varsites.all %>%
  mutate(maf = factor(maf,levels=c("0","1","2","3","4","5","10"))) %>%
  #group_by(height,maf,int) %>%
  #summarize(mean_snps = mean(nSNP,na.rm=TRUE)) %>%
  ggplot(aes(x=maf,y=nSNP,color=height,shape=int)) +
  stat_summary(fun=mean, 
               fun.min=function(z) { quantile(z,0.25) },
               fun.max=function(z) { quantile(z,0.75)},
               position = position_dodge(width = 0.75),
               size=1) +
  theme_custom() +
  xlab("Minor Allele Count Threshold") +
  ylab("SNPs per Locus")

varsites.all %>%
  pivot_wider(names_from=int, values_from=nSNP,names_prefix = "nSNP_") %>% 
  mutate(diff_nSNP_INT_EXT = nSNP_INT - nSNP_EXT) %>%
  ggplot() +
  geom_density_ridges(aes(y=height,fill=maf,x=diff_nSNP_INT_EXT), alpha=0.1) +
  #geom_point(aes(x=nSNP_INT/2000,y=nSNP_EXT/2000,color=maf),alpha=0.1) +
  #facet_grid(cols=vars(height),rows=vars(maf)) +
  theme_custom() +
  geom_abline(slope=1,intercept=0, lty=3)
  

## Putting the data together
varsites.all$int <- as.factor(varsites.all$int)

pre.post <- varsites.all %>%
  as_tibble() %>%
  filter(sim > 15 & sim < 26) %>%
  left_join(mut.all, by=c(
    "sim" = "sim_num",
    "gene" = "gene_name",
    "height" = "tree_height",
    "int" = "int_ext"
  )) %>%
  mutate(diff = variants - nSNP)

summary(pre.post)

ggdensity(pre.post,x="diff",fill=NA,alpha=0.3,color="maf",
          add="mean")

#plot(pre.post$variants,pre.post$diff,col=pre.post$maf)

pre.post %>%
  ggdensity(x="nSNP",fill="maf")

first_loss_df <- pre.post %>%
  filter(nSNP == 0 & miss == "0") %>%
  group_by(height,sim,int,gene,variants) %>%
  summarise(first_loss = min(as.numeric(maf),na.rm=T))  %>%
  ungroup() %>%
  group_by(first_loss,int) %>%
  mutate(mean = mean(variants/5000), # each locus is 5000bp
         median = median(variants/5000))%>%
  ungroup() %>%
  group_by(first_loss) %>%
  mutate(mean_fl = mean(variants/5000))
  # left_join(mut.all, by=c(
  #   "sim" = "sim_num",
  #   "gene" = "gene_name",
  #   "height" = "tree_height",
  #   "int" = "int_ext")) 
fl_maf <- first_loss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5,adjust=0.5) +
  geom_vline(aes(xintercept=mean_fl,group=height,col=first_loss),lty=2)  +
  theme_custom() +
  theme(legend.position="right",
        strip.background=element_rect(fill="white"),
        strip.text = element_text(size=rel(1.2),family="Open Sans Light"),
        legend.title = element_text(size=rel(1),family="Open Sans Light")) +
  scale_color_viridis(discrete=TRUE, option="viridis",end=0.9) +
  ylab("Density") +
  xlab("Mutation Rate") +
  guides(col=guide_legend(title="Minor\nAllele\nCount"))
fl_maf_facet <- first_loss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5) +
  geom_vline(aes(xintercept=mean,group=int,col=first_loss),lty=2) +
  theme_custom() +
  theme(legend.position="right",
        strip.background=element_rect(fill="white"),
        strip.text = element_text(size=rel(1.2),family="Open Sans Light"),
        legend.title = element_text(size=rel(1),family="Open Sans Light")) +
  scale_color_viridis(discrete=TRUE, option="viridis",end=0.9) +
  ylab("Density") +
  xlab("Mutation Rate") +
  facet_wrap(~int) +
  guides(col=guide_legend(title="Minor\nAllele\nCount"))

first_loss_miss_df <- pre.post %>%
  filter(nSNP == 0 & maf == "0") %>%
  group_by(height,sim,int,gene,variants) %>%
  summarise(first_loss = min(as.numeric(miss),na.rm=T))  %>%
  ungroup() %>%
  group_by(first_loss,int) %>%
  mutate(mean = mean(variants/5000)) %>%
  ungroup() %>%
  group_by(first_loss) %>%
  mutate(mean_fl = mean(variants/5000))
# left_join(mut.all, by=c(
#   "sim" = "sim_num",
#   "gene" = "gene_name",
#   "height" = "tree_height",
#   "int" = "int_ext")) 

fl_miss <- first_loss_miss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  filter(!first_loss %in% c(0.25,0.5)) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5) +
  geom_vline(aes(xintercept=mean_fl,col=first_loss),lty=2) +
  theme_custom() +
  theme(legend.position="right",
        strip.background=element_rect(fill="white"),
        strip.text = element_text(size=rel(1.2),family="Open Sans Light"),
        legend.title = element_text(size=rel(1),family="Open Sans Light")) +
  scale_color_viridis(discrete=TRUE, option="viridis",end=0.8) +
  ylab("Density") +
  xlab("Mutation Rate") +
  guides(col=guide_legend(title="Missing\nData\nThreshold"))
fl_miss_facet <- first_loss_miss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000)  %>%
  filter(!first_loss %in% c(0.25,0.5)) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5) +
  geom_vline(aes(xintercept=mean,group=int,col=first_loss),lty=2) +
  theme_custom() +
  theme(legend.position="right",
        strip.background=element_rect(fill="white"),
        strip.text = element_text(size=rel(1.2),family="Open Sans Light"),
        legend.title = element_text(size=rel(1),family="Open Sans Light")) +
  scale_color_viridis(discrete=TRUE, option="viridis",end=0.8) +
  ylab("Density") +
  xlab("Mutation Rate") +
  facet_wrap(~int) +
  guides(col=guide_legend(title="Missing\nData\nThreshold"))

## plotting together!
## exported at 1200x1200px
#library(patchwork)
(fl_maf | fl_miss) /
  fl_maf_facet /
  fl_miss_facet + plot_annotation(tag_levels = 'A') + 
  plot_layout(heights = c(1.5, 1, 1)) & 
  theme(plot.tag = element_text(size = rel(2),family="Open Sans")) 

  
first_loss_miss_df %>% 
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggplot() +
  #geom_point(aes(x=first_loss,y=mut_rate,col=int),alpha=0.5,height=0) +
  # stat_halfeye(data= . %>% filter(int=="EXT"),
  #              aes(y=mut_rate,x=first_loss),orientation="vertical",adjust=0.75,alpha=0.5,side="left",
  #              normalize="groups") +
  # stat_halfeye(data = . %>% filter(int=="INT"),
  #              aes(y=mut_rate,x=first_loss),orientation="vertical",adjust=0.75,alpha=0.5,side="right") +
  geom_boxplot(aes(x=first_loss,y=mut_rate,fill=int),draw_quantiles=c(0.5)) +
  #facet_wrap(~height) +
  theme_custom()
  
pre.post %>% 
  group_by(maf,height,sim) %>%
  summarise(mean_var = mean(variants)) %>%
  ggplot() +
  geom_boxplot(aes(x=height,y=mean_var,col=maf))+
  #geom_jitter(aes(x=height,y=mean_var,col=maf),alpha=0.5,width=0.1,height=0.1) +
  theme_custom()
  
  
pre.post %>%
  mutate(lost = case_when(nSNP == 0 ~ "Y",
                          TRUE ~ "N")) %>%
  group_by(maf,lost) %>%
  #summarize(n = n())
  #t.test(.$variants ~ .$lost)
  filter(lost == "Y") %>%
  ggdensity(x="variants",fill="maf",alpha=0.3,color="maf",
            add="median")
pre.post %>%
  mutate(lost = case_when(nSNP == 0 ~ "Y",
                          TRUE ~ "N")) %>%
  ggplot() +
  #geom_point(aes(x=maf,y=diff,col=maf)) +
  #geom_jitter(aes(col=height),width=0.1,height=0,alpha=0.2) +
  #geom_density(aes(x=variants,fill=maf),alpha=0.2) +
  geom_density(aes(x=diff,fill=maf), alpha=0.2, trim=TRUE) +
  facet_wrap(~lost) +
  theme_custom()

pre.post %>%
  ggplot(aes(x=variants,y=diff)) +
  stat_bin_hex(bins=100)

pre.post.wm <- pre.post %>%
  group_by(gene,height,int,sim) %>%
  mutate(max_snp=max(nSNP),
         lost_snp=max_snp-nSNP) %>% 
  ungroup() %>% 
  mutate(snps_lost = case_when(diff < 0 ~ 0,
                               TRUE ~ diff)) %>%
  filter(maf > 0) %>%
  group_by(height,int,sim,miss,maf) %>%
  summarize(mean = mean(variants),
            weighted_mean = weighted.mean(variants,snps_lost))

mutdist.maf <- pre.post.wm %>% 
  group_by(maf) %>%
  mutate(mean_mut = mean(weighted_mean/5000)) %>%
  ggplot() +
  geom_density(aes(x = weighted_mean/5000,
                   #weight = snps_lost/sum(snps_lost),
                   col=factor(maf,levels=c(1,2,3,4,5,10)),fill=factor(maf,levels=c(1,2,3,4,5,10))),
               stat = "density",size=0.5,adjust=0.5,alpha=0.3) +
  geom_vline(aes(xintercept=mean_mut, col=maf, group=maf), lty=2, alpha=1, size=1.3) +
#  geom_rect(aes(xmin = mean_mut - 0.0001, xmax = mean_mut+0.0001,
#                ymin = 1275 - 75, ymax = 1275 + 75), fill = "white") +
#  geom_text(aes(x=mean_mut,label=round(mean_mut,4)),y=1275,angle=90,
#            size=5, family="Open Sans Light", check_overlap=TRUE) +
  theme_custom() +
  #facet_wrap(~int) +
  xlab("Locus Mutation Rate") +
  ylab("Density of lost SNPs") +
  xlim(0.035,0.045) +
  scale_y_continuous(expand=c(0.01,100)) +
  scale_color_viridis_d(direction=-1,aesthetics=c("color","fill"))

mutdist.int <- pre.post.wm %>% 
  group_by(int,maf) %>%
  mutate(mean_mut = mean(weighted_mean/5000),
         med_mut = median(weighted_mean/5000)) %>%
  ggplot() +
  geom_density(aes(x = weighted_mean/5000,
                   #weight = snps_lost/sum(snps_lost),
                   col=factor(int),fill=factor(int)),
               stat = "density",size=0.5,adjust=0.5,alpha=0.3) +
  geom_vline(aes(xintercept=mean_mut, col=int, group=int), lty=2, alpha=1, size=1.3) +
  #  geom_rect(aes(xmin = mean_mut - 0.0001, xmax = mean_mut+0.0001,
  #                ymin = 1275 - 75, ymax = 1275 + 75), fill = "white") +
  #  geom_text(aes(x=mean_mut,label=round(mean_mut,4)),y=1275,angle=90,
  #            size=5, family="Open Sans Light", check_overlap=TRUE) +
  theme_custom() +
  facet_wrap(~factor(maf,levels=c(1,2,3,4,5,10)), ncol=1) +
  xlab("Locus Mutation Rate") +
  ylab("Density of lost SNPs") +
  #xlim(0.035,0.045) +
  #scale_y_continuous(expand=c(0.01,100)) +
  scale_color_viridis_d(direction=1,aesthetics=c("color","fill"),option="plasma")

ks.maf <- pairwise_ks_test(pre.post.wm$weighted_mean,pre.post.wm$maf,alternative="two.sided") %>%
  reshape2::melt() %>%
  mutate(Var1 = as.factor(Var1),
         Var2 = as.factor(Var2)) %>%
  mutate(value2 = case_when(value == 1 ~ FALSE,
                            value != 1 ~ TRUE),
         sig = case_when(value < 0.01 ~ TRUE,
                         value > 0.01 ~ FALSE)) %>%
  ggplot(aes(x=Var1,y=Var2)) +
  geom_tile(aes(fill=value)) +
  geom_tile(col="gray90",fill=NA,lwd=1.5) +
  theme_custom()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(2))) +
  #scale_fill_gradientn(colours = pal,values=seq(1,0,by=-0.01)) +
  scale_fill_distiller(type = "seq",
                        direction = -1,
                        palette = "Greys") +
  scale_x_discrete(position = "top") +
  xlab("Minor Allele Count") +
  ylab("Minor Allele Count") +
  geom_text(data = . %>% filter(sig), aes(Var2, Var1), color = "white", size = 12, label = "*") +
  geom_segment(aes(x=0.5,y=0.5,xend=6.5,yend=6.5),color="gray90",size=1.5)

ks.int <- pairwise_ks_test(pre.post.wm$weighted_mean,pre.post.wm$int,alternative="two.sided") %>%
  reshape2::melt() %>%
  mutate(Var1 = as.factor(Var1),
         Var2 = as.factor(Var2)) %>%
  mutate(value2 = case_when(value == 1 ~ FALSE,
                            value != 1 ~ TRUE),
         sig = case_when(value < 0.01 ~ TRUE,
                         value > 0.01 ~ FALSE)) %>%
  ggplot(aes(x=Var1,y=Var2)) +
  geom_tile(aes(fill=value)) +
  geom_tile(col="gray90",fill=NA,lwd=1.5) +
  theme_custom()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(2))) +
  #scale_fill_gradientn(colours = pal,values=seq(1,0,by=-0.01)) +
  scale_fill_distiller(type = "seq",
                       direction = -1,
                       palette = "Greys") +
  scale_x_discrete(position = "top") +
  xlab("Int") +
  ylab("Int") +
  geom_text(data = . %>% filter(sig), aes(Var2, Var1), color = "white", size = 12, label = "*") +
  geom_segment(aes(x=0.5,y=0.5,xend=6.5,yend=6.5),color="gray90",size=1.5)

pre.post.wm.miss <- pre.post %>%
  group_by(gene,height,int,sim) %>%
  mutate(max_snp=max(nSNP),
         lost_snp=max_snp-nSNP) %>% 
  ungroup() %>% 
  mutate(snps_lost = case_when(diff < 0 ~ 0,
                               TRUE ~ diff)) %>%
  #filter(miss > 0) %>%
  group_by(height,int,sim,miss,maf) %>%
  summarize(mean = mean(variants),
            weighted_mean = weighted.mean(variants,snps_lost)) 

mutdist.miss <- pre.post.wm.miss %>% 
  group_by(miss) %>%
  mutate(mean_mut = mean(weighted_mean/5000)) %>%
  ggplot() +
  geom_density(aes(x = weighted_mean/5000,
                   #weight = snps_lost/sum(snps_lost),
                   col=factor(miss,levels=c(0,0.25,0.5,0.75,0.9)),fill=factor(miss,levels=c(0,0.25,0.5,0.75,0.9))),
               stat = "density",size=0.5,adjust=0.5,alpha=0.3) +
  geom_vline(aes(xintercept=mean_mut, col=miss, group=miss), lty=2, alpha=1, size=1.3) +
#  geom_rect(aes(xmin = mean_mut - 0.0001, xmax = mean_mut+0.0001,
#                ymin = 400 - 25, ymax = 400 + 25), fill = "white") +
#  geom_text(data = . %>% group_by(miss) %>% slice(1), 
#            aes(x=mean_mut,label=round(mean_mut,4)),y=400,angle=90,
#            size=5, family="Open Sans Light") +
  theme_custom() +
  #facet_wrap(~sim) +
  xlab("Locus Mutation Rate") +
  ylab("Density of lost SNPs") +
  # xlim(0.035,0.045) +
  scale_y_continuous(expand=c(0.01,10)) +
  scale_color_viridis_d(direction=-1,aesthetics=c("color","fill")) +
  xlim(0.036,0.044)

ks.miss<-pairwise_ks_test(pre.post.wm$weighted_mean,pre.post.wm$miss,alternative="two.sided") %>%
  reshape2::melt() %>%
  mutate(Var1 = as.factor(Var1),
         Var2 = as.factor(Var2)) %>%
  mutate(value2 = case_when(value == 1 ~ FALSE,
                            value != 1 ~ TRUE),
         sig = case_when(value < 0.01 ~ TRUE,
                         value > 0.01 ~ FALSE)) %>%
  ggplot(aes(x=Var1,y=Var2)) +
  geom_tile(data = . %>% filter(value2),aes(fill=value)) +
  geom_tile(col="gray90",fill=NA,lwd=1.5) +
  theme_custom()+
  theme(panel.border = element_blank(),
        axis.text = element_text(size=rel(1.5)),
        axis.title = element_text(size=rel(2))) +
  #scale_fill_gradientn(colours = pal,values=seq(1,0,by=-0.01)) +
  scale_fill_distiller(type = "seq",
                       direction = -1,
                       palette = "Greys") +
  scale_x_discrete(position = "top") +
  xlab("Missing Data Threshold") +
  ylab("Missing Data Threshold") +
  geom_text(data = . %>% filter(sig), aes(Var2, Var1), color = "white", size = 12, label = "*") +
  geom_segment(aes(x=0.5,y=0.5,xend=5.5,yend=5.5),color="gray90",size=1.5)

ggarrange(mutdist.maf,mutdist.miss,
          font.label=list(family="Open Sans",size=24),
          label.x=-0.05)

pre.post.locus <- pre.post %>%
  group_by(gene,height,int,sim) %>%
  mutate(max_snp=max(nSNP),
         lost_snp=max_snp-nSNP) %>% 
  group_by(gene,iteration) %>%
  mutate(locus_lost = case_when(lost_snp == max_snp ~ "lost",
                                lost_snp < max_snp ~ "retained"),
         mut_rate = variants/5000) 

pre.post.locus %>%
  filter(locus_lost == "lost" & max_snp > 0) %>%
  group_by(maf,int) %>%
  mutate(med_mut = median(mut_rate),
         mean_mut = mean(mut_rate)) %>%
  ggplot(aes(x=mut_rate)) +
  geom_density(aes(x = mut_rate,
                   #weight = snps_lost/sum(snps_lost),
                   col=factor(maf,levels=c(0,1,2,3,4,5,10)),fill=factor(maf,levels=c(0,1,2,3,4,5,10))),
               stat = "density",size=0.5,adjust=0.8,alpha=0.05) +
  geom_vline(aes(xintercept=mean_mut, col=maf, group=maf), lty=2, alpha=1, size=1.3) +
  # geom_rect(aes(xmin = mut_rate - 0.0001, xmax = mut_rate+0.0001,
  #               ymin = 1275 - 75, ymax = 1275 + 75), fill = "white") +
  # geom_text(aes(x=mut_rate,label=round(mut_rate,4)),y=1275,angle=90,
  #           size=5, family="Open Sans Light", check_overlap=TRUE) +
  theme_custom() +
  facet_wrap(~int) +
  xlab("Locus Mutation Rate") +
  ylab("Density of lost Loci") +
  #xlim(0.035,0.045) +
  #scale_y_continuous(expand=c(0.01,100)) +
  scale_color_viridis_d(direction=-1,aesthetics=c("color","fill"))

pre.post.locus %>%
  filter(locus_lost == "lost" & max_snp > 0) %>%
  group_by(miss,int) %>%
  mutate(med_mut = median(mut_rate),
         mean_mut = mean(mut_rate)) %>%
  ggplot(aes(x=mut_rate)) +
  geom_density(aes(x = mut_rate,
                   #weight = snps_lost/sum(snps_lost),
                   col=miss,fill=miss),
               stat = "density",size=0.5,adjust=0.8,alpha=0.05) +
  geom_vline(aes(xintercept=mean_mut, col=miss, group=miss), lty=2, alpha=1, size=1.3) +
  # geom_rect(aes(xmin = mut_rate - 0.0001, xmax = mut_rate+0.0001,
  #               ymin = 1275 - 75, ymax = 1275 + 75), fill = "white") +
  # geom_text(aes(x=mut_rate,label=round(mut_rate,4)),y=1275,angle=90,
  #           size=5, family="Open Sans Light", check_overlap=TRUE) +
  theme_custom() +
  facet_wrap(~int) +
  xlab("Locus Mutation Rate") +
  ylab("Density of lost Loci") +
  #xlim(0.035,0.045) +
  #scale_y_continuous(expand=c(0.01,100)) +
  scale_color_viridis_d(direction=-1,aesthetics=c("color","fill"))

pairwise_ks_test(filter(pre.post.locus,locus_lost == "lost" & max_snp > 0)$mut_rate,
                   filter(pre.post.locus,locus_lost == "lost" & max_snp > 0)$maf,
                   alternative="less")

pairwise_ks_test(filter(pre.post.locus,locus_lost == "lost" & max_snp > 0 & int == "INT")$mut_rate,
                 filter(pre.post.locus,locus_lost == "lost" & max_snp > 0 & int == "INT")$maf,
                 alternative="less")
pairwise_ks_test(filter(pre.post.locus,locus_lost == "lost" & max_snp > 0 & int == "EXT")$mut_rate,
                 filter(pre.post.locus,locus_lost == "lost" & max_snp > 0 & int == "EXT")$maf,
                 alternative="less")
