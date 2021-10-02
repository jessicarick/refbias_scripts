library(tidyverse)
library(ggdist)

## Pulling info on simulated number of sites
mut.files <- list.files(path="output/new/",pattern="092321-varSites*")

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
ggpubr::ggdensity(data=mut.all[mut.all$sim_num > 15,],
                  x="variants",
                  color="tree_height",
                  fill="tree_height",
                  xlab="Simulated Variant Sites",
                  ylab="Frequency",
                  add="mean",rug=TRUE)

## Pulling information for post-filtering SNPs

var.list <- list.files(path="output/new/",pattern="092321-SNPs*")

varsites.all <- data.frame()

for (i in 1:length(var.list)){
  varsites <- read_csv(paste("output/new/",var.list[i],sep=""),
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

## Putting the data together
varsites.all$int <- as.factor(varsites.all$int)

pre.post <- varsites.all %>%
  as_tibble() %>%
  filter(sim > 15) %>%
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
  group_by(first_loss) %>%
  mutate(mean = mean(variants/5000))
  # left_join(mut.all, by=c(
  #   "sim" = "sim_num",
  #   "gene" = "gene_name",
  #   "height" = "tree_height",
  #   "int" = "int_ext")) 
fl_maf <- first_loss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5) +
  geom_vline(aes(xintercept=mean,group=height,col=first_loss),lty=2) +
  theme(axis.title = element_text(size=rel(1.5),family="Open Sans"),
        axis.text = element_text(size=rel(1.2),family="Open Sans"),
        strip.text = element_text(size=rel(1.5),family="Open Sans"),
        legend.text = element_text(size=rel(1.2),family="Open Sans"),
        legend.title = element_text(size=rel(1.5),family="Open Sans"),
        legend.position="right") +
  ylab("Density") +
  xlab("Mutation Rate") +
  guides(col=guide_legend(title="Minor\nAllele\nCount"))
fl_maf_facet <- first_loss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5) +
  geom_vline(aes(xintercept=mean,group=height,col=first_loss),lty=2) +
  theme(axis.title = element_text(size=rel(1.5),family="Open Sans"),
        axis.text = element_text(size=rel(1.2),family="Open Sans"),
        strip.text = element_text(size=rel(1.5),family="Open Sans"),
        legend.text = element_text(size=rel(1.2),family="Open Sans"),
        legend.title = element_text(size=rel(1.5),family="Open Sans"),
        legend.position="top") +
  ylab("Density") +
  xlab("Mutation Rate") +
  facet_wrap(~height) +
  guides(col=guide_legend(title="Minor Allele Count"))

first_loss_miss_df <- pre.post %>%
  filter(nSNP == 0 & maf == "0") %>%
  group_by(height,sim,int,gene,variants) %>%
  summarise(first_loss = min(as.numeric(miss),na.rm=T))  %>%
  ungroup() %>%
  group_by(first_loss) %>%
  mutate(mean = mean(variants/5000)) 
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
  geom_vline(aes(xintercept=mean,group=height,col=first_loss),lty=2) +
  theme(axis.title = element_text(size=rel(1.5),family="Open Sans"),
        axis.text = element_text(size=rel(1.2),family="Open Sans"),
        strip.text = element_text(size=rel(1.5),family="Open Sans"),
        legend.text = element_text(size=rel(1.2),family="Open Sans"),
        legend.title = element_text(size=rel(1.5),family="Open Sans"),
        legend.position="right") +
  ylab("Density") +
  xlab("Mutation Rate") +
  guides(col=guide_legend(title="Missing\nData\nThreshold"))
fl_miss_facet <- first_loss_miss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000)  %>%
  filter(!first_loss %in% c(0.25,0.5)) %>%
  ggdensity(x="mut_rate",col="first_loss",alpha=0.2,lwd=1.5) +
  geom_vline(aes(xintercept=mean,group=height,col=first_loss),lty=2) +
  theme(axis.title = element_text(size=rel(1.5),family="Open Sans"),
        axis.text = element_text(size=rel(1.2),family="Open Sans"),
        strip.text = element_text(size=rel(1.5),family="Open Sans"),
        legend.text = element_text(size=rel(1.2),family="Open Sans"),
        legend.title = element_text(size=rel(1.5),family="Open Sans"),
        legend.position="top") +
  ylab("Density") +
  xlab("Mutation Rate") +
  facet_wrap(~height) +
  guides(col=guide_legend(title="Missing Data Threshold"))

## plotting together!
## exported at 1200x1200px
library(patchwork)
(fl_maf | fl_miss) /
  fl_maf_facet /
  fl_miss_facet + plot_annotation(tag_levels = 'A') & 
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
  summarize(n = n())
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
