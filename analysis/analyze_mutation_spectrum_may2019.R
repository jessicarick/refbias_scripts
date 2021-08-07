library(tidyverse)

## Pulling info on simulated number of sites
mut.files <- list.files(path="output/new/",pattern="072821-varSites*")

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
  int <- as.character(regmatches(mut.files[i], regexec('-sim[0-9]+-([A-Z]+)-', mut.files[i]))[[1]][2])
  mut.df <- data.frame(gene=seq(1,1000),
                       variants=varsites[-1],
                       sim_num=sim,
                       tree_height=height,
                       int_ext=int)
  mut.all <- rbind(mut.all,mut.df)
}

mut.all$gene_name <- paste0("gene",sprintf("%04d",mut.all$gene))
mut.all$sim_num <- as.integer(as.character(mut.all$sim_num))

#mut.all$tree_height <- as.factor(mut.all$tree_height)
ggpubr::ggdensity(data=mut.all,x="variants",
                  color="tree_height",
                  fill="tree_height",
                  xlab="Simulated Variant Sites",
                  ylab="Frequency",
                  add="mean")

## Pulling information for post-filtering SNPs

var.list <- list.files(path="output/new/",pattern="072821-SNPs*")

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
           miss = gsub(".*miss(.*)\\_maf.*","\\1",iteration),
           maf = gsub(".*maf(.*)\\.REF.*","\\1",iteration),
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
  left_join(mut.all, by=c(
    "sim" = "sim_num",
    "gene" = "gene_name",
    "height" = "tree_height",
    "int" = "int_ext"
  ))

summary(pre.post)
pre.post$diff <- pre.post$variants - pre.post$nSNP
ggdensity(pre.post,x="diff",fill=NA,alpha=0.3,color="maf",
          add="mean")

#plot(pre.post$variants,pre.post$diff,col=pre.post$maf)

pre.post %>%
  ggdensity(x="nSNP",fill="maf")

first_loss_df <- pre.post %>%
  filter(nSNP == 0 & miss == "0") %>%
  group_by(height,sim,int,gene,variants) %>%
  summarise(first_loss = min(as.numeric(maf)))  %>%
  ungroup()
  # left_join(mut.all, by=c(
  #   "sim" = "sim_num",
  #   "gene" = "gene_name",
  #   "height" = "tree_height",
  #   "int" = "int_ext")) 
first_loss_df %>%
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggdensity(x="mut_rate",col="first_loss") +
  facet_wrap(~height)
  
first_loss_df %>% 
  mutate(first_loss = as.factor(first_loss),
         mut_rate = variants/5000) %>%
  ggplot() +
  geom_jitter(aes(x=first_loss,y=mut_rate,col=int),alpha=0.5) +
  #geom_violin(aes(x=first_loss,y=mut_rate,col=int),fill=NA,draw_quantiles=c(0.5)) +
  facet_wrap(~height) +
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
  #t.test(.$variants ~ .$lost)
  filter(miss == "0" & lost == "N") %>%
  ggdensity(x="nSNP",fill=NA,alpha=0.3,color="maf",
            add="mean")
  
  ggplot() +
  #geom_point(aes(x=maf,y=diff,col=maf)) +
  #geom_jitter(aes(col=height),width=0.1,height=0,alpha=0.2) +
  #geom_density(aes(x=variants,fill=maf),alpha=0.2) +
  geom_density(aes(x=diff,fill=maf), alpha=0.2, trim=TRUE) +
  facet_wrap(~lost) +
  theme_custom()
