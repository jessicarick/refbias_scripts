library(tidyverse)

## Pulling info on simulated number of sites
mut.files <- list.files(path="output/new/",pattern="072221-varSites*")

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

var.list <- list.files(path="output/new/",pattern="072221-SNPs*")

varsites.all <- data.frame()

for (i in 1:length(var.list)){
  varsites <- read_csv(paste("output/new/",var.list[i],sep=""),
                         #row.names=NULL,
                         col_names=c("gene","iteration","nSNP"),
                         #stringsAsFactors = FALSE,
                         na=c("","NA")) %>%
    mutate(height = regmatches(var.list[i], regexec('SNPs-([A-Z]+)\\-', var.list[i]))[[1]][2],
           int = regmatches(var.list[i], regexec('SNPs-[A-Z]+-([A-Z]+)', var.list[i]))[[1]][2],
           sim = as.integer(gsub("^s([0-9]+)\\_.*","\\1",iteration)),
           qual = gsub(".*q([0-9]+)\\_.*","\\1",iteration),
           miss = gsub(".*miss(.*)\\_maf.*","\\1",iteration),
           maf = gsub(".*maf([0-9]+\\.[0-9]+)\\..*","\\1",iteration),
           noref = gsub(".*\\.([A-Z]+)\\..*","\\1",iteration)) %>%
    mutate(sim = as.integer(sim))

  varsites.all <- varsites.all %>%
    bind_rows(varsites)
}

summary(varsites.all)

## Putting the data together
varsites.all$int <- as.factor(varsites.all$int)
varsites.all$sim <- as.integer(varsites.all$int)
pre.post <- varsites.all %>%
  left_join(mut.all, by=c(
    "sim" = "sim_num",
    "gene" = "gene_name",
    "height" = "tree_height",
    "int" = "int_ext"
  ))

summary(pre.post)
pre.post$diff <- pre.post$variants - pre.post$nSNP
hist(pre.post$diff)

plot(pre.post$variants,pre.post$diff,col=pre.post$maf)

pre.post %>%
  ggdensity(x="diff",fill="maf")

plot1 <- pre.post %>%
  mutate(lost = case_when(nSNP > 0 ~ "N",
                          nSNP == 0 ~ "Y")) %>%
  #filter(nSNP == 0) %>%
  ggplot() +
  #geom_point(aes(x=variants,y=nSNP,col=maf)) +
  #geom_jitter(aes(col=height),width=0.1,height=0,alpha=0.2) +
  geom_violin(aes(x=maf,y=variants,fill=lost)) +
  #geom_density(aes(x=variants,fill=maf)) +
  facet_wrap(~lost)
