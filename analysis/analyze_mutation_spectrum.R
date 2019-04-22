library(tidyverse)
library(ggpubr)

## Pulling info on simulated number of sites
mut.files <- list.files(path="output/",pattern="041119-varSites*")

mut.all <- data.frame(gene=integer(),
                      variants=integer(),
                      sim_num=integer(),
                      tree_height=integer(),
                      int_ext=character(),
                      stringsAsFactors = FALSE)

for (i in 1:length(mut.files)){
  varsites <- scan(paste("output/",mut.files[i],sep=""))
  sim <- as.integer(regmatches(mut.files[i], regexec('sim([0-9]+)', mut.files[i]))[[1]][2])
  height <- as.integer(regmatches(mut.files[i], regexec('-([0-9]+)\\-sim', mut.files[i]))[[1]][2])
  int <- as.character(regmatches(mut.files[i], regexec('-([A-Z]+)', mut.files[i]))[[1]][2])
  mut.df <- data.frame(gene=seq(1,1000),
                       variants=mut[-1],
                       sim_num=sim,
                       tree_height=height,
                       int_ext=int)
  mut.all <- rbind(mut.all,mut.df)
}
mut.all$tree_height <- as.factor(mut.all$tree_height)
ggpubr::ggdensity(data=mut.all,x="variants",
                    color="tree_height",
                    fill="tree_height",
                    xlab="Variant Sites",
                    ylab="Frequency",
                    add="mean")

## Pulling information from individual filtering combinations
mut.filtered <- read.csv("output/041119-mutations.txt",header=TRUE,stringsAsFactors = FALSE)
mut <- mut.filtered[mut.filtered$gene != "gene",]
mut$num_noRef[is.na(mut$num_noRef)] <- 0
mut$num_nonInv[is.na(mut$num_nonInv)] <- 0

for (i in 1:nrow(mut)){
  mut$tree_height[i] <- as.integer(regmatches(mut$gene[i], regexec('height([0-9]+)\\_s', mut$gene[i]))[[1]][2])
  mut$sim_num[i] <- as.integer(regmatches(mut$gene[i], regexec('sim([0-9]+)\\_g', mut$gene[i]))[[1]][2])
  mut$locus[i] <- as.integer(regmatches(mut$gene[i], regexec('gene([0-9]+)', mut$gene[i]))[[1]][2])
  mut$qual[i] <- as.integer(regmatches(mut$gene[i], regexec('q([0-9]+)\\_m', mut$gene[i]))[[1]][2])
  mut$miss[i] <- as.numeric(regmatches(mut$gene[i], regexec('miss(.*)\\_m', mut$gene[i]))[[1]][2])
  mut$maf[i] <- as.numeric(regmatches(mut$gene[i], regexec('maf(0.*)', mut$gene[i]))[[1]][2])
}
summary(mut)
mut$tree_height <- as.factor(mut$tree_height)
#mut$sim_num <- as.factor(mut$sim_num)

## combine the two tables together

mut.combined <- left_join(mut,mut.all[mut.all$int_ext == "EXT",],
                          by=c("locus" = "gene","sim_num","tree_height"))
summary(mut.combined)
mut.combined$diff <- mut.combined$variants - mut.combined$num_SNPs
mut.combined$retained <- case_when(mut.combined$num_SNPs == 0 ~ "no",
                                   TRUE ~ "yes")

hist(mut.combined$diff)
