library(tidyverse)
library(ggpubr)

mut <- read.csv("output/040319-mutations.txt", header=T, stringsAsFactors=F)

for (i in 1:nrow(mut)){
  mut$height[i] <- as.integer(regmatches(mut$gene[i], regexec('height([0-9]+)\\_s', mut$gene[i]))[[1]][2])
  mut$sim[i] <- as.integer(regmatches(mut$gene[i], regexec('sim([0-9]+)\\_g', mut$gene[i]))[[1]][2])
  mut$locus[i] <- as.integer(regmatches(mut$gene[i], regexec('gene([0-9]+)', mut$gene[i]))[[1]][2])
  mut$qual[i] <- as.integer(regmatches(mut$gene[i], regexec('q([0-9]+)\\_m', mut$gene[i]))[[1]][2])
  mut$miss[i] <- as.numeric(regmatches(mut$gene[i], regexec('miss(.*)\\_m', mut$gene[i]))[[1]][2])
  mut$maf[i] <- as.numeric(regmatches(mut$gene[i], regexec('maf(0.*)', mut$gene[i]))[[1]][2])
}
#mut$height <- as.factor(mut$height)
#mut$sim <- as.factor(mut$sim)

mut$orig <- case_when(is.na(mut$maf) ~ "orig",
                      !is.na(mut$maf) ~ "filtered")
mut$orig <- as.factor(mut$orig)

mut.orig <- mut[is.na(mut$maf),]
mut.filter <- mut[!is.na(mut$maf),]

hist(mut.orig$num_SNPs)

ggdensity(mut.orig[mut.orig$locus == 1,], x = "num_SNPs",
          add = "mean", rug = TRUE,
          color = "height", fill = "height",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))
ggdensity(mut[mut$locus == "1",], x = "num_SNPs",
          add = "mean", rug = TRUE,
          color = "orig", fill = "orig",
          palette = c("#00AFBB", "#E7B800"))

nloci <- 1000
diff <- data.frame(height=integer(length=3*nloci),
                   locus=integer(length=3*nloci),
                   mean_orig=numeric(length=3*nloci),
                   mean_filtered=numeric(length=3*nloci),
                   diff=numeric(length=3*nloci))
diff$locus <- rep(seq(1:nloci),3)
diff$height <- c(rep(500000,nloci),rep(2000000,nloci),rep(10000000,nloci))

for(i in 1:nrow(diff)){
  diff$mean_orig[i] <- mean(mut.orig$num_SNPs[mut.orig$height == diff$height[i] & 
                                                mut.orig$locus == diff$locus[i]])
  diff$mean_filtered[i] <- mean(mut.filter$num_SNPs[mut.filter$height == diff$height[i] & 
                                                    mut.filter$locus == diff$locus[i]])
}

diff$diff <- diff$mean_orig - diff$mean_filtered

diff$height <- as.factor(diff$height)
mean <- data.frame(mean=mean(diff$mean_filtered, na.rm=T))
ggdensity(diff, x = "mean_orig",
          add = "mean", rug = TRUE,
          color = "#00AFBB", fill = "#00AFBB",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  geom_density(aes(x=mean_filtered), color="#E7B800", fill=scales::alpha("#E7B800",0.5)) +
  geom_vline(data=mean,aes(xintercept=mean),
             linetype="dashed", color="#E7B800")

######

mut.mean <- mut[!is.na(mut$maf),] %>% 
  group_by(.dots=c("height","miss","qual","maf")) %>% 
  summarize(mean=mean(num_SNPs))
orig.mean <- mut[is.na(mut$maf),] %>% 
  group_by(height) %>% 
  summarize(mean=mean(num_SNPs))

mut.mean$orig_mean <- case_when(mut.mean$height == 500000 ~ orig.mean$mean[1],
                                  mut.mean$height == 2000000 ~ orig.mean$mean[2],
                                  mut.mean$height == 10000000 ~ orig.mean$mean[3])
mut.mean$diff <- mut.mean$orig_mean - mut.mean$mean

mod1 <- lm(diff ~ height + miss + maf + qual, data=mut.mean)
summary(mod1)

mut.mean$maf <- as.factor(mut.mean$maf)
mut.mean$miss <- as.factor(mut.mean$miss)
ggdensity(diff, x = c("mean_orig","mean_filtered"),
          add = "mean", rug = TRUE,
          color = "height", fill = "height",
          palette = "npg",combine=TRUE)
ggdensity(mut.mean, x = "diff",
          add = "mean", rug = TRUE,
          color = "miss", fill = "miss",
          palette = "npg")