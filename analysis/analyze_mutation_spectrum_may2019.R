library(tidyverse)

## Pulling info on simulated number of sites
mut.files <- list.files(path="output/",pattern="053019-varSites*")

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
                       variants=varsites[-1],
                       sim_num=sim,
                       tree_height=height,
                       int_ext=int)
  mut.all <- rbind(mut.all,mut.df)
}
#mut.all$tree_height <- as.factor(mut.all$tree_height)
ggpubr::ggdensity(data=mut.all,x="variants",
                  color="tree_height",
                  fill="tree_height",
                  xlab="Variant Sites",
                  ylab="Frequency",
                  add="mean")

## Pulling information for post-filtering SNPs

var.list <- list.files(path="output/",pattern="053019-variants*")

varsites.all <- data.frame(header=character(),
                           gene=character(),
                           start=integer(),
                           end=integer(),
                           SNPs=integer())

for (i in 1:length(var.list)){
  varsites <- read.table(paste("output/",var.list[i],sep=""),
                         fill=TRUE,
                         row.names=NULL,
                         col.names=c("v1","v2"),
                         stringsAsFactors = FALSE,
                         na.strings=c("","NA"))
  headers <- grep('height', varsites$v1)
  varsites$header <- character(nrow(varsites))
  
  for (k in 1:(length(headers)-1)){
    varsites$header[headers[k]:(headers[k+1]-1)] <- varsites$v1[headers[k]]
  }
  k <- length(headers)
  varsites$header[headers[k]:nrow(varsites)] <- varsites$v1[headers[k]]
  
  varsites$sim <- as.integer(regmatches(var.list[i], regexec('sim([0-9]+)', var.list[i]))[[1]][2])
  varsites$height <- as.integer(regmatches(var.list[i], regexec('-([0-9]+)\\-sim', var.list[i]))[[1]][2])
  varsites$int <- as.character(regmatches(var.list[i], regexec('-([A-Z]+)', var.list[i]))[[1]][2])
  varsites$gene <- substr(varsites[,2],1,7)
  varsites$locus <- as.integer(substr(varsites$gene,5,7))
  varsites$firstlast <- substr(varsites[,2],9,20)

  for (j in 1:nrow(varsites)){
    varsites$qual[j] <- as.integer(regmatches(varsites$header[j],regexec('q([0-9]+)',varsites$header[j]))[[1]][2])
    varsites$miss[j] <- as.numeric(regmatches(varsites$header[j],regexec('miss(.*[0-9]+)\\_',varsites$header[j]))[[1]][2])
    varsites$maf[j] <- as.numeric(regmatches(varsites$header[j],regexec('maf(.*[0-9]+)\\.[A-Z]+',varsites$header[j]))[[1]][2])
    varsites$noref[j] <- as.character(regmatches(varsites$header[j],regexec('\\.([A-Z]+)',varsites$header[j]))[[1]][2])
  }

  varsites.final <- varsites %>% separate(firstlast, "-",
                  into = c("start", "end"), 
                  remove = TRUE)
  varsites.final$SNPs <- as.integer(varsites.final$end) - as.integer(varsites.final$start)
  
  varsites.final <- varsites.final[!is.na(varsites.final$gene),-c(1,2)]
  
  varsites.all <- rbind(varsites.all,varsites.final)
}

summary(varsites.all)

## Putting the data together
varsites.all$int <- as.factor(varsites.all$int)
pre.post <- varsites.all %>%
  left_join(mut.all,by=c("sim" = "sim_num","height" = "tree_height", "int" = "int_ext", "locus" = "gene"))
summary(pre.post)
pre.post$diff <- pre.post$variants - pre.post$SNPs
hist(pre.post$diff)

plot(pre.post$variants,pre.post$diff,col=pre.post$maf)
