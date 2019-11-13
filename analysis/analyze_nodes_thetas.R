setwd("C:/Users/jrick/Dropbox/i/projects/reference_bias/scripts/refbias_scripts/")

nodes.500000 <- read.csv("output/101819_node_comp_500000.csv",header=TRUE,row.names=1)
nodes.2000000 <- read.csv("output/101819_node_comp_2000000.csv",header=TRUE,row.names=1)
nodes.10000000 <- read.csv("output/101819_node_comp_10000000.csv",header=TRUE,row.names=1)
nodes.all <- nodes.500000 %>%
  full_join(nodes.2000000) %>%
  full_join(nodes.10000000)

nodes.all$maf <- as.factor(nodes.all$maf)
nodes.all$height <- as.factor(nodes.all$height)

results.raxml <- read.csv("output/101819-output-raxml.csv",header=TRUE,row.names=1)
results.raxml$maf <- as.factor(results.raxml$maf)
results.raxml$height <- as.factor(results.raxml$height)

all.data <- nodes.all %>%
  left_join(results.raxml)

## ridgeline plots
plot1 <- ggplot(data = all.data, 
                aes(x=maf,
                    y=mean.support,
                    fill=height))

plot2 <- plot1 +
  #geom_density_ridges(scale = 1.5, rel_min_height = 0.1, alpha = 0.6,
  #                    stat="binline", bins=100)+
  geom_boxplot()+
  scale_fill_manual(values=viridis(3),name="",aesthetics = "fill")+
  theme(axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
        axis.text.y = element_text(size=rel(2)),
        axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
        axis.text.x = element_text(size=rel(2)),
        legend.text = element_text(size=rel(2)),
        legend.title = element_text(size=rel(2)),
        plot.margin = unit(c(6,5.5,20,10),"points"),
        line = element_line(size=1),
        panel.border = element_rect(color = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  scale_y_continuous(name="Number of SNPs")+
  scale_x_discrete(name="MAF")
  #facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
  #geom_hline(yintercept=0,cex=2,lty=2,col="gray")#+
  #theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)

print(plot2)


par(mfrow=c(1,3))
plot(mean_bs ~ sites, data=all.data)
plot(mean_bs ~ maf, data=all.data)
plot(sites ~ maf, data=all.data)

#################################

pre.thetas <- read.csv("output/101819_pre_thetas.csv",header=TRUE,row.names=1)
post.thetas <- read.csv("output/101819_post_thetas.csv",header=TRUE,row.names=1)

hist(pre.thetas$Theta)
hist(post.thetas$Theta)
pre.thetas$maf <- as.factor(pre.thetas$maf)

## ridgeline plots
cols <- viridis(8,alpha=0.7)
plot1 <- ggplot()

plot2 <- plot1 +
  #geom_density_ridges(scale = 1.5, rel_min_height = 0.1, alpha = 0.6)+
  geom_histogram(data = pre.thetas[pre.thetas$Height == 500000,], 
               aes(x=Theta,y=..density..), fill=cols[1])+
  geom_histogram(data = post.thetas[post.thetas$maf == "maf0.01" & pre.thetas$Height == 500000,], 
             aes(x=Theta,y=..density..), fill=cols[2])+
  geom_histogram(data = post.thetas[post.thetas$maf == "maf0.02" & pre.thetas$Height == 500000,], 
                 aes(x=Theta,y=..density..), fill=cols[3])+
  geom_histogram(data = post.thetas[post.thetas$maf == "maf0.03" & pre.thetas$Height == 500000,], 
                 aes(x=Theta,y=..density..), fill=cols[4])+
  geom_histogram(data = post.thetas[post.thetas$maf == "maf0.04" & pre.thetas$Height == 500000,], 
               aes(x=Theta,y=..density..), fill=cols[5])+
  geom_histogram(data = post.thetas[post.thetas$maf == "maf0.05" & pre.thetas$Height == 500000,], 
               aes(x=Theta,y=..density..), fill=cols[6])+
  geom_histogram(data = post.thetas[post.thetas$maf == "maf0.1" & pre.thetas$Height == 500000,], 
               aes(x=Theta,y=..density..), fill=cols[7])+
  #scale_fill_manual(values=viridis(3),name="",aesthetics = "fill")+
  theme(axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
        axis.text.y = element_text(size=rel(2)),
        axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
        axis.text.x = element_text(size=rel(2)),
        legend.text = element_text(size=rel(2)),
        legend.title = element_text(size=rel(2)),
        plot.margin = unit(c(6,5.5,20,10),"points"),
        line = element_line(size=1),
        panel.border = element_rect(color = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 16),
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(name="Theta")+
  scale_y_discrete(name="Density")+
#facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
#geom_hline(yintercept=0,cex=2,lty=2,col="gray")#+
  theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)

print(plot2)
