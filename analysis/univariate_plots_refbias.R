#######################
## univariate plots
######################
library(tidyverse)
library(ggridges)

output <- "072221-output"

results.raxml <- read.csv(paste("output/new/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

responses <- colnames(results.raxml[,9:31])

pdf(paste("output/",output,"-univariate-plots_refdist.pdf",sep=""),width=11,height=8)

for (i in 1:length(responses)){
  plot1 <- ggplot(data = results.raxml[results.raxml$noref == "REF",], 
                  aes(x=as.factor(results.raxml[results.raxml$noref == "REF",]$maf),
                      y=results.raxml[results.raxml$noref == "REF",responses[i]],
                      fill=results.raxml[results.raxml$noref == "REF",]$int))
  
  plot2 <- plot1 +
    geom_boxplot(alpha=0.7, notch=FALSE, varwidth=FALSE, weight=2)+
    #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
    scale_fill_manual(values=c("#009980", "#006699","turquoise","magenta","orange","green"),name="",aesthetics = "fill")+
    theme_classic()+
    theme(axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
          axis.text.y = element_text(size=rel(2)),
          axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
          axis.text.x = element_text(size=rel(2)),
          legend.text = element_text(size=rel(2)),
          legend.title = element_text(size=rel(2)),
          plot.margin = unit(c(6,5.5,20,10),"points"),
          line = element_line(size=1),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          strip.text.x = element_text(size = 16))+
    ylab(responses[i])+
    scale_x_discrete(name="MAF")+
    facet_wrap(vars(height),nrow=1,strip.position = "top")
    #geom_hline(yintercept=0,cex=2,lty=2,col="gray")#+
    #theme_ridges()
  
  # ylim1 = boxplot.stats(results.raxml[,responses[i]])$stats[c(1, 5)]
  # 
  # plot3 <- plot2 + coord_cartesian(ylim = ylim1 + c(-0.05, 0.05) * diff(ylim1) / 2)
  
  print(plot2)
  
  # ggviolin(results, x = "maf", y = "tree.height", fill = "int",
  #          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  #          add = "boxplot", add.params = list(fill = "white"))
}

dev.off()

############################

plot1 <- ggplot(data = results.raxml, 
                aes(y=as.factor(results.raxml$maf),
                    x=RF.Dist.ML,
                    fill=factor(int)))

plot2 <- plot1 +
  #geom_boxplot(alpha=0.7, notch=FALSE, varwidth=FALSE, weight=2)+
  geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  #scale_fill_manual(values=c("#009980", "#006699","turquoise","magenta","orange","green"),name="",aesthetics = "fill")+
  scale_fill_viridis_d(begin=0.2,end=0.8,alpha=0.5,name="",aesthetics = "fill") +
  theme_classic()+
  theme(axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
        axis.text.y = element_text(size=rel(2)),
        axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
        axis.text.x = element_text(size=rel(2)),
        legend.text = element_text(size=rel(2)),
        legend.title = element_text(size=rel(2)),
        #plot.margin = unit(c(6,5.5,20,10),"points"),
        line = element_line(size=1),
        panel.border = element_rect(color = "black", fill=NA, size=1),
        strip.text.x = element_text(size = 16))+
  #scale_x_continuous(name="Ingroup Gamma")+
  scale_y_discrete(name="MAF")+
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
  geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
  theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)

# ylim1 = boxplot.stats(results.raxml[,responses[i]])$stats[c(1, 5)]
# 
# plot3 <- plot2 + coord_cartesian(ylim = ylim1 + c(-0.05, 0.05) * diff(ylim1) / 2)

print(plot2)
