#######################
## univariate plots
######################
output <- "052419-output"

results.raxml <- read.csv(paste("output/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

responses <- colnames(results.raxml[,10:27])

pdf(paste("output/",output,"-univariate-plots.pdf",sep=""),width=11,height=8)

for (i in 1:length(responses)){
  plot1 <- ggplot(data = results.raxml, aes(factor(maf),results.raxml[,responses[i]],fill=factor(int)))
  
  plot2 <- plot1 +
    geom_boxplot(alpha=0.7)+
    scale_fill_manual(values=c("#009980", "#006699"),name="")+
    theme_classic()+
    theme(axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
          axis.text.y = element_text(size=rel(2)),
          axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
          axis.text.x = element_text(size=rel(2)),
          legend.text = element_text(size=rel(2)),
          legend.title = element_text(size=rel(2)),
          plot.margin = unit(c(6,5.5,20,10),"points"),
          line = element_line(size=1),
          panel.border = element_rect(color = "black", fill=NA, size=1))+
    scale_y_continuous(name=responses[i])+
    scale_x_discrete(name="Minor Allele Frequency Cutoff")+
    facet_wrap(vars(noref))
  
  # ylim1 = boxplot.stats(results.raxml[,responses[i]])$stats[c(1, 5)]
  # 
  # plot3 <- plot2 + coord_cartesian(ylim = ylim1 + c(-0.05, 0.05) * diff(ylim1) / 2)
  
  print(plot2)
  
  # ggviolin(results, x = "maf", y = "tree.height", fill = "int",
  #          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  #          add = "boxplot", add.params = list(fill = "white"))
}

dev.off()
