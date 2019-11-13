## ridgeline plots
plot1 <- ggplot(data = results.raxml, 
                aes(y=maf,
                    x=sites,
                    fill=height))

plot2 <- plot1 +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.1, alpha = 0.7)+
  scale_fill_manual(values=viridis(3),name="",aesthetics = "fill")+
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
  scale_x_continuous(name="Number of Sites")+
  scale_y_discrete(name="MAF")+
  #facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
  geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
  theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)

print(plot2)

