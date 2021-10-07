#######################
## univariate plots
######################
library(tidyverse)
library(ggridges)

output <- "082921-cichlids-emp-output"

#results.raxml <- read.csv(paste("output/new/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

responses <- colnames(results.raxml[,c(8,12:35)])
responses <- colnames(results.raxml[,c(9:20)])
#responses <- colnames(results.raxml.cichlids[,c(9:22)])

pdf(paste("output/new/",output,"-univariate-plots_refdist.pdf",sep=""),width=11,height=8)

results <- as.data.frame(results.mod)
for (i in 1:length(responses)){
  plot1 <- ggplot(data = results, 
                  aes(x=as.factor(maf),
                      y=results[,responses[i]],
                      fill=int))
  
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
    scale_x_discrete(name="MAF") +
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

# plot1 <- ggplot(data = results.raxml, 
#                 aes(y=as.factor(results.raxml$maf),
#                     x=RF.Dist.ML,
#                     fill=factor(simulation)))
# 
# plot2 <- plot1 +
#   #geom_boxplot(alpha=0.7, notch=FALSE, varwidth=FALSE, weight=2)+
#   geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
#   #scale_fill_manual(values=c("#009980", "#006699","turquoise","magenta","orange","green"),name="",aesthetics = "fill")+
#   scale_fill_viridis_d(begin=0.2,end=0.8,alpha=0.5,name="",aesthetics = "fill") +
#   theme_classic()+
#   theme(axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
#         axis.text.y = element_text(size=rel(2)),
#         axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
#         axis.text.x = element_text(size=rel(2)),
#         legend.text = element_text(size=rel(2)),
#         legend.title = element_text(size=rel(2)),
#         #plot.margin = unit(c(6,5.5,20,10),"points"),
#         line = element_line(size=1),
#         panel.border = element_rect(color = "black", fill=NA, size=1),
#         strip.text.x = element_text(size = 16))+
#   #scale_x_continuous(name="Ingroup Gamma")+
#   scale_y_discrete(name="MAF")+
#   #xlim(-50,10)+
#   facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
#   geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
#   theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)
# 
# # ylim1 = boxplot.stats(results.raxml[,responses[i]])$stats[c(1, 5)]
# # 
# # plot3 <- plot2 + coord_cartesian(ylim = ylim1 + c(-0.05, 0.05) * diff(ylim1) / 2)
# 
# print(plot2)


## correlation plots
cors <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df))
  # turn all three matrices (r, n, and P into a data frame)
  Mdf <- purrr::map(M, ~data.frame(.x))
  # return the three data frames in a list
  return(Mdf)
}
formatted_cors <- function(df){
  cors(df) %>%
    purrr::map(~rownames_to_column(.x, var="measure1")) %>%
    purrr::map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}

resp <- responses[c(3,4,7,15,16,17,22,23,24)]
pred <- c("height","missing","maf","sites","int","avg_dxy")
resp
formatted_cors(results.raxml[,resp]) %>%
  ggplot(aes(measure1, measure2, col=r)) + ## to get the rect filled
  geom_tile(col="black", fill="white") +
  geom_point(aes(size = abs(r)), shape=15) +
  labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation") +
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1))  +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_custom() +
  theme(text=element_text(family="Open Sans"),
        axis.text = element_text(size=rel(1.8)),
        axis.text.x = element_text(angle=60,hjust=1),
        legend.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.2)),
        axis.ticks = element_blank()) +
  scale_size(range=c(1,20), guide=NULL) 

formatted_cors(results.raxml[,pred]) %>%
  ggplot(aes(measure1, measure2, col=r)) + ## to get the rect filled
  geom_tile(col="black", fill="white") +
  geom_point(aes(size = abs(r)), shape=15) +
  labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation") +
  scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1))  +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme_custom() +
  theme(text=element_text(family="Open Sans"),
        axis.text = element_text(size=rel(1.8)),
        axis.text.x = element_text(angle=60,hjust=1),
        legend.title = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.2)),
        axis.ticks = element_blank()) +
  scale_size(range=c(1,20), guide=NULL) 
