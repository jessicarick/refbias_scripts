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
results <- results %>%
  mutate(height2 = case_when(height == "LONG" ~ "Low ILS",
                             height == "MED" ~ "Med ILS",
                             height == "SHORT" ~ "High ILS"))
for (i in 1:length(responses)){
  plot1 <- ggplot(data = results, 
                  aes(x=as.factor(maf),
                      y=results[,responses[i]]))
  
  plot2 <- plot1 +
    #geom_boxplot(alpha=0.7, notch=FALSE, varwidth=FALSE, weight=2)+
    geom_point(aes(color=int),alpha=0.2) +
    stat_summary(aes(color=int),fun.data="mean_cl_boot",position=position_nudge(x=c(0.1,-0.1)),size=1.2,alpha=0.9) +
    #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
    scale_color_manual(values=c("#009980", "#006699","turquoise","magenta","orange","green"),name="Reference")+
    theme_custom()+
    theme(#axis.title.y = element_text(angle=90, size=rel(2), face="plain"),
          #axis.text.y = element_text(size=rel(2)),
          #axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
          #axis.text.x = element_text(size=rel(2)),
          legend.text = element_text(size=rel(2)),
          legend.title = element_text(size=rel(2)),
          plot.margin = unit(c(6,5.5,20,10),"points"),
          line = element_line(size=1),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          strip.text.x = element_text(size = 16),
          strip.background = element_rect(color="black",fill=NA),
          legend.position="none")+
    ylab("Standardized ingroup Sackin")+
    #ylab(responses[i]) +
    scale_x_discrete(name="Minor Allele Count") +
    facet_wrap(vars(factor(height2,levels=c("Low ILS","Med ILS","High ILS"))),nrow=1,strip.position = "top") +
    geom_hline(yintercept=0,cex=1.5,lty=2,col="gray")#+
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

## plotting variables of interest together
resp1 <- responses[c(1,3,4)]
resp2 <- responses[c(5:7)]
resp3 <- responses[8:10]
resp4 <- responses[11:13]
results_long <- results %>% 
  pivot_longer(cols=c(sites,tree.height:gt_rf),
               names_to="stat",values_to="value") %>% 
  filter(stat %in% resp3)

results_long %>% 
  ggplot(aes(x=as.factor(maf),y=value)) + 
  geom_point(aes(color=int),alpha=0.2) + 
  stat_summary(aes(color=int), fun.data="mean_cl_boot",
               position=position_nudge(x=c(-0.1,0.1)),size=1.2,alpha=0.9) + 
  theme_custom() + 
  scale_color_manual(values=c("#009980", "#006699","turquoise","magenta","orange","green"),name="Reference")+
  facet_grid(stat ~ factor(height2,levels=c("Low ILS","Med ILS","High ILS")), scales="free_y") +
  theme(axis.title.y = element_blank(),
    #axis.text.y = element_text(size=rel(2)),
    #axis.title.x = element_text(size=rel(2),vjust=-2, face="plain"),
    #axis.text.x = element_text(size=rel(2)),
    legend.text = element_text(size=rel(1.1)),
    legend.title = element_text(size=rel(1.1)),
    plot.margin = unit(c(6,5.5,20,10),"points"),
    line = element_line(size=1),
    panel.border = element_rect(color = "black", fill=NA, size=1),
    strip.text = element_text(size = 16),
    strip.background = element_rect(color="black",fill=NA),
    legend.position="top")+
  scale_x_discrete(name="Minor Allele Count") +
  geom_hline(yintercept=0,cex=1.5,lty=2,col="gray")

#########################
## correlation plots ####
#########################
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

### for emp

resp <- c("ingroup.tree.height","Avg.BLs",
          "ingroup.gamma","ingroup.colless","ingroup.sackin","snps")
pred <- c("missing","maf","snps","int")
resp
cor.cich <- formatted_cors(results.emp.cichlids[,resp]) %>%
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
cor.lates <- formatted_cors(results.emp.lates[,resp]) %>%
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
ggarrange(cor.cich,cor.lates,labels=c("A. Tropheines","B. Lates"),
          font.label=list(size=20,font.family="Open Sans"))

formatted_cors(results.emp[,pred]) %>%
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
