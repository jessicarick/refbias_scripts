library(lme4)
library(tidyverse)
library(MuMIn)
library(car)
library(LMERConvenienceFunctions)
library(ggsci)
library(ggpubr)
library(ggridges)
library(mgcv)

cols <- c("#F2AD00","gray80","#00A08A")
here::i_am("analysis/model_refbias_trees_RF.R")

output <- "092321-subsamp-output"
output <- "092321-output"
results.raxml <- read.csv(here("output","new",paste0(output,"-raxml.csv")),header=TRUE,row.names=1,sep=",")

## preparing data object
results.mod <- results.raxml[results.raxml$simulation > 15 & results.raxml$simulation < 26 & results.raxml$RF.Dist.ML < 0.9,]
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.numeric(as.character(results.mod$missing))
results.mod$maf <- as.numeric(as.character(results.mod$maf))
results.mod$int <- as.factor(results.mod$int)

## rf distance
# all

m.rf.all <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                   avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod)
sum.rf.all <- summary(m.rf.all)
r.squaredGLMM(m.rf.all)

confint.rf.all<-data.frame(confint(m.rf.all)[-c(1:3),])
colnames(confint.rf.all) <- c("minCI","maxCI")
confint.rf.all$var <- rownames(confint.rf.all)
confint.rf.all$est <- sum.rf.all$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf.all$sig <- sum.rf.all$coefmat.full[-1,5]
confint.rf.all$sig <- case_when(confint.rf.all$minCI > 0 ~ "pos",
                                  confint.rf.all$maxCI < 0 ~ "neg",
                                  TRUE ~ "ns")
confint.rf.all$sig <- factor(confint.rf.all$sig, levels=c("pos","ns","neg"))
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.all <- ggplot(confint.rf.all, aes(x = var, y = est, color=sig))
vars.rf.all.bars <- vars.rf.all + geom_blank() +
  #color = "cyl",                                # Color by groups
  #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
  #sorting = "descending",                       # Sort value in descending order
  #add = "segments",                             # Add segments from y = 0 to dots
  #add.params = list(color = "lightgray", size = 2), # Change segment color and size
  #group = "cyl",                                # Order by groups
  #dot.size = 4,                                 # Large dot size
  #label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
  #font.label = list(color = "white", size = 9, 
  #                   vjust = 0.5),               # Adjust label parameters
  #ggtheme = theme_pubr(),                        # ggplot2 theme
xlab("")+
  ylab("Coefficient\nall Trees")+
  #xlim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) +
  theme(axis.title.x=element_blank())
vars.rf.all.bars

# short

m.rf.short <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                     avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod[results.mod$height == "SHORT",])
sum.rf.short <- summary(m.rf.short)
r.squaredGLMM(m.rf.short)

confint.rf.short<-data.frame(confint(m.rf.short)[-c(1:3),])
colnames(confint.rf.short) <- c("minCI","maxCI")
confint.rf.short$var <- rownames(confint.rf.short)
confint.rf.short$est <- sum.rf.short$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf.short$sig <- sum.rf.short$coefmat.full[-1,5]
confint.rf.short$sig <- case_when(confint.rf.short$minCI > 0 ~ "pos",
                            confint.rf.short$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.rf.short$sig <- factor(confint.rf.short$sig, levels=c("pos","ns","neg"))
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.short <- ggplot(confint.rf.short, aes(x = var, y = est, color=sig))
vars.rf.short.bars <- vars.rf.short + geom_blank() +
  #color = "cyl",                                # Color by groups
  #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
  #sorting = "descending",                       # Sort value in descending order
  #add = "segments",                             # Add segments from y = 0 to dots
  #add.params = list(color = "lightgray", size = 2), # Change segment color and size
  #group = "cyl",                                # Order by groups
  #dot.size = 4,                                 # Large dot size
  #label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
  #font.label = list(color = "white", size = 9, 
  #                   vjust = 0.5),               # Adjust label parameters
  #ggtheme = theme_pubr(),                        # ggplot2 theme
xlab("")+
  ylab("Coefficient\nSHORT Trees")+
  #xlim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) +
  theme(axis.title.x=element_blank())
vars.rf.short.bars

# med

m.rf.med <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                   avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                 data = results.mod[results.mod$height == "MED",])
sum.rf.med <- summary(m.rf.med)
r.squaredGLMM(m.rf.med)

confint.rf.med<-data.frame(confint(m.rf.med)[-c(1:3),])
colnames(confint.rf.med) <- c("minCI","maxCI")
confint.rf.med$var <- rownames(confint.rf.med)
confint.rf.med$est <- sum.rf.med$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf.med$sig <- sum.rf.med$coefmat.full[-1,5]
confint.rf.med$sig <- case_when(confint.rf.med$minCI > 0 ~ "pos",
                            confint.rf.med$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.rf.med$sig <- factor(confint.rf.med$sig, levels=c("pos","ns","neg"))

#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.med <- ggplot(confint.rf.med, aes(x = var, y = est, color=sig))
vars.rf.med.bars <- vars.rf.med + geom_blank() +
  #color = "cyl",                                # Color by groups
  #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
  #sorting = "descending",                       # Sort value in descending order
  #add = "segments",                             # Add segments from y = 0 to dots
  #add.params = list(color = "lightgray", size = 2), # Change segment color and size
  #group = "cyl",                                # Order by groups
  #dot.size = 4,                                 # Large dot size
  #label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
  #font.label = list(color = "white", size = 9, 
  #                   vjust = 0.5),               # Adjust label parameters
  #ggtheme = theme_pubr(),                        # ggplot2 theme
xlab("")+
  ylab("Coefficient\nMED Trees")+
  #xlim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) +
  theme(axis.title.x=element_blank())
vars.rf.med.bars

# long
m.rf.long <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                    avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                  data = results.mod[results.mod$height == "LONG",])
sum.rf.long <- summary(m.rf.long)
r.squaredGLMM(m.rf.long)

confint.rf.long<-data.frame(confint(m.rf.long)[-c(1:3),])
colnames(confint.rf.long) <- c("minCI","maxCI")
confint.rf.long$var <- rownames(confint.rf.long)
confint.rf.long$est <- sum.rf.long$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf.long$sig <- sum.rf.long$coefmat.full[-1,5]
confint.rf.long$sig <- case_when(confint.rf.long$minCI > 0 ~ "pos",
                            confint.rf.long$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.rf.long$sig <- factor(confint.rf.long$sig, levels=c("pos","ns","neg"))

#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.long <- ggplot(confint.rf.long, aes(x = var, y = est, color=sig))
vars.rf.long.bars <- vars.rf.long + geom_blank() +
  #color = "cyl",                                # Color by groups
  #palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
  #sorting = "descending",                       # Sort value in descending order
  #add = "segments",                             # Add segments from y = 0 to dots
  #add.params = list(color = "lightgray", size = 2), # Change segment color and size
  #group = "cyl",                                # Order by groups
  #dot.size = 4,                                 # Large dot size
  #label = round(dfm$mpg_z,1),                        # Add mpg values as dot labels
  #font.label = list(color = "white", size = 9, 
  #                   vjust = 0.5),               # Adjust label parameters
  #ggtheme = theme_pubr(),                        # ggplot2 theme
  xlab("")+
  ylab("Coefficient\nLONG Trees")+
  #ylim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) +
  theme(axis.title.x=element_blank())
vars.rf.long.bars

# put them all together
vars.plots.rf <- ggarrange(vars.rf.long.bars,
          vars.rf.med.bars,
          vars.rf.short.bars,
          #labels=c("A","B","C"),
          ncol=3)
print(vars.plots.rf)

## ridgeline plots
plot1 <- ggplot(data = results.mod, 
                aes(y=as.factor(maf),
                    x=RF.Dist.ML,
                    fill=int))

plot2 <- ggplot(data = results.mod, 
                aes(y=as.factor(maf),
                    x=RF.Dist.ML,
                    fill=int)) +
  #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  stat_halfeye(data= . %>% filter(int == "EXT"), aes(col=int), alpha=0.7, 
               adjust=1, side="top", slab_color = "black", slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=0.03),
               justification=-0.1, point_interval="median_qi", scale=0.6) +
  stat_halfeye(data= . %>% filter(int == "INT"), aes(col=int), alpha=0.7, 
               adjust=1, side="bottom", slab_color = "black",slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=-0.03),
               justification=1.1, point_interval="median_qi", scale=0.6) +
  scale_fill_manual(values=cols.int[c(1,3)],name="", aesthetics = "fill")+
  scale_color_manual(values=cols.int[c(1,3)],name="", aesthetics = "color")+
  #scale_fill_viridis_d(begin=0.2,end=0.8,alpha=0.5,name="",aesthetics = "fill")+
  theme_classic(base_size=12, base_family="Open Sans Light")+
  theme(axis.title.y = element_text(angle=90, size=rel(1.5), face="plain"),
        axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5),vjust=-2, face="plain"),
        axis.text.x = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(1.5)),
        legend.title = element_text(size=rel(1.5)),
        legend.position = "none",
        #plot.margin = unit(c(6,5.5,20,10),"points"),
        panel.border = element_rect(color = "black", fill=NA, size=1),
        axis.line = element_line(size=0.5),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_line(size=0.5,linetype="dashed")) +
  scale_x_continuous(name="RF Distance to True Tree")+
  scale_y_discrete(name="Minor Allele Count") +
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom") 
#geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
#theme_ridges(line_size = 0.5, grid = FALSE, center_axis_labels=TRUE)
plot3 <- ggplot(data = results.mod, 
                aes(y=as.factor(missing),
                    x=RF.Dist.ML,
                    fill=int)) +
  #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  stat_halfeye(data= . %>% filter(int == "EXT"), aes(col=int), alpha=0.7, 
               adjust=0.9, side="top", slab_color = "black", slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=0.03),
               justification=-0.1, point_interval="median_qi", scale=0.6,
               normalize="panels") +
  stat_halfeye(data= . %>% filter(int == "INT"), aes(col=int), alpha=0.7, 
               adjust=0.9, side="bottom", slab_color = "black",slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=-0.03),
               justification=1.1, point_interval="median_qi", scale=0.6,
               normalize="panels") +
  scale_fill_manual(values=cols.int[c(1,3)],name="", aesthetics = "fill")+
  scale_color_manual(values=cols.int[c(1,3)],name="", aesthetics = "color")+
  #scale_fill_viridis_d(begin=0.2,end=0.8,alpha=0.5,name="",aesthetics = "fill")+
  theme_classic(base_size=12, base_family="Open Sans Light")+
  theme(axis.title.y = element_text(angle=90, size=rel(1.5), face="plain"),
        axis.text.y = element_text(size=rel(1.5)),
        axis.title.x = element_text(size=rel(1.5),vjust=-2, face="plain"),
        axis.text.x = element_text(size=rel(1.5)),
        legend.text = element_text(size=rel(2)),
        legend.title = element_text(size=rel(2)),
        legend.position = "none",
        #plot.margin = unit(c(6,5.5,20,10),"points"),
        panel.border = element_rect(color = "black", fill=NA, size=1),
        axis.line = element_line(size=0.5),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_line(size=0.5,linetype="dashed")) +
  scale_x_continuous(name="RF Distance to True Tree")+
  scale_y_discrete(name="Missing Data") +
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom") 
#geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
#theme_ridges(line_size = 0.5, grid = FALSE, center_axis_labels=TRUE)


print(plot2)

# saved at 800x1400px
ggarrange(plot2, plot3, ncol=1,labels="AUTO",font.label=list(size=24),
          common.legend=TRUE,legend = "top")
ggarrange(vars.plots.rf,plot2,plot3,ncol=1,labels="AUTO")
