library(lme4)
library(tidyverse)
library(MuMIn)
library(car)
library(LMERConvenienceFunctions)
library(ggsci)
library(ggpubr)
library(ggridges)

cols <- c("#F2AD00","gray80","#00A08A")

output <- "072221-output"
results.raxml <- read.csv(paste("output/new/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

## preparing data object
results.mod <- results.raxml[results.raxml$simulation < 11,]
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.factor(results.mod$missing)
results.mod$maf <- as.factor(results.mod$maf)
results.mod$int <- as.factor(results.mod$int)

## rf distance
# short

m.rf.short <- lmer(RF.Dist.ML ~ int + maf + missing +
                     int:maf + int:missing + (1 | simulation),
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
  ylab("Coefficient")+
  #xlim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.rf.short.bars

# med

m.rf.med <- lmer(RF.Dist.ML ~ int + maf + missing +
                   int:maf + int:missing + (1 | simulation),
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
  ylab("Coefficient")+
  #xlim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.rf.med.bars

# long
m.rf.long <- lmer(RF.Dist.ML ~ int + maf + missing +
                    int:maf + int:missing + (1 | simulation),
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
  ylab("Coefficient")+
  #ylim(-20,120)+
  #scale_color_npg() +
  #scale_x_reverse() +
  #scale_y_continuous(limits=c(-25,150))+
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.rf.long.bars

# put them all together
ggarrange(vars.rf.short.bars,
          vars.rf.med.bars,
          vars.rf.long.bars,
          labels=c("A","B","C"),
          ncol=3)

## ridgeline plots
plot1 <- ggplot(data = results.mod, 
                aes(y=as.factor(maf),
                    x=RF.Dist.ML,
                    fill=int))

plot2 <- plot1 +
  geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  scale_fill_manual(values=cols[c(1,3)],name="",aesthetics = "fill")+
  #scale_fill_viridis_d(begin=0.2,end=0.8,alpha=0.5,name="",aesthetics = "fill")+
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
  scale_x_continuous(name="RF Distance to True Tree")+
  scale_y_discrete(name="MAF")+
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
  geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
  theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)

print(plot2)
