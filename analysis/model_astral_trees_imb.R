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
results.astral <- read.csv(here("output","new",paste0(output,"-astral.csv")),header=TRUE,row.names=1,sep=",")

## preparing data object
results.mod <- results.astral[results.astral$simulation < 26,]

results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.numeric(as.character(results.mod$missing))
results.mod$maf <- as.numeric(as.character(results.mod$maf))
results.mod$int <- as.factor(results.mod$int)
results.mod$method <- as.factor(results.mod$method)

## rf distance
# all

m.rf.all.astr <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                   avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod)
m.rf.all.astr <- lmer(RF.Dist.ML ~ int + maf + missing +
                        int:maf + int:missing + (1 | simulation),
                      data = results.mod)
sum.rf.all.astr <- summary(m.rf.all.astr)
r.squaredGLMM(m.rf.all.astr)

confint.rf.all.astr<-data.frame(confint(m.rf.all.astr)[-c(1:3),])
colnames(confint.rf.all.astr) <- c("minCI","maxCI")
confint.rf.all.astr$var <- rownames(confint.rf.all.astr)
confint.rf.all.astr$est <- sum.rf.all.astr$coefficients[-1,1]
confint.rf.all.astr$sig <- sum.rf.all.astr$coefmat.full[-1,5]
confint.rf.all.astr$sig <- case_when(confint.rf.all.astr$minCI > 0 ~ "pos",
                                  confint.rf.all.astr$maxCI < 0 ~ "neg",
                                  TRUE ~ "ns")
confint.rf.all.astr$sig <- factor(confint.rf.all.astr$sig, levels=c("pos","ns","neg"))


vars.rf.all.astr <- ggplot(confint.rf.all.astr, aes(x = var, y = est, color=sig))
vars.rf.all.bars.astr <- vars.rf.all.astr + geom_blank() +
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
vars.rf.all.bars.astr

# short

m.rf.short.astr <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                     avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod[results.mod$height == "SHORT",])
sum.rf.short.astr <- summary(m.rf.short.astr)
r.squaredGLMM(m.rf.short.astr)

confint.rf.short.astr<-data.frame(confint(m.rf.short.astr)[-c(1:3),])
colnames(confint.rf.short.astr) <- c("minCI","maxCI")
confint.rf.short.astr$var <- rownames(confint.rf.short.astr)
confint.rf.short.astr$est <- sum.rf.short.astr$coefficients[-1,1]
#confint.rf.astr$cond.est <- sum.rf.astr$coefficients[2,]
confint.rf.short.astr$sig <- sum.rf.short.astr$coefmat.full[-1,5]
confint.rf.short.astr$sig <- case_when(confint.rf.short.astr$minCI > 0 ~ "pos",
                            confint.rf.short.astr$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.rf.short.astr$sig <- factor(confint.rf.short.astr$sig, levels=c("pos","ns","neg"))
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.short.astr <- ggplot(confint.rf.short.astr, aes(x = var, y = est, color=sig))
vars.rf.short.bars.astr <- vars.rf.short.astr + geom_blank() +
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
vars.rf.short.bars.astr

# med

m.rf.med.astr <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                   avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                 data = results.mod[results.mod$height == "MED",])
sum.rf.med.astr <- summary(m.rf.med.astr)
r.squaredGLMM(m.rf.med.astr)

confint.rf.med.astr<-data.frame(confint(m.rf.med.astr)[-c(1:3),])
colnames(confint.rf.med.astr) <- c("minCI","maxCI")
confint.rf.med.astr$var <- rownames(confint.rf.med.astr)
confint.rf.med.astr$est <- sum.rf.med.astr$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf.med.astr$sig <- sum.rf.med.astr$coefmat.full[-1,5]
confint.rf.med.astr$sig <- case_when(confint.rf.med.astr$minCI > 0 ~ "pos",
                            confint.rf.med.astr$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.rf.med.astr$sig <- factor(confint.rf.med.astr$sig, levels=c("pos","ns","neg"))

#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.med.astr <- ggplot(confint.rf.med.astr, aes(x = var, y = est, color=sig))
vars.rf.med.bars.astr <- vars.rf.med.astr + geom_blank() +
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
vars.rf.med.bars.astr

# long
m.rf.long.astr <- lmer(RF.Dist.ML ~ avg_dxy + maf + missing +
                    avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                  data = results.mod[results.mod$height == "LONG",])
sum.rf.long.astr <- summary(m.rf.long.astr)
r.squaredGLMM(m.rf.long.astr)

confint.rf.long.astr<-data.frame(confint(m.rf.long.astr)[-c(1:3),])
colnames(confint.rf.long.astr) <- c("minCI","maxCI")
confint.rf.long.astr$var <- rownames(confint.rf.long.astr)
confint.rf.long.astr$est <- sum.rf.long.astr$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf.long.astr$sig <- sum.rf.long.astr$coefmat.full[-1,5]
confint.rf.long.astr$sig <- case_when(confint.rf.long.astr$minCI > 0 ~ "pos",
                            confint.rf.long.astr$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.rf.long$sig <- factor(confint.rf.long.astr$sig, levels=c("pos","ns","neg"))

#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf.long.astr <- ggplot(confint.rf.long.astr, aes(x = var, y = est, color=sig))
vars.rf.long.bars.astr <- vars.rf.long.astr + geom_blank() +
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
vars.rf.long.bars.astr

# put them all together
vars.plots.rf.astr <- ggarrange(vars.rf.long.bars.astr,
          vars.rf.med.bars.astr,
          vars.rf.short.bars.astr,
          #labels=c("A","B","C"),
          ncol=3)
print(vars.plots.rf.astr)
