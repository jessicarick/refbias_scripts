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

## imbalance
# all

m.imb.all.astr <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                   avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod)
# m.imb.all.astr <- lmer(std.ingroup.colless ~ int + maf + missing +
#                         int:maf + int:missing + (1 | simulation),
#                       data = results.mod)
sum.imb.all.astr <- summary(m.imb.all.astr)
r.squaredGLMM(m.imb.all.astr)

confint.imb.all.astr<-data.frame(confint(m.imb.all.astr)[-c(1:3),])
colnames(confint.imb.all.astr) <- c("minCI","maxCI")
confint.imb.all.astr$var <- rownames(confint.imb.all.astr)
confint.imb.all.astr$est <- sum.imb.all.astr$coefficients[-1,1]
confint.imb.all.astr$sig <- sum.imb.all.astr$coefmat.full[-1,5]
confint.imb.all.astr$sig <- case_when(confint.imb.all.astr$minCI > 0 ~ "pos",
                                  confint.imb.all.astr$maxCI < 0 ~ "neg",
                                  TRUE ~ "ns")
confint.imb.all.astr$sig <- factor(confint.imb.all.astr$sig, levels=c("pos","ns","neg"))


vars.imb.all.astr <- ggplot(confint.imb.all.astr, aes(x = var, y = est, color=sig))
vars.imb.all.bars.astr <- vars.imb.all.astr + geom_blank() +
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
vars.imb.all.bars.astr

# short

m.imb.short.astr <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                     avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod[results.mod$height == "SHORT",])
sum.imb.short.astr <- summary(m.imb.short.astr)
r.squaredGLMM(m.imb.short.astr)

confint.imb.short.astr<-data.frame(confint(m.imb.short.astr)[-c(1:3),])
colnames(confint.imb.short.astr) <- c("minCI","maxCI")
confint.imb.short.astr$var <- rownames(confint.imb.short.astr)
confint.imb.short.astr$est <- sum.imb.short.astr$coefficients[-1,1]
#confint.imb.astr$cond.est <- sum.imb.astr$coefficients[2,]
confint.imb.short.astr$sig <- sum.imb.short.astr$coefmat.full[-1,5]
confint.imb.short.astr$sig <- case_when(confint.imb.short.astr$minCI > 0 ~ "pos",
                            confint.imb.short.astr$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.imb.short.astr$sig <- factor(confint.imb.short.astr$sig, levels=c("pos","ns","neg"))
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.imb.short.astr <- ggplot(confint.imb.short.astr, aes(x = var, y = est, color=sig))
vars.imb.short.bars.astr <- vars.imb.short.astr + geom_blank() +
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
vars.imb.short.bars.astr

# med

m.imb.med.astr <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                   avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                 data = results.mod[results.mod$height == "MED",])
sum.imb.med.astr <- summary(m.imb.med.astr)
r.squaredGLMM(m.imb.med.astr)

confint.imb.med.astr<-data.frame(confint(m.imb.med.astr)[-c(1:3),])
colnames(confint.imb.med.astr) <- c("minCI","maxCI")
confint.imb.med.astr$var <- rownames(confint.imb.med.astr)
confint.imb.med.astr$est <- sum.imb.med.astr$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.imb.med.astr$sig <- sum.imb.med.astr$coefmat.full[-1,5]
confint.imb.med.astr$sig <- case_when(confint.imb.med.astr$minCI > 0 ~ "pos",
                            confint.imb.med.astr$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.imb.med.astr$sig <- factor(confint.imb.med.astr$sig, levels=c("pos","ns","neg"))

#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.imb.med.astr <- ggplot(confint.imb.med.astr, aes(x = var, y = est, color=sig))
vars.imb.med.bars.astr <- vars.imb.med.astr + geom_blank() +
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
vars.imb.med.bars.astr

# long
m.imb.long.astr <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                    avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                  data = results.mod[results.mod$height == "LONG",])
sum.imb.long.astr <- summary(m.imb.long.astr)
r.squaredGLMM(m.imb.long.astr)

confint.imb.long.astr<-data.frame(confint(m.imb.long.astr)[-c(1:3),])
colnames(confint.imb.long.astr) <- c("minCI","maxCI")
confint.imb.long.astr$var <- rownames(confint.imb.long.astr)
confint.imb.long.astr$est <- sum.imb.long.astr$coefficients[-1,1]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.imb.long.astr$sig <- sum.imb.long.astr$coefmat.full[-1,5]
confint.imb.long.astr$sig <- case_when(confint.imb.long.astr$minCI > 0 ~ "pos",
                            confint.imb.long.astr$maxCI < 0 ~ "neg",
                            TRUE ~ "ns")
confint.imb.long$sig <- factor(confint.imb.long.astr$sig, levels=c("pos","ns","neg"))

#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.imb.long.astr <- ggplot(confint.imb.long.astr, aes(x = var, y = est, color=sig))
vars.imb.long.bars.astr <- vars.imb.long.astr + geom_blank() +
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
vars.imb.long.bars.astr

# put them all together
vars.plots.imb.astr <- ggarrange(vars.imb.long.bars.astr,
          vars.imb.med.bars.astr,
          vars.imb.short.bars.astr,
          #labels=c("A","B","C"),
          ncol=3)
print(vars.plots.imb.astr)
