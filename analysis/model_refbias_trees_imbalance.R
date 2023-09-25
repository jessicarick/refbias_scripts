library(lme4)
library(dplyr)
library(MuMIn)
library(car)
library(LMERConvenienceFunctions)
library(ggsci)
library(ggpubr)
library(wesanderson)
library(ggridges)
library(ggdist)

cols.int <- c("#005F73","gray80","#5FB89D")
cols.sig <- c("#E6B749","gray80","#6A8D4E")

here::i_am("analysis/model_refbias_trees_RF.R")

output <- "092321-output"
results.raxml <- read.csv(here("output","new",paste0(output,"-raxml.csv")),header=TRUE,row.names=1,sep=",")

## preparing data object
results.mod <- results.raxml[results.raxml$simulation > 15 & results.raxml$simulation < 26,]
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.numeric(as.character(results.mod$missing))
results.mod$maf <- as.numeric(as.character(results.mod$maf))
results.mod$int <- as.factor(results.mod$int)

#results.mod$int <- as.factor(matrix(unlist(regmatches(results.mod$taxa_ref, regexec('([A-Z]+)-', results.mod$taxa_ref))),
#                                    nrow=5760,ncol=2,byrow=TRUE)[,2])

## imbalance
# all

m.imb.all <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                    avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                    data = results.mod)
sum.imb.all <- summary(m.imb.all)
r.squaredGLMM(m.imb.all)

confint.imb.all<-data.frame(confint(m.imb.all)[-c(1:3),])
colnames(confint.imb.all) <- c("minCI","maxCI")
confint.imb.all$var <- rownames(confint.imb.all)
confint.imb.all$est <- sum.imb.all$coefficients[-1,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb.all$sig <- sum.imb.all$coefmat.full[-1,5]
confint.imb.all$sig <- case_when(confint.imb.all$minCI > 0 ~ "pos",
                                   confint.imb.all$maxCI < 0 ~ "neg",
                                   TRUE ~ "ns")
confint.imb.all$sig <- factor(confint.imb.all$sig, levels=c("pos","ns","neg"))

#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.all <- ggplot(confint.imb.all, aes(x = var, y = est, color=sig))
vars.imb.all.bars <- vars.imb.all + geom_blank() +
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
  #ylim(-200,1200)+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols.sig[c(3,2,1)])+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
print(vars.imb.all.bars)

# short

m.imb.short <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                      avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                    data = results.mod[results.mod$height == "SHORT" & results.mod$noref == "REF",])
sum.imb.short <- summary(m.imb.short)
r.squaredGLMM(m.imb.short)

confint.imb.short<-data.frame(confint(m.imb.short)[-c(1:3),])
colnames(confint.imb.short) <- c("minCI","maxCI")
confint.imb.short$var <- rownames(confint.imb.short)
confint.imb.short$est <- sum.imb.short$coefficients[-1,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb.short$sig <- sum.imb.short$coefmat.full[-1,5]
confint.imb.short$sig <- case_when(confint.imb.short$minCI > 0 ~ "pos",
                             confint.imb.short$maxCI < 0 ~ "neg",
                             TRUE ~ "ns")
confint.imb.short$sig <- factor(confint.imb.short$sig, levels=c("pos","ns","neg"))

#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.short <- ggplot(confint.imb.short, aes(x = var, y = est, color=sig))
vars.imb.short.bars <- vars.imb.short + geom_blank() +
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
  ylab("Coefficient\nShort Trees")+
  #ylim(-200,1200)+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols.sig[c(3,2,1)])+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
print(vars.imb.short.bars)

# med

m.imb.med <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                    avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                  data = results.mod[results.mod$height == "MED" & results.mod$noref == "REF",])
sum.imb.med <- summary(m.imb.med)
r.squaredGLMM(m.imb.med)

confint.imb.med<-data.frame(confint(m.imb.med)[-c(1:3),])
colnames(confint.imb.med) <- c("minCI","maxCI")
confint.imb.med$var <- rownames(confint.imb.med)
confint.imb.med$est <- sum.imb.med$coefficients[-1,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb.med$sig <- sum.imb.med$coefmat.full[-1,5]
confint.imb.med$sig <- case_when(confint.imb.med$minCI > 0 ~ "pos",
                             confint.imb.med$maxCI < 0 ~ "neg",
                             TRUE ~ "ns")
confint.imb.med$sig <- factor(confint.imb.med$sig, levels=c("pos","ns","neg"))

#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.med <- ggplot(confint.imb.med, aes(x = var, y = est, color=sig))
vars.imb.med.bars <- vars.imb.med + geom_blank() +
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
  ylab("Coefficient\nMedium Trees")+
  #ylim(-200,1200)+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.imb.med.bars

# long

m.imb.long <- lmer(std.ingroup.colless ~ avg_dxy + maf + missing +
                     avg_dxy:maf + avg_dxy:missing + (1 | simulation),
                   data = results.mod[results.mod$height == "LONG" & results.mod$noref == "REF",])
sum.imb.long <- summary(m.imb.long)
r.squaredGLMM(m.imb.long)

confint.imb.long<-data.frame(confint(m.imb.long)[-c(1:3),])
colnames(confint.imb.long) <- c("minCI","maxCI")
confint.imb.long$var <- rownames(confint.imb.long)
confint.imb.long$est <- sum.imb.long$coefficients[-1,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb.long$sig <- sum.imb.long$coefmat.full[-1,5]
confint.imb.long$sig <- case_when(confint.imb.long$minCI > 0 ~ "pos",
                             confint.imb.long$maxCI < 0 ~ "neg",
                             TRUE ~ "ns")
confint.imb.long$sig <- factor(confint.imb.long$sig, levels=c("pos","ns","neg"))

#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.long <- ggplot(confint.imb.long, aes(x = var, y = est, color=sig))
vars.imb.long.bars <-vars.imb.long + geom_blank() +
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
  ylab("Coefficient\nLong Trees")+
  #ylim(-200,1200)+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.imb.long.bars

# put them all together
vars.plots.imb <- ggarrange(vars.imb.short.bars,
          vars.imb.med.bars,
          vars.imb.long.bars,
          #labels=c("A","B","C"),
          ncol=3)
vars.plots.imb

## ridgeline plots
plot1 <- ggplot(data = results.mod, 
                aes(y=as.factor(maf),
                    x=std.ingroup.colless,
                    fill=int))

plot2 <- ggplot(data = results.mod, 
                aes(y=as.factor(maf),
                    x=std.ingroup.colless,
                    fill=int)) +
  #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  stat_halfeye(data= . %>% filter(int == "EXT"), aes(col=int), alpha=0.7, 
               adjust=1, side="top", slab_color = "black", slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=0.03),
               justification=-0.1, point_interval="median_qi", scale=0.9,
               normalize="panels") +
  stat_halfeye(data= . %>% filter(int == "INT"), aes(col=int), alpha=0.7, 
               adjust=1, side="bottom", slab_color = "black",slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=-0.03),
               justification=1.1, point_interval="median_qi", scale=0.9,
               normalize="panels") +
  scale_fill_manual(values=cols.int[c(1,3)],name="", aesthetics = "fill")+
  scale_color_manual(values=cols.int[c(1,3)],name="", aesthetics = "color")+
  #scale_fill_viridis_d(begin=0.2,end=0.8,alpha=0.5,name="",aesthetics = "fill")+
  theme_classic(base_size=12, base_family="Open Sans Light")+
  theme(axis.title.y = element_text(angle=90, size=rel(1.5), face="plain",
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
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
  scale_x_continuous(name="Standardized Ingroup Colless")+
  scale_y_discrete(name="Minor Allele Count") +
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom") 
#geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
#theme_ridges(line_size = 0.5, grid = FALSE, center_axis_labels=TRUE)
plot3 <- ggplot(data = results.mod, 
                aes(y=as.factor(missing),
                    x=std.ingroup.colless,
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
  theme(axis.title.y = element_text(angle=90, size=rel(1.5), face="plain",
                                    margin = margin(t = 0, r = 7, b = 0, l = 0)),
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
  scale_x_continuous(name="Standardized Ingroup Colless")+
  scale_y_discrete(name="Missing Data") +
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom") 
#geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
#theme_ridges(line_size = 0.5, grid = FALSE, center_axis_labels=TRUE)

# saved at 800x1400px
ggarrange(plot2, plot3, ncol=1,labels="AUTO",font.label=list(size=24),
          common.legend=TRUE,legend = "top")
ggarrange(vars.plots.imb,plot2,plot3,ncol=1,labels="AUTO")
