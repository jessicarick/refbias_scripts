############################################
## Script for analyzing gamma empirial datasets
## Written by J. Rick
## Updated 18 Feb 2022 to use here package
############################################

library(lme4)
library(tidyverse)
library(MuMIn)
library(car)
library(LMERConvenienceFunctions)
library(ggsci)

cols <- c("#F2AD00","gray80","#00A08A")

date <- "010822"
results.emp.lates <- read.csv(here("output","new",paste(date,"-lates-emp-output-raxml.csv",sep="")),header=TRUE,row.names=1,sep=",")
results.emp.cichlids <- read.csv(here("output","new",paste(date,"-cichlids-emp-output-raxml.csv",sep="")),header=TRUE,row.names=1,sep=",")

## preparing data object
#results.mod <- results.raxml
# results.emp.lates$simulation <- as.factor(results.emp.lates$simulation)
# results.emp.lates$quality <- as.factor(results.emp.lates$quality)
# results.emp.lates$missing <- as.factor(results.emp.lates$missing)
# results.emp.lates$maf <- as.factor(results.emp.lates$maf)
# 
# results.emp.cichlids$simulation <- as.factor(results.emp.cichlids$simulation)
# results.emp.cichlids$quality <- as.factor(results.emp.cichlids$quality)
# results.emp.cichlids$missing <- as.factor(results.emp.cichlids$missing)
# results.emp.cichlids$maf <- as.factor(results.emp.cichlids$maf)

## Height
# lates

m.height.lates <- lmer(ingroup.tree.height ~ int + maf + missing +
                         int:maf + int:missing + (1 | simulation),
                    data = results.emp.lates)
sum.height.lates <- summary(m.height.lates)
r.squaredGLMM(m.height.lates)

confint.height.lates<-data.frame(confint(m.height.lates))[-c(1:3),]
colnames(confint.height.lates) <- c("minCI","maxCI")
confint.height.lates$var <- rownames(confint.height.lates)
confint.height.lates$est <- sum.height.lates$coefficients[-1,1]
#confint.height$cond.est <- sum.height$coefficients[2,]
confint.height.lates$sig <- sum.height.lates$coefmat.full[-1,5]
confint.height.lates$sig <- case_when(confint.height.lates$minCI > 0 ~ "pos",
                                   confint.height.lates$maxCI < 0 ~ "neg",
                                   TRUE ~ "ns")
confint.height.lates$sig <- factor(confint.height.lates$sig, levels=c("pos","ns","neg"))


vars.height.lates <- ggplot(confint.height.lates, aes(x = var, y = est, color=sig))
vars.height.lates.bars <- vars.height.lates + geom_blank() +
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
  ylab("Coefficient\nlates Trees (High ILS)")+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1,shape=15) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) +
  theme(axis.title.x=element_blank())
vars.height.lates.bars

## Height
# cichlids

m.height.cichlids <- lmer(ingroup.tree.height ~ int + maf + missing +
                         int:maf + int:missing + (1 | simulation),
                       data = results.emp.cichlids)
sum.height.cichlids <- summary(m.height.cichlids)
r.squaredGLMM(m.height.cichlids)

confint.height.cichlids<-data.frame(confint(m.height.cichlids))[-c(1:3),]
colnames(confint.height.cichlids) <- c("minCI","maxCI")
confint.height.cichlids$var <- rownames(confint.height.cichlids)
confint.height.cichlids$est <- sum.height.cichlids$coefficients[-1,1]
#confint.height$cond.est <- sum.height$coefficients[2,]
confint.height.cichlids$sig <- sum.height.cichlids$coefmat.full[-1,5]
confint.height.cichlids$sig <- case_when(confint.height.cichlids$minCI > 0 ~ "pos",
                                      confint.height.cichlids$maxCI < 0 ~ "neg",
                                      TRUE ~ "ns")
confint.height.cichlids$sig <- factor(confint.height.cichlids$sig, levels=c("pos","ns","neg"))


vars.height.cichlids <- ggplot(confint.height.cichlids, aes(x = var, y = est, color=sig))
vars.height.cichlids.bars <- vars.height.cichlids + geom_blank() +
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
  ylab("Coefficient\ncichlids Trees (High ILS)")+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1,shape=15) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) +
  theme(axis.title.x=element_blank())
vars.height.cichlids.bars

## ridgeline plots
plot1 <- ggplot(data = results.emp.lates, 
                aes(y=as.factor(maf),
                    x=ingroup.tree.height,
                    fill=int))

plot2 <- ggplot(data = results.emp.lates, 
                aes(y=as.factor(maf),
                    x=ingroup.tree.height,
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
  scale_x_continuous(name="Ingroup Height")+
  scale_y_discrete(name="Minor Allele Count") +
  #xlim(-50,10)+
  #facet_wrap(vars(height),nrow=1,strip.position = "bottom") +
  geom_vline(xintercept=0,cex=0.8,lty=2,col="gray")
#theme_ridges(line_size = 0.5, grid = FALSE, center_axis_labels=TRUE)
plot3 <- ggplot(data = results.mod, 
                aes(y=as.factor(missing),
                    x=ingroup.tree.height,
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
        legend.text = element_text(size=rel(2)),
        legend.title = element_text(size=rel(2)),
        legend.position = "none",
        #plot.margin = unit(c(6,5.5,20,10),"points"),
        panel.border = element_rect(color = "black", fill=NA, size=1),
        axis.line = element_line(size=0.5),
        strip.text.x = element_text(size = 16),
        panel.grid.major = element_line(size=0.5,linetype="dashed")) +
  scale_x_continuous(name="Ingroup Height")+
  scale_y_discrete(name="Missing Data") +
  #xlim(-50,10)+
  #facet_wrap(vars(height),nrow=1,strip.position = "bottom") +
  geom_vline(xintercept=0,cex=0.8,lty=2,col="gray")
#theme_ridges(line_size = 0.5, grid = FALSE, center_axis_labels=TRUE)

print(plot2)

height.plot <- ggarrange(plot2, plot3, ncol=1,labels="AUTO",font.label=list(size=24))
ggarrange(vars.plots.height,plot2,plot3,ncol=1,labels="AUTO")

