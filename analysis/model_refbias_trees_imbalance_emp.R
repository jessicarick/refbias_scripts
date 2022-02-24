############################################
## Script for analyzing gamma empirial datasets
## Written by J. Rick
## Updated 18 Feb 2022 to use here package
############################################

library(lme4)
library(dplyr)
library(MuMIn)
library(car)
library(LMERConvenienceFunctions)
library(ggsci)
library(ggpubr)
library(wesanderson)
library(ggridges)
library(here)

cols <- c("#F2AD00","gray80","#00A08A")

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


## imbalance
# lates

m.imb.lates <- lmer(ingroup.colless ~ int + maf + missing +
                      int:maf + int:missing + (1 | simulation),
                       data = results.emp.lates)
sum.imb.lates <- summary(m.imb.lates)
r.squaredGLMM(m.imb.lates)

confint.imb.lates<-data.frame(confint(m.imb.lates))[-c(1:3),]
colnames(confint.imb.lates) <- c("minCI","maxCI")
confint.imb.lates$var <- rownames(confint.imb.lates)
confint.imb.lates$est <- sum.imb.lates$coefficients[-1,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb.lates$sig <- sum.imb.lates$coefmat.full[-1,5]
confint.imb.lates$sig <- case_when(confint.imb.lates$minCI > 0 ~ "pos",
                                      confint.imb.lates$maxCI < 0 ~ "neg",
                                      TRUE ~ "ns")
confint.imb.lates$sig <- factor(confint.imb.lates$sig, levels=c("pos","ns","neg"))


vars.imb.lates <- ggplot(confint.imb.lates, aes(x = var, y = est, color=sig))
vars.imb.lates.bars <- vars.imb.lates + geom_blank() +
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
vars.imb.lates.bars

## imbalance
# cichlids

m.imb.cichlids <- lmer(ingroup.colless ~ int + maf + missing +
                            int:maf + int:missing + (1 | simulation),
                          data = results.emp.cichlids)
sum.imb.cichlids <- summary(m.imb.cichlids)
r.squaredGLMM(m.imb.cichlids)

confint.imb.cichlids<-data.frame(confint(m.imb.cichlids))[-c(1:3),]
colnames(confint.imb.cichlids) <- c("minCI","maxCI")
confint.imb.cichlids$var <- rownames(confint.imb.cichlids)
confint.imb.cichlids$est <- sum.imb.cichlids$coefficients[-1,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb.cichlids$sig <- sum.imb.cichlids$coefmat.full[-1,5]
confint.imb.cichlids$sig <- case_when(confint.imb.cichlids$minCI > 0 ~ "pos",
                                         confint.imb.cichlids$maxCI < 0 ~ "neg",
                                         TRUE ~ "ns")
confint.imb.cichlids$sig <- factor(confint.imb.cichlids$sig, levels=c("pos","ns","neg"))


vars.imb.cichlids <- ggplot(confint.imb.cichlids, aes(x = var, y = est, color=sig))
vars.imb.cichlids.bars <- vars.imb.cichlids + geom_blank() +
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
vars.imb.cichlids.bars

## ridgeline plots
plot1 <- ggplot(data = results.emp.lates, 
                aes(y=as.factor(maf),
                    x=ingroup.colless,
                    fill=int))

plot2 <- ggplot(data = results.emp.lates, 
                aes(y=as.factor(maf),
                    x=ingroup.colless,
                    fill=int)) +
  #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  stat_halfeye(data= . %>% filter(int == "EXT"), aes(col=int), alpha=0.7, 
               adjust=1, side="top", slab_color = "black", slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=0.03),
               justification=-0.1, point_interval="median_qi", scale=0.9) +
  stat_halfeye(data= . %>% filter(int == "INT"), aes(col=int), alpha=0.7, 
               adjust=1, side="bottom", slab_color = "black",slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=-0.03),
               justification=1.1, point_interval="median_qi", scale=0.9) +
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
  scale_x_continuous(name="Ingroup Colless")+
  scale_y_discrete(name="Minor Allele Count") 

plot3 <- ggplot(data = results.emp.lates, 
                aes(y=as.factor(missing),
                    x=ingroup.colless,
                    fill=int)) +
  #geom_density_ridges(scale = 0.95, rel_min_height = 0.1, alpha = 0.5)+
  stat_halfeye(data= . %>% filter(int == "EXT"), aes(col=int), alpha=0.7, 
               adjust=1, side="top", slab_color = "black", slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=0.03),
               justification=-0.1, point_interval="median_qi", scale=0.9) +
  stat_halfeye(data= . %>% filter(int == "INT"), aes(col=int), alpha=0.7, 
               adjust=1, side="bottom", slab_color = "black",slab_size=0.2,
               interval_alpha=0.7, point_alpha=0.8, position=position_nudge(y=-0.03),
               justification=1.1, point_interval="median_qi", scale=0.9) +
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
  scale_x_continuous(name="Ingroup Colless")+
  scale_y_discrete(name="Missing Data")

imb.plot <- ggarrange(plot2, plot3, ncol=1,labels=c("C","D"),font.label=list(size=24))
ggarrange(vars.plots.imb,plot2,plot3,ncol=1,labels="AUTO")

# exported at 1000x800px
ggarrange(gam.plot,imb.plot,nrow=1)
