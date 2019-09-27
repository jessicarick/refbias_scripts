library(lme4)
library(dplyr)
library(MuMIn)
library(car)
library(LMERConvenienceFunctions)
library(ggsci)
library(ggpubr)
library(wesanderson)
library(ggridges)

cols <- c("#F2AD00","gray80","#00A08A")

output <- "091819-output"
results.raxml <- read.csv(paste("output/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

## preparing data object
results.mod <- results.raxml
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.factor(results.mod$missing)
results.mod$maf <- as.factor(results.mod$maf)

#results.mod$int <- as.factor(matrix(unlist(regmatches(results.mod$taxa_ref, regexec('([A-Z]+)-', results.mod$taxa_ref))),
#                                    nrow=5760,ncol=2,byrow=TRUE)[,2])

## imbalance
# short

m.imb.short <- lmer(ingroup.sackin ~ int + maf + missing + quality +
                      int:maf + int:missing + int:quality + (1 | simulation),
                    data = results.mod[results.mod$height == "500000" & results.mod$noref == "REF",])
sum.imb.short <- summary(m.imb.short)
r.squaredGLMM(m.imb.short)

confint.imb<-data.frame(confint(m.imb.short)[-c(1:2),])
colnames(confint.imb) <- c("minCI","maxCI")
confint.imb$var <- rownames(confint.imb)
confint.imb$est <- sum.imb.short$coefficients[,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb$sig <- sum.imb.short$coefmat.full[,5]
confint.imb$sig <- case_when(confint.imb$minCI > 0 ~ "pos",
                             confint.imb$maxCI < 0 ~ "neg",
                             TRUE ~ "ns")
#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.short <- ggplot(confint.imb, aes(x = var, y = est, color=sig))
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
  ylab("Model Averaged Parameter Estimate\nHeight = 500,000")+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
print(vars.imb.short.bars)

# med

m.imb.med <- lmer(ingroup.sackin ~ int + maf + missing + quality +
                    int:maf + int:missing + int:quality + (1 | simulation),
                  data = results.mod[results.mod$height == "2000000"  & results.mod$noref == "REF",])
sum.imb.med <- summary(m.imb.med)
r.squaredGLMM(m.imb.med)

confint.imb<-data.frame(confint(m.imb.med)[-c(1:2),])
colnames(confint.imb) <- c("minCI","maxCI")
confint.imb$var <- rownames(confint.imb)
confint.imb$est <- sum.imb.med$coefficients[,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb$sig <- sum.imb.med$coefmat.full[,5]
confint.imb$sig <- case_when(confint.imb$minCI > 0 ~ "pos",
                             confint.imb$maxCI < 0 ~ "neg",
                             TRUE ~ "ns")
#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.med <- ggplot(confint.imb, aes(x = var, y = est, color=sig))
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
  ylab("Model Averaged Parameter Estimate\nHeight = 2,000,000")+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.imb.med.bars

# long

m.imb.long <- lmer(ingroup.sackin ~ int + maf + missing + quality +
                     int:maf + int:missing + int:quality + (1 | simulation),
                   data = results.mod[results.mod$height == "10000000" & results.mod$noref == "REF",])
sum.imb.long <- summary(m.imb.long)
r.squaredGLMM(m.imb.long)

confint.imb<-data.frame(confint(m.imb.long)[-c(1:2),])
colnames(confint.imb) <- c("minCI","maxCI")
confint.imb$var <- rownames(confint.imb)
confint.imb$est <- sum.imb.long$coefficients[,1]
#confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb$sig <- sum.imb.long$coefmat.full[,5]
confint.imb$sig <- case_when(confint.imb$minCI > 0 ~ "pos",
                             confint.imb$maxCI < 0 ~ "neg",
                             TRUE ~ "ns")
#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65


vars.imb.long <- ggplot(confint.imb, aes(x = var, y = est, color=sig))
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
  ylab("Model Averaged Parameter Estimate\nHeight = 10,000,000")+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.imb.long.bars

# put them all together
ggarrange(vars.imb.short.bars,
          vars.imb.med.bars,
          vars.imb.long.bars,
          labels=c("A","B","C"),
          ncol=3)

## ridgeline plots
plot1 <- ggplot(data = results.raxml[results.raxml$noref == "REF",], 
                aes(y=as.factor(results.raxml$maf[results.raxml$noref == "REF"]),
                    x=results.raxml[results.raxml$noref == "REF","ingroup.sackin"],
                    fill=results.raxml$int[results.raxml$noref == "REF"]))

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
  scale_x_continuous(name="Std Ingroup Sackin")+
  scale_y_discrete(name="MAF")+
  #xlim(-50,10)+
  facet_wrap(vars(height),nrow=1,strip.position = "bottom")+
  geom_hline(yintercept=0,cex=2,lty=2,col="gray")+
  theme_ridges(line_size = 1, grid = TRUE, center_axis_labels=TRUE)

print(plot2)