############################################
## Script for modeling relationships of refbias trees
## Copied by J. Rick, 17 sept 2019 from other script
## Last update 17 sept 2019
############################################

library(CCA)
library(CCP)
library(yacca)
library(sem)
library(dplyr)
library(MuMIn)
library(car)
library(lme4)
library(LMERConvenienceFunctions)
library(argparse)
library(tidyverse)
library(ggsci)
library(ggpubr)

## for debugging
output <- "091119-output"
results.raxml <- read.csv(paste("output/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

## for actual run
# parse.args <- function() {
#   parser <- ArgumentParser()
#   parser$add_argument('--ml.trees', help='file with ML (true) trees')
#   parser$add_argument('--ml.tree.names', help='file with names of ML trees')
#   parser$add_argument('--raxml.trees', help='file with raxml trees')
#   parser$add_argument('--raxml.tree.names', help='file with raxml tree names')
#   parser$add_argument('--astral.trees', help='file with astral trees')
#   parser$add_argument('--astral.tree.names', help='file with astral tree names')
#   parser$add_argument('-o', '--output', help='output file prefix')
#   return(parser$parse_args())
# }
# args <- parse.args()

##################################
## modeling! RF distance
#################################
mod.results <- data.frame(call=character(),r2m=numeric(),r2c=numeric())

## rf dist
results.mod <- results.raxml
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.factor(results.mod$missing)
results.mod$maf <- as.factor(results.mod$maf)

## univariate modeling
m.rf.height <- lm(RF.Dist.ML ~ height * simulation, data=results.mod)
summary(m.rf.height)
r.squaredGLMM(m.rf.height)

m.rf.quality <- lm(RF.Dist.ML ~ quality * simulation, data=results.mod)
summary(m.rf.quality)
r.squaredGLMM(m.rf.quality)

m.rf.missing <- lm(RF.Dist.ML ~ missing * simulation, data=results.mod)
summary(m.rf.missing)
r.squaredGLMM(m.rf.missing)

m.rf.maf <- lm(RF.Dist.ML ~ maf * simulation, data=results.mod)
summary(m.rf.maf)
r.squaredGLMM(m.rf.maf)

m.rf.sites <- lm(RF.Dist.ML ~ std.sites * simulation, data=results.mod)
summary(m.rf.sites)
r.squaredGLMM(m.rf.sites)

## multivariate, no interactions
# using LMER convienence functions #######################
m.rf <- lmer(RF.Dist.ML ~ (quality + missing + maf + refdist + height)
             + (1 | simulation), data = results.mod)

m.rf1 <- bfFixefLMER_F.fnc(m.rf, method="AIC",threshold=10)
pamer.fnc(m.rf1)
sum.rf <- summary(m.rf1)
r.squaredGLMM(m.rf1)

car::vif(m.rf)

# using MuMIn
options(na.action="na.fail")
m_list <- dredge(m.rf, rank="AIC",trace=1)
head(m_list, 10)

m.rf_list <- get.models(m_list, subset = delta < 10)

m.rf_avg <- model.avg(m.rf_list, revised.var = TRUE)
sum.rf <- summary(m.rf_avg)

confint.rf<-data.frame(confint(m.rf_avg,type="FULL"))
colnames(confint.rf) <- c("minCI","maxCI")
confint.rf$var <- rownames(confint.rf)
confint.rf$est <- sum.rf$coefficients[1,]
confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf$sig <- sum.rf$coefmat.full[,5]
confint.rf$sig2 <- case_when(confint.rf$sig < 0.05 & confint.rf$est > 0 ~ "pos",
                             confint.rf$sig < 0.05 & confint.rf$est < 0 ~ "neg",
                             TRUE ~ "ns")
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf <- ggplot(confint.rf, aes(x = var, y = est, color=sig2))
vars.rf + geom_blank() +
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
  ylab("Model Averaged Parameter Estimate")+
  scale_color_npg() +
  #scale_x_reverse() +
  #scale_color_manual(values=c("orchid","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))

vars <- (importance(m.rf_avg))
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))

imp.rf <- ggplot(vars2, aes(x = var, y = importance))
imp.rf + geom_point(size=4)+
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
  ylab("")+
  scale_color_npg()+
  #scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none",axis.title = element_text(size=15))+
  ylab("Summed AIC Weight")

## multivariate with interactions
# using LMER convienence functions #######################
m.rf.int <- lmer(RF.Dist.ML ~ (quality + missing + maf + refdist + height)*
               (quality + missing + maf + refdist + height) 
             + (1 | simulation), data = results.mod)

m.rf.int1 <- bfFixefLMER_F.fnc(m.rf.int, method="AIC",threshold=10)
pamer.fnc(m.rf.int1)
sum.rf <- summary(m.rf.int1)
r.squaredGLMM(m.rf.int1)

options(na.action="na.fail")
m_list_rf_int <- dredge(m.rf.int, rank="AIC",trace=TRUE)
head(m_list_rf_int, 20)

m.rf_avg_int <- model.avg(m_list_rf_int, revised.var = TRUE)
sum.rf_int <- summary(m.rf_avg_int)

confint.rf_int<-data.frame(confint(m.rf_avg_int,type="FULL"))
colnames(confint.rf_int) <- c("minCI","maxCI")
confint.rf_int$var <- rownames(confint.rf_int)
confint.rf_int$est <- sum.rf_int$coefficients[1,]
confint.rf_int$cond.est <- sum.rf_int$coefficients[2,]
confint.rf_int$sig <- sum.rf_int$coefmat.full[,5]
confint.rf_int$sig2 <- case_when(confint.rf_int$sig < 0.05 & confint.rf_int$est > 0 ~ "pos",
                             confint.rf_int$sig < 0.05 & confint.rf_int$est < 0 ~ "neg",
                             TRUE ~ "ns")
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65


vars.rf_int <- ggplot(confint.rf_int[confint.rf_int$sig != "ns",], aes(x = var, y = est, color=sig2))
vars.rf_int + geom_blank() +
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
  ylab("Model Averaged Parameter Estimate")+
  scale_color_npg() +
  #scale_x_reverse() +
  #scale_color_manual(values=c("orchid","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))

vars_int <- (importance(m.rf_avg_int))
vars2_int <- data.frame(importance=as.numeric(vars_int),var=names(vars_int))

imp.rf_int <- ggplot(vars2_int, aes(x = var, y = importance))
imp.rf_int + geom_point(size=4)+
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
  ylab("")+
  scale_color_npg()+
  #scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none",axis.title = element_text(size=15))+
  ylab("Summed AIC Weight")
