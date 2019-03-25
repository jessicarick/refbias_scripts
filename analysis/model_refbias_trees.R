############################################
## Script for modeling relationships of refbias trees
## Copied by J. Rick, 11 March 2019 from other script
############################################

##############################
## starting working with canonical correlation
## (JAR on my own)
## added some of chad's script
##############################
results <- read.table()

#check for the argparse package installed
if (!('argparse' %in% installed.packages()[, 'Package'])) {
  install.packages('argparse', repos='http://cran.rstudio.com/')
}

suppressPackageStartupMessages(library('argparse'))

parse.args <- function() {
  parser <- ArgumentParser()
  parser$add_argument('--ml.trees', help='file with ML (true) trees')
  parser$add_argument('--ml.tree.names', help='file with names of ML trees')
  parser$add_argument('--raxml.trees', help='file with raxml trees')
  parser$add_argument('--raxml.tree.names', help='file with raxml tree names')
  parser$add_argument('--astral.trees', help='file with astral trees')
  parser$add_argument('--astral.tree.names', help='file with astral tree names')
  parser$add_argument('-o', '--output', help='output file prefix')
  return(parser$parse_args())
}
args <- parse.args()

library(CCA)
library(CCP)
library(yacca)
library(sem)
library(dplyr)
library(MuMIn)
library(car)

results.num <- results
results.num$int <- as.numeric(results.num$int)
results.num$noref <- as.numeric(results.num$noref)
results.num$method <- as.numeric(results.num$method)
results.num$missing <- as.numeric(results.num$missing)

preds <- as.matrix(results.num[,c(2:8,10)])
resp <- as.matrix(results.num[,c(12:21,23:27)])

cca.out <- cca(preds,  resp[,-c(12)],  xcenter = TRUE, ycenter = TRUE, xscale = TRUE, yscale = TRUE)
cca.out
p.perm(preds, resp[,-c(10:11)])

par(mar = c(2,2,2,2), mfrow = c(1,2))
yacca::helio.plot(cca.out, cv = 1, x.name = "Predictors", y.name = "Response")
yacca::helio.plot(cca.out, cv = 2, x.name = "Predictors", y.name = "Response")
yacca::helio.plot(cca.out, cv = 3, x.name = "Predictors", y.name = "Response")
yacca::helio.plot(cca.out, cv = 4, x.name = "Predictors", y.name = "Response")

## a different CCA
# CCA <- CCA::cc(preds,resp)
# summary(CCA)
# plt.cc(CCA, var.label = TRUE, type="v")

##################################
## modeling! 
#################################

# ## rf dist
results.mod <- results
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$maf <- as.factor(results.mod$maf)
# 
# m.rf <- lmer(RF.Dist.ML ~ quality * missing * maf * int * noref * std.sites * height + (1 | simulation), data = results.mod)
# 
# m.rf1 <- bfFixefLMER_F.fnc(m.rf, method="AIC",threshold=5)
# pamer.fnc(m.rf1)
# summary(m.rf1)
# 
# options(na.action="na.fail")
# m_list <- dredge(m.rf, rank="AIC", fixed=(1 | simulation))
# head(m_list, 10)
# 
# ## gamma
# m.gam <- lmer(std.gamma ~ maf * int + (1 | simulation), data = results[results$noref == "REF",])
# 
# m.gam1 <- bfFixefLMER_F.fnc(m.gam, method="AIC",threshold=5)        
# pamer.fnc(m.gam1)
# summary(m.gam1)

# JMA Modeling ###########################################
# RF Distance ############################################
dat_rf <- subset(results.mod, select=c(RF.Dist.ML, quality, missing, maf, int, noref, sites, simulation, height, method))
dat_rf$missing <- as.factor(dat_rf$missing)

dat_astral <- dat_rf[dat_rf$method == "astral",]
dat_raxml <- dat_rf[dat_rf$method == "raxml",]

mod <- lm(RF.Dist.ML ~ .*., data=dat_rf)
formula(mod)

m.rf <- lm(RF.Dist.ML ~ quality + missing + maf + int + missing + noref, data = dat_raxml)
summary(m.rf)
r.squaredGLMM(m.rf)
car::vif(m.rf)

m.rf1 <- lmer(RF.Dist.ML ~ (quality + missing + maf + int + noref)#*(quality + missing + maf + int + noref + std.sites + height) 
              + (1 | simulation), data = dat_astral, REML=FALSE)

summary(m.rf1)
vif(m.rf1)
r.squaredGLMM(m.rf1)

m.rf2 <- lm(RF.Dist.ML ~ (quality + missing + maf + int + noref + height)*(quality + missing + maf + int + noref + height), 
            data = dat_raxml)
summary(m.rf2)
r.squaredGLMM(m.rf2)

options(na.action="na.fail")
m_list <- dredge(m.rf2, rank="AIC",trace=1)
head(m_list, 10)

m.rf_list <- get.models(m_list, subset = delta < 5)

m.rf_avg <- model.avg(m.rf_list, revised.var = TRUE)
sum.rf <- summary(m.rf_avg)

confint.rf<-data.frame(confint(m.rf_avg,method="full"))
colnames(confint.rf) <- c("minCI","maxCI")
confint.rf$var <- rownames(confint.rf)
confint.rf$est <- sum.rf$coefficients[1,]
confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf$sig <- sum.rf$coefmat.full[,5]
confint.rf$sig2 <- case_when(confint.rf$sig < 0.05 & confint.rf$est > 0 ~ "pos",
                             confint.rf$sig < 0.05 & confint.rf$est < 0 ~ "neg",
                             confint.rf$sig >= 0.05 ~ "ns")
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65
vars <- data.frame(sum.rf$importance)
vars$var <- rownames(vars)

vars.rf <- ggplot(confint.rf[-c(1),], aes(x = var, y = est, color=sig2))
vars.rf + geom_blank()+
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
  #scale_color_npg()+
  scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none")

imp.rf <- ggplot(vars, aes(x = var, y = sum.rf.importance))
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

## testing models

m.rf.1 <- lm(RF.Dist.ML ~ (quality + missing + maf + int + noref + std.sites + height)*(quality + missing + maf + int + noref + std.sites + height), 
             data = dat_rf)
summary(m.rf.1)
m.rf.2 <- lm(RF.Dist.ML ~ (missing + maf + int + noref + std.sites + height)*(missing + maf + int + noref + std.sites + height), 
             data = dat_rf)
summary(m.rf.2)
m.rf.3 <- lm(RF.Dist.ML ~ (maf + int + noref + std.sites + height)*(maf + int + noref + std.sites + height), 
             data = dat_rf)
summary(m.rf.3)
m.rf.4 <- lm(RF.Dist.ML ~ (maf + int + noref + height)*(maf + int + noref + height), 
             data = dat_rf)
summary(m.rf.4)
m.rf.5 <- lm(RF.Dist.ML ~ (maf + int + height)*(maf + int + height), 
             data = dat_rf)
summary(m.rf.5)
m.rf.6 <- lm(RF.Dist.ML ~ (maf + height)*(maf + height), 
             data = dat_rf)
summary(m.rf.6)
m.rf.7 <- lm(RF.Dist.ML ~ maf, 
             data = dat_rf)
summary(m.rf.7)
m.rf.8 <- lm(RF.Dist.ML ~ maf + height, 
             data = dat_rf)
summary(m.rf.8)
m.rf.9 <- lm(RF.Dist.ML ~ maf + height + int, 
             data = dat_rf)
summary(m.rf.9)
m.rf.10 <- lm(RF.Dist.ML ~ maf + height + int + int*maf, 
              data = dat_rf)
summary(m.rf.10)

my.models<-model.sel(m.rf.1,m.rf.2,m.rf.3,m.rf.4,m.rf.5,m.rf.6,m.rf.7,m.rf.8,m.rf.9,m.rf.10,rank=AIC)
my.models

options(na.action="na.fail")
m_list <- dredge(m.rf.3, rank="AIC")
head(m_list, 10)

m.rf_list <- get.models(m_list, subset = TRUE)

m.rf_avg <- model.avg(m.rf_list, revised.var = TRUE)
sum.rf <- summary(m.rf_avg)

sum.rf <- summary(m.rf.3)

confint.rf<-data.frame(confint(m.rf_avg))
colnames(confint.rf) <- c("minCI","maxCI")
confint.rf$var <- rownames(confint.rf)
confint.rf$est <- sum.rf$coefficients[1,]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf$sig <- sum.rf$coefmat.full[,5]
confint.rf$sig2 <- case_when(confint.rf$sig < 0.05 & confint.rf$est > 0 ~ "pos",
                             confint.rf$sig < 0.05 & confint.rf$est < 0 ~ "neg",
                             confint.rf$sig >= 0.05 ~ "ns")

vars <- data.frame(sum.rf$importance)
vars$var <- rownames(vars)

vars.rf <- ggplot(confint.rf[-c(1),], aes(x = var, y = est, color=sig2))
vars.rf + geom_blank()+
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
  scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none")

imp.rf <- ggplot(vars, aes(x = var, y = sum.rf.importance))
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

## playing around with the random effect in the model
beta <- vector(length = 20)
intercept <- vector(length=20)
for (i in 1:20){
  tmp <- summary(lm(RF.Dist.ML ~ maf, subset=(simulation == i),data=dat_rf))
  beta[i] <- tmp$coefficients[2,1]
  intercept[i] <- tmp$coefficients[1,1]
}
plot(beta,intercept)

m.rf1 <- lmer(RF.Dist.ML ~ quality + missing + maf + int + noref + std.sites + height
              + (1 | simulation), data = dat_rf, REML=FALSE)

summary(m.rf1)
vif(m.rf1)
r.squaredGLMM(m.rf1)

options(na.action="na.fail")
m_list <- dredge(m.rf1, rank="AIC",trace=1)
head(m_list, 10)

m.rf_list <- get.models(m_list, subset = TRUE)
m.rf_avg <- model.avg(m.rf_list, revised.var = TRUE)
sum.rf <- summary(m.rf_avg)

confint.rf<-data.frame(confint(m.rf_avg))
colnames(confint.rf) <- c("minCI","maxCI")
confint.rf$var <- rownames(confint.rf)
confint.rf$est <- sum.rf$coefficients[1,]
confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf$sig <- sum.rf$coefmat.full[,5]
confint.rf$sig2 <- case_when(confint.rf$sig < 0.05 & confint.rf$est > 0 ~ "pos",
                             confint.rf$sig < 0.05 & confint.rf$est < 0 ~ "neg",
                             confint.rf$sig >= 0.05 ~ "ns")
#confint.rf[6,c(1:2,4:5)] <- confint.rf[6,c(1:2,4:5)] - 65
vars <- data.frame(sum.rf$importance)
vars$var <- rownames(vars)

vars.rf <- ggplot(confint.rf[-c(1),], aes(x = var, y = est, color=sig2))
vars.rf + geom_blank()+
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
  #scale_color_npg()+
  scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none")

imp.rf <- ggplot(vars, aes(x = var, y = sum.rf.importance))
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