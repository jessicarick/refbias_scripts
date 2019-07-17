############################################
## Script for modeling relationships of refbias trees
## Copied by J. Rick, 11 March 2019 from other script
############################################

##############################
## starting working with canonical correlation
## (JAR on my own)
## added some of chad's script
##############################

output <- "053019-output"
results.raxml <- read.csv(paste("output/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

# #check for the argparse package installed
# if (!('argparse' %in% installed.packages()[, 'Package'])) {
#   install.packages('argparse', repos='http://cran.rstudio.com/')
# }
# 
# suppressPackageStartupMessages(library('argparse'))
# 
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

library(CCA)
library(CCP)
library(yacca)
library(sem)
library(dplyr)
library(MuMIn)
library(car)
library(lme4)
library(LMERConvenienceFunctions)

results.num <- results.raxml
results.num$int <- as.numeric(as.factor(results.num$int))
results.num$noref <- as.numeric(as.factor(results.num$noref))
results.num$method <- as.numeric(as.factor(results.num$method))
results.num$missing <- as.numeric(results.num$missing)

preds <- as.matrix(results.num[,c(2:8,10,12)])
resp <- as.matrix(results.num[,c(13:28)])

cca.out <- cca(preds[,-c(2)],  resp[,-c(10:13)],  xcenter = TRUE, ycenter = TRUE, xscale = TRUE, yscale = TRUE)
cca.out
p.perm(preds[,-2], resp[,-c(10:12)])

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
## modeling! RF distance
#################################

# ## rf dist
results.mod <- results.raxml
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.factor(results.mod$missing)
results.mod$maf <- as.factor(results.mod$maf)
 

# CDB Modeling ###########################################
# using LMER convienence functions #######################
m.rf <- lmer(RF.Dist.ML ~ (quality + missing + maf + refdist + noref + height)*
               (quality + missing + maf + refdist + noref + height) 
             + (1 | simulation), data = results.mod)

m.rf1 <- bfFixefLMER_F.fnc(m.rf, method="AIC",threshold=5)
pamer.fnc(m.rf1)
sum.rf <- summary(m.rf1)
r.squaredGLMM(m.rf1)

options(na.action="na.fail")
#clust <- try(makeCluster(getOption("cl.cores", 4), type = "SOCK"))
m_list_rf <- dredge(m.rf, rank="AIC",trace=TRUE)
head(m_list_rf, 10)
subset(m_list_rf, delta < 10)
plot(m_list_rf)

## gamma
results.ref <- results.mod[results.mod$noref == "REF",]
m.gam <- lmer(std.gamma ~ (quality + missing + maf + refdist + height)*
                (quality + missing + maf + refdist + height) + (1 | simulation), 
              data = results.ref)

m.gam1 <- bfFixefLMER_F.fnc(m.gam, method="AIC",threshold=5)
pamer.fnc(m.gam1)
summary(m.gam1)
r.squaredGLMM(m.gam1)

options(na.action="na.fail")
m_list_gam <- dredge(m.gam, rank="AIC")
head(m_list_gam, 10)

## colless
m.imb <- lmer(ingroup.colless ~ (quality + missing + maf + refdist + height)*
                (quality + missing + maf + refdist + height) + (1 | simulation), 
              data = results.ref)

m.imb1 <- bfFixefLMER_F.fnc(m.imb, method="AIC",threshold=5)
pamer.fnc(m.imb1)
summary(m.imb1)
r.squaredGLMM(m.imb1)

options(na.action="na.fail")
m_list_imb <- dredge(m.imb, rank="AIC")
head(m_list_imb, 10)
subset(m_list_imb, delta < 10)
plot(m_list_imb)

##########################################################
# JMA Modeling ###########################################
# RF Distance ############################################
##########################################################
dat_rf <- subset(results.mod, select=c(RF.Dist.ML, quality, missing, maf, int, noref, std.sites, simulation, height, refdist))
dat_rf$missing <- as.factor(dat_rf$missing)
dat_rf$maf <- as.factor(dat_rf$maf)
dat_rf_ref <- dat_rf[dat_rf$noref == "REF",]

dat_astral <- dat_rf[dat_rf$method == "astral",]
dat_raxml <- dat_rf[dat_rf$method == "raxml",]

mod <- lm(RF.Dist.ML ~ .*., data=dat_rf)
formula(mod)
summary(mod)
r.squaredGLMM(mod)

m.rf <- lm(RF.Dist.ML ~ quality + maf + missing + noref + std.sites + height + refdist, data = dat_rf)
summary(m.rf)
r.squaredGLMM(m.rf)
car::vif(m.rf)

m.rf1 <- lmer(RF.Dist.ML ~ (quality + missing + maf + refdist + noref + std.sites + height)#*(quality + missing + maf + int + noref + std.sites + height) 
              + (1 | simulation), data = dat_rf, REML=FALSE)

summary(m.rf1)
vif(m.rf1)
r.squaredGLMM(m.rf1)

m.rf2 <- lmer(RF.Dist.ML ~ (quality + missing + maf + refdist + height + missing)*
                (quality + missing + maf + refdist + height + missing) + (1 | simulation), 
            data = dat_rf_ref)
summary(m.rf2)
print(summary(m.rf2))
vif(m.rf2)
r.squaredGLMM(m.rf2)

options(na.action="na.fail")
m_list <- dredge(m.rf2, rank="AIC",trace=1)
head(m_list, 10)

m.rf_list <- get.models(m_list, subset = cumsum(weight) < 0.999)
m.rf_list_short <- get.models(m_list, subset = delta < 10)

m.rf_avg <- model.avg(m.rf_list_short, revised.var = TRUE)
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

## testing models

m.rf.1 <- lm(RF.Dist.ML ~ (quality + missing + maf + refdist + noref + std.sites + height)*(quality + missing + maf + refdist + noref + std.sites + height), 
             data = dat_rf)
summary(m.rf.1)
m.rf.2 <- lm(RF.Dist.ML ~ (missing + maf + refdist + noref + std.sites + height)*(missing + maf + refdist + noref + std.sites + height), 
             data = dat_rf)
summary(m.rf.2)
m.rf.3 <- lm(RF.Dist.ML ~ (maf + refdist + noref + std.sites + height)*(maf + refdist + noref + std.sites + height), 
             data = dat_rf)
summary(m.rf.3)
m.rf.4 <- lm(RF.Dist.ML ~ (maf + refdist + noref + height)*(maf + refdist + noref + height), 
             data = dat_rf)
summary(m.rf.4)
m.rf.5 <- lm(RF.Dist.ML ~ (maf + refdist + height)*(maf + refdist + height), 
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
m.rf.9 <- lm(RF.Dist.ML ~ maf + height + refdist, 
             data = dat_rf)
summary(m.rf.9)
m.rf.10 <- lm(RF.Dist.ML ~ maf + height + refdist + refdist*maf, 
              data = dat_rf)
summary(m.rf.10)
m.rf.11 <- lmer(RF.Dist.ML ~ (maf + refdist + noref + std.sites + height)*(maf + refdist + noref + std.sites + height) +
                  (1 | simulation), data = dat_rf)
summary(m.rf.11)
r.squaredGLMM(m.rf.11)

m.rf.12 <- lmer(RF.Dist.ML ~ quality + missing + maf + refdist + noref + height + (1| simulation),
                data = dat_rf)
summary(m.rf.12)
r.squaredGLMM(m.rf.12)

my.models<-model.sel(m.rf.1,m.rf.2,m.rf.3,m.rf.4,m.rf.5,m.rf.6,m.rf.7,m.rf.8,m.rf.9,m.rf.10,rank=AIC)
my.models

options(na.action="na.fail")
m_list <- dredge(m.rf.12, rank="AIC")
head(m_list, 10)

m.rf_list <- get.models(m_list, subset = delta < 10)

m.rf_avg <- model.avg(m.rf_list, revised.var = TRUE)
sum.rf <- summary(m.rf_avg)

#sum.rf <- summary(m.rf.3)

confint.rf<-data.frame(confint(m.rf_avg))
colnames(confint.rf) <- c("minCI","maxCI")
confint.rf$var <- rownames(confint.rf)
confint.rf$est <- sum.rf$coefficients[1,]
#confint.rf$cond.est <- sum.rf$coefficients[2,]
confint.rf$sig <- sum.rf$coefmat.full[,5]
confint.rf$sig2 <- case_when(confint.rf$sig < 0.05 & confint.rf$est > 0 ~ "pos",
                             confint.rf$sig < 0.05 & confint.rf$est < 0 ~ "neg",
                             confint.rf$sig >= 0.05 ~ "ns")

vars <- importance(sum.rf)
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))
vars2 <- vars2 %>% mutate(var = factor(var, levels = rev(levels(var))))

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
  ylab("Model Averaged Coefficient")+
  #scale_color_manual(values=c("orchid","black","turquoise"))+
  scale_color_npg()+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=8)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) #+
  ylim(-200,200)

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

m.rf_list <- get.models(m_list, subset = delta < 10)
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

## no interactions

m.rf2 <- lmer(RF.Dist.ML ~ (quality + missing + maf + refdist + height) + (1 | simulation), 
               data = dat_rf)
summary(m.rf2,correlation=TRUE)
vif(m.rf2)
r.squaredGLMM(m.rf2)

options(na.action="na.fail")
m_list2 <- dredge(m.rf2, rank="AIC",trace=1)
head(m_list2, 10)

m.rf2_list <- get.models(m_list2, subset = delta < 10)

m.rf2_avg <- model.avg(m.rf2_list, revised.var = TRUE)
sum.rf2 <- summary(m.rf2_avg)

confint.rf2<-data.frame(confint(m.rf2_avg,method="full"))
colnames(confint.rf2) <- c("minCI","maxCI")
confint.rf2$var <- rownames(confint.rf2)
confint.rf2$est <- sum.rf2$coefficients[1,]
confint.rf2$cond.est <- sum.rf2$coefficients[2,]
confint.rf2$sig <- sum.rf2$coefmat.full[,5]
confint.rf2$sig2 <- case_when(confint.rf2$sig < 0.05 & confint.rf2$est > 0 ~ "pos",
                               confint.rf2$sig < 0.05 & confint.rf2$est < 0 ~ "neg",
                               confint.rf2$sig >= 0.05 ~ "ns")
#confint.rf2[6,c(1:2,4:5)] <- confint.rf2[6,c(1:2,4:5)] - 65
vars <- importance(m.rf2_avg)
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))
vars2 <- vars2 %>% mutate(var = factor(var, levels = rev(levels(var))))

vars.rf2 <- ggplot(confint.rf2[-c(1),], aes(x = var, y = est, color=sig2))
vars.rf2 + geom_blank()+
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
  ylab("Model Averaged Coefficient")+
  scale_color_npg()+
  #scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=8)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_text(size=18))

imp.rf2 <- ggplot(vars2, aes(x = var, y = importance))
imp.rf2 + geom_point(size=4)+
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
  ylab("Summed AIC Weight - RF Distance")


####################################################
# Gamma ############################################
####################################################

dat_gam <- subset(results.mod, select=c(ingroup.gamma, quality, missing, maf, int, noref, std.sites, simulation, height, refdist))
dat_gam$missing <- as.factor(dat_gam$missing)
dat_gam$maf <- as.factor(dat_gam$maf)

dat_gam_ref <- dat_gam[dat_gam$noref == "REF",]

dat_astral <- dat_gam[dat_gam$method == "astral",]
dat_raxml <- dat_gam[dat_gam$method == "raxml",]

mod <- lm(ingroup.gamma ~ .*., data=dat_gam)
formula(mod)
r.squaredGLMM(mod)

m.imb <- lm(ingroup.gamma ~ quality + maf + missing + noref + height + refdist, data = dat_gam)
summary(m.imb)
r.squaredGLMM(m.imb)
car::vif(m.imb)

m.imb1 <- lmer(ingroup.gamma ~ (quality + missing + maf + refdist + height)*(quality + missing + maf + refdist + height) 
              + (1 | simulation), data = dat_gam[dat_gam$noref == "REF",], REML=FALSE)

summary(m.imb1)
vif(m.imb1)
r.squaredGLMM(m.imb1)

m.imb2 <- lmer(ingroup.gamma ~ (quality + missing + maf + refdist + height) + (1 | simulation), 
              data = dat_gam[dat_gam$noref == "REF",])
summary(m.imb2,correlation=TRUE)
vif(m.imb2)
r.squaredGLMM(m.imb2)

options(na.action="na.fail")
m_list <- dredge(m.imb1, rank="AIC",trace=1)
head(m_list, 10)

m.imb_list <- get.models(m_list, subset = cumsum(weight) < 0.95)

m.imb_avg <- model.avg(m.imb_list, revised.var = TRUE)
sum.imb <- summary(m.imb_avg)

confint.gam<-data.frame(confint(m.imb_avg,method="full"))
colnames(confint.gam) <- c("minCI","maxCI")
confint.gam$var <- rownames(confint.gam)
confint.gam$est <- sum.imb$coefficients[1,]
confint.gam$cond.est <- sum.imb$coefficients[2,]
confint.gam$sig <- sum.imb$coefmat.full[,5]
confint.gam$sig2 <- case_when(confint.gam$sig < 0.05 & confint.gam$est > 0 ~ "pos",
                             confint.gam$sig < 0.05 & confint.gam$est < 0 ~ "neg",
                             confint.gam$sig >= 0.05 ~ "ns")
#confint.gam[6,c(1:2,4:5)] <- confint.gam[6,c(1:2,4:5)] - 65
vars <- importance(m.imb_avg)
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))
vars2 <- vars2 %>% mutate(var = factor(var, levels = rev(levels(var))))

vars.gam <- ggplot(confint.gam[confint.gam$sig2 != "ns",], aes(x = var, y = est, color=sig2))
vars.gam + geom_blank()+
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
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none")

imp.gam <- ggplot(vars, aes(x = var, y = sum.imb.importance))
imp.gam + geom_point(size=4)+
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
  ylab("Summed AIC Weight - ingroup.gamma")

## no interactions

m.gam2 <- lmer(ingroup.gamma ~ (quality + missing + maf + refdist + height + noref) + (1 | simulation), 
               data = dat_gam)
summary(m.gam2,correlation=TRUE)
vif(m.gam2)
r.squaredGLMM(m.gam2)

options(na.action="na.fail")
m_list2 <- dredge(m.gam2, rank="AIC",trace=1)
head(m_list2, 10)

m.gam2_list <- get.models(m_list2, subset = delta < 10)

m.gam2_avg <- model.avg(m.gam2_list, revised.var = TRUE)
sum.gam2 <- summary(m.gam2_avg)

confint.gam2<-data.frame(confint(m.gam2_avg,method="full"))
colnames(confint.gam2) <- c("minCI","maxCI")
confint.gam2$var <- rownames(confint.gam2)
confint.gam2$est <- sum.gam2$coefficients[1,]
confint.gam2$cond.est <- sum.gam2$coefficients[2,]
confint.gam2$sig <- sum.gam2$coefmat.full[,5]
confint.gam2$sig2 <- case_when(confint.gam2$sig < 0.05 & confint.gam2$est > 0 ~ "pos",
                              confint.gam2$sig < 0.05 & confint.gam2$est < 0 ~ "neg",
                              confint.gam2$sig >= 0.05 ~ "ns")
#confint.gam2[6,c(1:2,4:5)] <- confint.gam2[6,c(1:2,4:5)] - 65
vars <- importance(m.gam2_avg)
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))
vars2 <- vars2 %>% mutate(var = factor(var, levels = rev(levels(var))))

vars.gam2 <- ggplot(confint.gam2[-c(1),], aes(x = var, y = est, color=sig2))
vars.gam2 + geom_blank()+
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
  ylab("Model Averaged Coefficient")+
  scale_color_npg()+
  #scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=8)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_text(size=18))

imp.gam2 <- ggplot(vars2, aes(x = var, y = importance))
imp.gam2 + geom_point(size=4)+
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
  ylab("Summed AIC Weight - ingroup.gamma")

####################################################
# Colless ##########################################
####################################################

dat_imb <- subset(results.mod, select=c(ingroup.colless, quality, missing, maf, int, noref, std.sites, simulation, height, refdist))
dat_imb$missing <- as.factor(dat_imb$missing)
dat_imb$maf <- as.factor(dat_imb$maf)

dat_imb_ref <- dat_imb[dat_imb$noref == "REF",]

mod <- lm(ingroup.colless ~ .*., data=dat_imb)
formula(mod)
r.squaredGLMM(mod)

m.imb <- lm(ingroup.colless ~ quality + maf + missing + noref + height + refdist, data = dat_imb)
summary(m.imb)
r.squaredGLMM(m.imb)
car::vif(m.imb)

m.imb1 <- lmer(ingroup.colless ~ (quality + missing + maf + refdist + height)*(quality + missing + maf + refdist + height) 
               + (1 | simulation), data = dat_imb_ref, REML=FALSE)

summary(m.imb1)
r.squaredGLMM(m.imb1)

m.imb2 <- lmer(ingroup.colless ~ (quality + missing + maf + refdist + height) + (1 | simulation), 
               data = dat_imb_ref)
summary(m.imb2,correlation=TRUE)
vif(m.imb2)
r.squaredGLMM(m.imb2)

options(na.action="na.fail")
m_list <- dredge(m.imb1, rank="AIC",trace=1)
head(m_list, 10)

m.imb_list <- get.models(m_list, subset = delta < 10)

m.imb_avg <- model.avg(m.imb_list, revised.var = TRUE)
sum.imb <- summary(m.imb_avg)

confint.imb<-data.frame(confint(m.imb_avg,method="full"))
colnames(confint.imb) <- c("minCI","maxCI")
confint.imb$var <- rownames(confint.imb)
confint.imb$est <- sum.imb$coefficients[1,]
confint.imb$cond.est <- sum.imb$coefficients[2,]
confint.imb$sig <- sum.imb$coefmat.full[,5]
confint.imb$sig2 <- case_when(confint.imb$sig < 0.05 & confint.imb$est > 0 ~ "pos",
                              confint.imb$sig < 0.05 & confint.imb$est < 0 ~ "neg",
                              confint.imb$sig >= 0.05 ~ "ns")
#confint.imb[6,c(1:2,4:5)] <- confint.imb[6,c(1:2,4:5)] - 65
vars <- importance(m.imb_avg)
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))
vars2 <- vars2 %>% mutate(var = factor(var, levels = rev(levels(var))))

vars.imb <- ggplot(confint.imb[confint.imb$sig2 != "ns",], aes(x = var, y = est, color=sig2))
vars.imb + geom_blank()+
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
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=5)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),legend.position="none")

imp.imb <- ggplot(vars2, aes(x = var, y = importance))
imp.imb + geom_point(size=4)+
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
  ylab("Summed AIC Weight - ingroup.colless")

## no interactions

m.imb2 <- lmer(ingroup.colless ~ (quality + missing + maf + refdist + height) + (1 | simulation), 
               data = dat_imb_ref)
summary(m.imb2,correlation=TRUE)
vif(m.imb2)
r.squaredGLMM(m.imb2)

options(na.action="na.fail")
m_list2 <- dredge(m.imb2, rank="AIC",trace=1)
head(m_list2, 10)

m.imb2_list <- get.models(m_list2, subset = delta < 15)

m.imb2_avg <- model.avg(m.imb2_list, revised.var = TRUE)
sum.imb2 <- summary(m.imb2_avg)

confint.imb2<-data.frame(confint(m.imb2_avg,method="full"))
colnames(confint.imb2) <- c("minCI","maxCI")
confint.imb2$var <- rownames(confint.imb2)
confint.imb2$est <- sum.imb2$coefficients[1,]
confint.imb2$cond.est <- sum.imb2$coefficients[2,]
confint.imb2$sig <- sum.imb2$coefmat.full[,5]
confint.imb2$sig2 <- case_when(confint.imb2$sig < 0.05 & confint.imb2$est > 0 ~ "pos",
                               confint.imb2$sig < 0.05 & confint.imb2$est < 0 ~ "neg",
                               confint.imb2$sig >= 0.05 ~ "ns")
#confint.imb2[6,c(1:2,4:5)] <- confint.imb2[6,c(1:2,4:5)] - 65
vars <- importance(m.imb2_avg)
vars2 <- data.frame(importance=as.numeric(vars),var=names(vars))
vars2 <- vars2 %>% mutate(var = factor(var, levels = rev(levels(var))))

vars.imb2 <- ggplot(confint.imb2[-c(1),], aes(x = var, y = est, color=sig2))
vars.imb2 + geom_blank()+
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
  ylab("Model Averaged Coefficient")+
  scale_color_npg()+
  #scale_color_manual(values=c("orchid","black","turquoise"))+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=8)+
  coord_flip()+
  theme_minimal()+
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_text(size=18))

imp.imb2 <- ggplot(vars2, aes(x = var, y = importance))
imp.imb2 + geom_point(size=4)+
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
  ylab("Summed AIC Weight - ingroup.colless")
