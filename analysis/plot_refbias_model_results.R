########################
### Summary plots for refbias analyses
### and plots exploring the data
### 
### J. Rick, Sept 2021
### Updated Feb 2022
########################

library(ggplot2)
library(viridis)
library(ggpubr)
library(ggdist)
library(kableExtra)
library(here)

here::i_am("analysis/plot_refbias_model_results.R")
source(here("analysis","theme_custom.R"))

output <- "092321-output"
output <- "092321-subsamp-output"
results.raxml <- read_csv(here("output","new",paste0("092321-output","-raxml.csv")))[,-1]
results.astral <- read_csv(here(here("output","new",paste0("092321-output","-astral.csv"))))[,-1]
results.raxml.lates <- read_csv(here("output","new",paste("010822-lates-emp-output","-raxml.csv",sep="")))[,-1]
results.astral.lates <- read_csv(here("output","new",paste("031623-lates-emp-output","-astral.csv",sep="")))[,-1]
results.raxml.cichlids <- read_csv(here("output","new",paste("010822-cichlids-emp-output","-raxml.csv",sep="")))[,-1]
results.astral.cichlids <- read_csv(here("output","new",paste("031623-cichlids-emp-output","-astral.csv",sep="")))[,-1]
figdir <- "~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/figures/"

#cols.int <- c("#005F73","gray80","#5FB89D")
cols.int <- c("#e2d200","#2A9D8F","#e76f51") # colors for tree heights
cols.sig <- c("black","gray80","black")

## preparing data objects
results.emp <- results.raxml.lates %>%
  add_row(results.raxml.cichlids)
results.mod <- results.raxml[results.raxml$simulation > 15 & results.raxml$simulation < 26,]
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.numeric(as.character(results.mod$missing))
results.mod$maf <- as.numeric(as.character(results.mod$maf))
results.mod$int <- as.factor(results.mod$int)

results.mod.astr <- results.astral[results.astral$simulation > 15 & results.astral$simulation < 26,]
results.mod.astr$simulation <- as.factor(results.mod.astr$simulation)
results.mod.astr$height <- as.factor(results.mod.astr$height)
results.mod.astr$quality <- as.factor(results.mod.astr$quality)
results.mod.astr$missing <- as.factor(as.character(results.mod.astr$missing))
results.mod.astr$maf <- as.numeric(as.character(results.mod.astr$maf))
results.mod.astr$int <- as.factor(results.mod.astr$int)

#----------------------------------------------#
### plot for maf & miss : sites correlation ####
#----------------------------------------------#
maf.sites <- results.mod %>% 
  ggplot(aes(x=maf,y=sites,fill=height)) + 
  geom_boxplot() + theme_custom() + 
  scale_fill_manual(values=cols.int,label=c("Low ILS","Med ILS","High ILS")) + 
  xlab("Minor Alelle Count") + 
  ylab("Number of SNPs")
miss.sites <- results.mod %>% 
  ggplot(aes(x=missing,y=sites,fill=height)) + 
  geom_boxplot() + theme_custom() + 
  scale_fill_manual(values=cols.int,label=c("Low ILS","Med ILS","High ILS")) + 
  xlab("Missing Data") + 
  ylab("Number of SNPs")
ggarrange(maf.sites,miss.sites,
          labels="AUTO",common.legend=TRUE,
          font.label=list(family="Open Sans", size=24))

## for empirical datasets
maf.sites.emp <- results.emp %>% 
  #filter(snps < 1000000) %>% 
  ggplot(aes(x=as.factor(maf),y=snps,fill=int)) + 
  geom_boxplot() + theme_custom() + 
  scale_fill_manual(values=cols.int) + 
  xlab("Minor Alelle Count") + 
  ylab("Number of SNPs") +
  facet_wrap(~height)
miss.sites.emp <- results.emp %>% 
  #filter(snps < 1000000) %>%
  ggplot(aes(x=as.factor(missing),y=snps,fill=int)) + 
  geom_boxplot() + theme_custom() + 
  scale_fill_manual(values=cols.int,label=c("Low ILS","Med ILS","High ILS")) + 
  xlab("Missing Data") + 
  ylab("Number of SNPs") +
  facet_wrap(~height)
ggarrange(maf.sites.emp,miss.sites.emp,
          labels="AUTO",common.legend=TRUE,
          font.label=list(family="Open Sans", size=24))

#----------------------------------------------------------------------------#
## heatmaps to show relationship of maf and miss to variables of interest ####
#----------------------------------------------------------------------------#

## for sims
results_mean <- results.raxml %>%
  select(int,height,maf,missing,avg_dxy,RF.Dist.ML,ingroup.gamma,ingroup.colless,Q.Dist.ML,CI.Dist.ML,ingroup.tree.height) %>%
  #pivot_wider(id_cols="maf",names_from="missing",values_from="RF.Dist.ML",values_fn=mean) %>%
  group_by(maf,missing,int,height) %>%
  mutate(mean_dxy = mean(avg_dxy),
         mean_rf = mean(RF.Dist.ML),
         mean_gam = mean(ingroup.gamma),
         mean_colless = mean(ingroup.colless),
         mean_qdist = mean(Q.Dist.ML),
         mean_cidist = mean(CI.Dist.ML),
         mean_height = mean(ingroup.tree.height)) %>%
  select(maf,mean_dxy:mean_height) %>%
  pivot_longer(cols=mean_rf:mean_height,names_to="param",values_to="mean") 

## for emp
results_mean.lates <- results.raxml.lates %>%
  select(int,height,maf,missing,ingroup.gamma,ingroup.colless,ingroup.tree.height) %>%
  #pivot_wider(id_cols="maf",names_from="missing",values_from="RF.Dist.ML",values_fn=mean) %>%
  group_by(maf,missing,int,height) %>%
  mutate(mean_gam = mean(ingroup.gamma),
         mean_colless = mean(ingroup.colless),
         mean_height = mean(ingroup.tree.height)) %>%
  select(maf,mean_gam:mean_height) %>%
  pivot_longer(cols=mean_gam:mean_height,names_to="param",values_to="mean") %>%
  filter(height=="Lates")
results_mean.cichlids <- results.raxml.cichlids %>%
  select(int,height,maf,missing,ingroup.gamma,ingroup.colless,ingroup.tree.height) %>%
  #pivot_wider(id_cols="maf",names_from="missing",values_from="RF.Dist.ML",values_fn=mean) %>%
  group_by(maf,missing,int,height) %>%
  mutate(mean_gam = mean(ingroup.gamma),
         mean_colless = mean(ingroup.colless),
         mean_height = mean(ingroup.tree.height)) %>%
  select(maf,mean_gam:mean_height) %>%
  pivot_longer(cols=mean_gam:mean_height,names_to="param",values_to="mean") %>%
  filter(height=="Cichlids")

# rf_interp <- interp::interp(x=results_mean$maf[results_mean$param == "mean_rf"], y=results_mean$missing[results_mean$param == "mean_rf"], z=results_mean$mean[results_mean$param == "mean_rf"],
#               output="grid",duplicate = "mean")
categ <- "fullsim"
categ <- "subsamp"
categ <- "emp"
for (h in unique(results.raxml$height)) {
  results_mean_h <- results_mean.lates%>%
    filter(height == h)
  
  contour_rf <- results_mean_h %>%
    filter(param == "mean_rf") %>%
    ggplot() +
      geom_contour_filled(aes(x=maf,y=missing,z=mean),bins=15) +
      theme_custom() +
    scale_fill_viridis(discrete=TRUE,option="C", name="RF Distance") +
    theme(legend.title = element_text(size=rel(1.8)),
          axis.title = element_text(size=rel(1.8)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(fill="white")) +
    ylab("Missing Data") +
    xlab("Minor Allele Count") +
    #facet_grid(vars(int),vars(height))
    facet_grid(rows=vars(int))

  contour_gam <- results_mean_h %>%
    filter(param == "mean_gam") %>%
    ggplot() +
    geom_contour_filled(aes(x=maf,y=missing,z=mean),bins=15) +
    theme_custom() +
    scale_fill_viridis(discrete=TRUE,option="C", name="Ingroup\nGamma") +
    theme(legend.title = element_text(size=rel(1.8)),
          axis.title = element_text(size=rel(1.8)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(fill="white")) +
    ylab("Missing Data") +
    xlab("Minor Allele Count") +
    #facet_grid(vars(int),vars(height))
    facet_grid(rows=vars(int))

  contour_imb <- results_mean_h %>%
    filter(param == "mean_colless") %>%
    ggplot() +
    geom_contour_filled(aes(x=maf,y=missing,z=mean),bins=15) +
    theme_custom() +
    scale_fill_viridis(discrete=TRUE,option="C", name="Ingroup\nColless\nImbalance") +
    theme(legend.title = element_text(size=rel(1.8)),
          axis.title = element_text(size=rel(1.8)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(fill="white")) +
    ylab("Missing Data") +
    xlab("Minor Allele Count") +
    #facet_grid(vars(int),vars(height))
    facet_grid(rows=vars(int))

  contour_qdist <- results_mean_h %>%
    filter(param == "mean_qdist") %>%
    ggplot() +
    geom_contour_filled(aes(x=maf,y=missing,z=mean),bins=15) +
    theme_custom() +
    scale_fill_viridis(discrete=TRUE,option="C", name="Quartet\nDistance") +
    theme(legend.title = element_text(size=rel(1.8)),
          axis.title = element_text(size=rel(1.8)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(fill="white")) +
    ylab("Missing Data") +
    xlab("Minor Allele Count") +
    #facet_grid(vars(int),vars(height))
    facet_grid(rows=vars(int))

  contour_ci <- results_mean_h %>%
    filter(param == "mean_cidist") %>%
    ggplot() +
    geom_contour_filled(aes(x=maf,y=missing,z=mean),bins=15) +
    theme_custom() +
    scale_fill_viridis(discrete=TRUE,option="C", name="CI Distance") +
    theme(legend.title = element_text(size=rel(1.8)),
          axis.title = element_text(size=rel(1.8)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(fill="white")) +
    ylab("Missing Data") +
    xlab("Minor Allele Count") +
    #facet_grid(vars(int),vars(height))
    facet_grid(rows=vars(int))

  contour_height <- results_mean_h %>%
    filter(param == "mean_height") %>%
    ggplot() +
    geom_contour_filled(aes(x=maf,y=missing,z=mean),bins=15) +
    theme_custom() +
    scale_fill_viridis(discrete=TRUE,option="C", name="Ingroup\nTree Height") +
    theme(legend.title = element_text(size=rel(1.8)),
          axis.title = element_text(size=rel(1.8)),
          strip.text = element_text(size=rel(1.3)),
          strip.background = element_rect(fill="white")) +
    ylab("Missing Data") +
    xlab("Minor Allele Count") +
    #facet_grid(vars(int),vars(height))
    facet_grid(rows=vars(int))

  if (categ == "emp") {
    # saved at 2000x800px for each clade
    heatmap <- ggarrange(contour_height, contour_imb, contour_gam, ncol=3,
              labels="AUTO",font.label=list(size=24))
    heatmap %>%
      ggexport(filename = paste0(figdir,"empirical/refbias_emp_heatmap_",h,"_byINT.png"), width=2000,height=800)
               
  } else if (categ == "subsamp") {
    # exported at 2000 x 1400 for each height
    heatmap <- ggarrange(contour_rf, contour_qdist, contour_ci,
              contour_gam, contour_imb, contour_height, nrow=2, ncol=3,
              labels="AUTO",font.label=list(size=24))
    heatmap %>%
      ggexport(filename = paste0(figdir,"subsamp/refbias_subsamp_heatmap_",h,"_byINT.png"), width=2000,height=1400)
  } else {
    # exported at 2000 x 1400 for each height
    heatmap <- ggarrange(contour_rf, contour_qdist, contour_ci,
                         contour_gam, contour_imb, contour_height, nrow=2, ncol=3,
                         labels="AUTO",font.label=list(size=24))
    #heatmap %>%
    #  ggexport(filename = paste0(figdir,"full_sims/refbias_heatmap_",h,"_byINT.png"), width=2000,height=1400)
  }
  print(heatmap)
}

#----------------------------------#
### plots of model coefficients ####
#----------------------------------#
ggarrange(vars.plots.rf,vars.plot.gam,vars.plots.imb,labels="AUTO",font.label=list(size=24),ncol=1)
ggarrange(vars.rf.all.bars,vars.gam.all.bars,vars.imb.all.bars,labels="AUTO",font.label=list(size=24),nrow=1)

# raxml, full sims
all.confint <- confint.gam.long %>% 
  as_tibble() %>% 
  mutate(resp = "std.ingroup.gamma") %>% 
  add_row(confint.rf.long %>% 
            as_tibble() %>% 
            mutate(resp = "RF.Dist.ML")) %>%
  add_row(confint.imb.long %>% 
            as_tibble() %>% 
            mutate(resp = "std.ingroup.colless")) %>%
  mutate(height = "long") %>%
  add_row(confint.gam.med %>% 
            as_tibble() %>% 
            mutate(resp = "std.ingroup.gamma") %>% 
            add_row(confint.rf.med %>% 
                      as_tibble() %>% 
                      mutate(resp = "RF.Dist.ML")) %>%
            add_row(confint.imb.med %>% 
                      as_tibble() %>% 
                      mutate(resp = "std.ingroup.colless")) %>%
            mutate(height = "med")) %>%
  add_row(confint.gam.short %>% 
            as_tibble() %>% 
            mutate(resp = "std.ingroup.gamma") %>% 
            add_row(confint.rf.short %>% 
                      as_tibble() %>% 
                      mutate(resp = "RF.Dist.ML")) %>%
            add_row(confint.imb.short %>% 
                      as_tibble() %>% 
                      mutate(resp = "std.ingroup.colless")) %>%
            mutate(height = "short")) %>%
  add_row(confint.gam.all %>% 
             as_tibble() %>% 
             mutate(resp = "std.ingroup.gamma") %>% 
             add_row(confint.rf.all %>% 
                       as_tibble() %>% 
                       mutate(resp = "RF.Dist.ML")) %>%
             add_row(confint.imb.all %>% 
                       as_tibble() %>% 
                       mutate(resp = "std.ingroup.colless")) %>%
             mutate(height = "all")) %>%
  mutate(categ = case_when(height == "all" ~ "all",
                           TRUE ~ "by.height"),
         nudge = case_when(height == "short" ~ -0.075,
                           height == "long" ~ 0.225,
                           height == "med" ~ 0.075,
                           height == "all" ~ -0.225)) %>%
  mutate(col_sig = case_when(sig == "ns" ~ "ns",
                             TRUE ~ height))
cols.int.sig <- c("#e2d200","#2A9D8F","gray","#e76f51")
cols.int.plusAll <- c("black","#e2d200","#2A9D8F","#e76f51") # colors for tree heights

## FIG 3, PANEL A ######
confint.plot.bygroup <- all.confint  %>%
  filter(categ == "by.height") %>%
  mutate(var = factor(var,levels=c("maf","missing","avg_dxy","avg_dxy:maf","avg_dxy:missing"))) %>%
  ggplot(aes(x=var,y=est,col=height)) +
  xlab("")+
  ylab("Mixed Model Coefficient")+
  scale_color_manual(values=c(cols.int),guide="none")+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI, shape=col_sig, alpha=col_sig),
                  fill="white",fatten=8,lwd=2, 
                  position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height"])) +
  # geom_point(data=. %>% filter(col_sig == "ns"), aes(x=var,y=est,group=height),
  #            color="black",shape=16, fill="black",size=4,
  #            position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height"])) +
  scale_shape_manual(values=c(16,16,21,16,16)) +
  scale_alpha_manual(values=c(1,1,0.4,1,1)) +
  #coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=c("MAC Threshold",
                            "Missing data\nthreshold",
                            "Reference\ngenome",
                            "Reference\ngenome x MAC",
                            "Reference genome\nx Missing data")) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(angle=45,hjust=1),
        #strip.text = element_text(size=15),
        plot.margin=margin(t=45),
        strip.background = element_rect(),
        #strip.text.y.right = element_text(size=15,angle=90),
        strip.text = element_blank(),
        panel.spacing = unit(1, "lines")) +
  facet_grid(rows=vars(resp),scales="free_y",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Distance to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))
  #facet_wrap(vars(resp,height),scales="free",as.table=TRUE)

confint.plot.all <- all.confint %>%
  filter(height == "all") %>%
  ggplot(aes(x=var,y=est,col=sig)) +
  #scale_color_manual(values=cols.sig[c(2,3)]) +
  scale_color_manual(values=cols.sig) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge), alpha=1, shape=4) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

## saved at 1000 x 1000px
p1 <- ggarrange(confint.plot.all,confint.plot.bygroup,ncol=1,
                hjust=-0.2,
          heights=c(1,1.5),labels=c("a) Models for all trees","b) Models by tree height"),
          font.label=list(size=24,family="Open Sans"),
          common.legend=TRUE,legend.grob=get_legend(confint.plot.bygroup),legend="none")

# astral, full sims
all.confint.astr <- confint.imb.long.astr %>% 
  as_tibble() %>% 
  mutate(resp = "std.ingroup.colless") %>% 
  add_row(confint.rf.long.astr %>% 
            as_tibble() %>% 
            mutate(resp = "RF.Dist.ML")) %>%
  mutate(height = "long") %>%
  add_row(confint.imb.med.astr %>% 
            as_tibble() %>% 
            mutate(resp = "std.ingroup.colless") %>% 
            add_row(confint.rf.med.astr %>% 
                      as_tibble() %>% 
                      mutate(resp = "RF.Dist.ML")) %>%
            mutate(height = "med")) %>%
  add_row(confint.imb.short.astr %>% 
            as_tibble() %>% 
            mutate(resp = "std.ingroup.colless") %>% 
            add_row(confint.rf.short.astr %>% 
                      as_tibble() %>% 
                      mutate(resp = "RF.Dist.ML")) %>%
            mutate(height = "short")) %>%
  add_row(confint.imb.all.astr %>% 
            as_tibble() %>% 
            mutate(resp = "std.ingroup.colless") %>% 
            add_row(confint.rf.all.astr %>% 
                      as_tibble() %>% 
                      mutate(resp = "RF.Dist.ML")) %>%
            mutate(height = "all")) %>%
  mutate(categ = case_when(height == "all" ~ "all",
                           TRUE ~ "by.height"),
         nudge = case_when(height == "short" ~ -0.2,
                           height == "long" ~ 0.2,
                           height == "med" ~ 0,
                           height == "all" ~ 0)) %>%
  mutate(col_sig = case_when(sig == "ns" ~ "ns",
                             TRUE ~ height))
cols.int.sig <- c("#e2d200","#2A9D8F","gray","#e76f51")
confint.plot.bygroup.astr <- all.confint.astr  %>%
  filter(categ == "by.height") %>%
  ggplot(aes(x=var,y=est,col=height)) +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols.int,guide="none")+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI, shape=col_sig, alpha=col_sig),
                  fill="white",fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint.astr$nudge[all.confint.astr$categ == "by.height"])) +
  # geom_point(data=. %>% filter(col_sig == "ns"), aes(x=var,y=est,group=height),
  #            color="black",shape=16, fill="black",size=4,
  #            position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height"])) +
  scale_shape_manual(values=c(16,16,21,16,16)) +
  scale_alpha_manual(values=c(1,1,0.4,1,1)) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance")
             ))
#facet_wrap(vars(resp,height),scales="free",as.table=TRUE)

confint.plot.all.astr <- all.confint.astr %>%
  filter(height == "all") %>%
  ggplot(aes(x=var,y=est,col=sig)) +
  #scale_color_manual(values=cols.sig[c(2,3)]) +
  scale_color_manual(values=cols.sig) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint.astr$nudge), alpha=1, shape=4) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

## saved at 1000 x 1000px
p1 <- ggarrange(confint.plot.all.astr,confint.plot.bygroup.astr,ncol=1,
                hjust=-0.2,
                heights=c(1,1.5),labels=c("a) Models for all trees","b) Models by tree height"),
                font.label=list(size=24,family="Open Sans"),
                common.legend=TRUE,legend.grob=get_legend(confint.plot.bygroup.astr),legend="none")

#----------------------------------#
### table with param estimates #####
#----------------------------------#
confint.table <- all.confint %>%
  relocate(resp,height,var,est,minCI,maxCI) %>%
  select(resp:sig) %>%
  mutate(sig2 = case_when(sig == "pos" ~ "*",
                         sig == "neg" ~ "*",
                         sig == "ns" ~ "")) %>%
  pivot_wider(id_cols=c(resp,var),names_from=height,values_from=c(est:maxCI,sig2),names_glue="{height}_{.value}",
              names_sort=TRUE) %>%
  relocate(resp,var,starts_with("all"),starts_with("long"),starts_with("med"),starts_with("short")) %>%
  select(!resp) %>%
  kbl(booktabs = T,digits=3,
      format="latex",
      col.names=c("Variable","Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig")) %>%
  add_header_above(c(" ","All Trees" = 4, "Low ILS" = 4, "Medium ILS" = 4,  "High ILS" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Standardized ingroup gamma", 1, 5) %>%
  pack_rows("Robinson-Foulds distance to true tree", 6, 10) %>%
  pack_rows("Standardized ingroup Colless imbalance", 11, 15)
#writeLines(confint.table, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_fullsims.tex')
#writeLines(confint.table, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_subsamp.tex')

#---------------------------------------------------#
### table with param estimates for astral trees #####
#---------------------------------------------------#
confint.table.astr <- all.confint.astr %>%
  relocate(resp,height,var,est,minCI,maxCI) %>%
  select(resp:sig) %>%
  mutate(sig2 = case_when(sig == "pos" ~ "*",
                          sig == "neg" ~ "*",
                          sig == "ns" ~ "")) %>%
  pivot_wider(id_cols=c(resp,var),names_from=height,values_from=c(est:maxCI,sig2),names_glue="{height}_{.value}",
              names_sort=TRUE) %>%
  relocate(resp,var,starts_with("all"),starts_with("long"),starts_with("med"),starts_with("short")) %>%
  select(!resp) %>%
  kbl(booktabs = T,digits=3,
      format="latex",
      col.names=c("Variable","Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig")) %>%
  add_header_above(c(" ","All Trees" = 4, "Low ILS" = 4, "Medium ILS" = 4,  "High ILS" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Robinson-Foulds distance to true tree", 1, 5) %>%
  pack_rows("Standardized ingroup Colless imbalance", 6, 10)
#writeLines(confint.table.astr, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_fullsims_astral.tex')


#-----------------------------------------------#
### table with param estimates for emp runs #####
#-----------------------------------------------#

all.confint.emp <- confint.gam.cichlids %>% 
  as_tibble() %>% 
  mutate(resp = "ingroup.gamma") %>% 
  add_row(confint.height.cichlids %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.tree.height")) %>%
  add_row(confint.imb.cichlids %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.colless")) %>%
  mutate(height = "Tropheines") %>%
  add_row(confint.gam.lates %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.gamma") %>% 
            add_row(confint.height.lates %>% 
                      as_tibble() %>% 
                      mutate(resp = "ingroup.tree.height")) %>%
            add_row(confint.imb.lates %>% 
                      as_tibble() %>% 
                      mutate(resp = "ingroup.colless")) %>%
            mutate(height = "Lates")) %>%
  mutate(col_sig = case_when(sig == "ns" ~ "ns",
                             TRUE ~ height),
         newvar = case_when(var == "intINT" ~ "ref",
                            var == "intINT:maf" ~ "ref:mac",
                            var == "intINT:missing" ~ "ref:missing",
                            var == "maf" ~ "mac",
                            TRUE ~ var))

confint.table.emp <- all.confint.emp %>%
  relocate(resp,height,newvar,est,minCI,maxCI) %>%
  select(resp:sig) %>%
  mutate(sig2 = case_when(sig == "pos" ~ "*",
                          sig == "neg" ~ "*",
                          sig == "ns" ~ "")) %>%
  pivot_wider(id_cols=c(resp,newvar),names_from=height,values_from=c(est:maxCI,sig2),names_glue="{height}_{.value}",
              names_sort=TRUE) %>%
  relocate(resp,newvar,starts_with("Tropheines"),starts_with("Lates")) %>%
  select(!resp) %>%
  kbl(booktabs = T,digits=3,
      format="latex",
      col.names=c("Variable","Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig")) %>%
  add_header_above(c(" ","Tropheine Topologies" = 4, "Lates Topologies" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Ingroup gamma", 1, 5) %>%
  pack_rows("Ingroup tree height", 6, 10) %>%
  pack_rows("Ingroup Colless imbalance", 11, 15)
#writeLines(confint.table.emp, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_emp.tex')

#-------------------------------------------------------#
### table with param estimates for emp astral trees #####
#-------------------------------------------------------#

all.confint.emp.astr <- confint.imb.cichlids.astr %>% 
  as_tibble() %>% 
  mutate(resp = "ingroup.colless") %>%
  mutate(height = "Tropheines") %>%
  add_row(confint.imb.lates.astr %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.colless") %>%
            mutate(height = "Lates")) %>%
  mutate(col_sig = case_when(sig == "ns" ~ "ns",
                             TRUE ~ height),
         newvar = case_when(var == "intINT" ~ "ref",
                            var == "intINT:maf" ~ "ref:mac",
                            var == "intINT:missing" ~ "ref:missing",
                            var == "maf" ~ "mac",
                            TRUE ~ var))
confint.table.emp.astr <- all.confint.emp.astr %>%
  relocate(resp,height,var,est,minCI,maxCI) %>%
  select(resp:sig) %>%
  mutate(sig2 = case_when(sig == "pos" ~ "*",
                          sig == "neg" ~ "*",
                          sig == "ns" ~ "")) %>%
  pivot_wider(id_cols=c(resp,var),names_from=height,values_from=c(est:maxCI,sig2),names_glue="{height}_{.value}",
              names_sort=TRUE) %>%
  relocate(resp,var,starts_with("Tropheines"),starts_with("Lates")) %>%
  select(!resp) %>%
  kbl(booktabs = T,digits=3,
      format="latex",
      col.names=c("Variable","Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig")) %>%
  add_header_above(c(" ","Tropheine Topologies" = 4, "Lates Topologies" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Ingroup Colless imbalance", 1, 5) 
#writeLines(confint.table.emp.astr, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_emp_astral.tex')


#----------------------------------------------------------#
### plot of distances to true tree by maf & int ############
### because this seems to be most interesting from models ##
#----------------------------------------------------------#
rf.plot <- results.mod %>%
  ggplot(aes(x=maf,y=RF.Dist.ML,col=height,shape=int)) +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15)) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Robinson-Foulds\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
        legend.text = element_text(size=rel(1.5)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  coord_cartesian(ylim=c(0,0.35),clip="off") +
  geom_hline(yintercept=0)
q.plot <- results.mod %>%
  ggplot(aes(x=maf,y=Q.Dist.ML,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Quartet\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim=c(0,0.045),clip="off") +
  geom_hline(yintercept=0)
ci.plot <- results.mod %>%
  ggplot(aes(x=maf,y=CI.Dist.ML,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15)) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Clustering Information\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim=c(0,0.28),clip="off") +
  geom_hline(yintercept=0)

results.mod.astr$maf <- as.numeric(as.character(results.mod.astr$maf))
rf.plot.astr <- results.mod.astr %>%
  ggplot(aes(x=maf,y=RF.Dist.ML,col=height,shape=int)) +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15)) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Robinson-Foulds\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
        legend.text = element_text(size=rel(1.5)))  +
  coord_cartesian(ylim=c(0,0.35),clip="off") +
  geom_hline(yintercept=0)
q.plot.astr <- results.mod.astr %>%
  ggplot(aes(x=maf,y=Q.Dist.ML,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Quartet\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim=c(0,0.045),clip="off") +
  geom_hline(yintercept=0)
ci.plot.astr <- results.mod.astr %>%
  ggplot(aes(x=maf,y=CI.Dist.ML,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15)) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Clustering Information\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  coord_cartesian(ylim=c(0,0.28),clip="off") +
  geom_hline(yintercept=0)

## saved at 800px x 1000px
p2 <- ggarrange(rf.plot,
                q.plot,
                ci.plot,
                ncol=2,nrow=3,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)
p2.both <- ggarrange(rf.plot,rf.plot.astr,
                q.plot,q.plot.astr,
                ci.plot,ci.plot.astr,
                
                ncol=2,nrow=3,labels=c("A","D","B","E","C","F"),
                font.label=list(size=24,font.family="Open Sans"),
                common.legend=TRUE,legend="right",
                label.x=0.3,label.y=0.95)

gam.plot <- results.mod %>%
  ggplot(aes(x=maf,y=std.ingroup.gamma,col=height,shape=int)) +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized\ningroup gamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  geom_hline(yintercept=0)
ic.plot <- results.mod %>%
  ggplot(aes(x=maf,y=std.ingroup.colless,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.03,0.09),clip="off")
ic.plot.astr <- results.mod.astr %>%
  ggplot(aes(x=maf,y=std.ingroup.colless,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.03,0.09),clip="off")
is.plot <- results.mod %>%
  ggplot(aes(x=maf,y=std.ingroup.sackin,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.025,0.06),clip="off")
is.plot.astr <- results.mod.astr %>%
  ggplot(aes(x=maf,y=std.ingroup.sackin,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0) +
  coord_cartesian(ylim=c(-0.025,0.06), clip="off")
th.plot.maf <- results.mod %>%
  ggplot(aes(x=maf,y=ingroup.tree.height,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup Tree Height") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggarrange(rf.plot,ic.plot,
          q.plot,is.plot,
          ci.plot,gam.plot,
          ncol=2,nrow=3,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.1,label.y=0.95)
ggarrange(rf.plot,ic.plot,
          gam.plot,
          ncol=1,nrow=3,labels=NULL,
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.1,label.y=0.95) 

## FIG 3, PANEL B ######
rf.ic.gam.plot <- results.mod %>%
  select(height,missing,maf,int,RF.Dist.ML,std.ingroup.colless,std.ingroup.gamma) %>%
  pivot_longer(cols=RF.Dist.ML:std.ingroup.gamma,names_to="statistic",values_to="estimate") %>%
  ggplot(aes(x=maf,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("") +
  ylab("\n") +
  theme(axis.title.y = element_text(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines")) +
  geom_hline(yintercept=0)

## FIGURE 3 ######
# edited in AI to add annotations
fig3 <- confint.plot.bygroup + rf.ic.gam.plot +
  plot_layout(widths = c(1.5,1), guides = "collect")
ggsave(file="Fig3.svg", plot=fig3, height=9.5, width=12)

## FIGURE 5 ######
# edited in AI to add annotations
# same as Fig 3 except with subsampled data instead of full data set
fig5 <- confint.plot.bygroup + rf.ic.gam.plot +
  plot_layout(widths = c(1.5,1), guides = "collect")
ggsave(file="Fig5.svg", plot=fig5, height=9.5, width=12)

## saved at 800px x 1000px
p3 <- ggarrange(gam.plot,ic.plot,is.plot,ncol=1,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)

#library(patchwork)
p4 <- ggarrange(rf.plot, rf.plot.astr, ic.plot, gam.plot, ncol=1, 
                labels=c("c)","d)","e)","f)"), font.label=list(size=24,font.family="Open Sans"),
                common.legend=TRUE, legend="right",
                label.x=0.2, label.y=0.95) + theme(plot.margin=margin(l=20))
## saved at 1800 x 1200 as refbias_coeffs_maf_multipanel_120421.png and refbias_coefficients_all_multipanel_subsamp_120421.png
p5 <- p1 + p4 + 
  plot_layout(widths = c(1.2, 1))
p5
#ggexport(p5,filename=paste0(figdir,"full_sims/refbias_coeffs_maf_multipanel_101322.png"),width=1300,height=1000)

## saved at 1600x1200 as refbias_raxml_astral_comparison_byMAF.png
ggarrange(rf.plot,rf.plot.astr,
          q.plot,q.plot.astr,
          ci.plot,ci.plot.astr,
          ic.plot,ic.plot.astr,
          is.plot,is.plot.astr,
          ncol=2,nrow=5,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.28,label.y=0.95)

#-------------------------------------------------------#
### plot of distances to true tree by missing & int #####
### for the supplement ##################################
#-------------------------------------------------------#
rf.plot.miss <- results.mod %>%#
  ggplot(aes(x=missing,y=RF.Dist.ML,col=height,shape=int)) +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15)) +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Robinson-Foulds\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 25, b = 0, l = 0)),
        legend.text = element_text(size=rel(1.5)))
q.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=Q.Dist.ML,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Quartet\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ci.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=CI.Dist.ML,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15)) +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Clustering Information\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

gam.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=std.ingroup.gamma,col=height,shape=int)) +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Standardized\ningroup gamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 45, b = 0, l = 0))) +
  geom_hline(yintercept=0)
ic.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=std.ingroup.colless,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Standardized ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0)
is.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=std.ingroup.sackin,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Standardized ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0)
th.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=ingroup.tree.height,col=height,shape=int))  +
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(labels=c("Low ILS","Med ILS","High ILS"),
                     values=cols.int) +
  scale_shape_manual(values=c(17,15))  +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup Tree Height") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

## saved at 1200 x 1000px
ggarrange(rf.plot.miss,ic.plot.miss,
          q.plot.miss,is.plot.miss,
          ci.plot.miss,gam.plot.miss,
          ncol=2,nrow=3,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.1,label.y=0.95)

#-----------------------------------------#
### confint plots for empirical trees #####
#-----------------------------------------#
all.confint.lates <- confint.gam.lates %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.gamma") %>% 
            add_row(confint.height.lates %>% 
                      as_tibble() %>% 
                      mutate(resp = "ingroup.height")) %>%
            add_row(confint.imb.lates %>% 
                      as_tibble() %>% 
                      mutate(resp = "ingroup.colless")) %>%
            mutate(height = "lates",
                   resp = factor(resp, levels = c("ingroup.height","ingroup.colless","ingroup.gamma")),
                   sig_shape = case_when(sig == "ns" ~ "ns",
                                         TRUE ~ "sig")) %>%
  mutate(var = factor(var,levels=c("maf","missing","intINT","intINT:maf","intINT:missing")))
all.confint.cich <- confint.gam.cichlids %>% 
  as_tibble() %>% 
  mutate(resp = "ingroup.gamma") %>% 
  add_row(confint.height.cichlids %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.height")) %>%
  add_row(confint.imb.cichlids %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.colless")) %>%
  mutate(height = "cichlid",
         resp = factor(resp, levels = c("ingroup.height","ingroup.colless","ingroup.gamma")),
         sig_shape = case_when(sig == "ns" ~ "ns",
                               TRUE ~ "sig")) %>%
  mutate(var = factor(var,levels=c("maf","missing","intINT","intINT:maf","intINT:missing")))

## FIG 8 PANEL A ########
confint.plot.emp.all <- all.confint.lates %>%
  
  ggplot(aes(x=var,y=est)) +
  scale_color_manual(values=c(cols.int[2],cols.sig[2],cols.int[2])) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI,alpha=sig,shape=sig_shape),fatten=8,lwd=1,
                  col=cols.int[2],position=position_nudge(x=0.1)) +
  geom_pointrange(data=all.confint.cich, aes(ymin = minCI, ymax = maxCI, alpha=sig,shape=sig_shape),
                 fatten=8,lwd=1, col=cols.int[3],position=position_nudge(x=-0.1)) +
  scale_alpha_manual(values=c(1,0.3,1)) +
  scale_shape_manual(values=c(21,16)) +
  scale_x_discrete(labels=c("MAC Threshold",
                            "Missing data\nthreshold",
                            "Reference\ngenome",
                            "Reference\ngenome x MAC",
                            "Reference genome\nx Missing data")) +
  #coord_flip() +
  theme_custom() +
  ylab("Mixed Model Coefficient") +
  xlab("") +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        #axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1),
        panel.spacing = unit(1, "lines")) +
  facet_grid(rows=vars(resp),scales="free_y",
             labeller = labeller(
               resp = c("ingroup.height" = "Ingroup Height",
                        "ingroup.colless" = "Ingroup Imbalance",
                        "ingroup.gamma" = "Ingroup Gamma")
             ))
#confint.plot.all.lates <- confint.plot.all
#confint.plot.all.cich <- confint.plot.all

# p.confint <- ggarrange(confint.plot.all.cich,confint.plot.all.lates,ncol=1,
#           labels="AUTO",
#           font.label=list(size=24,font.family="Open Sans"),
#           label.x=0,label.y=1,
#           common.legend=TRUE,legend="right")

#------------------------------------------------#
### confint plots for astral empirical trees #####
#------------------------------------------------#
all.confint.astr.lates <- confint.imb.lates.astr %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.colless") %>%
  mutate(height = "lates",
         resp = factor(resp))
all.confint.astr.cich <- confint.imb.cichlids.astr %>% 
            as_tibble() %>% 
            mutate(resp = "ingroup.colless") %>%
  mutate(height = "cichlid",
         resp = factor(resp))
confint.plot.astr.emp.all <- all.confint.astr.lates %>%
  ggplot(aes(x=var,y=est)) +
  scale_color_manual(values=c(cols.int[2],cols.sig[2],cols.int[2])) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI,alpha=sig),fatten=5,lwd=1,
                  shape=4, col=cols.int[2],position=position_nudge(x=0.1)) +
  geom_pointrange(data=all.confint.astr.cich, aes(ymin = minCI, ymax = maxCI, alpha=sig),
                  fatten=5,lwd=1,shape=4, col=cols.int[3], position=position_nudge(x=-0.1)) +
  scale_alpha_manual(values=c(1,0.3,1)) +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  coord_flip() +
  theme_custom() +
  ylab("Coefficient Estimate") +
  theme(axis.text = element_text(size=15),
        legend.position="right",
        axis.title.y = element_blank(),
        strip.text = element_text(size=15)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("ingroup.height" = "Ingroup Height",
                        "ingroup.colless" = "Ingroup Imbalance",
                        "ingroup.gamma" = "Ingroup Gamma")
             ))

# raxml plot only for imbalance for comparison
confint.plot.emp.imb <- all.confint.lates %>%
  filter(resp == "ingroup.colless") %>%
  ggplot(aes(x=var,y=est)) +
  scale_color_manual(values=c(cols.int[2],cols.sig[2],cols.int[2])) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI,alpha=sig),fatten=5,lwd=1,
                  shape=4, col=cols.int[2],position=position_nudge(x=0.1)) +
  geom_pointrange(data=all.confint.cich %>% filter(resp == "ingroup.colless"), 
                  aes(ymin = minCI, ymax = maxCI, alpha=sig),
                  fatten=5,lwd=1,shape=4, col=cols.int[3],position=position_nudge(x=-0.1)) +
  scale_alpha_manual(values=c(1,0.3,1)) +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  coord_flip() +
  theme_custom() +
  ylab("Coefficient Estimate") +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title.y = element_blank(),
        strip.text = element_text(size=15)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("ingroup.height" = "Ingroup Height",
                        "ingroup.colless" = "Ingroup Imbalance",
                        "ingroup.gamma" = "Ingroup Gamma")))

# export at 600x800
# as refbias_emp_model_coeff_raxml_astral.png
p.confint.emp.astr <- ggarrange(confint.plot.emp.imb,confint.plot.astr.emp.all,ncol=2,
                       labels=c("a)","b)"),
                       font.label=list(size=24,font.family="Open Sans"),
                       label.x=0,label.y=1,
                       common.legend=TRUE,legend="none")

#--------------------------------------------------------#
### plots across maf/miss for empirical trees #####
#--------------------------------------------------------#
ic.plot.maf <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.colless,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.colless,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.colless,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  ylim(0.3,1) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 
is.plot.maf <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.sackin,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.sackin,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.sackin,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  ylim(c(0.6,1.4)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
gam.plot.maf <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.gamma,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.gamma,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.gamma,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nGamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
height.plot.maf <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.tree.height,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.tree.height,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.tree.height,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nTree Height") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# across missing data
ic.plot.miss <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.colless,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.colless,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.colless,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  ylim(0.3,1) +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 21, b = 0, l = 0))) 
is.plot.miss <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.sackin,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.sackin,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.sackin,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  ylim(c(0.6,1.4)) +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 
gam.plot.miss <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.gamma,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.gamma,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.gamma,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nGamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
height.plot.miss <- results.raxml.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.tree.height,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.tree.height,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.raxml.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.tree.height,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nTree Height") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

#--------------------------------------------------------#
### plots across maf/miss for astral empirical trees #####
#--------------------------------------------------------#
ic.plot.maf.astr <- results.astral.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.colless,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.colless,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.colless,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  ylim(0.3,1) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  annotate(geom="text",x=7,y=0.4,label="ASTRAL",
           size=5)
is.plot.maf.astr <- results.astral.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.sackin,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.sackin,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.sackin,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  ylim(c(0.6,1.4)) +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# across missing data
ic.plot.miss.astr <- results.astral.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.colless,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.colless,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.colless,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  ylim(0.3,1) +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 21, b = 0, l = 0))) +
  annotate(geom="text",x=0.5,y=0.4,label="ASTRAL",
           size=5)
is.plot.miss.astr <- results.astral.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.sackin,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.sackin,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.astral.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.sackin,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  ylim(c(0.6,1.4)) +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 

### plot combining raxml & astral emp results ####
## plot of distribution of imbalance
cols.int.alpha <- c(scales::alpha("#EE9882",0.8),"#B03618",
                    scales::alpha("#96E3DA",0.8),"#1E7067")

## FIG 8 PANEL x ######
p.density.emp.imb <- results.raxml.lates %>%
  add_row(results.raxml.cichlids %>% mutate(height = "Cichlids")) %>%
  left_join(results.astral.cichlids %>% ungroup() %>% add_row(results.astral.lates) %>%
              group_by(simulation,height,method,quality,missing,maf,int,noref) %>%
              slice_head(n=1) %>% ungroup(), by=c("height","simulation","quality","missing","maf","int","noref"),
            suffix=c(".raxml",".astral")) %>%
  select(height,missing,maf,int,ingroup.colless.raxml,ingroup.colless.astral) %>%
  pivot_longer(cols=c(ingroup.colless.raxml,ingroup.colless.astral),
               names_to="method",names_prefix="ingroup.colless.",values_to="ingroup.colless") %>%
  mutate(method_group = paste(height,method,sep="_")) %>%
  ggplot() +
  geom_density(aes(x=ingroup.colless,fill=method_group,group=method_group), alpha=0.8) +
  theme_custom() +
  facet_grid(cols=vars(height),
             labeller=as_labeller(c("Cichlids" = "Tropheines","Lates" = "Lates")),
             drop=TRUE) +
  xlab("Ingroup Colless Imbalance") +
  scale_fill_manual(values=rev(cols.int.alpha)) +
  theme(strip.text.x = element_text(size=rel(1.2)),
        #legend.position = c(0.9,0.9),
        legend.position = "none",
        plot.margin = unit(c(0,0,0.5,1.5),"cm"))

library(patchwork)
p3.emp <- ggarrange(height.plot.maf,height.plot.miss,
                    gam.plot.maf,gam.plot.miss,
                    ic.plot.maf,ic.plot.miss,
                    ic.plot.maf.astr,ic.plot.miss.astr,
                    
                    #is.plot.maf,is.plot.miss,
                    ncol=2,nrow=4,labels=c("c)","d)","e)","f)","g)","h)","i)","j)"),
                    font.label=list(size=24,font.family="Open Sans"),
                    label.x=0.85,label.y=0.95,
                    common.legend=TRUE,legend="right")

## FIG 8, PANEL C ######
max <- results.raxml.lates %>%
  add_row(results.raxml.cichlids) %>%
  select(height,missing,maf,int,tree.height,ingroup.colless,ingroup.gamma) %>%
  pivot_longer(cols=tree.height:ingroup.gamma,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("tree.height",
                                                "ingroup.colless",
                                                "ingroup.gamma"))) %>%
  group_by(height,maf,missing,int,statistic) %>%
  summarize(mean = mean(estimate)) %>%
  mutate(mean = case_when(statistic == "tree.height" & mean > 0.75 ~ 0.75,
                          statistic == "ingroup.colless" & mean > 0.9 ~ 0.9,
                                 TRUE ~ mean))
height.ic.gam.plot.emp.maf <- results.raxml.lates %>%
  add_row(results.raxml.cichlids) %>%
  select(height,missing,maf,int,tree.height,ingroup.colless,ingroup.gamma) %>%
  pivot_longer(cols=tree.height:ingroup.gamma,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("tree.height",
                                                "ingroup.colless",
                                                "ingroup.gamma"))) %>%
  ggplot(aes(x=maf,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(values=cols.int[c(2,3)]) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  #scale_y_continuous(expand=c(0,0)) +
  xlab("") +
  ylab("\n") +
  theme(axis.title.y = element_text(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines")) +
  geom_blank(data=max,aes(x=maf,y=mean))
height.ic.gam.plot.emp.miss <- results.raxml.lates %>%
  add_row(results.raxml.cichlids) %>%
  select(height,missing,maf,int,tree.height,ingroup.colless,ingroup.gamma) %>%
  pivot_longer(cols=tree.height:ingroup.gamma,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("tree.height",
                                                "ingroup.colless",
                                                "ingroup.gamma"))) %>%
  ggplot(aes(x=missing,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(values=cols.int[c(2,3)]) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=seq(0,1,0.25)) +
  #scale_y_continuous(expand=c(0,0)) +
  xlab("") +
  ylab("\n") +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  expand_limits(results.raxml.cichlids) +
  geom_blank(data=max,aes(x=maf,y=mean))

fig8c <- ggarrange(height.ic.gam.plot.emp.maf,height.ic.gam.plot.emp.miss,
          widths=c(1.1,1), common.legend = TRUE, legend = "right")

fig8 <- confint.plot.emp.all + height.ic.gam.plot.emp.maf + plot_spacer() + height.ic.gam.plot.emp.miss + 
  plot_layout(widths=c(1.5,1,-0.05,1),  nrow=1, guides = "collect")
ggsave(file="Fig8.svg", plot=fig8, height=9.5, width=14)


#p.confint + p3.emp
## saved at 1900 x 1000px
## fig 8 -- refbias_emp_coefficients_all_multipanel
ggarrange(confint.plot.emp.all, p.density.emp.imb, labels=c("a)","b)"),
          font.label=list(size=24,font.family="Open Sans"),ncol=1,heights=c(1.5,1)) + p3.emp

## saved at 1900 x 900px
## 
library(patchwork)
p3.emp.astr <- ggarrange(ic.plot.maf.astr,ic.plot.miss.astr,
                    is.plot.maf.astr,is.plot.miss.astr,
                    
                    #is.plot.maf,is.plot.miss,
                    ncol=2,nrow=3,labels=c("b)","c)","d)","e)","f)","g)","h)","i)"),
                    font.label=list(size=24,font.family="Open Sans"),
                    label.x=0.85,label.y=0.95,
                    common.legend=TRUE,legend="right")

#p.confint + p3.emp

ggarrange(confint.plot.all, labels=c("a)"),
          font.label=list(size=24,font.family="Open Sans")) + p3.emp

## SUPP FIG COMPARING RAXML-ASTRAL EMP ########
max.astr.emp <- results.raxml.lates %>%
  add_row(results.raxml.cichlids) %>%
  select(height,missing,maf,int,ingroup.colless,ingroup.sackin) %>%
  add_row(results.astral.lates %>% select(height,missing,maf,int,ingroup.colless,ingroup.sackin)) %>%
  add_row(results.astral.cichlids %>% select(height,missing,maf,int,ingroup.colless,ingroup.sackin)) %>%
  pivot_longer(cols=ingroup.colless:ingroup.sackin,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("ingroup.colless",
                                                "ingroup.sackin"))) %>%
  group_by(height,maf,missing,int,statistic) %>%
  summarize(mean = mean(estimate)) %>%
  mutate(mean = case_when(statistic == "ingroup.colless" & mean > 1 ~ 1,
                          statistic == "ingroup.sackin" & mean > 1.5 ~ 1.5,
                          TRUE ~ mean))

p3.emp.maf.raxml.astr <- ggarrange(ic.plot.maf,ic.plot.maf.astr,
                                   ic.plot.miss,ic.plot.miss.astr,
                                   is.plot.maf,is.plot.maf.astr,
                                   is.plot.miss,is.plot.miss.astr,
                         #is.plot.maf,is.plot.miss,
                         ncol=4,nrow=2,labels=c("a)","b)","c)","d)","e)","f)","g)","h)","i)"),
                         font.label=list(size=24,font.family="Open Sans"),
                         label.x=0.85,label.y=0.95,
                         common.legend=TRUE,legend="right")

(p.density.emp.imb + plot_spacer()) / p3.emp.maf.raxml.astr +
  plot_layout(heights=c(1,3),widths=c(2,1))

height.ic.gam.plot.emp.maf <- results.raxml.lates %>%
  add_row(results.raxml.cichlids) %>%
  select(height,missing,maf,int,ingroup.colless,ingroup.sackin) %>%
  pivot_longer(cols=ingroup.colless:ingroup.sackin,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("ingroup.colless",
                                                "ingroup.sackin"))) %>%
  ggplot(aes(x=maf,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(values=cols.int[c(2,3)]) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  #scale_y_continuous(expand=c(0,0)) +
  xlab("") +
  ylab("\n") +
  theme(axis.title.y = element_text(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines")) +
  geom_blank(data=max.astr.emp %>% filter(statistic %in% c("ingroup.colless","ingroup.sackin")),aes(x=maf,y=mean))
height.ic.gam.plot.emp.miss <- results.raxml.lates %>%
  add_row(results.raxml.cichlids) %>%
  select(height,missing,maf,int,ingroup.colless,ingroup.sackin) %>%
  pivot_longer(cols=ingroup.colless:ingroup.sackin,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("ingroup.colless",
                                                "ingroup.sackin"))) %>%
  ggplot(aes(x=missing,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(values=cols.int[c(2,3)]) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=seq(0,1,0.25)) +
  #scale_y_continuous(expand=c(0,0)) +
  xlab("") +
  #ylab("\n") +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  expand_limits(results.raxml.cichlids) +
  geom_blank(data=max.astr.emp %>% filter(statistic %in% c("ingroup.colless","ingroup.sackin")),aes(x=maf,y=mean))
height.ic.gam.plot.emp.maf.astr <- results.astral.lates %>%
  add_row(results.astral.cichlids) %>%
  select(height,missing,maf,int,ingroup.colless,ingroup.sackin) %>%
  pivot_longer(cols=ingroup.colless:ingroup.sackin,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("ingroup.colless",
                                                "ingroup.sackin"))) %>%
  ggplot(aes(x=maf,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(values=cols.int[c(2,3)]) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  scale_y_continuous(expand=c(0,0)) +
  xlab("") +
  #ylab("\n") +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.ticks.y = element_blank()) +
  geom_blank(data=max.astr.emp %>% filter(statistic %in% c("ingroup.colless","ingroup.sackin")),aes(x=maf,y=mean))
height.ic.gam.plot.emp.miss.astr <- results.astral.lates %>%
  add_row(results.astral.cichlids) %>%
  select(height,missing,maf,int,ingroup.colless,ingroup.sackin) %>%
  pivot_longer(cols=ingroup.colless:ingroup.sackin,names_to="statistic",values_to="estimate") %>%
  mutate(statistic = factor(statistic, levels=c("ingroup.colless",
                                                "ingroup.sackin"))) %>%
  ggplot(aes(x=missing,y=estimate,col=height,shape=int)) + 
  stat_summary(alpha=0.8,size=4,geom="point") +
  stat_summary(geom="line",lty=2) +
  scale_color_manual(values=cols.int[c(2,3)]) +
  scale_shape_manual(values=c(17,15),
                     labels=c("Outgroup\nReference","Ingroup\nReference"))  +
  theme_custom() +
  facet_grid(rows=vars(statistic), scales="free_y") +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25),labels=seq(0,1,0.25)) +
  #scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.25)) +
  xlab("") +
  ylab("\n") +
  theme(axis.title.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  expand_limits(results.raxml.cichlids) +
  geom_blank(data=max.astr.emp %>% filter(statistic %in% c("ingroup.colless","ingroup.sackin")),aes(x=maf,y=mean))

height.ic.gam.plot.emp.maf + height.ic.gam.plot.emp.maf.astr +
  height.ic.gam.plot.emp.miss + height.ic.gam.plot.emp.miss.astr +
  plot_layout(guides="collect", widths=c(1,1,1,1)) &
  theme(legend.position='bottom')

##------------------------------------##
## PLOTS FOR ASTRAL TREES ##############
##------------------------------------##

all.confint.astr <- confint.rf.long.astr %>% 
            as_tibble() %>% 
            mutate(resp = "RF.Dist.ML") %>%
  mutate(height = "long") %>%
  add_row(confint.rf.med.astr %>%
            as_tibble() %>%
            mutate(resp = "RF.Dist.ML") %>%
            mutate(height = "med")) %>%
  add_row(confint.rf.short.astr %>% 
            as_tibble() %>%
            mutate(resp = "RF.Dist.ML") %>%
            mutate(height = "short")) %>%
  add_row(confint.rf.all.astr %>% 
            as_tibble() %>%
            mutate(resp = "RF.Dist.ML") %>%
            mutate(height = "all")) %>%
  mutate(categ = case_when(height == "all" ~ "all",
                           TRUE ~ "by.height"),
         nudge = case_when(height == "short" ~ -0.2,
                           height == "long" ~ 0.2,
                           height == "med" ~ 0,
                           height == "all" ~ 0)) %>%
  mutate(col_sig = case_when(sig == "ns" ~ "ns",
                             TRUE ~ height))
cols.int.sig <- c("#e2d200","#2A9D8F","gray","#e76f51")
confint.plot.bygroup.astr <- all.confint.astr  %>%
  filter(categ == "by.height") %>%
  filter(resp %in% c("RF.Dist.ML")) %>%
  ggplot(aes(x=var,y=est,col=height)) +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols.int,guide="none")+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI, shape=col_sig, alpha=col_sig),
                  fill="white",fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint.astr$nudge[all.confint.astr$categ == "by.height" & all.confint.astr$resp == "RF.Dist.ML"])) +
  # geom_point(data=. %>% filter(col_sig == "ns"), aes(x=var,y=est,group=height),
  #            color="black",shape=16, fill="black",size=4,
  #            position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height"])) +
  scale_shape_manual(values=c(16,16,21,16,16)) +
  scale_alpha_manual(values=c(1,1,0.4,1,1)) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

confint.plot.all.astr <- all.confint.astr %>%
  filter(height == "all") %>%
  filter(resp %in% c("RF.Dist.ML")) %>%
  ggplot(aes(x=var,y=est,col=sig)) +
  scale_color_manual(values=cols.sig[c(2,3)]) +
  #scale_color_manual(values=cols.sig) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint.astr$nudge), alpha=1, shape=4) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

## subset raxml plots to only rfdist
confint.plot.bygroup.rf <- all.confint  %>%
  filter(categ == "by.height") %>%
  filter(resp %in% c("RF.Dist.ML")) %>%
  ggplot(aes(x=var,y=est,col=height)) +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols.int,guide="none")+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI, shape=col_sig, alpha=col_sig),
                  fill="white",fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height" & all.confint$resp == "RF.Dist.ML"])) +
  # geom_point(data=. %>% filter(col_sig == "ns"), aes(x=var,y=est,group=height),
  #            color="black",shape=16, fill="black",size=4,
  #            position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height"])) +
  scale_shape_manual(values=c(16,16,21,16,16)) +
  scale_alpha_manual(values=c(1,1,0.4,1,1)) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

confint.plot.all.rf <- all.confint %>%
  filter(height == "all") %>%
  filter(resp %in% c("RF.Dist.ML")) %>%
  ggplot(aes(x=var,y=est,col=sig)) +
  scale_color_manual(values=cols.sig[c(2,3)]) +
  #scale_color_manual(values=cols.sig) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge[all.confint$resp == "RF.Dist.ML"]), alpha=1, shape=4) +
  coord_flip() +
  theme_custom() +
  scale_x_discrete(labels=rev(c("miss","mac","ref:miss","ref:mac","ref"))) +
  theme(axis.text = element_text(size=15),
        legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=15),
        plot.margin=margin(t=45)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

## saved at 1000 x 1000px
p1.astr <- ggarrange(confint.plot.all.astr,confint.plot.bygroup.astr,ncol=1,
                hjust=-0.2,
                heights=c(1,1.5),labels=c("a) Models for all trees","b) Models by tree height"),
                font.label=list(size=24,family="Open Sans"),
                common.legend=TRUE,legend.grob=get_legend(confint.plot.bygroup),legend="none")
p1.top <- ggarrange(confint.plot.all.rf,confint.plot.all.astr)
p1.bottom <- ggarrange(confint.plot.bygroup.rf,confint.plot.bygroup.astr)

p1.raxml.astr.comp <- ggarrange(p1.top,p1.bottom,ncol=1,
                                hjust=-0.2,
                                heights=c(1,1.5),labels=c("a) Models for all trees","b) Models by tree height"),
                                font.label=list(size=24,family="Open Sans"),
                                common.legend=TRUE,legend.grob=get_legend(confint.plot.bygroup),legend="none")

#----------------------------------#
### table with param estimates #####
#----------------------------------#
confint.table.astr <- all.confint.astr %>%
  relocate(resp,height,var,est,minCI,maxCI) %>%
  select(resp:sig) %>%
  mutate(sig2 = case_when(sig == "pos" ~ "*",
                          sig == "neg" ~ "*",
                          sig == "ns" ~ "")) %>%
  pivot_wider(id_cols=c(resp,var),names_from=height,values_from=c(est:maxCI,sig2),names_glue="{height}_{.value}",
              names_sort=TRUE) %>%
  relocate(resp,var,starts_with("all"),starts_with("long"),starts_with("med"),starts_with("short")) %>%
  select(!resp) %>%
  kbl(booktabs = T,digits=3,
      format="latex",
      col.names=c("Variable","Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig")) %>%
  add_header_above(c(" ","All Trees" = 4, "Low ILS" = 4, "Medium ILS" = 4,  "High ILS" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Robinson-Foulds distance to true tree", 1, 5)
#writeLines(confint.table, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_fullsims_astral.tex')
#writeLines(confint.table, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table_subsamp.tex')


###############################
### plot of dxy comparison ####
###############################
# saved as refbias_dxy_comparison_sims_emp.png at 800x600px
refdist %>% 
  select(simulation,height,int,avg_dxy) %>% 
  rename(sim=simulation, dxy=avg_dxy) %>% 
  add_row(dxy) %>% 
  ggplot(aes(x=height,y=dxy)) + 
    stat_summary(aes(color=int, shape=int), position = position_nudge(x=0.1),
                 size=1) + 
    gghalves::geom_half_point(aes(color=int, shape=int),
                              side="l",alpha=0.3,range_scale=0.3,
                              position=position_nudge(x=0.1)) + 
    theme_custom() + 
    geom_vline(xintercept=2.5) + 
    ylim(0,0.8) + 
    geom_rect(xmin=0,xmax=6,ymin=0.76,ymax=0.9) + 
    annotate(geom="text",x=1.5,y=0.8,label="Empirical",
             color="white", family="Open Sans", size=5) + 
    annotate(geom="text",x=4,y=0.8,label="Simulations",
             color="white",family="Open Sans",size=5) + 
    theme(axis.title.x = element_blank()) + 
    ylab("Average dxy") + 
    scale_x_discrete(labels=c("Tropheines","Lates","Low ILS","Med ILS","High ILS")) +
    scale_shape_manual(values=c(17,15)) +
    scale_color_manual(values=c("gray60","black"))
