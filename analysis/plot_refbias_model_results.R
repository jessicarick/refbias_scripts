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
results.raxml <- read.csv(here("output","new",paste0(output,"-raxml.csv")),header=TRUE,row.names=1,sep=",")
figdir <- "~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/figures/"

#cols.int <- c("#005F73","gray80","#5FB89D")
cols.int <- c("#e2d200","#2A9D8F","#e76f51") # colors for tree heights
cols.sig <- c("black","gray80","black")

## preparing data object
results.mod <- results.raxml[results.raxml$simulation > 15 & results.raxml$simulation < 26 & results.raxml$RF.Dist.ML < 0.9,]
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.numeric(as.character(results.mod$missing))
results.mod$maf <- as.numeric(as.character(results.mod$maf))
results.mod$int <- as.factor(results.mod$int)

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
  filter(snps < 1000000) %>% 
  ggplot(aes(x=as.factor(maf),y=snps,fill=int)) + 
  geom_boxplot() + theme_custom() + 
  scale_fill_manual(values=cols.int) + 
  xlab("Minor Alelle Count") + 
  ylab("Number of SNPs") +
  facet_wrap(~height)
miss.sites.emp <- results.emp %>% 
  filter(snps < 1000000) %>%
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
results_mean <- results.mod %>%
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
results_mean <- results.emp.lates %>%
  select(int,height,maf,missing,ingroup.gamma,ingroup.colless,ingroup.tree.height) %>%
  #pivot_wider(id_cols="maf",names_from="missing",values_from="RF.Dist.ML",values_fn=mean) %>%
  group_by(maf,missing,int,height) %>%
  mutate(mean_gam = mean(ingroup.gamma),
         mean_colless = mean(ingroup.colless),
         mean_height = mean(ingroup.tree.height)) %>%
  select(maf,mean_gam:mean_height) %>%
  pivot_longer(cols=mean_gam:mean_height,names_to="param",values_to="mean") %>%
  filter(height=="Lates")

# rf_interp <- interp::interp(x=results_mean$maf[results_mean$param == "mean_rf"], y=results_mean$missing[results_mean$param == "mean_rf"], z=results_mean$mean[results_mean$param == "mean_rf"],
#               output="grid",duplicate = "mean")
categ <- "fullsim"
categ <- "subsamp"
categ <- "emp"
for (h in unique(results.mod$height)) {
  results_mean_h <- results_mean %>%
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
    heatmap %>%
      ggexport(filename = paste0(figdir,"full_sims/refbias_heatmap_",h,"_byINT.png"), width=2000,height=1400)
  }
  print(heatmap)
}

#----------------------------------#
### plots of model coefficients ####
#----------------------------------#
ggarrange(vars.plots.rf,vars.plot.gam,vars.plots.imb,labels="AUTO",font.label=list(size=24),ncol=1)
ggarrange(vars.rf.all.bars,vars.gam.all.bars,vars.imb.all.bars,labels="AUTO",font.label=list(size=24),nrow=1)

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
         nudge = case_when(height == "short" ~ -0.2,
                           height == "long" ~ 0.2,
                           height == "med" ~ 0,
                           height == "all" ~ 0)) %>%
  mutate(col_sig = case_when(sig == "ns" ~ "ns",
                             TRUE ~ height))
cols.int.sig <- c("#e2d200","#2A9D8F","gray","#e76f51")
confint.plot.bygroup <- all.confint  %>%
  filter(categ == "by.height") %>%
  ggplot(aes(x=var,y=est,col=height)) +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols.int,guide="none")+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI, shape=col_sig, alpha=col_sig),
                  fill="white",fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge[all.confint$categ == "by.height"])) +
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
        strip.text = element_text(size=15)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
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
        strip.text = element_text(size=15)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

## saved at 1000 x 1000px
p1 <- ggarrange(confint.plot.all,confint.plot.bygroup,ncol=1,
          heights=c(1,1.5),labels=c("a)","b)"),font.label=list(size=24,family="Open Sans"),
          common.legend=TRUE,legend.grob=get_legend(confint.plot.bygroup),legend="none")

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
        legend.text = element_text(size=rel(1.5)))
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
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
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
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

## saved at 800px x 1000px
p2 <- ggarrange(rf.plot,q.plot,ci.plot,ncol=1,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)

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
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0)
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
  ylab("Standardized ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))) +
  geom_hline(yintercept=0)
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
## saved at 800px x 1000px
p3 <- ggarrange(gam.plot,ic.plot,is.plot,ncol=1,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)

#library(patchwork)
p4 <- ggarrange(rf.plot, ic.plot, gam.plot, ncol=1, 
                labels=c("c)","d)","e)"), font.label=list(size=24,font.family="Open Sans"),
                common.legend=TRUE, legend="right",
                label.x=0.3, label.y=0.95)
p1 + p4 + 
  plot_layout(widths = c(1.2, 1))

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
                   resp = factor(resp, levels = c("ingroup.height","ingroup.colless","ingroup.gamma")))
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
         resp = factor(resp, levels = c("ingroup.height","ingroup.colless","ingroup.gamma")))
confint.plot.all <- all.confint.lates %>%
  ggplot(aes(x=var,y=est)) +
  scale_color_manual(values=c(cols.int[2],cols.sig[2],cols.int[2])) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI,alpha=sig),fatten=5,lwd=1,
                  shape=4, col=cols.int[2],position=position_nudge(x=0.1)) +
  geom_pointrange(data=all.confint.cich, aes(ymin = minCI, ymax = maxCI, alpha=sig),
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
                        "ingroup.gamma" = "Ingroup Gamma")
             ))
#confint.plot.all.lates <- confint.plot.all
#confint.plot.all.cich <- confint.plot.all

p.confint <- ggarrange(confint.plot.all.cich,confint.plot.all.lates,ncol=1,
          labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          label.x=0,label.y=1,
          common.legend=TRUE,legend="right")

# across maf
th.plot.maf <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.tree.height,shape=as.factor(int))) +
  stat_summary(alpha=0.8,aes(col=cols.int[3])) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates %>% mutate(height = "lates"),
               aes(x=as.numeric(as.character(maf)),y=ingroup.tree.height,group=int),
               geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates %>% mutate(height = "lates"),
               aes(x=as.numeric(as.character(maf)),y=ingroup.tree.height,shape=int,col=cols.int[2]),
               alpha=0.8,inherit.aes=FALSE) +
  theme_custom() +
  scale_shape_manual(values=c(17,15), labels=c("Outgroup","Ingroup"), name="Reference\nGenome")  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nTree Height") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        legend.position="right",
        legend.text = element_text(size=rel(1.5))) +
  scale_color_identity("Dataset", labels=c("Lates","Tropheines"), guide="legend") +
  guides(shape = guide_legend(override.aes = list(size = 1)),
         color = guide_legend(override.aes = list(size = 1)))
gam.plot.maf <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.gamma,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.gamma,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.gamma,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\ngamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 17, b = 0, l = 0)))
ic.plot.maf <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.colless,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.colless,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.colless,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 
is.plot.maf <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=maf,y=ingroup.sackin,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.sackin,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(maf)),y=ingroup.sackin,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# across missing data
th.plot.miss <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.tree.height,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.tree.height,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.tree.height,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nTree Height") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
gam.plot.miss <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.gamma,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.gamma,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.gamma,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\ngamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0)))
ic.plot.miss <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.colless,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.colless,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.colless,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 21, b = 0, l = 0))) 
is.plot.miss <- results.emp.cichlids %>% 
  mutate(height = "cichlids") %>%
  ggplot(aes(x=missing,y=ingroup.sackin,shape=as.factor(int))) +
  stat_summary(alpha=0.8,col=cols.int[3]) +
  stat_summary(geom="line",lty=2,col=cols.int[3]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.sackin,group=int),geom="line",lty=2,inherit.aes=FALSE,col=cols.int[2]) +
  stat_summary(data=results.emp.lates,aes(x=as.numeric(as.character(missing)),y=ingroup.sackin,shape=int),alpha=0.8,inherit.aes=FALSE,col=cols.int[2]) +
  theme_custom() +
  scale_shape_manual(values=c(17,15))  +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 

## saved at 1900 x 900px
library(patchwork)
p3.emp <- ggarrange(th.plot.maf,th.plot.miss,
                    ic.plot.maf,ic.plot.miss,
                    gam.plot.maf,gam.plot.miss,
                    
                    #is.plot.maf,is.plot.miss,
                    ncol=2,nrow=3,labels=c("b)","c)","d)","e)","f)","g)","h)","i)"),
                    font.label=list(size=24,font.family="Open Sans"),
                    label.x=0.85,label.y=0.95,
                    common.legend=TRUE,legend="right")

#p.confint + p3.emp

ggarrange(confint.plot.all, labels=c("a)"),
          font.label=list(size=24,font.family="Open Sans")) + p3.emp
