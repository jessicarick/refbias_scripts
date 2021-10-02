########################
### Summary plots for refbias analyses
### and plots exploring the data
### 
### J. Rick, Sept 2021
########################

output <- "092321-output"
results.raxml <- read.csv(paste("output/new/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

########################
## heatmaps to show relationship of maf and miss to variables of interest
cols.int <- c("#005F73","gray80","#5FB89D")
cols.sig <- c("#E6B749","gray80","#6A8D4E")

results_mean <- results.mod %>%
  select(int,height,maf,missing,RF.Dist.ML,std.ingroup.gamma,std.ingroup.colless,Q.Dist.ML,CI.Dist.ML) %>%
  #pivot_wider(id_cols="maf",names_from="missing",values_from="RF.Dist.ML",values_fn=mean) %>%
  group_by(maf,missing,int,height) %>%
  mutate(mean_rf = mean(RF.Dist.ML),
         mean_gam = mean(std.ingroup.gamma),
         mean_colless = mean(std.ingroup.colless),
         mean_qdist = mean(Q.Dist.ML),
         mean_cidist = mean(CI.Dist.ML)) %>%
  select(maf,mean_rf:mean_cidist) %>%
  pivot_longer(cols=mean_rf:mean_cidist,names_to="param",values_to="mean") 

# rf_interp <- interp::interp(x=results_mean$maf[results_mean$param == "mean_rf"], y=results_mean$missing[results_mean$param == "mean_rf"], z=results_mean$mean[results_mean$param == "mean_rf"],
#               output="grid",duplicate = "mean")

contour_rf <- results_mean %>%
  filter(param == "mean_rf") %>%
  ggplot() +
    geom_contour_filled(aes(x=maf,y=missing,z=mean)) +
    theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="RF Distance") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count") +
  facet_grid(vars(int),vars(height))

contour_gam <- results_mean %>%
  filter(param == "mean_gam") %>%
  ggplot() +
  geom_contour_filled(aes(x=maf,y=missing,z=mean)) +
  theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="Gamma") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count") +
  facet_grid(vars(int),vars(height))

contour_imb <- results_mean %>%
  filter(param == "mean_colless") %>%
  ggplot() +
  geom_contour_filled(aes(x=maf,y=missing,z=mean)) +
  theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="Colless\nImbalance") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count") +
  facet_grid(vars(int),vars(height))

contour_qdist <- results_mean %>%
  filter(param == "mean_qdist") %>%
  ggplot() +
  geom_contour_filled(aes(x=maf,y=missing,z=mean)) +
  theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="Quartet\nDistance") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count") +
  facet_grid(vars(int),vars(height))

contour_ci <- results_mean %>%
  filter(param == "mean_cidist") %>%
  ggplot() +
  geom_contour_filled(aes(x=maf,y=missing,z=mean)) +
  theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="CI Distance") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count") +
  facet_grid(vars(int),vars(height))

# exported at 1800 x 800px
ggarrange(contour_rf, contour_qdist, contour_ci,
          contour_gam, contour_imb, nrow=2, ncol=3,
          labels="AUTO",font.label=list(size=24))

########################
### and plots of model coefficients
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
                           height == "all" ~ 0))
confint.plot.bygroup <- all.confint  %>%
  filter(categ == "by.height") %>%
  ggplot(aes(x=var,y=est,col=sig,shape=height)) +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols.sig)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge), alpha=1) +
  coord_flip() +
  theme_custom() +
  theme(axis.text = element_text(size=15),
        legend.position="right",
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
  scale_color_manual(values=cols.sig) +
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge), alpha=1, shape=4) +
  coord_flip() +
  theme_custom() +
  theme(axis.text = element_text(size=15),
        legend.position="right",
        axis.title = element_blank(),
        strip.text = element_text(size=15)) +
  facet_grid(cols=vars(resp),scales="free_x",
             labeller = labeller(
               resp = c("RF.Dist.ML" = "RF Dist to True Tree",
                        "std.ingroup.colless" = "Ingroup Imbalance",
                        "std.ingroup.gamma" = "Ingroup Gamma")
             ))

## saved at 1000 x 1000px
ggarrange(confint.plot.all,confint.plot.bygroup,ncol=1,
          heights=c(1,2),labels="AUTO",font.label=list(size=24,family="Open Sans"),
          common.legend=TRUE,legend.grob=get_legend(confint.plot.bygroup),legend="right")

## table with param estimates
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
  kbl(booktabs = T,digits=3,format="latex",
      col.names=c("Variable","Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig",
                  "Est","Min CI","Max CI","Sig")) %>%
  add_header_above(c(" ","All Trees" = 4, "Low ILS" = 4, "Medium ILS" = 4,  "High ILS" = 4)) %>%
  kable_styling(latex_options = c("repeat_header")) %>%
  pack_rows("Standardized ingroup gamma", 1, 5) %>%
  pack_rows("Robinson-Foulds distance to true tree", 6, 10) %>%
  pack_rows("Standardized ingroup Colless imbalance", 11, 15) %>%
  landscape()
writeLines(confint.table, '~/Dropbox/Apps/Overleaf/Reference Bias in SNP-based Phylogenetics/confint_table.tex')


########################
## plot of distances to true tree by maf & int
## because this seems to be most interesting from models
rf.plot <- results.mod %>%
  ggplot(aes(x=maf,y=RF.Dist.ML,shape=int,col=height)) +
  stat_summary(alpha=0.8) +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Robinson-Foulds\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
q.plot <- results.mod %>%
  ggplot(aes(x=maf,y=Q.Dist.ML,shape=int,col=height)) +
  stat_summary(alpha=0.8) +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Quartet\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
ci.plot <- results.mod %>%
  ggplot(aes(x=maf,y=CI.Dist.ML,shape=int,col=height)) +
  stat_summary(alpha=0.8) +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Clustering Information\nDistance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

## saved at 800px x 1000px
ggarrange(rf.plot,q.plot,ci.plot,ncol=1,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)

gam.plot <- results.mod %>%
  ggplot(aes(x=maf,y=std.ingroup.gamma,shape=int,col=height)) +
  stat_summary(alpha=0.8) +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized\ningroup gamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  geom_hline(yintercept=0)
ic.plot <- results.mod %>%
  ggplot(aes(x=maf,y=std.ingroup.colless,shape=int,col=height)) +
  stat_summary(alpha=0.8) +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_hline(yintercept=0)
is.plot <- results.mod %>%
  ggplot(aes(x=maf,y=std.ingroup.sackin,shape=int,col=height)) +
  stat_summary(alpha=0.8) +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Count") +
  ylab("Standardized ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  geom_hline(yintercept=0)
## saved at 800px x 1000px
ggarrange(gam.plot,ic.plot,is.plot,ncol=1,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)
