########################
### Summary plots for refbias analyses with empirical data
### and plots exploring the data
### 
### J. Rick, Sept 2021
########################

output <- "082021-lates-emp-output"
results.raxml <- read.csv(paste("output/new/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")

########################
## heatmaps to show relationship of maf and miss to variables of interest
cols <- c("#F2AD00","gray80","#00A08A")

results_mean <- results.mod %>%
  #select(maf,missing,ingroup.gamma,ingroup.colless) %>%
  #pivot_wider(id_cols="maf",names_from="missing",values_from="RF.Dist.ML",values_fn=mean) %>%
  group_by(maf,missing,int) %>%
  mutate(mean_th = mean(ingroup.tree.height),
         mean_gam = mean(ingroup.gamma),
         mean_colless = mean(ingroup.colless)) %>%
  select(maf,mean_th:mean_colless) %>%
  pivot_longer(cols=mean_th:mean_colless,names_to="param",values_to="mean") 

# rf_interp <- interp::interp(x=results_mean$maf[results_mean$param == "mean_rf"], y=results_mean$missing[results_mean$param == "mean_rf"], z=results_mean$mean[results_mean$param == "mean_rf"],
#               output="grid",duplicate = "mean")

contour_gam <- results_mean %>%
  filter(param == "mean_gam") %>%
  ggplot() +
  geom_contour_filled(aes(x=maf,y=missing,z=mean),binwidth=0.5) +
  theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="Gamma") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count")# +
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
  xlab("Minor Allele Count") #+
  facet_grid(vars(int),vars(height))

contour_th <- results_mean %>%
  filter(param == "mean_th") %>%
  ggplot() +
  geom_contour_filled(aes(x=maf,y=missing,z=mean)) +
  theme_custom() +
  scale_fill_viridis(discrete=TRUE,option="C", name="Tree Height") +
  theme(legend.title = element_text(size=rel(1.8)),
        axis.title = element_text(size=rel(1.8)),
        strip.text = element_text(size=rel(1.3)),
        strip.background = element_rect(fill="white")) +
  ylab("Missing Data") +
  xlab("Minor Allele Count") #+
  facet_grid(vars(int),vars(height))

# exported at 1800 x 400px
ggarrange(contour_th,
          contour_gam, contour_imb, nrow=1, ncol=3,
          labels="AUTO",font.label=list(size=24))

########################
### and plots of model coefficients
#ggarrange(vars.plot.gam,vars.plots.imb,labels="AUTO",font.label=list(size=24),ncol=1)
#ggarrange(vars.gam.all.bars,vars.imb.all.bars,labels="AUTO",font.label=list(size=24),nrow=1)

all.confint <- confint.gam.all %>% 
             as_tibble() %>% 
             mutate(resp = "std.ingroup.gamma") %>% 
             add_row(confint.imb.all %>% 
                       as_tibble() %>% 
                       mutate(resp = "std.ingroup.colless")) %>%
             mutate(height = "all")

confint.plot.all <- all.confint %>%
  ggplot(aes(x=var,y=est,col=sig)) +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=5,lwd=1, 
                  position=position_nudge(x=all.confint$nudge), alpha=1) +
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

########################
## plot of distances to true tree by maf & int
## because this seems to be most interesting from models
gam.plot <- results.mod %>%
  ggplot(aes(x=maf,y=ingroup.gamma,shape=int,col=int)) +
  stat_summary() +
  stat_summary(geom="line",lty=2) +
  geom_point() +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Frequency") +
  ylab("Ingroup\ngamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  geom_hline(yintercept=0)
ic.plot <- results.mod %>%
  ggplot(aes(x=maf,y=ingroup.colless,shape=int,col=int)) +
  stat_summary() +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Frequency") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 
is.plot <- results.mod %>%
  ggplot(aes(x=maf,y=ingroup.sackin,shape=int,col=int)) +
  stat_summary() +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Minor Allele Frequency") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
## saved at 800px x 1000px
ggarrange(gam.plot,ic.plot,is.plot,ncol=1,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)

gam.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=ingroup.gamma,shape=int,col=int)) +
  stat_summary() +
  stat_summary(geom="line",lty=2) +
  geom_point() +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data Threshold") +
  ylab("Ingroup\ngamma") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  geom_hline(yintercept=0)
ic.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=ingroup.colless,shape=int,col=int)) +
  stat_summary() +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data Threshold") +
  ylab("Ingroup\nColless imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) 
is.plot.miss <- results.mod %>%
  ggplot(aes(x=missing,y=ingroup.sackin,shape=int,col=int)) +
  stat_summary() +
  stat_summary(geom="line",lty=2) +
  theme_custom() +
  #scale_x_continuous(limits=c(0,10),breaks=seq(0,10,1),labels=seq(0,10,1)) +
  xlab("Missing Data Threshold") +
  ylab("Ingroup\nSackin imbalance") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
## saved at 800px x 1000px
ggarrange(gam.plot.miss,ic.plot.miss,is.plot.miss,ncol=1,labels=c("D","E","F"),
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.x=0.17,label.y=0.95)

## saved at 1000x1000px
ggarrange(gam.plot,gam.plot.miss,
          ic.plot,ic.plot.miss,
          is.plot,is.plot.miss,
          ncol=2,nrow=3,
          byrow=FALSE,labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          common.legend=TRUE,legend="right",
          label.y=0.3,label.x=0.9)
