## correlation plots for refbias analyses
## updated apr 2023, j. rick

pal <- wesanderson::wes_palette("Zissou1", 200, type = "continuous")

corr.astral <- results.astral %>% 
  select(maf,missing,avg_dxy,std.ingroup.sackin,std.ingroup.colless,RF.Dist.ML,Q.Dist.ML,CI.Dist.ML,mean.support) %>% 
  cor(method="pearson",use="pairwise.complete.obs") %>% 
  reshape2::melt() 

corr.raxml <- results.raxml %>% 
  select(maf,missing,avg_dxy,
         tree.height,std.ingroup.sackin,std.ingroup.gamma,std.ingroup.colless,
         RF.Dist.ML,Q.Dist.ML,CI.Dist.ML,mean.support,ingroup.tree.height) %>% 
  cor(method="pearson",use="pairwise.complete.obs") %>% 
  reshape2::melt()

corr.raxml %>%
  full_join(corr.astral,by=c("Var1","Var2"),suffix=c(".raxml",".astral")) %>%
  ggplot(aes(x=as.integer(Var1),y=as.integer(Var2))) + 
  geom_tile(data= . %>% filter(as.integer(Var1) < as.integer(Var2)), 
            aes(width=abs(value.raxml),height=abs(value.raxml),fill=value.raxml)) +
  geom_tile(data= . %>% filter(as.integer(Var1) > as.integer(Var2)), 
            aes(width=abs(value.astral),height=abs(value.astral),fill=value.astral)) +
  geom_text(data= . %>% filter(is.na(value.astral)) %>% filter(as.integer(Var1) > as.integer(Var2)),
            label="X") +
  geom_tile(color="gray50",fill=NA) +
  theme_custom() + 
  scale_fill_gradientn(limits = c(-1,1),colours = rev(pal),name="Pearson\nCorrelation\nCoefficient") +
  scale_size_continuous(guide="none",limits=c(0,1),range=c(0.5,3)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text())

corr.plot.raxml <- results.raxml %>% 
  select(maf,missing,avg_dxy,
         tree.height,std.ingroup.sackin,std.ingroup.gamma,std.ingroup.colless,
         RF.Dist.ML,Q.Dist.ML,CI.Dist.ML,mean.support,ingroup.tree.height) %>% 
  cor() %>% 
  reshape2::melt() %>% 
  ggplot(aes(x=Var1,y=Var2,fill=value)) + 
  geom_tile(aes(width=value,height=value,size=value)) +
  geom_tile(color="gray50",fill=NA) +
  theme_custom() + 
  scale_fill_gradientn(limits = c(-1,1),colours = rev(pal),name="Pearson\nCorrelation\nCoefficient") +
  scale_size_continuous(guide="none",limits=c(-1,1)) +
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text())

corr.plot.astral <- results.astral %>% 
  select(maf,missing,avg_dxy,std.ingroup.sackin,std.ingroup.colless,RF.Dist.ML,Q.Dist.ML,CI.Dist.ML,mean.support) %>% 
  cor() %>% 
  reshape2::melt() %>% 
  ggplot(aes(x=Var1,y=Var2,fill=value)) + 
    geom_tile(aes(width=value,height=value,size=value)) +
    geom_tile(color="gray50",fill=NA) +
    theme_custom() + 
    scale_fill_gradientn(limits = c(-1,1),colours = rev(pal),name="Pearson\nCorrelation\nCoefficient") +
    scale_size_continuous(guide="none",limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 60, hjust=1),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text())

# saved at 1600x800pt as refbias_raxml_astral_param_corr.png
ggarrange(corr.plot.raxml, corr.plot.astral, common.legend = TRUE,
          labels=c("a) RAxML","b) ASTRAL"))


################
## heatmap for correlation between raxml and astral stats

library(tidyverse)
library(here)

corr.raxml.astral <- read_csv(here("analysis","interaction_matrix_refbias")) %>%
  pivot_longer(cols=mean.support:maf,names_to="var2",values_to = "value")
corr.raxml.astral %>%
  ggplot(aes(x=`...1`,y=var2)) +
  geom_tile(aes(fill=value)) +
  geom_text(aes(label=round(value,3))) +
  scale_fill_viridis()
