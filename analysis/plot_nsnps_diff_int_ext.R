## plots of difference in number of sites INT-EXT
## for refbias paper

results.raxml.all <- results.raxml %>%
  bind_rows(results.raxml.cichlids %>% mutate(sites = snps),
            results.raxml.lates %>% mutate(sites = snps))

p.emp.diff.sites <- results.raxml.all %>%
  filter(sites < 1e6 & height %in% c("Cichlids","Lates") & maf > 0) %>%
  group_by(height,simulation,method,missing,maf) %>%
  pivot_wider(names_from = int, names_prefix="sites",values_from=sites) %>% 
  select(sitesINT,sitesEXT) %>%
  summarize(diff_int_ext = max(sitesINT,na.rm=T)-max(sitesEXT,na.rm=T)) %>%
  ggplot() +
    geom_density(aes(x=diff_int_ext,fill=as.factor(maf)),alpha=0.3) +
    geom_blank(aes(x=-diff_int_ext)) +
    facet_wrap(vars(height),scales="fixed") +
    geom_vline(xintercept=0,lty=2,lwd=0.5,alpha=0.2) +
    #annotate(geom="text",x=max(diff_int_ext),y=)
    theme_custom() +
    theme(legend.title = element_text()) +
    xlab("SNP Difference INT-EXT") +
    labs(fill="Minor Allele\nCount Threshold")
p.sim.diff.sites <- results.raxml.all %>%
  filter(!(height %in% c("Cichlids","Lates")) & maf > 0) %>%
  group_by(height,simulation,method,missing,maf) %>%
  pivot_wider(names_from = int, names_prefix="sites",values_from=sites) %>% 
  select(sitesINT,sitesEXT) %>%
  summarize(diff_int_ext = max(sitesINT,na.rm=T)-max(sitesEXT,na.rm=T)) %>%
  ggplot() +
  geom_density(aes(x=diff_int_ext,fill=as.factor(maf)),alpha=0.3) +
  geom_blank(aes(x=-diff_int_ext)) +
  facet_wrap(vars(height),scales="fixed") +
  geom_vline(xintercept=0,lty=2,lwd=0.5,alpha=0.2) +
  #annotate(geom="text",x=max(diff_int_ext),y=)
  theme_custom() +
  xlab("SNP Difference INT-EXT")

ggarrange(p.emp.diff.sites,p.sim.diff.sites,
          nrow=2,common.legend=TRUE,legend="right")


p.emp.diff.sites <- results.raxml.all %>%
  filter(sites < 1e6 & height %in% c("Cichlids","Lates")) %>%
  group_by(height,simulation,method,missing,maf) %>%
  pivot_wider(names_from = int, names_prefix="sites",values_from=sites) %>% 
  select(sitesINT,sitesEXT) %>%
  summarize(diff_int_ext = max(sitesINT,na.rm=T)-max(sitesEXT,na.rm=T)) %>%
  ggplot() +
  geom_density(aes(x=diff_int_ext,fill=as.factor(missing)),alpha=0.3) +
  geom_blank(aes(x=-diff_int_ext)) +
  facet_wrap(vars(height),scales="fixed") +
  geom_vline(xintercept=0,lty=2,lwd=0.5,alpha=0.2) +
  #annotate(geom="text",x=max(diff_int_ext),y=)
  theme_custom() +
  theme(legend.title = element_text()) +
  xlab("SNP Difference INT-EXT") +
  labs(fill="Missing Data\nThreshold")
p.sim.diff.sites <- results.raxml.all %>%
  filter(!(height %in% c("Cichlids","Lates"))) %>%
  group_by(height,simulation,method,missing,maf) %>%
  pivot_wider(names_from = int, names_prefix="sites",values_from=sites) %>% 
  select(sitesINT,sitesEXT) %>%
  summarize(diff_int_ext = max(sitesINT,na.rm=T)-max(sitesEXT,na.rm=T)) %>%
  ggplot() +
  geom_density(aes(x=diff_int_ext,fill=as.factor(missing)),alpha=0.3) +
  geom_blank(aes(x=-diff_int_ext)) +
  facet_wrap(vars(height),scales="fixed") +
  geom_vline(xintercept=0,lty=2,lwd=0.5,alpha=0.2) +
  #annotate(geom="text",x=max(diff_int_ext),y=)
  theme_custom() +
  xlab("SNP Difference INT-EXT")

results.raxml.all %>%
  filter(sites < 1e6) %>%
  group_by(height,simulation,method,missing,maf) %>%
  pivot_wider(names_from = int, names_prefix="sites",values_from=sites) %>% 
  select(sitesINT,sitesEXT) %>%
  summarize(diff_int_ext = max(sitesINT,na.rm=T)-max(sitesEXT,na.rm=T)) %>%
  ggplot() +
  geom_jitter(aes(x=maf,y=missing,color=diff_int_ext)) +
  facet_grid(rows=vars(height)) +
  scale_color_viridis()


results.raxml.all %>%
  filter(maf == 0 & sites < 1e7) %>%
  group_by(height,simulation,method,missing,int) %>%
  mutate(sites=sites, 
         diff_sites_maf = case_when(maf == 0 ~ NA,
                                    maf > 0 ~ sites-lag(sites)),
            .groups="keep",order_by=maf) %>%
  ggplot() +
  geom_line(aes(x=missing,y=sites,group=interaction(height,simulation,method,int),color=int),
            alpha=0.3) +
  facet_grid(rows=vars(height),scales="free_y") +
  theme_custom()
  
