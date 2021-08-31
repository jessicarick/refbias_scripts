require(MuMIn)
require(brms)
require(lme4)
require(tidybayes)

cols <- c("#F2AD00","gray80","#00A08A")

## preparing data object
results.mod <- results.raxml
results.mod$simulation <- as.factor(results.mod$simulation)
results.mod$height <- as.factor(results.mod$height)
results.mod$quality <- as.factor(results.mod$quality)
results.mod$missing <- as.factor(results.mod$missing)
results.mod$maf <- as.factor(results.mod$maf)
results.mod$int <- as.factor(results.mod$int)

for (param in c("gamma","colless","sackin","ingroup.gamma","ingroup.colless","ingroup.sackin")) {
    formula <- paste0(param," ~ int + maf + missing + int:maf + int:missing + (1|simulation)")
    m <- lmer(formula,
                 data = results.mod)
    m.sum <- summary(m)
    m.brms <- brms::brm(formula,
                           data = results.mod,
                           iter = 55000,
                           chains = 3,
                           cores = 3,
                           thin = 5,
                           warmup = 5000,
                           control = list(adapt_delta = 0.99, max_treedepth = 15))
    m.brms.sum <- summary(m.brms)
    #plot(m.brms)
    #brms::marginal_effects(m.brms)
    
    plot.brms <- as_tibble(m.brms.sum$fixed,rownames="var")%>%
      rename(minCI = `l-95% CI`, maxCI = `u-95% CI`) %>%
      mutate(sig = case_when(minCI > 0 ~ "pos",
                             maxCI < 0 ~ "neg",
                             TRUE ~ "ns")) %>%
      filter(var != "Intercept") %>%
      ggplot(aes(x = var, y = Estimate, col=sig)) +
      xlab("")+
      ylab(paste0("Effect Size, ",param)) +
      #scale_color_npg() +
      #scale_x_reverse() +
      scale_color_manual(values=cols)+
      geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
      #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
      geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1,shape=16) +
      coord_flip() +
      theme_custom() +
      theme(axis.text = element_text(size=15),legend.position="none",
            axis.title = element_text(size=18))
      
    assign(paste0("plot.",param,".brms.lates"),plot.brms)
    assign(paste0("brms.",param,".mod.lates"),m.brms)
    assign(paste0("lmer.",param,".mod.lates"),m)
  }

posterior <- as.array(brms.RF.Dist.ML.mod.LONG)
dim(posterior)
dimnames(posterior)

#bayesplot::mcmc_pairs(posterior,pars=dimnames(posterior)$parameters[2:10])

pdf("brms_plots_082021.pdf")
for (mod in mget(ls(pattern="brms.*.mod*"))) {
  p <- mod %>%
    #recover_types(results.sub) %>%
    tidybayes::spread_draws(b_intINT,b_maf0.01,b_maf0.02,b_maf0.03,b_maf0.04,b_maf0.05,b_maf0.1,
                            `b_intINT:maf0.01`,`b_intINT:maf0.02`,`b_intINT:maf0.03`,`b_intINT:maf0.04`,`b_intINT:maf0.05`,`b_intINT:maf0.1`) %>%
    pivot_longer(cols=starts_with("b_"),names_to="variable") %>%
    ggplot(aes(y = variable, x = value)) +
    stat_slab(position=position_nudge(y=0.1),height=0.8,aes(fill=stat(x))) +
    geom_boxplot(position=position_nudge(y=-0.1),width=0.2,outlier.shape=NA) +
    geom_vline(xintercept = 0, lty=2, alpha=0.5) +
    scale_fill_viridis_c(option = "C") +
    theme(legend.position="none") +
    ggtitle(colnames(mod$data)[1]) 
  print(p)
}
dev.off()

model_refbias_brms <- function(data,height,param,iter=12000,burnin=2000,chains=3) {
  results.sub <- data[data$height == height,]
  formula <- paste0(param," ~ avg_dxy + maf + missing + avg_dxy:maf + avg_dxy:missing + (1|simulation)")
  # m <- lmer(formula,
  #           data = results.sub)
  # m.sum <- summary(m)
  
  m.brms <- brms::brm(formula,
                      data = results.sub,
                      iter = iter,
                      chains = chains,
                      cores = chains,
                      warmup = burnin,
                      control = list(adapt_delta = 0.99))
  m.brms.sum <- summary(m.brms)
  #plot(m.gam.med.brms)
  #brms::marginal_effects(m.gam.med.brms)
  
  plot.brms <- as_tibble(m.brms.sum$fixed,rownames="var")%>%
    rename(minCI = `l-95% CI`, maxCI = `u-95% CI`) %>%
    mutate(sig = case_when(minCI > 0 ~ "pos",
                           maxCI < 0 ~ "neg",
                           TRUE ~ "ns")) %>%
    filter(var != "Intercept") %>%
    ggplot(aes(x = var, y = Estimate, col=sig)) +
    xlab("")+
    ylab(paste0("Coefficient\n",height," Trees")) +
    #scale_color_npg() +
    #scale_x_reverse() +
    scale_color_manual(values=cols)+
    geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
    #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
    geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1,shape=16) +
    coord_flip() +
    theme_custom() +
    theme(axis.text = element_text(size=15),legend.position="none",
          axis.title = element_text(size=18))
  print(plot.brms)
  
  results <- list(m.brms,m.brms.sum)
  names(results) <- c(paste0("m.brms.",height,".","param_model"),paste0("m.brms.",height,".","param_summary"))
  return(results)
}

plot_refbias_brms <- function(model) {
  plot <- m.brms %>%
    tidybayes::gather_draws(b_Intercept,b_avg_dxy,b_maf,b_missing) %>%
    #pivot_longer(cols=starts_with("b_"),names_to="variable") %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_slab(position=position_nudge(y=0.1),height=0.8) +
    geom_boxplot(position=position_nudge(y=-0.1),width=0.2,outlier.shape=NA) +
    geom_vline(xintercept = 0, lty=2, alpha=0.5) +
    theme_custom()
  return(plot)
}

##############
d=dredge(m.gam.med, beta="sd")
print(d)
coefficients(d)

sum <- summary(model.avg(d))

confint.avg.med <- as_tibble(confint(model.avg(d))[-c(1),],rownames="var") %>%
  rename(minCI = `2.5 %`, maxCI = `97.5 %`) %>%
  mutate(est = t(sum$coefficients)[-1,1],
         pval = sum$coefmat.full[-1,5],
         sig = case_when(minCI > 0 ~ "pos",
                         maxCI < 0 ~ "neg",
                         TRUE ~ "ns")) 
  


vars.gam.med <- confint.avg.med %>%
  arrange(pval) %>%
  filter(pval < 0.2) %>%
  ggplot(aes(x = var, y = est, color=sig))
vars.gam.med.bars <- vars.gam.med + geom_blank() +
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
  ylab("Coefficient\nMED Trees (Med ILS)")+
  #scale_color_npg() +
  #scale_x_reverse() +
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_linerange(aes(ymin = minCI, ymax = maxCI),lwd=7) +
  geom_pointrange(aes(ymin = minCI, ymax = maxCI),fatten=4,lwd=1,shape=16) +
  coord_flip() +
  theme_custom() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18)) 
vars.gam.med.bars


as_tibble(sum$coefmat.full, col_names=c("Estimate","StdErr","AdjSE","zvalue","pvalue")) %>%
  mutate(Parameter = row.names(sum$coefmat.full)) %>%
  filter(Parameter != "(Intercept)") %>%
  ggplot(aes(x=Parameter,y=Estimate)) +
    geom_point() +
    coord_flip() +
    theme_custom()
