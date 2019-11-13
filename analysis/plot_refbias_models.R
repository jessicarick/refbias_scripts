## Summary figures for models

# rf distance to true tree
confint.rf.short$ils <- "short"
confint.rf.med$ils <- "med"
confint.rf.long$ils <- "long"

confint.rf <- confint.rf.short %>%
  full_join(confint.rf.med) %>%
  full_join(confint.rf.long)

vars.rf <- ggplot(confint.rf, aes(x = var, y = est, color=sig, shape=ils))
vars.rf.plot <- vars.rf + geom_blank() +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_point(aes(y=est),size=6,alpha=0.7,shape=15) +
  #geom_point(data=confint.rf.med,aes(y=est),size=6,alpha=0.7,shape=16) +
  #geom_point(data=confint.rf.long,aes(y=est),size=6,alpha=0.7,shape=17) +
  geom_pointrange(aes(ymin=minCI,ymax=maxCI),alpha=0.7,size=1.5)+
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.rf.plot

# gamma
confint.gam.short$ils <- "short"
confint.gam.med$ils <- "med"
confint.gam.long$ils <- "long"

confint.gam <- confint.gam.short %>%
  full_join(confint.gam.med) %>%
  full_join(confint.gam.long)

vars.gam <- ggplot(confint.gam, aes(x = var, y = est, color=sig, shape=ils))
vars.gam.plot <- vars.gam + geom_blank() +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_point(aes(y=est),size=6,alpha=0.7,shape=15) +
  geom_pointrange(aes(ymin=minCI,ymax=maxCI),alpha=0.7,size=1.5) +
  #geom_point(data=confint.gam.med,aes(y=est),size=6,alpha=0.7,shape=16) +
  #geom_pointrange(data=confint.gam.med,aes(ymin=minCI,ymax=maxCI),alpha=0.7,shape=16,size=1.5) +
  #geom_point(data=confint.gam.long,aes(y=est),size=6,alpha=0.7,shape=17) +
  #geom_pointrange(data=confint.gam.long,aes(ymin=minCI,ymax=maxCI),alpha=0.7,shape=17,size=1.5) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.gam.plot

# imbalance
confint.imb.short$ils <- "short"
confint.imb.med$ils <- "med"
confint.imb.long$ils <- "long"

confint.imb <- confint.imb.short %>%
  full_join(confint.imb.med) %>%
  full_join(confint.imb.long)

vars.imb <- ggplot(confint.imb.short, aes(x = var, y = est, color=sig))
vars.imb.plot <- vars.imb + geom_blank() +
  xlab("")+
  ylab("Coefficient")+
  scale_color_manual(values=cols)+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
  #geom_point(aes(y=est),size=6,alpha=0.7,shape=15) +
  geom_pointrange(aes(ymin=minCI,ymax=maxCI),alpha=0.7,size=1.5,shape=15) +
  #geom_point(data=confint.imb.med,aes(y=est),size=6,alpha=0.7,shape=16) +
  geom_pointrange(data=confint.imb.med,aes(x=var, y=est,ymin=minCI,ymax=maxCI),alpha=0.7,shape=16,size=1.5) +
  #geom_point(data=confint.imb.long,aes(y=est),size=6,alpha=0.7,shape=17) +
  geom_pointrange(data=confint.imb.long,aes(x=var,y=est,ymin=minCI,ymax=maxCI),alpha=0.7,shape=17,size=1.5) +
  coord_flip() +
  #geom_jitter()+
  theme_minimal() +
  theme(axis.text = element_text(size=15),legend.position="none",
        axis.title = element_text(size=18))
vars.imb.plot
