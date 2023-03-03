## node descendants analysis
rf.test <- InfoRobinsonFoulds(trees.subset[[1]],trees.subset[[71]],normalize=TRUE,reportMatching=TRUE)
VisualizeMatching(Func=RobinsonFouldsMatching,root(raxml.trees[[7]],"sim_0_0_0",resolve.root=TRUE),root(raxml.trees[[210]],"sim_0_0_0",resolve.root=TRUE))

output <- "092321-output"
results.raxml <- read.csv(paste("output/new/",output,"-raxml.csv",sep=""),header=TRUE,row.names=1,sep=",")
raxml.trees <- read.tree(paste("output/new/",gsub("-output","-all-raxml.trees",output),sep=""))
raxml.trees2 <- raxml.trees[results.raxml$tree.num]

# function to match nodes and calculate number of descendants/support for the nodes lost or retained
lostnodes <- function(i) {
  tree1 <- raxml.trees2[[i]]
  tree1.info <- results.raxml[i,]
  orig.info <- results.raxml[results.raxml$simulation == tree1.info$simulation &
                               results.raxml$height == tree1.info$height &
                               results.raxml$int == tree1.info$int &
                               results.raxml$maf == 0 & results.raxml$missing == 0,] %>%
    slice(1)
  orig.tree <- raxml.trees2[[which(results.raxml$simulation == tree1.info$simulation &
                                    results.raxml$height == tree1.info$height &
                                    results.raxml$int == tree1.info$int &
                                    results.raxml$maf == 0 & results.raxml$missing == 0)[1]]]
  sp.tree <- unroot(ml.tree[[which(ml.tree.info$simulation == tree1.info$simulation &
                       ml.tree.info$height == tree1.info$height)]])
  
  
  match <- matchNodes(sp.tree,tree1,method="descendants")
  match  %>%
    as_tibble() %>%
    mutate(lost = case_when(is.na(tr2) ~ TRUE,
                            TRUE ~ FALSE),
           descendants = sapply(tr1, function(x) length(Descendants(orig.tree,x,"tips")[[1]])),
           bs_support = case_when(tr1 == 52 ~ as.integer(100),
                                    TRUE ~ as.integer(orig.tree$node.label)),
           desc_category = case_when(descendants < 3 ~ "veryfew",
                                     descendants > 2 & descendants < 10 ~ "few",
                                     descendants > 40 ~ "many",
                                     TRUE ~ "moderate")) %>%
    # group_by(lost,desc_category) %>%
    # summarize(mean_desc = mean(descendants,na.rm=T),
    #           mean_supp = mean(mean_support,na.rm=T),
    #           #max_desc = max(descendants, na.rm=T),
    #           num = n(),
    #           mean_nodedepth = mean(tr1),
    #           median_nodedepth= median(tr1)) %>%
    add_column(tree_num = tree1.info$tree.num,.before = TRUE)
}

lostnodes.emp <- function(i) {
  tree1 <- raxml.trees2[[i]]
  tree1.info <- results.raxml[i,]
  orig.info <- results.raxml[results.raxml$simulation == tree1.info$simulation &
                               results.raxml$int == tree1.info$int &
                               results.raxml$maf == 0 & results.raxml$missing == 0,] %>%
    slice(1)
  orig.tree <- raxml.trees2[[which(results.raxml$simulation == tree1.info$simulation &
                                     results.raxml$int == tree1.info$int &
                                     results.raxml$maf == 0 & results.raxml$missing == 0)[1]]]
  
  match <- matchNodes(orig.tree,tree1,method="descendants")
  match  %>%
    as_tibble() %>%
    mutate(lost = case_when(is.na(tr2) ~ TRUE,
                            TRUE ~ FALSE),
           descendants = sapply(tr1, function(x) length(Descendants(orig.tree,x,"tips")[[1]])),
           bs_support = case_when(tr1 == 52 ~ as.integer(100),
                                  TRUE ~ as.integer(orig.tree$node.label)),
           desc_category = case_when(descendants < 3 ~ "veryfew",
                                     descendants > 2 & descendants < 10 ~ "few",
                                     descendants > 40 ~ "many",
                                     TRUE ~ "moderate")) %>%
    # group_by(lost,desc_category) %>%
    # summarize(mean_desc = mean(descendants,na.rm=T),
    #           mean_supp = mean(mean_support,na.rm=T),
    #           #max_desc = max(descendants, na.rm=T),
    #           num = n(),
    #           mean_nodedepth = mean(tr1),
    #           median_nodedepth= median(tr1)) %>%
    add_column(tree_num = tree1.info$tree.num,.before = TRUE)
}

# couldn't get it to work with purrr::map_df for some reason
# trees.lostnodes <- results.raxml %>%
#   filter(simulation > 15) %>%
#   select(tree.num) %>%
#   map_df(~lostnodes(.$tree.num))

lostnode.info <- tibble(tree_num=integer(), desc_category=character(),lost = logical(), mean_desc=numeric(), mean_supp=numeric(),
                        num=integer(), mean_nodedepth=numeric(), median_nodedepth=numeric())
lostnode.info <- tibble(tree_num=integer(), tr1=integer(), tr2=integer(),
                        lost=logical(), descendants=integer(), bs_support=numeric(),
                        desc_category=character())
for (i in 1:nrow(results.raxml)){
  sum <- lostnodes(i)
  lostnode.info <- lostnode.info %>%
    add_row(sum)
}

## simulation trees
lostnode.summary <- lostnode.info %>%
  drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num")) %>%
  #filter(maf !=0) %>%
 # mutate(missing=as.factor(missing)) %>%
  group_by(lost,maf,height,int) %>%
  summarize(mean_supp = mean(bs_support,na.rm=T),
            med_supp = median(bs_support,na.rm=T),
            mean_desc = mean(descendants,na.rm=T),
            med_desc = median(descendants,na.rm=T),
            mean_nodedepth = mean(tr1,na.rm=T),
            med_nodedepth = median(tr1,na.rm=T),
            #mean_med_nodedepth = mean(median_nodedepth,na.rm=T),
            #med_med_nodedepth = median(median_nodedepth, na.rm=T)
            ) %>%
  ungroup()
lostnode.info2 <- lostnode.info %>%
  #drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num"))

p1 <- lostnode.info2 %>% 
  #filter(bs_support < 100) %>%
  #group_by(tree_num,maf,lost) %>%
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=bs_support)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2,fun="median") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=mean_supp,col=as.factor(maf),group=height),lty=2) +
  theme_custom() +
  facet_grid(~int,scales="free") +
  #coord_cartesian(ylim=c(99,100)) +
  xlab("Minor Allele Count") +
  ylab("Node Bootstrap Support") +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
p2 <- lostnode.info2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=descendants)) +
  geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  facet_grid(~int) +
  #coord_cartesian(xlim=c(0,100)) +
  ylab("Mean Node Descendants") +
  xlab("Minor Allele Count") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
p2.2 <- lostnode.info2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=tr1)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  facet_grid(~int) +
  #coord_cartesian(xlim=c(0,100)) +
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
ggarrange(p1,p2,p2.2,common.legend=TRUE, ncol=1)

p3 <- lostnode.info2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_supp,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  geom_vline(data=lostnode.summary %>% filter(lost == FALSE),
             aes(xintercept=mean_mean_supp,col=missing,group=height),lty=2) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(50,100)) +
  xlab("Mean Node BS Support") +
  theme(strip.text = element_text(size=rel(1.2)))
p4 <- lostnode.info2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_desc,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(0,100)) +
  xlab("Mean Node Descendants") +
  geom_vline(data=lostnode.summary %>% filter(lost == FALSE),
             aes(xintercept=med_mean_desc,col=missing,group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))
p.density <- ggarrange(p1,p2,p2.2,nrow=1,
          common.legend=TRUE,
          labels="AUTO",
          font.label=list(size=24,font.family="Open Sans"),
          legend="right")

lostnode.info2 %>%
  filter(maf !=0 ) %>%
  ggplot(aes(x=maf,y=descendants)) +
  stat_halfeye(data= . %>% filter(lost==TRUE),aes(fill=height,col=height),
               slab_alpha=0.6,interval_alpha=0,adjust=0.6,trim=TRUE,
               position=position_nudge(x=0.03)) +
  #geom_jitter(data= . %>% filter(lost==TRUE),aes(fill=height,col=height)) +
  # stat_halfeye(data= . %>% filter(lost==FALSE),aes(fill=height,col=height),
  #              slab_alpha=0.6,interval_alpha=0,side="left",adjust=0.6,trim=TRUE,
  #              position=position_nudge(x=-0.03)) +
  theme_custom() +
  coord_flip()

p.bars <- lostnode.info2 %>%
  #filter(height != "LONG") %>%
  filter(maf != 0) %>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  group_by(lost,maf,desc_category) %>%
  summarize(sum_num = n()) %>% 
  #group_by(desc_category) %>% 
  #mutate(total_categ = sum(sum_num)) %>% ungroup() %>% mutate(prop=sum_num/total_categ) %>%
  ggplot(aes(x=maf)) +
  #stat_halfeye(aes(fill=desc_category,col=desc_category),slab_alpha=0.6,interval_alpha=0,adjust=1.2,trim=FALSE) +
  geom_bar(data=. %>% filter(lost == FALSE), aes(y=sum_num, fill=desc_category), position="dodge", stat="identity") +
  geom_bar(data=. %>% filter(lost == TRUE), aes(y=-sum_num, fill=desc_category), position="dodge", stat="identity") +
  geom_hline(yintercept=0) +
  theme_custom() +
  #ylim(-2800,11000) +
  ylab("Number of Nodes") +
  xlab("Minor Alelle Count") +
  annotate(geom="text",x=0.5,y=-4000,label="LOST",hjust=0,size=6,family="Open Sans Light") +
  annotate(geom="text",x=0.5,y=10500,label="RETAINED",hjust=0,size=6, family="Open Sans Light") +
  theme(strip.text = element_text(size=rel(1.2))) +
  #facet_grid(~lost) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8)

p.bars.prior <- lostnode.info2 %>%
  #filter(height != "LONG")%>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  filter(maf == 0 & missing == 0 & lost == FALSE) %>%
  ggplot() +
  geom_histogram(aes(x=descendants,y=..density..,fill=desc_category),bins=50) +
  #geom_histogram(aes(x=mean_supp,y=..density..),fill="black",alpha=0.5,bins=50) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8) +
  xlab("Mean Node Descendants") +
  theme_custom() +
  theme(legend.position="none")

p.bars.lab <- ggarrange(p.bars,p.bars.prior,labels=c("D","E"), 
          common.legend=TRUE,
          font.label=list(size=24,font.family="Open Sans"),
          legend="right",widths=c(3,1)) 

p.bars.lab <- p.bars +
  annotation_custom(
    ggplotGrob(p.bars.prior), 
    xmin = 5.5, xmax = 9, ymin=3000,ymax=11000
  )
  

#library(patchwork)
p.density / p.bars.lab +
  plot_layout(heights=c(1,1.5))
  
## modeling?
lostnode.summary2 <- lostnode.info2 %>%
  group_by(tree_num,lost,maf,missing,height,simulation,int) %>%
  summarize(mean_bs = mean(bs_support,na.rm=T),
         mean_desc = mean(descendants, na.rm=T),
         count = n(),
         mean_depth = mean(tr1))
model1<-lm(mean_bs ~ maf, data=lostnode.summary2)
model2<-lm(mean_desc ~ maf + int, data=lostnode.summary2)
summary(model1)
summary(model2)
r.squaredGLMM(model1)
r.squaredGLMM(model2)
plot(model2)

psum1 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_desc,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +    stat_smooth(method="lm") +
    theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Descendants")
psum2 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_depth,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +
  stat_smooth(method="lm") +
  theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth")
psum3 <- lostnode.summary2 %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=mean_desc,y=mean_bs,col=as.factor(maf))) +
  geom_jitter() +
  theme_custom() 
psum4 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=count,fill=lost)) +
  geom_bar(stat="identity") +
  theme_custom() +
  ylab("Number of Nodes") +
  xlab("Minor Allele Count") +
  scale_fill_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
  
ggarrange(psum1,psum2,psum4,nrow=1,common.legend=TRUE)

lostnode.summary2 %>%
  drop_na() %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=as.integer(maf),y=mean_desc,col=height)) +
  geom_jitter(width=0.001) +
  stat_summary(aes(x=maf,y=mean_desc,group=int),col="black",geom="line",lty=2) +
  stat_summary(aes(x=maf,y=mean_desc,group=int),col="black",geom="point",size=5,pch=4) +
  #geom_point(aes(x=maf,y=mean_retained_bs),col="black")+
  #facet_grid(vars(int)) +
  theme_custom() +
  ylab("Mean descendants of lost nodes") +
  xlab("Minor Allele Count") +
  theme(legend.position="none")

##########################
## empirical trees
##########################
# lostnode.info <- tibble(tree_num=integer(), desc_category=character(),lost = logical(), mean_desc=numeric(), mean_supp=numeric(),
#                         num=integer(), mean_nodedepth=numeric(), median_nodedepth=numeric())
lostnode.info.emp <- tibble(tree_num=integer(), tr1=integer(), tr2=integer(),
                        lost=logical(), descendants=integer(), bs_support=numeric(),
                        desc_category=character())
for (i in 1:nrow(results.raxml)){
  sum <- lostnodes.emp(i)
  lostnode.info.emp <- lostnode.info.emp %>%
    add_row(sum)
}

lostnode.summary.emp <- lostnode.info.emp %>%
  drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num")) %>%
  #filter(maf !=0) %>%
  # mutate(missing=as.factor(missing)) %>%
  group_by(lost,maf) %>%
  summarize(mean_supp = mean(bs_support,na.rm=T),
            med_supp = median(bs_support,na.rm=T),
            mean_desc = mean(descendants,na.rm=T),
            med_desc = median(descendants,na.rm=T),
            #mean_mean_nodedepth = mean(mean_nodedepth,na.rm=T),
            #med_mean_nodedepth = median(mean_nodedepth,na.rm=T),
            #mean_med_nodedepth = mean(median_nodedepth,na.rm=T),
            #med_med_nodedepth = median(median_nodedepth, na.rm=T)
  ) %>%
  ungroup()
lostnode.info.emp2 <- lostnode.info.emp %>%
  #drop_na() %>%
  left_join(results.raxml, by=c("tree_num" = "tree.num"))

p1.emp <- lostnode.info.emp2 %>% 
  #group_by(tree_num,maf,lost) %>%
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=bs_support)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=mean_supp,col=as.factor(maf),group=height),lty=2) +
  theme_custom() +
  #facet_grid(~height,scales="free") +
  #coord_cartesian(ylim=c(99,100)) +
  xlab("Minor Allele Count") +
  ylab("Node Bootstrap Support") +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
p2.emp <- lostnode.info.emp2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=descendants)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  #facet_grid(~height) +
  #coord_cartesian(xlim=c(0,100)) +
  ylab("Mean Node Descendants") +
  xlab("Minor Allele Count") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
p2.2.emp <- lostnode.info.emp2 %>%  
  #mutate(missing=as.factor(missing)) %>%
  #filter(lost == TRUE) %>% 
  ggplot(aes(x=maf,y=tr1)) +
  #geom_jitter(width=0,alpha=0.2,aes(color=lost)) +
  stat_summary(aes(color=lost),alpha=0.8,size=2) +
  #geom_density(alpha=0.5) +
  theme_custom() +
  #facet_grid(~height) +
  #coord_cartesian(xlim=c(0,100)) +
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth") +
  # geom_vline(data=lostnode.summary %>% filter(lost == TRUE),
  #            aes(xintercept=med_mean_desc,col=as.factor(maf),group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))
ggarrange(p1.emp,p2.emp,p2.2.emp,common.legend=TRUE, nrow=1)

p3.emp <- lostnode.info.emp2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_supp,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  geom_vline(data=lostnode.summary.emp %>% filter(lost == FALSE),
             aes(xintercept=mean_mean_supp,col=missing,group=height),lty=2) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(50,100)) +
  xlab("Mean Node BS Support") +
  theme(strip.text = element_text(size=rel(1.2)))
p4.emp <- lostnode.info.emp2 %>%  
  mutate(missing=as.factor(missing)) %>%
  filter(lost == FALSE) %>% 
  ggplot(aes(x=mean_desc,fill=missing,col=missing)) +
  geom_density(alpha=0.5) +
  theme_custom() +
  facet_grid(~height) +
  coord_cartesian(xlim=c(0,100)) +
  xlab("Mean Node Descendants") +
  geom_vline(data=lostnode.summary.emp %>% filter(lost == FALSE),
             aes(xintercept=med_mean_desc,col=missing,group=height),lty=2) +
  theme(strip.text = element_text(size=rel(1.2)))
p.density.emp <- ggarrange(p1.emp,p2.emp,p2.2.emp,nrow=1,
                       common.legend=TRUE,
                       labels="AUTO",
                       font.label=list(size=24,font.family="Open Sans"),
                       legend="right")

lostnode.info.emp2 %>%
  filter(maf !=0 ) %>%
  ggplot(aes(x=maf,y=descendants)) +
  stat_halfeye(data= . %>% filter(lost==TRUE),aes(fill=height,col=height),
               slab_alpha=0.6,interval_alpha=0,adjust=0.6,trim=TRUE,
               position=position_nudge(x=0.03)) +
  #geom_jitter(data= . %>% filter(lost==TRUE),aes(fill=height,col=height)) +
  # stat_halfeye(data= . %>% filter(lost==FALSE),aes(fill=height,col=height),
  #              slab_alpha=0.6,interval_alpha=0,side="left",adjust=0.6,trim=TRUE,
  #              position=position_nudge(x=-0.03)) +
  theme_custom() +
  coord_flip()

p.bars.emp <- lostnode.info.emp2 %>%
  #filter(height != "LONG") %>%
  filter(maf != 0) %>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  group_by(lost,maf,desc_category) %>%
  summarize(sum_num = n()) %>% 
  #group_by(desc_category) %>% 
  #mutate(total_categ = sum(sum_num)) %>% ungroup() %>% mutate(prop=sum_num/total_categ) %>%
  ggplot(aes(x=maf)) +
  #stat_halfeye(aes(fill=desc_category,col=desc_category),slab_alpha=0.6,interval_alpha=0,adjust=1.2,trim=FALSE) +
  geom_bar(data=. %>% filter(lost == FALSE), aes(y=sum_num, fill=desc_category), position="dodge", stat="identity") +
  geom_bar(data=. %>% filter(lost == TRUE), aes(y=-sum_num, fill=desc_category), position="dodge", stat="identity") +
  geom_hline(yintercept=0) +
  theme_custom() +
  #ylim(-2800,11000) +
  ylab("Number of Nodes") +
  xlab("Minor Alelle Count") +
  annotate(geom="text",x=0.5,y=-1300,label="LOST",hjust=0,size=6,family="Open Sans Light") +
  annotate(geom="text",x=0.5,y=1200,label="RETAINED",hjust=0,size=6, family="Open Sans Light") +
  theme(strip.text = element_text(size=rel(1.2))) +
  #facet_grid(~lost) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8)

p.bars.prior.emp <- lostnode.info.emp2 %>%
  #filter(height != "LONG")%>%
  mutate(desc_category=factor(desc_category,levels=c("veryfew","few","moderate","many"))) %>%
  filter(maf == 0 & missing == 0 & lost == FALSE) %>%
  ggplot() +
  geom_histogram(aes(x=descendants,y=..density..,fill=desc_category),bins=50) +
  #geom_histogram(aes(x=mean_supp,y=..density..),fill="black",alpha=0.5,bins=50) +
  scale_fill_viridis(discrete=TRUE,alpha=0.8) +
  xlab("Mean Node Descendants") +
  theme_custom() +
  theme(legend.position="none")

p.bars.lab.emp <- ggarrange(p.bars.emp,p.bars.prior.emp,labels=c("D","E"), 
                        common.legend=TRUE,
                        font.label=list(size=24,font.family="Open Sans"),
                        legend="right",widths=c(3,1)) 

p.bars.lab.emp <- p.bars +
  annotation_custom(
    ggplotGrob(p.bars.prior), 
    xmin = 5.5, xmax = 9, ymin=3000,ymax=11000
  )


#library(patchwork)
p.density / p.bars.lab +
  plot_layout(heights=c(1,1.5))

## modeling?
lostnode.summary2 <- lostnode.info2 %>%
  group_by(tree_num,lost,maf,missing,height,simulation,int) %>%
  summarize(mean_bs = mean(bs_support,na.rm=T),
            mean_desc = mean(descendants, na.rm=T),
            count = n(),
            mean_depth = mean(tr1))
model1<-lm(mean_bs ~ maf, data=lostnode.summary2)
model2<-lm(mean_desc ~ maf + int, data=lostnode.summary2)
summary(model1)
summary(model2)
r.squaredGLMM(model1)
r.squaredGLMM(model2)
plot(model2)

psum1 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_desc,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +    stat_smooth(method="lm") +
  theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Descendants")
psum2 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=mean_depth,col=lost)) +
  geom_point(data=. %>% filter(lost == TRUE), alpha=0.6, position=position_nudge(x=0.1)) +
  geom_point(data=. %>% filter(lost == FALSE), alpha=0.6, position=position_nudge(x=-0.1)) +
  stat_smooth(method="lm") +
  theme_custom()+
  scale_color_manual(labels=c("Retained","Lost"),
                     values=c("#0A9396","#9B2226"))+
  xlab("Minor Allele Count") +
  ylab("Mean Node Depth")
psum3 <- lostnode.summary2 %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=mean_desc,y=mean_bs,col=as.factor(maf))) +
  geom_jitter() +
  theme_custom() 
psum4 <- lostnode.summary2 %>%
  ggplot(aes(x=maf,y=count,fill=lost)) +
  geom_bar(stat="identity") +
  theme_custom() +
  ylab("Number of Nodes") +
  xlab("Minor Allele Count") +
  scale_fill_manual(labels=c("Retained","Lost"),
                    values=c("#0A9396","#9B2226"))

ggarrange(psum1,psum2,psum4,nrow=1,common.legend=TRUE)

lostnode.summary2 %>%
  drop_na() %>%
  filter(lost == TRUE) %>%
  ggplot(aes(x=as.integer(maf),y=mean_desc,col=height)) +
  geom_jitter(width=0.001) +
  stat_summary(aes(x=maf,y=mean_desc),col="black",geom="line",lty=2) +
  stat_summary(aes(x=maf,y=mean_desc),col="black",geom="point",size=5,pch=4) +
  #geom_point(aes(x=maf,y=mean_retained_bs),col="black")+
  facet_grid(vars(int)) +
  theme_custom() +
  ylab("Mean descendants of lost nodes") +
  xlab("Minor Allele Count") +
  theme(legend.position="none")
