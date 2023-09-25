# library(extrafont)
# font_import()
# loadfonts(device = "all")
library(tidyverse)

theme_custom <- function () { 
  theme_bw(base_size=12, base_family="Open Sans") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", color=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="Open Sans"),
      axis.text = element_text(size=14),
      legend.text = element_text(size=14),
      legend.title = element_blank()
    )
}

cols.int <- c("#e2d200","#2A9D8F","#e76f51") # colors for tree heights
cols.sig <- c("black","gray80","black")

loadfonts(device = "postscript")
