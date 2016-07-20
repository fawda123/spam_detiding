######
# checking metabolism estimates, July 2016
# uses 'spam1_input.RData' from Mike

# # install dev version of SWMPr
# install.packages('devtools')
# devtools::install_github('fawda123/SWMPr', ref = 'development')

library(SWMPr)
library(dplyr)
library(tidyr)
library(ggplot2)

load(file = 'spam1_input.RData')

# metadata
tz <- 'America/Belize'
lat <- 30.355
long <- -87.202  

# format sites for ecometab
sp01 <- select(sp01, datetime, do, depth, atemp, sal, temp, wsp, bp) %>% 
  rename(
    datetimestamp = datetime, 
    do_mgl = do, 
    wspd = wsp
  )
sp02 <- select(sp02, datetime, do, depth, atemp, sal, temp, wsp, bp) %>% 
  rename(
    datetimestamp = datetime, 
    do_mgl = do, 
    wspd = wsp
  )

# ecometab
met_sp01 <- ecometab(sp01, tz = tz, lat = lat, long = long, depth_val = 1.5)
met_sp02 <- ecometab(sp02, tz = tz, lat = lat, long = long, depth_val = 5.9)

save(met_sp01, met_sp02, file = 'spam1_output.RData')

# a plot
to_plo <- rbind(
  data.frame(met_sp01, site = 'sp01'), 
  data.frame(met_sp02, site = 'sp02')
  ) %>% 
  select(date, Pg, Rt, NEM, site) %>% 
  gather('var', 'val', Pg:NEM)

pdf('spam1_output.pdf', height = 6, width = 9, family = 'serif')
ggplot(to_plo, aes(x = date, y = val, group = var, colour = var)) + 
  geom_line() + 
  facet_wrap(~site) + 
  theme_minimal()
dev.off()

  
