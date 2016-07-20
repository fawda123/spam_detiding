library(WtRegDO)
# devtools::load_all('M:/docs/swmp/WtRegDO')
library(dplyr)
library(tidyr)
library(ggplot2)

# wqm data
load(file = 'sp02.RData')

# bottle data
bot <- read.csv('M:/docs/SPAM/sp02bottles.csv', header = T) %>% 
  mutate(Date = as.Date(as.character(Date), format = '%m/%d/%Y'))

# format data
sp02 <- select(sp02, DateTime, Temp, Sal, DO, Depth, Atemp, BP, Wsp, PAR.t) %>% 
  rename(
    DateTimeStamp = DateTime, 
    DO_obs = DO,
    ATemp = Atemp,
    WSpd = Wsp, 
    TotPAR = PAR.t
  )

# get ecosystem metabolism

tz <- 'CST6CDT'
lat <- 30.355
long <- -87.202

# setup depth variable, 6 meter before aug 1, 3 m after
depth <- rep(6, nrow(sp02))
depth[sp02$DateTimeStamp >= as.POSIXct('2013-08-01 0:0', tz = tz)] <- 3

# estimate ecosystem metabolism using observed DO time series
# put in areal rates
metab <- ecometab(sp02, DO_var = 'DO_obs', tz = tz,
 lat = lat, long = long, depth_val = depth)


# estimate ecosystem metabolism using observed DO time series
# put in areal rates
metab2 <- ecometab(sp02, DO_var = 'DO_obs', tz = tz,
 lat = lat, long = long, depth_val = 6) 

##
# subset metab ests by dates of bottle measurements
uni_dts <- unique(bot$Date)

tocomb <- vector('list', length = length(uni_dts))
names(tocomb) <- as.character(uni_dts)
wins <- for(dt in uni_dts){
  
  dts <- c(dt - 7, dt + 7)
  sel <- with(metab, Date <= dts[2] & Date >= dts[1])
  sel <- metab[sel, ]
  sel$dtcomb <- dt
  
  tocomb[[as.character(dt)]] <- sel

}
  
comp <- do.call('rbind', tocomb) %>% 
  group_by(dtcomb) %>% 
  summarise(
    Pg_pla = mean(Pg, na.rm = T),
    Rt_pla = mean(Rt, na.rm = T)
  ) %>% 
  ungroup %>% 
  mutate(Date = as.Date(dtcomb, origin = c('1970-01-01'))) %>%
  select(-dtcomb) %>% 
  left_join(., bot, by = 'Date') %>% 
  mutate(
    Pg_bot = Pg * 6, 
    Rt_bot = Rt * 6
  ) %>% 
  select(-Pg, -Rt) %>% 
  gather(.,'var', 'val', -Date) %>% 
  separate(., var, c('est', 'type'), sep = '_')

tocomb2 <- vector('list', length = length(uni_dts))
names(tocomb2) <- as.character(uni_dts)
wins <- for(dt in uni_dts){
  
  dts <- c(dt - 7, dt + 7)
  sel <- with(metab2, Date <= dts[2] & Date >= dts[1])
  sel <- metab2[sel, ]
  sel$dtcomb <- dt
  
  tocomb2[[as.character(dt)]] <- sel

}

comp2 <- do.call('rbind', tocomb2) %>% 
  group_by(dtcomb) %>% 
  summarise(
    Pg_pla = mean(Pg, na.rm = T),
    Rt_pla = mean(Rt, na.rm = T)
  ) %>% 
  ungroup %>% 
  mutate(Date = as.Date(dtcomb, origin = c('1970-01-01'))) %>%
  select(-dtcomb) %>% 
  left_join(., bot, by = 'Date') %>% 
  mutate(
    Pg_bot = Pg * 6, 
    Rt_bot = Rt * 6
  ) %>% 
  select(-Pg, -Rt) %>% 
  gather(.,'var', 'val', -Date) %>% 
  separate(., var, c('est', 'type'), sep = '_')

ylabs <- expression(paste('mmol ', O [2], ' ', m^-2, d^-1))

lims <- c(-700, 700)

p1 <- plot(metab, by = 'days') + scale_y_continuous(ylabs, limits = lims) +
  ggtitle('Depth change') +
  theme_classic()

p2 <- plot(metab2, by = 'days') + scale_y_continuous(ylabs, limits = lims) +
  ggtitle('Depth constant') +
  theme_classic()

p3 <- ggplot(comp, aes(x = Date, y = val, colour = est, linetype = type)) + 
  stat_hline(yintercept = 0, colour = 'grey') +
  geom_line(size = 1) + 
  geom_point(size = 4) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(9, 'Set1')[1:2]) + 
  theme_classic() + 
  scale_y_continuous(ylabs, limits = lims)

p4 <- ggplot(comp2, aes(x = Date, y = val, colour = est, linetype = type)) + 
  stat_hline(yintercept = 0, colour = 'grey') +
  geom_line(size = 1) + 
  geom_point(size = 4) + 
  scale_colour_manual(values = RColorBrewer::brewer.pal(9, 'Set1')[1:2]) + 
  theme_classic() + 
  scale_y_continuous(ylabs, limits = lims)

grid.arrange(p1, p2, p3, p4)


