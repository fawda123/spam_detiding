######
# created Aug. 18, 2014, M. Beck

######
# load prepped data for wtreg

site <- 'sp01'

load(paste0(site, '_prep.RData'))

dat <- get(paste0(site, '_prep'))

######
# detide DO signal

cl <- makeCluster(8)
registerDoParallel(cl)

dtd <- wtreg_fun(dat, wins = c(8, 12, 1), parallel = T, progress = T)

stopCluster(cl)

to_plo <- dtd[2000:3000,]
ggpoly <- poly.fun(to_plo$solar, to_plo)

ylab<-expression(paste('DO (mg ',L^-1,')'))
p2 <- ggplot(to_plo, aes(x = DateTimeStamp)) + 
  ggpoly +
  geom_line(aes(y = DO_obs, colour = 'Observed')) +
  geom_line(aes(y = DO_prd, colour = 'Predicted')) +
  coord_cartesian(ylim = rng.fun(to_plo$DO_obs)) +
  scale_fill_manual(values='orange',labels='Day') +
  theme_bw() +
  scale_y_continuous(ylab)  +
  my_theme
p2
    
to_nem <- dtd[, !names(dtd) %in% c('dTide', 'met.date', 'solar', 'value', 'day.hrs')]
to_nem$DO_nrm2 <- with(to_nem, DO_obs - DO_prd + DO_nrm)

nem_obs <- nem.fun(to_nem, stat = site, DO_var = 'DO_obs')
nem_obs <- nem_obs[, c('Date', 'Pg', 'Rt', 'NEM')]
nem_obs$var <- 'Observed'
nem_dtd <- nem.fun(to_nem, stat = site,  DO_var = 'DO_nrm')
nem_dtd <- nem_dtd[, c('Date', 'Pg', 'Rt', 'NEM')]
nem_dtd$var <- 'Detided'
nem_dtd2 <- nem.fun(to_nem, stat = site,  DO_var = 'DO_nrm2')
nem_dtd2 <- nem_dtd2[, c('Date', 'Pg', 'Rt', 'NEM')]
nem_dtd2$var <- 'Detided2'

to_plo <- rbind(nem_obs, nem_dtd, nem_dtd2)
to_plo <- melt(to_plo, id.var = c('Date', 'var'), 
  measure.var = c('Pg', 'Rt', 'NEM'))

p1 <- ggplot(to_plo, aes(x = Date, y = 0.032 * value, group = variable, 
  colour = variable)) +
  geom_line() + 
  geom_point() +
  facet_wrap(~ var, ncol = 1)
  
p1

to_plo <- dtd
ggpoly <- poly.fun(to_plo$solar, to_plo)

ylab<-expression(paste('DO (mg ',L^-1,')'))
p2 <- ggplot(to_plo, aes(x = DateTimeStamp)) + 
  ggpoly +
  geom_line(aes(y = DO_obs, colour = 'Observed')) +
  geom_line(aes(y = DO_prd, colour = 'Predicted')) +
  coord_cartesian(ylim = rng.fun(to_plo$DO_obs)) +
  scale_fill_manual(values='orange',labels='Day') +
  theme_bw() +
  scale_y_continuous(ylab)  +
  my_theme
p2

# some diagnostic plots
names(dat)[names(dat) %in% 'variable'] <- 'solar'

to_plo <- dat
ggpoly <- poly.fun(to_plo$solar, to_plo)

ylab<-expression(paste('DO (mg ',L^-1,')'))
p3 <- ggplot(to_plo, aes(x = DateTimeStamp)) + 
  ggpoly +
  geom_line(aes(y = DO_obs, colour = 'Observed')) +
  coord_cartesian(ylim = rng.fun(to_plo$DO_obs)) +
  scale_fill_manual(values='orange',labels='Day') +
  theme_bw() +
  scale_y_continuous(ylab)  +
  my_theme

ylab<-expression(paste('Tide (m)'))
p4 <- ggplot(to_plo, aes(x = DateTimeStamp)) + 
  ggpoly +
  geom_line(aes(y = Tide), size = 1.2) +
  coord_cartesian(ylim = rng.fun(to_plo$Tide)) +
  scale_fill_manual('', values='orange',labels='Day') +
  theme_bw() +  
  theme(
    legend.position = 'top',
    axis.title.x = element_blank(),legend.box= 'horizontal',
    plot.margin= unit(c(0, 1, 0, 1), "lines"),
#         axis.text = element_text(size = 14), 
    text = element_text(size = 16)
    ) +
  scale_y_continuous(ylab)

grid.arrange(p3, p4, ncol = 1)

    
