# README
Marcus W. Beck, beck.marcus@epa.gov  


```r
library(WtRegDO)
# devtools::load_all('M:/docs/SWMP/WtRegDO')
library(tidyverse)
library(doParallel)
```

Prep data for weighted regression:

```r
# load, prep data
load('data/sp01.RData')
load('data/sp02.RData')

sp01 <- sp01 %>% 
  mutate(site = 'sp01')
sp02 <- sp02 %>% 
  mutate(site = 'sp02')

# format for wtregdo
spdat <- rbind(sp01, sp02) %>% 
  select(site, DateTime, Temp, Sal, DO, Atemp, BP, Wsp, Depth) %>%
  rename(
    DateTimeStamp = DateTime, 
    DO_obs = DO, 
    ATemp = Atemp, 
    WSpd = Wsp,
    Tide = Depth
  ) %>% 
  mutate(
    DateTimeStamp = as.POSIXct(DateTimeStamp, tz = 'America/Regina')
  ) %>% 
  filter(!is.na(Tide))

sp01 <- filter(spdat, site %in% 'sp01') %>% 
  data.frame %>% 
  select(-site)
sp02 <- filter(spdat, site %in% 'sp02') %>% 
  data.frame %>% 
  select(-site)
```

Check correlation of tidal height with sun angle:

```r
# run weighted regression in parallel
# requires parallel backend
registerDoParallel(cores = 7)

# metadata for the location
tz <- 'America/Regina'
lat <- 30.36
long <- -87.20

# check correlations
sp01eval <- evalcor(sp01, tz, lat, long, progress = TRUE, plot = FALSE)
sp02eval <- evalcor(sp02, tz, lat, long, progress = TRUE, plot = FALSE)

# format for plot
sp01cor <- data.frame(DateTimeStamp = sp01$DateTimeStamp, Tide = sp01$Tide, corval = sp01eval, site = 'sp01')
sp02cor <- data.frame(DateTimeStamp = sp02$DateTimeStamp, Tide = sp02$Tide, corval = sp02eval, site = 'sp02')
corval <- rbind(sp01cor, sp02cor) %>%
  gather('var', 'val', Tide:corval)
save(corval, file = 'corval.RData', compress = 'xz')
```

Evaluate corelation results:

```r
load(file = 'data/corval.RData')

ggplot(corval, aes(x = DateTimeStamp, y = val, colour = site)) +
  geom_line() +
  facet_wrap(~var, ncol = 1, scales = 'free_y') +
  theme_bw()
```

![](README_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

Run weighted regression for each site, estimate metabolism:

```r
# run weighted regression in parallel
# requires parallel backend
registerDoParallel(cores = 7)

# metadata for the location
tz <- 'America/Regina'
lat <- 30.36
long <- -87.20

# weighted regression
sp01wtreg <- wtreg(sp01, parallel = F, wins = list(3, 1, 0.6), progress = TRUE,
  tz = tz, lat = lat, long = long)
sp02wtreg <- wtreg(sp02, parallel = TRUE, wins = list(3, 1, 0.6), progress = TRUE,
  tz = tz, lat = lat, long = long)

save(sp01wtreg, file = 'data/sp01wtreg.RData', compress = 'xz')
save(sp02wtreg, file = 'data/sp02wtreg.RData', compress = 'xz')

# estimate ecosystem metabolism using observed and detided DO time series
sp01metobs <- ecometab(sp01wtreg, DO_var = 'DO_obs', tz = tz,
  lat = lat, long = long)
sp01metdtd <- ecometab(sp01wtreg, DO_var = 'DO_nrm', tz = tz,
  lat = lat, long = long)
sp02metobs <- ecometab(sp02wtreg, DO_var = 'DO_obs', tz = tz,
  lat = lat, long = long)
sp02metdtd <- ecometab(sp02wtreg, DO_var = 'DO_nrm', tz = tz,
  lat = lat, long = long)

save(sp01metobs, sp01metdtd, sp02metobs, sp02metdtd, file = 'data/metabs.RData', compress = 'xz')
```

Check results:

```r
load('data/sp01wtreg.RData')
load('data/sp02wtreg.RData')

# Observed and detided DO
sp01wtreg <- mutate(sp01wtreg, site = 'sp01')
sp02wtreg <- mutate(sp02wtreg, site = 'sp02')
toplo1 <- rbind(sp01wtreg, sp02wtreg) %>% 
  select(site, DateTimeStamp, DO_obs, DO_nrm) %>% 
  gather('var', 'val', DO_obs:DO_nrm)

ggplot(toplo1, aes(x = DateTimeStamp, y = val, colour = var)) +
  geom_line() +
  facet_wrap(~site, ncol = 1, scales = 'free_y') +
  theme_bw() +
  scale_y_continuous('DO_mgl')
```

![](README_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
# observed and detided metabolism
load(file = 'data/metabs.RData')

plot(sp01metobs, by = 'days') + 
  ggtitle('sp01, observed DO')
```

![](README_files/figure-html/unnamed-chunk-6-2.png)<!-- -->

```r
plot(sp01metdtd, by = 'days') + 
  ggtitle('sp01, detided DO')
```

![](README_files/figure-html/unnamed-chunk-6-3.png)<!-- -->

```r
plot(sp02metobs, by = 'days') + 
  ggtitle('sp02, observed DO')
```

![](README_files/figure-html/unnamed-chunk-6-4.png)<!-- -->

```r
plot(sp02metdtd, by = 'days') + 
  ggtitle('sp02, detided DO')
```

![](README_files/figure-html/unnamed-chunk-6-5.png)<!-- -->

```r
tokp <- c('meanPg', 'sdPg', 'anomPg', 'meanRt', 'sdRt', 'anomRt')
metcomp <- list(
  sp01metobs = meteval(sp01metobs)[tokp],
  sp01metdtd = meteval(sp01metdtd)[tokp],
  sp02metobs = meteval(sp02metobs)[tokp],
  sp02metdtd = meteval(sp02metdtd)[tokp]
  ) %>% 
  do.call('cbind', .)
metcomp
```

```
##        sp01metobs sp01metdtd sp02metobs sp02metdtd
## meanPg 187.036    195.1479   51.11921   49.01884  
## sdPg   112.9164   70.38351   39.98914   39.60008  
## anomPg 1.273885   0.6369427  3.821656   8.280255  
## meanRt -186.8889  -192.2115  -51.25198  -48.08603 
## sdRt   117.9725   66.18281   34.66544   36.2747   
## anomRt 1.273885   0          2.547771   8.280255
```
