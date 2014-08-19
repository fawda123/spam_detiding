######
# load the data and prep

site <- 'sp01'
load(paste0(site, '.RData'))

dat <- get(site)

# names to match those from SWMP for processing
names(dat) <- c('DateTimeStamp', 'Temp', 'Depth', 'Sal', 'DO_mgl', 'Turb',
  'Chl', 'CDOM', 'WSpd', 'Wdir', 'ATemp', 'RH', 'TotPAR', 'BP', 'PAR.t',
  'PAR.b')

# tidal components to predict
tide.comp <- c('M2', 'S2', 'N2', 'K1', 'O1', 'P1', 'Q1', 'MF', 'MM', 'SSA',
  'M4', 'M6', 'S4', 'MS4')

#tidal mod and predictions using all comps
mod.all <- tidem(dat[, 'Depth'], dat[, 'DateTimeStamp'], 
  constituents = tide.comp)
mod.all <- predict(mod.all, newdata = dat[,'DateTimeStamp'])

dat <- data.frame(dat, Tide = mod.all)

# prep for weighted regression
dat <- prep_wtreg(dat)

# save output
nm <- paste0(site, '_prep')
assign(nm, dat)

save(list = nm, file = paste0(nm, '.RData'))
