#functions for NEM processing of NERRS data
#created Dec. 2013 by M. Beck, adapted from 'spam_NEM_fun.r' and M. Murrell

#funcion that splits dataset into 24hr days based on sunrise
#merge with original data
met.day.fun <- function(dat.in){

  require(StreamMetabolism)  #for sunrise.set function
  
  #get sunrise/sunset times using sunrise.st function from StreamMetabolism
  lat<-30.355
  long<--87.202
  tz<-'America/Chicago' #central standard
  start.day<-format(dat.in$DateTimeStamp[which.min(dat.in$DateTimeStamp)]-(60*60*24),format='%Y/%m/%d')
  tot.days<-1+length(unique(as.Date(dat.in$DateTimeStamp)))
  
  #ss.dat is matrix of sunrise/set times for each days  within period of obs
  ss.dat<-suppressWarnings(sunrise.set(lat,long,start.day,tz,tot.days))
  
  #remove duplicates, sometimes sunrise.set screws up
  ss.dat<-ss.dat[!duplicated(strftime(ss.dat[,1],format='%Y-%m_%d')),]
  ss.dat<-data.frame(
    ss.dat,
    met.date=as.Date(ss.dat$sunrise,tz=tz)
    )
  ss.dat<-melt(ss.dat,id.vars='met.date')
  if(!"POSIXct" %in% class(ss.dat$value))
    ss.dat$value<-as.POSIXct(ss.dat$value, origin='1970-01-01',tz=tz)
  ss.dat<-ss.dat[order(ss.dat$value),]
  ss.dat$day.hrs<-unlist(lapply(
    split(ss.dat,ss.dat$met.date),
    function(x) rep(as.numeric(x[2,'value']-x[1,'value']),2) 
    ))
  
  #matches is vector of row numbers indicating starting value that each
  #unique DateTimeStamp is within in ss.dat
  #output is meteorological day matches appended to dat.in
  matches<-findInterval(dat.in$DateTimeStamp,ss.dat$value)
  data.frame(dat.in,ss.dat[matches,])
      
  }

#calculates oxygen mass transfer coefficient, from Thiebault et al. 2008
#output from this can be used  to get volumetric rearation coefficient
#input is water temp, salinity, air temp, wind speed, barometric press, height of anemometer
f_calcKL<-function(Temp,Sal,ATemp,WSpd,BP,Height=10){

  require(oce) #for swSigmaT
  
  #celsius to kelvin conversion
  CtoK<-function(val) val+273.15 
  sig.fun<-Vectorize(swSigmaT)
  
  to.vect<-function(Temp,Sal,ATemp,WSpd,BP,Height=10){
    
    Patm<-BP*100; # convert from millibars to Pascals
    zo<-1e-5; # assumed surface roughness length (m) for smooth water surface
    U10<-WSpd*log(10/zo)/log(Height/zo)
    TempK<-CtoK(Temp)
    ATempK<-CtoK(ATemp)
    sigT<-sig.fun(Sal,Temp,10) # set for 10 decibars = 1000mbar = 1 bar = 1atm
    rho_w<-1000+sigT #density of SW (kg m-3)
    Upw<-1.002e-3*10^((1.1709*(20-Temp)-(1.827*10^-3*(Temp-20)^2))/(Temp+89.93)) #dynamic viscosity of pure water (Sal=0);
    Uw<-Upw*(1+(5.185e-5*Temp+1.0675e-4)*(rho_w*Sal/1806.55)^0.5+(3.3e-5*Temp+2.591e-3)*(rho_w*Sal/1806.55))  # dynamic viscosity of SW
    Vw<-Uw/rho_w  #kinematic viscosity
    Ew<-6.112*exp(17.65*ATemp/(243.12+ATemp))  # water vapor pressure (hectoPascals)
    Pv<-Ew*100 # Water vapor pressure in Pascals
    Rd<-287.05  # gas constant for dry air ( kg-1 K-1)
    Rv<-461.495  # gas constant for water vapor ( kg-1 K-1)
    rho_a<-(Patm-Pv)/(Rd*ATempK) +Pv/(Rv*TempK)
    kB<-1.3806503e-23 # Boltzman constant (m2 kg s-2 K-1)
    Ro<-1.72e-10     #radius of the O2 molecule (m)
    Dw<-kB*TempK/(4*pi*Uw*Ro)  #diffusivity of O2 in water 
    KL<-0.24*170.6*(Dw/Vw)^0.5*(rho_a/rho_w)^0.5*U10^1.81  #mass xfer coef (m d-1)
   
    return(KL)
    
    }
  
  out.fun<-Vectorize(to.vect)
  
  out.fun(Temp,Sal,ATemp,WSpd,BP,Height=10)
  
  }

######
# NEM function, uses all functions above
# estimates daily integrated rates, gross production, total respiraiton
# 'dat.in' input is station data frame 
# 'stat' is character string for station, five letters
# 'DO_var' is chr string of column for DO 
# 'depth.val' is value for station depth if needed to add manually
# 'bott.stat' is logical indicating if station is below pycnocline, default to surface (T) accounts for air-sea exchange
nem.fun<-function(dat.in, stat, DO_var = 'DO_mgl', depth.val = NULL, 
  bott.stat = F){
  
  ##dependent packages
  require(reshape) #data reshape
  require(wq) #for oxySol function
  require(oce) #for swSigma function

  ##begin calculations

  cat(stat,'\n')
  flush.console()
  strt<-Sys.time()
  
  #columns to be removed prior to processing
  to.rem<-c('flag', 'dTide', 'met.date', 'variable', 'value', 'day.hrs', 
    'dec_time', 'hour')
  dat.in<-dat.in[, !names(dat.in) %in% to.rem]
  
  #convert DO from mg/L to mmol/m3
  dat.in$DO<-dat.in[, DO_var]/32*1000
  
  # get change in DO per hour, as mmol m^-3 hr^-1
  # scaled to time interval to equal hourly rates
  # otherwise, mmol m^-3 0.5hr^-1
  dDO_scl <- as.numeric(diff(dat.in$DateTimeStamp)/60)
  dDO<-diff(dat.in$DO)/dDO_scl
  
  #take diff of each column, divide by 2, add original value
  DateTimeStamp<-diff(dat.in$DateTimeStamp)/2 + dat.in$DateTimeStamp[-c(nrow(dat.in))]
  dat.in<-apply(
    dat.in[,2:ncol(dat.in)],
    2,
    function(x) diff(x)/2 + x[1:(length(x) -1)]
    )
  dat.in<-data.frame(DateTimeStamp,dat.in)
  DO <- dat.in$DO
  
  ##
  # replace missing wx values with climatological means
  # only ATemp, WSpd, and BP
  
  # monthly and hourly averages
  months <- format(dat.in$DateTimeStamp, '%m')
  hours <- format(dat.in$DateTimeStamp, '%H')
  clim_means <- ddply(data.frame(dat.in, months, hours),
    .variables=c('months', 'hours'),
    .fun = function(x){
      data.frame(
        ATemp = mean(x$ATemp, na.rm = T),
        WSpd = mean(x$WSpd, na.rm = T), 
        BP = mean(x$BP, na.rm = T)
      )   
    }
  )
  clim_means <- merge(
    data.frame(DateTimeStamp = dat.in$DateTimeStamp, months,hours),
    clim_means, by = c('months','hours'),
    all.x = T
  )
  clim_means <- clim_means[order(clim_means$DateTimeStamp),]

  # DateTimeStamp order in dat.in must be ascending to match
  if(is.unsorted(dat.in$DateTimeStamp))
    stop('DateTimeStamp is unsorted')
  
  # reassign empty values to means, objects are removed later
  ATemp_mix <- dat.in$ATemp
  WSpd_mix <- dat.in$WSpd
  BP_mix <- dat.in$BP
  ATemp_mix[is.na(ATemp_mix)] <- clim_means$ATemp[is.na(ATemp_mix)]
  WSpd_mix[is.na(WSpd_mix)] <- clim_means$WSpd[is.na(WSpd_mix)]
  BP_mix[is.na(BP_mix)] <- clim_means$BP[is.na(BP_mix)]

  ##
  # get sigma_t estimates
  SigT<-with(dat.in,swSigmaT(Sal,Temp,mean(dat.in$BP/100,na.rm=T)))
  
  #DOsat is DO at saturation given temp (C), salinity (st. unit), and press (atm)
  #DOsat converted to mmol/m3
  #used to get loss of O2 from diffusion
  DOsat<-with(dat.in,get(DO_var)/(oxySol(Temp*(1000+SigT)/1000,Sal)))
  
  #station depth, defaults to mean depth value plus 0.5 in case not on bottom
  #uses 'depth.val' if provided
  if(is.null(depth.val))
    H<-rep(0.5+mean(pmax(1,dat.in$Depth),na.rm=T),nrow(dat.in))
  else H<-rep(depth.val,nrow(dat.in))
  
  #use met.day.fun to add columns indicating light/day, date, and hours of sunlight
  dat.in <- met.day.fun(dat.in)
  
  #get air sea gas-exchange using wx data with climate means
  KL<-with(dat.in, f_calcKL(Temp,Sal,ATemp_mix,WSpd_mix,BP_mix))
  rm(list = c('ATemp_mix', 'WSpd_mix', 'BP_mix'))
  
  #get volumetric reaeration coefficient from KL
  Ka<-KL/24/H
  
  #get exchange at air water interface
  D=Ka*(DO/DOsat-DO)
  
  #combine all data for processing
  proc.dat<-dat.in[,!names(dat.in) %in% c('DateTimeStamp','cDepth','Wdir',
    'SDWDir','ChlFluor','Turb','pH','RH',DO_var,'DO_pct','SpCond','TotPrcp',
    'CumPrcp','TotSoRad','Depth')]
  proc.dat<-data.frame(proc.dat,DOsat,dDO,SigT,H,D)

  #get daily/nightly flux estimates for Pg, Rt, NEM estimates
  out<-lapply(
    split(proc.dat,proc.dat$met.date),
    function(x){
      
      #filter for minimum no. of records 
      if(length(with(x[x$variable=='sunrise',],na.omit(dDO))) < 3 |
         length(with(x[x$variable=='sunset',],na.omit(dDO))) < 3 ){
        DOF_d<-NA; D_d<-NA; DOF_n<-NA; D_n<-NA
        }
      
      else{
        #day
        DOF_d<-mean(with(x[x$variable=='sunrise',],dDO*H),na.rm=T)
        D_d<-mean(with(x[x$variable=='sunrise',],D),na.rm=T)
        
        #night
        DOF_n<-mean(with(x[x$variable=='sunset',],dDO*H),na.rm=T)
        D_n<-mean(with(x[x$variable=='sunset',],D),na.rm=T)
        }
      
      #metabolism
      #account for air-sea exchange if surface station
      #else do not
      if(!bott.stat){
        Pg<-((DOF_d-D_d) - (DOF_n-D_n))*unique(x$day.hrs)
        Rt<-(DOF_n-D_n)*24
      } else {
        Pg<-(DOF_d - DOF_n)*unique(x$day.hrs)
        Rt<-DOF_n*24
        }
      NEM<-Pg+Rt
      Pg_vol<-Pg/mean(x$H,na.rm=T)
      Rt_vol<-Rt/mean(x$H,na.rm=T)
      
      #dep vars to take mean
      var.out<-x[!names(x) %in% c('variable','value','met.date',
        'day.hrs')] 
      var.out<-data.frame(rbind(apply(var.out,2,function(x) mean(x,na.rm=T))))
      data.frame(Station=stat,Date=unique(x$met.date),var.out,DOF_d,D_d,DOF_n,D_n,Pg,Rt,NEM,
        Pg_vol,Rt_vol,numrecs=length(na.omit(x$dDO)))
      
      }
    )
  out<-do.call('rbind',out)
  
  return(out)
  
  }

######
#calculate number of anomolous estimates from nem output
#'nem.in' is output from 'nem.fun'
anoms.fun<-function(nem.in){
  Pg<-nem.in$Pg
  Pg<-sum(Pg<=0,na.rm=T)/length(na.omit(Pg))
  Rt<-nem.in$Rt
  Rt<-sum(Rt>=0,na.rm=T)/length(na.omit(Rt))
  return(data.frame(Pg,Rt))
  }

##
# create dec time using day on 24 hour scale
# 'dat_in' is data frame input with time vector as posix
# output is same data frame including new columns column
dec_fun <- function(dat_in){
  
  # get decimal value by metabolic date for hour/min
  by_met <- dlply(dat_in,
    .variable = 'met.date',
    .fun = function(x){
    
      strt <- (48 - nrow(x))/48
      out <- seq(strt, 1, length = 1 + nrow(x)) 
      out <- out[1:(length(out) - 1)]
      
      out
      
      }
    )
  
  # get continuous day value
  days <- as.character(seq(1:(length(by_met))) - 1)
  names(by_met) <- days
  by_met <- melt(by_met)
  by_met$L1 <- as.numeric(by_met$L1)
  
  # add continuous day value to decimal value
  out <- rowSums(by_met)

  # add to dat_in
  dat_in$dec_time <- out
  
  return(dat_in)
  
  }

######
#function for getting regression weights
# note that this subsets the input data frame for faster wt selection
# subset is by limiting window for product of weights (dec_time)
# subsetted weights are recombined to equal vector of length = nrow(dat_in)
#'wt_vars' is name of three variables to weight
#'ref_in' is row of dat.in that is used as reference
#'dat_in' is data to get weights from
#'wins' are the windows for the three wt.vars, values represent halves
#'all' will return all weights, rather than the product of all three
#'slice' is logical for subsetting 'dat_in' for faster wt selection
#'subs_only' is logical for returning only wt vectors that are non-zero
wt_fun <- function(ref_in, dat_in,
  wt_vars = c('dec_time', 'hour', 'Tide'),
  wins = list(4, 12, NULL),
  all = F, 
  slice = T, 
  subs_only = F){
  
  # sanity check
  if(sum(wt_vars %in% names(dat_in)) != length(wt_vars))
    stop('Weighting variables must be named in "dat_in"')
  
  # windows for each of three variables
  wins_1<-wins[[1]]
  wins_2<-wins[[2]]
  wins_3<-wins[[3]]
  
  # default window width for third variable is half its range
  if(is.null(wins[[3]])) wins_3 <- diff(range(dat_in[, wt_vars[3]]))/2
  
  # weighting tri-cube function
  # mirror extends weighting function if vector repeats, e.g. monthly
  # 'dat_cal' is observation for weight assignment
  # 'ref' is reference observation for fitting the model
  # 'win' is window width from above (divided by two)
  # 'mirr' is logical indicating if distance accounts for repeating variables (e.g., month)
  # 'scl_val' is range for the ref vector of obs, used to get correct distance for mirrored obs
  wt_fun_sub <- function(dat_cal, ref, win, mirr = F, scl_val = 1){
    
    # dist_val is distance of value from the ref
    dist_val <- sapply(ref, function(x) abs(dat_cal - x))
    
    # repeat if distance is checked on non-continuous number line
    if(mirr){
      
        dist_val <- pmin(
          sapply(ref, function(x)
            abs(x + scl_val - dat_cal)),
          sapply(ref, function(x) abs(dat_cal + scl_val - x)),
          dist_val
          )
      
      }
    
    # get wts within window, otherwise zero
    win_out <- dist_val > win
    dist_val <- (1 - (dist_val/win)^3)^3
    dist_val[win_out] <- 0
      
    return(dist_val)
      
    }

  #reference (starting) data
  ref_1 <- as.numeric(ref_in[, wt_vars[1]])
  ref_2 <- as.numeric(ref_in[, wt_vars[2]])
  ref_3 <- as.numeric(ref_in[, wt_vars[3]])

  ##
  # subset 'dat_in' by max window size for faster calc
  # this is repeated if min number of wts > 0 is not met
  # subset vector is all T if not using subset
  dec_rng <- range(dat_in$dec_time)
  ref_time <- unique(ref_in$dec_time)
  dec_sub <- with(dat_in, 
    dec_time > 
      ref_time - wins_1 * 5 & dec_time < ref_time + wins_1 * 5
    )
  if(!slice) dec_sub <- rep(T, length = nrow(dat_in))
  dat_sub <- dat_in[dec_sub, ]

  ##
  # weights for each observation in relation to reference
  # see comments for 'wt_fun_sub' for 'scl_val' argument
  
  # jday
  wts_1 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[1]]), 
    ref = ref_1, win = wins_1, mirr = F) 
  # hour
  wts_2 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[2]]), 
    ref = ref_2, win = wins_2, mirr = T, scl_val = 24)
  # tide
  wts_3 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[3]]), 
    ref = ref_3, win = wins_3, mirr = F)
  # all as product 
  out <- sapply(1:nrow(ref_in), function(x) wts_1[, x] * wts_2[, x] * wts_3[, x])
  
  gr_zero <- colSums(out > 0)
  #cat('   Number of weights greater than zero =',gr.zero,'\n')
  
  # extend window widths of weight vector is less than 100
  while(any(gr_zero < 100)){
    
    # increase window size by 10%
    wins_1 <- 1.1 * wins_1
    wins_2 <- 1.1 * wins_2
    wins_3 <- 1.1 * wins_3 
    
    # subset again
    dec_sub <- with(dat_in, 
      dec_time > ref_time - wins_1 * 5 & dec_time < ref_time + wins_1 * 5
      )
    if(!slice) dec_sub <- rep(T, length = nrow(dat_in))
    dat_sub <- dat_in[dec_sub, ]
    
    #weights for each observation in relation to reference
    wts_1 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[1]]), 
      ref = ref_1, win = wins_1, mirr = F)
    wts_2 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[2]]), 
      ref = ref_2, win = wins_2, mirr = T, scl_val = 24)
    wts_3 <- wt_fun_sub(as.numeric(dat_sub[, wt_vars[3]]), 
      ref = ref_3, win = wins_3, mirr = F)
    
    out <- sapply(1:nrow(ref_in), 
      function(x) wts_1[, x] * wts_2[, x] * wts_3[, x])
    
    gr_zero <- colSums(out > 0)
    
    }
  
  if(subs_only){
    
    nms <- which(dec_sub)
    out <- alply(out, 2, function(x) {
    
      to_sel <- x > 0
      tmp <- x[to_sel]
      names(tmp) <- which(dec_sub)[to_sel]
      tmp
    
      })
    
    return(out)
    
    }
  
  # extend weight vectors to length of dat_in
  empty_mat <- matrix(0, ncol = nrow(ref_in), nrow = nrow(dat_in))
  empty_fill <- function(wts_in) {
    out <- empty_mat
    out[dec_sub,] <- wts_in
    out
    }
  wts_1 <- empty_fill(wts_1)
  wts_2 <- empty_fill(wts_2)
  wts_3 <- empty_fill(wts_3)  
  out <- empty_fill(out)

  #return all weights if T
  if(all){
    out <- data.frame(dat_in$DateTimeStamp, 
      wts_1, wts_2, wts_3, out)
    names(out) <- c('DateTimeStamp', wt_vars, 'final')
    return(out)
    }
  
  #final weights are product of all three
  out
  
  }
  
######
# preps SWMP data in 'proc5' and 'tide_pred' for weighted regression
# 'site_in' input is character string of five letter site name
# 'wq_path' is character of path for wq files
# 'tide_path' is character of path for tidal prediction files
# 'DO_var' is character of name of column with DO, renamed as 'DO_obs'
# 'interp' is logical indicating if missing DO values are interpolated, max gap is four hours
prep_wtreg <- function(dat_in, 
  DO_var = 'DO_mgl',
  interp = T){
  
  to_proc <- dat_in
  
  # get dTide
  to_proc$dTide <- with(to_proc, c(diff(Tide)[1], diff(Tide)))
  
  # get metabolic day info
  to_proc <- met.day.fun(to_proc)
  
  # setup as decimal time
  to_proc <- dec_fun(to_proc)
  
  to_proc$hour <- as.numeric(format(to_proc$DateTimeStamp, '%H')) +
    as.numeric(format(to_proc$DateTimeStamp, '%M'))/60
  
  # remove extra cols
  to_rm <- c('SpCond', 'DO_pct', 'cDepth', 'pH', 'Turb', 'ChlFluor',
    'RH', 'Wdir', 'SDWDir', 'TotPrcp', 'CumPrcp', 'TotSoRad', 
    'PO4H', 'NH4F', 'NO2F', 'NO3F', 'NO23F', 'CHLA_N', 'flag')
  to_proc <- to_proc[, !names(to_proc) %in% to_rm]
  
  # reassign name for DO_var
  names(to_proc)[names(to_proc) %in% DO_var] <- 'DO_obs'
  
  # interp missing values using max gap of four hours (8 obs)
  # na.rm = F keeps leading and trailing NA vals
  if(interp)
    to_proc$DO_obs <- with(to_proc, 
      na.approx(DO_obs, x = DateTimeStamp, maxgap = 8, na.rm = F)
      )
  
  return(to_proc)
    
  }
  
######
# get predicted, normalized values not using interp grid, tide as predictor
# 'dat_in' is raw data used to create 'grd_in' and used to get predictions
# 'DO_obs' is string indicating name of col for observed DO values from 'dat_in'
# output is data frame same as 'dat_in' but includes predicted and norm columns
wtreg_fun <- function(dat_in, DO_obs = 'DO_obs', wins = list(4, 12, NULL),
  parallel = F, progress = F){

  # get mean tidal height from empirical data
  mean_tide <- mean(dat_in$Tide)

  #for counter
  strt <- Sys.time()
  
  out <- ddply(dat_in, 
    .variable = 'DateTimeStamp',
    .parallel = parallel, 
    .fun = function(row){
      
      # row for prediction
      ref_in <- row
      ref_in <- ref_in[rep(1, 2),]
      ref_in$Tide <- c(unique(ref_in$Tide), mean_tide)
      
      # progress
      if(progress){
        prog <- which(row$DateTimeStamp == dat_in$DateTimeStamp)
        sink('log.txt')
        cat('Log entry time', as.character(Sys.time()), '\n')
        cat(prog, ' of ', nrow(dat_in), '\n')
        print(Sys.time() - strt)
        sink()
        }
      
      # get wts
      ref_wts <- wt_fun(ref_in, dat_in, wins = wins, slice = T, 
        subs_only = T, wt_vars = c('dec_time', 'hour', 'Tide'))
  
      #OLS wtd model
      out <- lapply(1:length(ref_wts),
        function(x){
          
          # subset data for weights > 0
          dat_proc <- dat_in[as.numeric(names(ref_wts[[x]])),]
          
          # if no DO values after subset, return NA
          # or if observed DO for the row is NA, return NA
          if(sum(is.na(dat_proc$DO_obs)) == nrow(dat_proc)|
              any(is.na((ref_in$DO_obs)))){
            
            DO_pred <- NA
            beta <- NA
            Tide <- ref_in$Tide[x]
            
            } else {
            
              # subset weigths > 0, rescale weights average
              ref_wts <- ref_wts[[x]]/mean(ref_wts[[x]])
            
              # get model
              mod_md <- lm(
                DO_obs ~ dec_time + Tide, # + sin(2*pi*dec_time) + cos(2*pi*dec_time),
                weights = ref_wts,
                data = dat_proc
                )
            
              # get prediction from model
              Tide <- ref_in$Tide[x]
              DO_pred <- predict(
                mod_md, 
                newdata = data.frame(dec_time = ref_in$dec_time[x], Tide = Tide)
                )
            
              # get beta from model
              beta <- mod_md$coefficients['Tide']
            
            }
          
          # output
          DO_pred
          
          }
        
        )

      out <- unlist(out)
      names(out) <- c('DO_prd', 'DO_nrm')
      out
      
      })
  
  out$DateTimeStamp <- NULL
  out <- cbind(dat_in, out)

  return(out)
  
  } 
  

