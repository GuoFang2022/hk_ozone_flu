########################################################
## "HK Ozone influenza activity-CCM Analysis"
##     Analysis Code File
##
##             BY Fang
########################################################

################################## Load packages and data ###################################
packages=c('tidyverse','knitr','lubridate','rEDM',
           'doParallel','foreach','tseries', 'ggbeeswarm')
lapply(packages, require, character.only=T)

cores_all=detectCores()
cores=ifelse(cores_all<9,4,cores_all-2)
core_type='PSOCK'

load("dat_HK_fluP.rda")
dfA=dat_HK_fluP_wk_C
dfA$fluP <- dfA$logitFluP
str(dfA)
summary(dfA)

################################## Stationarity assumption test ###################################
df=dfA %>% select(-c(WeekID, WeekEnd))
statationary_test=function(x){
  oldw <- getOption("warn")
  options(warn = -1)
  A=kpss.test(unlist(df[,x]), null="Trend") ## a low p-value tells not stationary
  B=adf.test(unlist(df[,x])) # a large p-value tell not stationarity
  options(warn = oldw)
  
  AB=data.frame(kpss=A$p.value,adf=B$p.value,plt=x)
  AB
}

norm_out=names(df) %>% map_df(statationary_test) %>% print()

################################## Data normalization ###################################
nomz=function(x, normalization=T, dseason=T, season_sd=T, sea=365, dtrend=T, dTtype="linear"){
  x=as.numeric(x)
  xt=x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt=c(NA,diff(xt))} else if (dtrend==T & dTtype=="linear"){
    lm.t=lm(xt~c(1:length(xt)))
    xt=xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt=xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt=xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}

fn_nomz=function(x){x=nomz(x,
                              normalization=T, dtrend=F, 
                              dseason=F, season_sd=F,
                              sea=52.18, dTtype="linear")}
df_smapc=dfA %>%
  rename(date=WeekEnd) %>% 
  mutate_at(vars(-c(date, WeekID)),fn_nomz) 

df_smapc$state <- "all"

str(df_smapc)
summary(df_smapc)

################################## System property test ###################################
######## Determine optimal E for the system ######## 
fn_E_smapc=function(data,ST,y){
  
  df=data %>%
    filter(state==ST) %>%
    select(date,y) %>%
    na.omit()
  M <- nrow(df)
  lib <- c(1,M)
  pred <- c(1,M)
  E=EmbedDimension(dataFrame=df,
                   lib=lib,
                   pred=pred,
                   columns=y, target=y,
                   Tp = 1, #default
                   maxE=20,
                   showPlot=F)
  temp=data.frame(dis=y,E,ST)
  temp
}
plist=list(data=list(df_smapc),y=c('fluP'),ST=unique(df_smapc$state)) %>%
  cross_df()
E_smapc_out=plist %>% pmap_df(fn_E_smapc)

find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pk <- unlist(pks)
  pks <- ifelse(is.null(pk)==T, which.max(x), pk)
  pks
}
E_smapc =E_smapc_out %>% filter(E %in% 1:15) %>%
  group_by(dis) %>%
  mutate(row_number=1:n()) %>%
  filter(row_number==find_peaks(rho,m=0)[1])  %>%
  as.data.frame() %>%
  select(-row_number)
E_smapc

######## Evidence of nonlinearity in the system/Determine optimal theta for S-map ######## 
fn_theta_justY=function(data, ST, dis, theta){
  E=E_smapc[E_smapc[,'ST']==ST & E_smapc[,'dis']==dis,'E']
  
  df=data %>%
    filter(state==ST) %>%
    select(date,dis)
  M <- nrow(df)
  lib <- c(1,M)
  pred <- c(1,M)
  rho_theta = PredictNonlinear(dataFrame = df,
                               embedded = FALSE,
                               columns = dis,
                               target = dis,
                               Tp=1,
                               theta=theta,
                               lib=lib,
                               pred=pred,
                               showPlot = FALSE,
                               E = E)
  best_theta_df=rho_theta %>%
    mutate(dis=dis, state=ST, theta=theta)
}

plist=list(data=list(df_smapc),ST=unique(df_smapc$state),
           dis=c('fluP'),
           theta=c(0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4,
                   5, 6, 7, 8, 9)) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
theta_out=foreach(i = 1:nrow(plist),
                  .packages = c("rEDM","tidyverse",'lubridate'),
                  .combine=rbind,
                  .inorder=FALSE)  %dopar% {
                    theta_out=plist[i,] %>% pmap_df(fn_theta_justY)
                  }
stopCluster(cl)

best_theta=theta_out %>%
  group_by(dis) %>%
  mutate(row_number=1:n()) %>%
  filter(row_number==find_peaks(rho,m=0)[1]) %>%
  as.data.frame()
best_theta

######## Evidence of determinism in the system ######## 
rho_Tp_fluP <- PredictInterval(dataFrame = df_smapc %>% select(date, fluP), lib = "1 626",#"1 200",
                               pred = "1 626", target = 'fluP', columns = 'fluP', E = 5)


##################################  1. Surrogate CCM causality test ###################################
######## Determine E for pair-wise cross mapping ######## 
fn_E_ccm=function(data,plt,dis,tp_value,E){
  df=data %>%
    mutate(plt_tp=lag(.data[[plt]],-tp_value)) %>%
    filter(!(is.na(plt_tp))) %>%
    select(date,dis,plt_tp)  
  names(df) = c("date", dis, plt)
  M=nrow(df)
  E=E 
  libSize =c(E+2,M-E-2)
  ccm_out = CCM(dataFrame = df, E = E, 
                Tp = 0, 
                columns = dis,
                target = plt,
                libSizes = libSize,
                sample = 100,
                random=T,
                seed=2019)
  outB=ccm_out %>% 
    select(LibSize,rho=paste0(dis, ":", plt, sep = "")) %>% 
    pivot_wider(names_from='LibSize',values_from ='rho') %>% 
    rename_at(names(.),~c('minRho','maxRho')) %>% 
    mutate(dis=dis,plt=plt,E,tp_value)
  outB
}

plist=list(data=list(df_smapc),
           dis=c('fluP'),
           plt=c('o3h8max','ah','temp'),
           E=2:10,
           tp_value=-2:0) %>% cross_df()
cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
ccm_out=foreach(i = 1:nrow(plist),
                .packages = c("rEDM","tidyverse","lubridate"),
                .combine=rbind,
                .inorder=FALSE)  %dopar% {
                  plist[i,] %>% pmap_df(fn_E_ccm)
                }
stopCluster(cl)

E_ccm=ccm_out %>%
  mutate(vergeRho=maxRho-minRho) %>%
  group_by(dis,plt,tp_value) %>%
  mutate(row_number=1:n()) %>%
  filter(row_number==find_peaks(vergeRho,m=0)[1])  %>%
  as.data.frame() %>% 
  filter(tp_value==0)
E_ccm

######## Surrogate data ######## 
df_flu=df_smapc %>% select(date, o3h8max, ah,temp) %>%
  gather(key,value,-c(date)) %>% mutate(state="all")
str(df_flu)
# calculate spar values for independent variables
fn_state_spar=function(states,keys){
  splineres <- function(spar){
    res <- rep(0, length(x))
    for (i in 1:length(x)){
      mod <- smooth.spline(x[-i], y[-i], spar = spar)
      res[i] <- predict(mod, x[i])$y - y[i]
    }
    return(sum(res^2))
  }
  
  x=df_flu %>% ungroup() %>%
    filter(state==states & key==keys) %>%
    select(date) %>% pull() %>% as.numeric()
  y=df_flu %>% ungroup() %>%
    filter(state==states & key==keys) %>%
    select(value)  %>% pull() %>% as.numeric()
  
  spars <- seq(0, 1.5, by = 0.1)
  ss <- rep(0, length(spars))
  ss=foreach(i = 1:length(spars),
             .combine=rbind,
             .inorder=FALSE)  %dopar% {
               targetCol = paste("T", i, sep = "")
               ss[i] <- splineres(spars[i])
             }
  spar=spars[which.min(ss)]
  data.frame(state=states,plt=keys,spar)
}

plist=list(states=unique(df_flu$state),
           keys=c('o3h8max', 'ah','temp')) %>%
  cross_df()
cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
flu_spar=plist %>% pmap_df(fn_state_spar)
stopCluster(cl)

flu_spar

# calculate sd values (alpha) for independent variables (additive noise factor to produce surrogate data)
yearday_anom <- function(t,x,spars){
  # t: date formatted with POSIXt
  # x: time-series values to compute seasonal mean and anomaly
  doy <- as.numeric(strftime(t, format = "%j"))
  I_use <- which(!is.na(x))
  # create time indices to use for smoothing, replicating data to "wrap around"
  doy_sm <- rep(doy[I_use],3) + rep(c(-366,0,366),each=length(I_use))
  x_sm <- rep(x[I_use],3)
  xsp <- smooth.spline(doy_sm, y = x_sm, w = NULL,
                       spar = spars, cv = NA,
                       all.knots = TRUE,keep.data = TRUE, df.offset = 0)
  xbar <- data.frame(t=t,doy=doy) %>%
    left_join(data.frame(doy=xsp$x,xbar=xsp$y),by='doy') %>%
    select(xbar)
  out = data.frame(t=t,mean=xbar,anomaly=(x - xbar))
  names(out) <- c('date','mean','anomaly')
  return(out)
}

fn_anomaly_PNAS=function(states,plts){
  vec_t=df_flu %>%
    ungroup() %>%
    filter(state==states & key==plts) %>%
    select(date) %>% pull()
  vec_x=df_flu %>%
    ungroup() %>%
    filter(state==states & key==plts) %>%
    select(value) %>% pull()
  spars=flu_spar %>% filter(state==states & plt==plts) %>%
    select(spar) %>% pull()
  df_9=yearday_anom(vec_t,vec_x,spars)
  sd=sd(df_9$anomaly,na.rm=TRUE)
  data.frame(states,plt=plts,spar=spars,sd_PNAS=sd)
}

plist=list(states=unique(df_flu$state),
           plts=c('o3h8max','ah','temp')) %>% cross_df()
sd_data=plist %>% pmap_df(fn_anomaly_PNAS)

sd_data

# Compilation of function for producing surrogate data
fn_surr_data=function(data,ST,plts, tp_value){
  df=data %>%
    filter(state==ST) %>%
    mutate(plt=lag(.data[[plts]],-tp_value)) %>%
    filter(!(is.na(plt))) %>%
    select(date,"plt")
  alpha=sd_data %>% filter(states==ST & plt==plts) %>%
    select(sd_PNAS) %>% pull()
  set.seed(2019)
  surr_data=
    SurrogateData(unlist(df[,"plt"]), method = "ebisuzaki", # T_period = 52.18,
                  num_surr = num_surr,
                  alpha=alpha) %>%
    as.data.frame()
  df=df %>% select(-"plt")
  surrA=bind_cols(df,surr_data)
} 

######## CCM causality test with 1000 surrogate data ######## 
num_sample=100
num_surr=100

fn_season_ccm=function(data,ST,x,y,tp_value){
  df=data %>% 
    filter(state==ST) %>%
    select(date,y,x)
  
  E=E_ccm[E_ccm[,'plt']==x & E_ccm[,'dis']==y,'E']
  
  surr_data <- fn_surr_data(df_smapc,ST,x,0)
  all_data <- df %>% left_join(surr_data,by="date")
  names(all_data) = c("date", y, 'T1',paste0("T", 2:(num_surr+1)))
  
  m=nrow(all_data) %>% as.data.frame()
  libSize =c(E+2,m-E-2)
  
  rho_surr <- NULL
  for (i in 1:(num_surr+1)) {
    targetCol = paste("T", i, sep = "")
    ccm_out = CCM(dataFrame = all_data, E = E, Tp = tp_value,
                  columns = y,
                  target = targetCol,
                  libSizes = libSize,
                  random=T,
                  sample = num_sample,
                  seed=2019)
    col = paste(y, ":", targetCol, sep = "")
    dat=ccm_out %>% select(LibSize,col)
    names(dat)=c("lib","rho")
    test1=mutate(dat,i=i,dis=y,plt=x,E=E,
                 tp_value=tp_value,state=ST)
    rho_surr <- rbind(rho_surr,test1)
    
  }
  rho_surr
}

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           y=c('fluP'),
           x=c('o3h8max','ah','temp'),
           tp_value=-2:0) %>% cross_df()
cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
ccm_out=foreach(j = 1:nrow(plist),
                .packages = c("rEDM","tidyverse"),
                .combine=rbind,
                .inorder=FALSE)  %dopar% {
                  ccm_out=plist[j,] %>% pmap_df(fn_season_ccm)
                }
stopCluster(cl)

# Calculate the difference in cross-mapping skills obtained by the maximum and the minimum library to test convergence property
dat_min <- ccm_out %>% filter(lib<50) 
dat_min <- dat_min[order(dat_min$state,dat_min$tp_value,dat_min$plt,dat_min$dis,dat_min$i),]
dat_max <- ccm_out %>% filter(lib>50) 
dat_max <- dat_max[order(dat_max$state,dat_max$tp_value,dat_max$plt,dat_max$dis,dat_max$i),]

ccm_dat <- cbind(dat_max,dat_min[,"rho"]) 
names(ccm_dat) <- c("lib","rho_max","i","dis","plt","E","tp_value","ST","rho_min")

ccm_dat1 <- ccm_dat %>% mutate(rho=rho_max-rho_min)

ccm_out_raw = ccm_dat1  %>%
  filter(i==1) %>%
  select(plt, dis, ST,tp_value, rho)

ccm_p=ccm_dat1 %>%
  group_by(dis,plt,E,ST,tp_value) %>%
  summarise(p=1-ecdf(rho[i != 1])(rho[i == 1])) %>% 
  left_join(ccm_out_raw, by=c("plt", "dis", "ST", "tp_value")) %>% 
  rename(rho_raw=rho) %>%
  arrange(dis,plt) %>% 
  print(n=Inf)

ccm_causal <- ccm_dat1 %>%
  full_join(ccm_p, by = c("ST","plt","dis",'E','tp_value')) %>% 
  mutate(i=i-1,grp=ifelse(i==0,'raw','surr'),
         sig=ifelse(p<0.05 & rho_raw>0, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig"))) 

##################################################
p_function_short=function(dat, x,pH,pL){
  dfg2_Q1=dat %>% filter(grp=='surr') %>%
    group_by(dis,plt,tp_value) %>%
    summarise(ph=quantile(.data[[x]], pH,na.rm=T),
              pl=quantile(.data[[x]], pL,na.rm=T))
  dfg2_Q0=dat %>% filter(grp=='raw') %>% select(-grp) %>%
    group_by(dis,plt,tp_value) %>%
    select(p0=x)
  dfg2_Q=left_join(dfg2_Q1,dfg2_Q0,by=c('tp_value', 'dis','plt'))
  print(dfg2_Q)
}

ccm_plot_data <- p_function_short(dat=ccm_causal, x='rho',pH=0.95,pL=0.05) %>%
  mutate(plt= factor(plt,levels = c("o3h8max", "ah", "temp"),
                     labels = c("O[3]","AH","T")))  

ccm_causal %>% group_by(dis, tp_value, plt) %>%
  filter(grp=='raw') %>%
  select(dis, plt, tp_value, p) %>% 
  arrange(plt) %>% 
  print(Inf)

data_anno = ccm_plot_data %>% group_by(dis, tp_value, plt) %>%
  filter(p0>ph) %>%
  arrange(plt)

dat_text = data.frame(
  label = c('*'),
  plt   = data_anno$plt,
  x     = data_anno$tp_value,
  y     = data_anno$p0
)

myTheme <- theme_bw() +
  theme(
    axis.title = element_text(size = 12,family = 'serif'),
    axis.title.y = element_text(size = 17, vjust= 0.01),
    axis.title.x = element_text(margin = margin(t = 1)),
    axis.text=element_text(size = 11,family = 'serif'),
    plot.margin = margin(0, 4, 2, 5, "pt"),
    legend.position = 'none',
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(size = 15, hjust=0.5),
    strip.text = element_text(size = rel(1), face = "bold",family = 'serif'),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.switch.pad.grid = unit(-0.18, "cm")) 

p_ccm_surr <- ggplot(ccm_plot_data) +
  geom_hline(yintercept = 0,color='grey50', linetype='dashed') +
  geom_line(aes(tp_value, p0), color='gray20',   size=1) +
  geom_point(aes(tp_value, p0), color='black', size=2) +
  geom_ribbon(aes(x=tp_value, ymin =pl, ymax = ph), fill='gray20', color=NA,alpha=0.2) +
  facet_grid(.~plt,
             labeller = label_parsed) +
  labs(x=expression(paste("Lag (", italic("week"), ")")),
       y=expression(atop(NA,atop(paste("CCM Causality Test"), paste(Delta, rho["lib"]))))) +
  scale_x_continuous(breaks=-2:0,labels = c(-2,-1,0), expand = c(0.15, 0.15)) +
  scale_y_continuous(limits=c(-0.10, 0.30), breaks=seq(-0.10, 0.30, 0.10),
                     labels = c('-0.10', '0','0.10','0.20','0.30')) +
  myTheme +
  geom_text(data = dat_text,
            mapping = aes(x = x, y = y, label = label), parse = F, 
            vjust = 0.66, color='red', size=10) 


p_ccm_surr

##################################  2. S-map: effect strength ###################################
lead_lag_custom <- function(plt, n) {
  if (tp_value > 0) {
    plt <- lead(plt, n = n) # Use lead if n is positive
  } else {
    plt <- lag(plt, n = abs(n)) # Use lag if n is 0 or negative
  }
}

fn_smapc=function(data,ST,plt,dis,tp_value){
  df=data %>% filter(state==ST) %>%
    mutate(plt_tp=lead_lag_custom(.data[[plt]],tp_value)) %>%
    filter(!(is.na(plt_tp))) %>%
    select(date,all_of(dis),plt_tp) %>%
    na.omit()
  
  M <- nrow(df)
  E=E_smapc[E_smapc$ST==ST & E_smapc$dis==dis,"E"]
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis)
  
  dataFrame = cbind(df[E:M, 'date'],df[E:M, dis],
                    embed_1[E:M, 1:(E-1)], df[E:M, 'plt_tp']) %>%
    as.data.frame()
  names(dataFrame)=c('date',dis,letters[1:(E-1)],plt)
  
  columns = paste(paste(letters[1:(E-1)],collapse =' '), plt, sep=' ')

  smap = SMap(dataFrame = dataFrame,
              embedded = TRUE,
              columns = columns,
              target = dis,
              lib = c(1,nrow(dataFrame)),
              pred = c(1,nrow(dataFrame)),
              theta=best_theta[best_theta$state==ST & best_theta$dis==dis, "theta"],
              Tp = 1,  # one-week forward forecast
              E = E)
  smapc_df=smap$coefficients[c(1,2+E)]
  names(smapc_df)=c('date','effect')
  smapc_df=smapc_df %>%
    mutate(date=lubridate::as_date(date, origin = lubridate::origin)) %>%
    mutate(dis=dis,ST=ST, plt=plt,E=E,tp_value=tp_value)
}

plist=list(data=list(df_smapc),
           ST=unique(df_smapc$state),
           dis=c('fluP'),
           plt=c('o3h8max','ah','temp'),
           tp_value=-1:1) %>% 
  cross_df() #tp_value=1 here studies lag-0 effect given setting of "Tp=1",so on and  so forth

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
C_out=foreach(i = 1:nrow(plist),
              .packages = c("rEDM","tidyverse","lubridate"),
              .combine=rbind,
              .export='best_theta',
              .inorder=FALSE)  %dopar% {
                C_out=plist[i,] %>% pmap_df(fn_smapc)
              }

stopCluster(cl)

dat_plt = C_out %>% 
  group_by(ST, tp_value,plt,dis) %>%
  filter(effect < quantile(effect, probs=.95, na.rm = T),
         effect > quantile(effect, probs=.05, na.rm = T)) %>% 
  mutate(plt= factor(plt,levels = c("o3h8max", "ah", "temp"),
                     labels = c("O[3]","AH","T")),
         tp_value=tp_value-1) # define lag-0 as lag 0 properly

p_effect_lag <- ggplot()+
  geom_violin(data=dat_plt, aes(x= factor(tp_value), y=effect),width = 0.8,
              trim = FALSE, color='gray30') +
  geom_jitter(data=dat_plt, aes(x=factor(tp_value), y=effect, color=effect),
              size = 0.01, width = 0.3, alpha=0.5)+
  facet_grid(.~plt, labeller = label_parsed) +
  labs(x=expression(paste("Lag (", italic("week"), ")")),
       y=expression(atop(NA,atop(paste("S-map Effect Estimate"), paste(partialdiff, 'Flu/', partialdiff, 'Env')))))+
  scale_x_discrete(expand = c(0.1, 0.25)) +
  scale_y_continuous(limits=c(-0.13, 0.02), breaks=seq(-0.12, 0, 0.03),
                     labels = c('-0.12','-0.09', '-0.06','-0.03', '0')) +
  geom_hline(yintercept = 0,color='grey50', linetype='dashed') +
  scale_colour_gradientn(
    colors = c("#053061",  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white", "#FDDBC7", "#F4A582"),
    values = scales::rescale(c(quantile(dat_plt$effect, seq(0, 0.95, 0.19)), 0, quantile(dat_plt$effect, 0.985), max(dat_plt$effect)))
  ) +
  myTheme

p_effect_lag
