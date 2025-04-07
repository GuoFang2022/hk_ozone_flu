########################################################
## "HK Ozone influenza activity-GLM Analysis"
##     Analysis Code File
##
##             BY Fang
########################################################

################################## Load packages and data ###################################
library(tidyverse); library(lubridate); library(dlnm); library(splines);
library(tsModel); library(gnm);library(ggpubr); library(mgcv);library(tseries)

load("dat_HK_fluP.rda")

dfA=dat_HK_fluP_wk_C %>% 
  mutate(date=WeekEnd) %>%
  mutate(time = as.numeric(date), 
         year = year(date),
         month = month(date))

str(dfA)
summary(dfA)

################################## Diagnosis of core model ###################################
#To examine the residuals of the core model
formula_str=paste('fluP', "~as.factor(year) + as.factor(month) + loglag1+ loglag2+ loglag3+ loglag4", sep='')
dfAcore=dfA%>%
  mutate(loglag1 = log(Lag(fluP, 1)),
         loglag2 = Lag(loglag1, 1),
         loglag3 = Lag(loglag2, 1),
         loglag4 = Lag(loglag3, 1),
         loglag5 = Lag(loglag4, 1))
core= gam(as.formula(formula_str), data=dfAcore, family=quasibinomial, na.action=na.omit)
plot(residuals(core),xlab="Time", ylim=c(-1,1),main='Residuals Plot')
pacf(residuals(core),ylim=c(-0.15, 0.15), xlim=c(0,20), main = "PACF of Model Residuals")

# Ljung-Box test p-value ≥ 0.05, no significant autocorrelation exists.
Box.test(residuals(core), lag = 20, type = "Ljung-Box")

################################### Function for GLM analysis ###################################
glm_function <- function(dat, Y='fluP', x1='value', xlag=0, variable, datname) {
  
  df <- dat %>%
    mutate(loglag1 = log(Lag(fluP, 1)),
           loglag2 = Lag(loglag1, 1),
           loglag3 = Lag(loglag2, 1),
           loglag4 = Lag(loglag3, 1))
  
  formula_str <- switch(variable,
                        o3h8max= paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(ah,', xlag, ')+ as.factor(year) + as.factor(month) + loglag1+ loglag2+ loglag3+ loglag4', sep=''),
                        ah   = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(temp,', xlag,  ")+ as.factor(year) + as.factor(month) + loglag1+ loglag2+ loglag3+ loglag4", sep=''),
                        temp = paste(Y, '~Lag(', x1, ',', xlag, ')+ Lag(ah,', xlag, ")+ as.factor(year) + as.factor(month) + loglag1+ loglag2+ loglag3+ loglag4", sep='')
  )
  
  fit <- gam(as.formula(formula_str), data=df, family=quasibinomial, na.action=na.omit)
  summ <- summary(fit)
  
  outA <- data.frame(
    beta = summ$p.coeff[2],
    se   = summ$se[2],
    t    = summ$p.t[2],
    p    = summ$p.pv[2],
    lag   = xlag
  ) 
  
  return(outA)
}

variables <- c("o3h8max", "ah", "temp")
results <- list()

for (var in variables) {
  plist <- list(
    dat      = list(dfA),
    Y        = 'fluP',
    x1       = var,
    xlag     = 0:2,
    variable = var
  ) %>% cross_df()
  
  results[[var]] <- plist %>% pmap_df(glm_function)
}


#Calculate SD of each environmental predictor
df_SD=dfA %>%
  dplyr::summarize(o3h8max=sd(o3h8max), ah=sd(ah), temp=sd(temp)) %>%
  gather(plt, SD)

# Store the results for easy access
out_o3   <- results$o3h8max %>% as.data.frame() %>% mutate(plt='o3h8max')
out_ah   <- results$ah %>% as.data.frame() %>% mutate(plt='ah')
out_temp <- results$temp %>% as.data.frame() %>% mutate(plt='temp')

df_gam_final=rbind(out_o3,out_ah,out_temp) %>%
  left_join(df_SD,by='plt') %>%
  mutate(betalow=beta-1.96*se, 
         betahigh=beta+1.96*se) %>% # 95% CI; alpha=0.05
  mutate(Size=beta*SD,SizeL=betalow*SD,SizeH=betahigh*SD) %>%
  select(plt, lag, beta, betalow, betahigh, SD, Size, SizeL, SizeH, glm.p=p) %>%
  mutate(sig=ifelse(glm.p<0.05, 'sig', 'non_sig'),
         sig=factor(sig, levels=c("non_sig","sig"))) %>% 
  mutate(plt=factor(plt, levels=c( "o3h8max","ah", "temp"),
                    labels=c("O[3]", "AH", "T"))) %>% 
  arrange(plt)

df_gam_final

p.glm=ggplot(df_gam_final) +
  geom_hline(yintercept = 0,color='grey50', linetype='dashed') +
  geom_errorbar(aes(x=-lag, ymin =SizeL, ymax = SizeH), 
                color='gray10', width = 0.1, linewidth=0.8, alpha=0.5) +
  geom_point(aes(x=-lag,  y=Size, fill=Size), shape = 21, size=1.5, alpha=1) +
  facet_grid(.~plt,labeller = label_parsed) +
  scale_x_continuous(name=expression(paste("Lag (", italic("week"), ")")),
                     breaks=-2:0,labels = c(-2,-1,0), expand = c(0.15, 0.15)) +
  scale_y_continuous(name = expression(atop(NA, atop(paste("GLM Effect Estimate"), beta))),
                     labels = scales::number_format(accuracy = 0.01)) +
  scale_fill_gradientn(
    colors = c("#053061", "#2166AC", "#4393C3", "#D1E5F0", "white", "#FDDBC7", "#F4A582", "#D6604D","#B2182B"),
    values = scales::rescale(c(quantile(df_gam_final$Size, seq(0, 0.5, 0.16)), 0, quantile(df_gam_final$Size, seq(0.6, 1, 0.13))))) +
  theme_bw() +
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

p.glm
