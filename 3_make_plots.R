######################################################
# INSTALL AND LOAD DEPENDENCIES ----

#CRAN packages --- install each with install.packages()
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggsci)
# # # library(scales) #alpha
library(RColorBrewer)
# library(gridExtra)
library(viridisLite)
library(lme4)
library(splines)

######################################################
# SETUP ----

main_dir <- '/my/project/dir' #path to repo contents
data_dir <- file.path(main_dir,'data')
result_dir <- file.path(main_dir,'results') 
setwd(data_dir)
load(file='subjects.Rdata') #subjects
task <- 'RRC'

######################################################
# VISIT TIMING FIGURE ----

groups <- ifelse(grepl('A', subjects), 'ALS', 'HC')
n_ALS <- sum(groups=='ALS')
n_HC <- sum(groups=='HC')
colors_ALS <- viridis(n_ALS, begin = 0, end=0.8)
colors_HC <- viridis(n_HC, option='magma', begin = 0, end = 0.8)

# SCAN TIMING DATA
scanDaysALS <- read.csv('scanDaysALS.csv', header=FALSE)
scanDaysHC <- read.csv('scanDaysHC.csv', header=FALSE)
names(scanDaysALS) <- names(scanDaysHC) <- c('subject','visit','days_since_visit1','days_since_onset')
scanDaysHC$days_since_onset <- NA
scanDays <- rbind(scanDaysALS, scanDaysHC)
scanDays$group <- substr(scanDays$subject, 1, 1) #A or C
scanDays$group <- factor(scanDays$group, levels=c('A','C'), labels=c('ALS','HC'))

# EXCLUDE VISIT 1 FROM C09 (occurred 1055 days before visit 2) & A14
scanDays <- scanDays[!(scanDays$subject=='C09' & scanDays$visit == 1),]
scanDays$visit[scanDays$subject=='C09'] <- scanDays$visit[scanDays$subject=='C09'] - 1
scanDays$days_since_visit1[scanDays$subject=='C09'] <- scanDays$days_since_visit1[scanDays$subject=='C09'] - min(scanDays$days_since_visit1[scanDays$subject=='C09'])

ggplot(scanDays, aes(x=days_since_visit1, y=subject)) +
  geom_line(aes(group=subject), color='gray') + geom_point() +
  #scale_color_manual(values=c(colors_ALS, colors_HC)) +
  facet_grid(group ~ ., scales='free', space='free') +
  theme_few() + theme(legend.position='none') + xlab('Days Since Enrollment')

# GROUP ALS SUBJECTS BASED ON FOLLOW-UP TIME (RELATED TO DISEASE PROGRESSION AND SURVIVAL)
max_days <- scanDaysALS %>%
  group_by(subject) %>%
  summarize(length=max(days_since_onset)) %>%
  arrange(length)
max_days$years <- ceiling(max_days$length/365)
max_days$years[max_days$length < 365/2] <- 0.5
max_days$ALS_group <- 'group1' #up to 2 years following onset
max_days$ALS_group[max_days$years > 2] <- 'group2' #2-5 years following onset
max_days$ALS_group[max_days$years > 5] <- 'group3' #over 5 years following onset
max_days$ALS_group[max_days$subject=='A04'] <- 'group4'
max_days$ALS_group[max_days$subject=='A11'] <- 'group2'

######################################################
# ALSFRS-R FIGURE ----

ALSFRS <- read.csv('scanALSFRS.csv', skip=3)[,-1]
cols_ALSFRS <- c(3:7,9:15) #ALFRS score columns -- exclude Total and Cut_with (always missing)

#compute composite motor ALSFRS scores
ALSFRS$ALSFRS <- rowSums(as.matrix(ALSFRS[,cols_ALSFRS]))
ALSFRS$HandMotor <- ALSFRS$Handwriting + ALSFRS$Cut_without + ALSFRS$Dressing
ALSFRS$Other <- ALSFRS$ALSFRS - ALSFRS$HandMotor

#Plot ALSFRS over time
ALSFRS2 <- left_join(ALSFRS, scanDaysALS) #bring in days since first visit and since onset
pal_ALS <- c('black','grey',brewer.pal(9, 'Set1')[-6],brewer.pal(7, 'Set3')[-2]) #exclude yellows
pal_grp <- c(brewer.pal(3,'Paired'),'black')

ggplot(ALSFRS2, aes(x=days_since_onset/365*12, y=ALSFRS)) +
  geom_line(aes(color=subject)) + geom_point(aes(color=subject), fill='white', pch=21) +
  theme_bw() + guides(color=guide_legend(nrow=2)) +
  xlab('Months since Onset') + ylab('ALSFRS-R Score') +
  theme(legend.position='bottom', legend.title=element_blank()) +
  scale_color_manual(values=pal_ALS) + scale_x_continuous(breaks=seq(0,72,12)) +
  ggtitle('ALS Functional Ratings over Time')

# Supplementary Table A1

#first visit ALSFRS
ALSFRS2 %>% filter(days_since_visit1==0) %>% select(subject, visit, ALSFRS)
#last visit ALSFRS
tmp <- ALSFRS2 %>% group_by(subject) %>% summarize(lastvisit = max(visit))
ALSFRS2 <- left_join(ALSFRS2, tmp)
ALSFRS2 %>% filter(visit == lastvisit) %>% select(subject, visit, ALSFRS)

######################################################
# PROGRESSION RATE PLOTS ----

ALSFRS2_lastvisit <- ALSFRS2 %>% dplyr::filter(visit == lastvisit)
ALSFRS2_lastvisit$progression_rate <- (48 - ALSFRS2_lastvisit$ALSFRS)/(ALSFRS2_lastvisit$days_since_onset/365*12)
ALSFRS2_lastvisit <- select(ALSFRS2_lastvisit, subject, progression_rate)
ALSFRS2 <- left_join(ALSFRS2, ALSFRS2_lastvisit)
ALSFRS2$progressors <- 'moderate'
ALSFRS2$progressors[ALSFRS2$subject=='A04'] <- 'slow'
ALSFRS2$progressors[ALSFRS2$progression_rate > 0.7] <- 'fast'

ALSFRS3 <- filter(ALSFRS2, visit==lastvisit)
ALSFRS3 <- arrange(ALSFRS3, progression_rate)
ALSFRS3$subject <- factor(ALSFRS3$subject, levels=ALSFRS3$subject) #order subjects by progression rate

ggplot(ALSFRS2, aes(x=days_since_onset/365*12, y=ALSFRS, color=progressors)) +
  geom_line(aes(group=subject)) + geom_point(fill='white', pch=21) + theme_bw() +
  xlab('Months since Onset') + ylab('ALSFRS-R Score') +
  theme(legend.position='bottom') + scale_color_jco() +
  scale_x_continuous(breaks=seq(0,72,12)) +
  ggtitle('ALSFRS over Time')

ggplot(ALSFRS3, aes(x=subject, y=progression_rate, fill=progressors)) +
  #geom_hline(yintercept=c(0.1, 0.7), color='gray', linetype=2) +
  geom_bar(stat = 'identity') + theme_few() +
  xlab('Subject ID') + ylab('Progression Rate') +
  theme(legend.position='bottom') + scale_fill_jco(name='') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('')

range(ALSFRS2$progression_rate[ALSFRS2$progressors=='slow'], na.rm=TRUE)
range(ALSFRS2$progression_rate[ALSFRS2$progressors=='moderate'], na.rm=TRUE)
range(ALSFRS2$progression_rate[ALSFRS2$progressors=='fast'], na.rm=TRUE)

######################################################
# VALIDATION STUDY FIGURES (FIGURE E7) ----

load(file=file.path(result_dir,'activeareas_df.Rdata')) #active_areas_df, comptime_act_df
active_areas_df <- unique(active_areas_df)
active_areas_df$group <- substr(active_areas_df$subject, 1, 1) #A or C
active_areas_df$group <- factor(active_areas_df$group, levels=c('A','C'), labels=c('ALS','HC'))
active_areas_df$visit <- as.numeric(as.character(active_areas_df$visit))

active_areas_df <- dplyr::filter(active_areas_df, method %in% c('excur','FWER','FDR'))
active_areas_df$threshold2 <- active_areas_df$threshold
active_areas_df$threshold2[active_areas_df$method=='FWER'] <- '0_FWER'
active_areas_df$threshold2[active_areas_df$method=='FDR'] <- '0_FDR'
active_areas_df$threshold2 <- factor(active_areas_df$threshold2,
                                  levels=c('0_FWER','0_FDR','0','1','2'),
                                  labels=c('FWER','FDR','Bayesian (0%)','Bayesian (1%)','Bayesian (2%)'))

#longitudinal variance within HC subjects
active_areas_HC_var <- active_areas_df %>%
  filter(group == 'HC') %>%
  group_by(subject, hemisphere, method, threshold, threshold2) %>%
  summarize(variance_over_visits = var(area),
            mean_over_visits = mean(area),
            cv_over_visits = sqrt(variance_over_visits)/mean_over_visits,
            num_visits=n())
active_areas_HC_var_lh <- filter(active_areas_HC_var, hemisphere=='lh')

#variance BETWEEN HC subjects
var_between <- active_areas_HC_var_lh %>%
  group_by(threshold2) %>%
  summarize(var_between = var(mean_over_visits),
            mean_between = mean(mean_over_visits),
            cv_between = sqrt(var_between)/mean_between)

active_areas_HC_var_lh <- left_join(active_areas_HC_var_lh, var_between)

pal_thr <- c(rgb(128,0,38,maxColorValue=255), rgb(227,26,28,maxColorValue=255), rgb(254,178,76,maxColorValue=255))
pal_corr <- c('gray5','darkgray')
pal <- c(pal_corr, pal_thr)

ggplot(active_areas_HC_var_lh, aes(x=threshold2, y=cv_over_visits)) +
  geom_boxplot(aes(group=threshold2, fill=threshold2), alpha=0.8) +
  geom_point(data = var_between, aes(x=threshold2, y=cv_between, fill=threshold2), color='black', size=4, shape=23) +
  scale_fill_manual(values=pal) +
  ylab('CV of Activation Size') + xlab('') +
  ggtitle('CV of Activation Size in HC Subjects') +
  theme_few() + theme(legend.position='none', legend.title=element_blank())

ggplot(filter(active_areas_HC_var_lh, hemisphere=='lh'), aes(x=mean_over_visits, y=sqrt(variance_over_visits), group=threshold2)) +
  geom_point(aes(fill=threshold2), size=3, alpha=0.8, colour="black", pch=21) +
  geom_smooth(aes(color=threshold2), method='lm', se=FALSE, size=1.5) +
  geom_smooth(color='black', method='lm', se=FALSE, size=0.5) +
  scale_color_manual(values=pal) + scale_fill_manual(values=pal) +
  theme_few() + theme(legend.position='bottom', legend.title=element_blank()) +
  ylab('SD of Activation Size over Visits')  + xlab('Mean of Activation Size over Visits') +
  ggtitle('Variation in Activation Size within HC Subjects')


######################################################
# FIT LONGITUDINAL RANDOM EFFECTS MODELS ----

#Define BurdenOverall, BurdenHandMotor & BurdenOther
ALSFRS$BurdenOverall <- (12*4 - ALSFRS$ALSFRS)/(12*4)
ALSFRS$BurdenHandMotor <- (3*4 - ALSFRS$HandMotor)/(3*4)
ALSFRS$BurdenOther <- (9*4 - ALSFRS$Other)/(9*4) #everything but hand
quantile(ALSFRS$BurdenOverall, 0.9) #0.375
quantile(ALSFRS$BurdenHandMotor, 0.9) #0.5
quantile(ALSFRS$BurdenOther, 0.9) #0.361

active_areas_df2 <- dplyr::filter(active_areas_df, method=='excur', threshold %in% c(0,1,2))
active_areas_df2_classical <- dplyr::filter(active_areas_df, method %in% c('FWER','FDR'))
active_areas_df2 <- rbind(active_areas_df2, active_areas_df2_classical)
for(hem in c('lh','rh')){

  active_areas <- dplyr::filter(active_areas_df2, hemisphere==hem)#, method==mthd)
  active_areas <- left_join(active_areas, scanDays) #merge in days since onset/enrollment
  active_areas <- left_join(active_areas, select(max_days, subject, ALS_group)) #merge in ALS group ID
  active_areas$ALS_group[active_areas$group=='HC'] <- 'HC'
  active_areas_ALSFRS <- left_join(active_areas, select(ALSFRS, subject, visit, BurdenOverall, BurdenHandMotor, BurdenOther))
  
  ### MODEL SELECTION USING LRT ----
  
  #limit data to a 2-year window for each subject
  active_areas_ALSFRS <- filter(active_areas_ALSFRS, !(group=='HC' & days_since_visit1 > 365*2))
  active_areas_ALSFRS <- filter(active_areas_ALSFRS, !(subject=='A14' & visit == 1)) #filter out earlier visits to see greater change in ALSFRS

  #use left hemisphere, 0% threshold
  if(hem=='lh'){

    ### ALS
    dat_thr_ALS <- dplyr::filter(active_areas_ALSFRS, threshold==0, group=='ALS', method=='excur')
    lmerHandMotor_thr <- lmer(area ~ 1 + ns(BurdenHandMotor, df=3) + ns(BurdenOther, df=3) + (1 | subject ), data = dat_thr_ALS)
    lmerHandMotor2_thr <- lmer(area ~ 1 + BurdenHandMotor + BurdenOther + (1 | subject ), data = dat_thr_ALS)
    lmerHandMotor3_thr <- lmer(area ~ 1 + ns(BurdenHandMotor, df=3) + BurdenOther + (1 | subject ), data = dat_thr_ALS)
    lmerHandMotor4_thr <- lmer(area ~ 1 + ns(BurdenHandMotor, df=3) + (1 | subject ), data = dat_thr_ALS)
    #lmerHandMotor5_thr <- lmer(area ~ 1 + BurdenOther + (1 | subject ), data = dat_thr_ALS)
    anova(lmerHandMotor2_thr, lmerHandMotor3_thr) #test spline on BurdenHandMotor -- p=0.02211
    anova(lmerHandMotor_thr, lmerHandMotor3_thr) #test spline on BurdenOther -- p=0.7989
    anova(lmerHandMotor3_thr, lmerHandMotor4_thr) #test inclusion of BurdenOther -- p=0.09264

    ### ALS with Days Since Onset
    lmerHandMotor6_thr <- lmer(area ~ 1 + ns(days_since_onset, df=3) + (1 | subject ), data = dat_thr_ALS)
    lmerHandMotor7_thr <- lmer(area ~ 1 + days_since_onset + (1 | subject ), data = dat_thr_ALS)
    AIC(lmerHandMotor6_thr) #1208.877
    AIC(lmerHandMotor7_thr) #1250.991
    AIC(lmerHandMotor3_thr) #1185.976
    sqrt(mean(residuals(lmerHandMotor6_thr)^2)) #708
    sqrt(mean(residuals(lmerHandMotor7_thr)^2)) #711
    sqrt(mean(residuals(lmerHandMotor3_thr)^2)) #663

    ### HC
    dat_thr_HC <- dplyr::filter(active_areas_ALSFRS, threshold==0, group=='HC', method=='excur')
    lmerHC_thr <- lmer(area ~ 1 + days_since_visit1 + (1 | subject), data = dat_thr_HC)
    lmerHC2_thr <- lmer(area ~ 1 + ns(days_since_visit1, df=3) + (1 | subject), data = dat_thr_HC)
    anova(lmerHC2_thr, lmerHC_thr) #test spline on days_since_visit1 -- p=0.7816
    lmerHC3_thr <- lmer(area ~ 1 + (1 | subject), data = dat_thr_HC)
    anova(lmerHC3_thr, lmerHC_thr) #test slope on days_since_visit1 -- p=0.5275

  }

  ### DO ROBUSTNESS CHECKS

  groups <- list(
    'all' = 1:4, #all subjects together (including A04) -- for appendix Figure F8
    'excludeA04' = 1:3, #all subjects except A04 -- for Figure 4 in main text
    'fast' = 1, #fast progressors -- for Figure 5 and appendix
    'slow' = 2:3 #moderate progressors -- for Figure 5 and appendix
  )

  dat_pred_slowfast <- NULL
  for(g in 1:4){

    grp_name <- names(groups)[g]
    grp_inds <- groups[[g]]

    dat_pred_all <- dat_pred0_all <- NULL #for concatenating over thresholds
    thresholds <- unique(active_areas_ALSFRS$threshold2)
    HC_est <- HC_SE <- NULL
    for(thr in thresholds){

      ### ALS
      dat_thr_ALS <- dplyr::filter(active_areas_ALSFRS, threshold2==thr, group=='ALS')
      dat_thr_ALS <- dplyr::filter(dat_thr_ALS, ALS_group %in% paste0('group',grp_inds))
      lmerHandMotor_thr <- lmer(area ~ 1 + ns(BurdenHandMotor, df=3) + BurdenOther + (1 | subject), data = dat_thr_ALS)
      xvals <- seq(0,0.5,0.01) #c(0,0.1,0.2,0.3,0.4,0.5)
      dat_pred <- expand.grid(BurdenHandMotor = xvals, BurdenOther = xvals, days_since_visit1 = NA, threshold2=thr, group = 'ALS')
      mm <- model.matrix( ~ 1 + ns(BurdenHandMotor, df=3) + BurdenOther, dat_pred)
      dat_pred$Area_pred <- mm %*% fixef(lmerHandMotor_thr) # X * beta-hat
      dat_pred$SE_pred <- sqrt(diag(mm %*% tcrossprod(vcov(lmerHandMotor_thr),mm))) #Var(y-hat) = Var(X*beta-hat) = X*Var(beta-hat)*X'
      dat_pred$SE_lo <- dat_pred$Area_pred - dat_pred$SE_pred
      dat_pred$SE_hi <- dat_pred$Area_pred + dat_pred$SE_pred
      dat_pred$group <- 'ALS'
      dat_pred_all <- rbind(dat_pred_all, dat_pred)

      ### ALS with Overall Burden
      lmerOverall_thr <- lmer(area ~ 1 + ns(BurdenOverall, df=3) + (1 | subject), data = dat_thr_ALS)
      dat_pred0 <- expand.grid(BurdenOverall = xvals, days_since_visit1 = NA, threshold2=thr, group = 'ALS')
      mm <- model.matrix( ~ 1 + ns(BurdenOverall, df=3), dat_pred0)
      dat_pred0$Area_pred <- mm %*% fixef(lmerOverall_thr) # X * beta-hat
      dat_pred0$SE_pred <- sqrt(diag(mm %*% tcrossprod(vcov(lmerOverall_thr),mm))) #Var(y-hat) = Var(X*beta-hat) = X*Var(beta-hat)*X'
      dat_pred0$SE_lo <- dat_pred0$Area_pred - dat_pred0$SE_pred
      dat_pred0$SE_hi <- dat_pred0$Area_pred + dat_pred0$SE_pred
      dat_pred0$group <- 'ALS'
      dat_pred0_all <- rbind(dat_pred0_all, dat_pred0)

      ### HC model (intercept-only model)
      dat_thr_HC <- dplyr::filter(active_areas_ALSFRS, threshold2==thr, group=='HC')
      lmerHC0_thr <- lmer(area ~ 1 + (1 | subject), data = dat_thr_HC)
      HC_est = c(HC_est, fixef(lmerHC0_thr)[1])
      HC_SE = c(HC_SE, sqrt(vcov(lmerHC0_thr)[1,1]))
    } #end loop over thresholds

    df_HC <- data.frame(estimate = HC_est, SE = HC_SE, threshold = thresholds)

    #coefs_df$hem <- hem
    #coefs_df_all <- rbind(coefs_df_all, coefs_df)

    thr_long <- c('Bayesian (0%)', 'Bayesian (1%)', 'Bayesian (2%)', 'FWER', 'FDR')
    thr_short <- c('0%', '1%', '2%', 'FWER', 'FDR')
    dat_pred_all$group <- factor(dat_pred_all$group, levels=c('ALS','HC'))
    dat_pred0_all$group <- factor(dat_pred0_all$group, levels=c('ALS','HC'))
    dat_pred_all$threshold2 <- factor(dat_pred_all$threshold2, levels=thr_long, labels=thr_short)
    dat_pred0_all$threshold2 <- factor(dat_pred0_all$threshold2, levels=thr_long, labels=thr_short)
    df_HC$threshold <- factor(df_HC$threshold, levels=thr_long, labels=thr_short)

    dat_pred_all <- filter(dat_pred_all, (BurdenHandMotor==0 & BurdenOther <= 0.35) | BurdenOther==0)
    ymin <- min(dat_pred_all$SE_lo)
    ymax <- max(dat_pred_all$SE_hi)

    dat_pred_all$model_group <- dat_pred0_all$model_group <- grp_name
    if(g >= 3) dat_pred_slowfast <- rbind(dat_pred_slowfast, dat_pred_all)
    fname_suff <- ifelse(grp_name=='all', '', paste0('_',grp_name))

    title_hem <- ifelse(hem=='lh', "Contralateral (Left) Hemisphere", "Ipsilateral (Right) Hemisphere")

    # BAYESIAN GLM
    # xmax=0.5 is approximately equal to 90th quantile of Burden Hand Motor
    df_tmp <- filter(dat_pred_all, BurdenOther==0, group=='ALS', !(threshold2 %in% c('FWER','FDR')))
    #df_tmp <- filter(dat_pred_all, group=='ALS', !(threshold2 %in% c('FWER','FDR')))
    print(ggplot(df_tmp,  aes(x=BurdenHandMotor)) +
            #geom_point(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.02, y = estimate), size=3) +
            geom_text(data = NULL, aes(x=-0.05, y=ymax*0.95, label='HC'), size=4) +
            geom_text(data = NULL, aes(x=0.07, y=ymax*0.95, label='ALS'), size=4) +
            geom_point(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.05, y = estimate, fill=threshold), size=3, pch=21) +
            geom_errorbar(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.05, ymin = estimate - SE, ymax = estimate + SE), width=0) +
            #geom_line(aes(y=Area_pred, group=threshold2), color='black', size=1.2) +
            geom_line(aes(y=Area_pred, color=threshold2), size=1) +
            geom_ribbon(aes(ymin=(SE_lo), ymax=(SE_hi), group=threshold2), alpha=0.2) +
            #ylim(ymin,ymax) +
            scale_y_continuous(breaks=seq(0,4000,1000)) +
            scale_x_continuous(limits=c(-0.07, 0.5), breaks=seq(0,0.5,0.1)) +
            geom_vline(xintercept=0) +
            scale_color_manual('Effect Size', values=pal_thr) + scale_fill_manual('Effect Size', values=pal_thr) +
            ylab('Size of Activation') + xlab('Hand Motor Disability') + ggtitle(title_hem) +
            guides(fill=FALSE) + theme_few() + theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)))
    # xmax=0.35 is approximately equal to 90th quantile of Burden Other
    df_tmp <- filter(dat_pred_all, BurdenHandMotor==0, BurdenOther <= 0.35, group=='ALS', !(threshold2 %in% c('FWER','FDR')))
    print(ggplot(df_tmp,  aes(x=BurdenOther)) +
            #geom_point(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.02, y = estimate), size=3) +
            geom_text(data = NULL, aes(x=-0.035, y=ymax*0.95, label='HC'), size=4) +
            geom_text(data = NULL, aes(x=0.05, y=ymax*0.95, label='ALS'), size=4) +
            geom_point(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.035, y = estimate, fill=threshold), size=3, pch=21) +
            geom_errorbar(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.035, ymin = estimate - SE, ymax = estimate + SE), width=0) +
            #geom_line(aes(y=Area_pred, group=threshold2), color='black', size=1.2) +
            geom_line(aes(y=Area_pred, color=threshold2), size=1) +
            geom_ribbon(aes(ymin=(SE_lo), ymax=(SE_hi), group=threshold2), alpha=0.2) +
            #ylim(ymin,ymax) +
            scale_y_continuous(breaks=seq(0,4000,1000)) +
            scale_x_continuous(limits=c(-0.05, 0.35), breaks=seq(0,0.3,0.1)) +
            geom_vline(xintercept=0) +
            scale_color_manual('Effect Size', values=pal_thr) + scale_fill_manual('Effect Size', values=pal_thr) +
            ylab('Size of Activation') + xlab('Other Disability') + ggtitle(title_hem) +
            guides(fill=FALSE) +
            theme_few() + theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)))
    df_tmp <- filter(dat_pred0_all, !(threshold2 %in% c('FWER','FDR')))
    print(ggplot(df_tmp,  aes(x=BurdenOverall)) +
            #geom_point(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.02, y = estimate), size=3) +
            geom_text(data = NULL, aes(x=-0.05, y=ymax*0.95, label='HC'), size=4) +
            geom_text(data = NULL, aes(x=0.07, y=ymax*0.95, label='ALS'), size=4) +
            geom_point(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.05, y = estimate, fill=threshold), size=3, pch=21) +
            geom_errorbar(data = filter(df_HC, !(threshold %in% c('FWER','FDR'))), aes(x = -0.05, ymin = estimate - SE, ymax = estimate + SE), width=0) +
            #geom_line(aes(y=Area_pred, group=threshold2), color='black', size=1.2) +
            geom_line(aes(y=Area_pred, color=threshold2), size=1) +
            geom_ribbon(aes(ymin=(SE_lo), ymax=(SE_hi), group=threshold2), alpha=0.2) +
            scale_y_continuous(breaks=seq(0,4000,1000)) +
            scale_x_continuous(limits=c(-0.07, 0.5), breaks=seq(0,0.5,0.1)) +
            geom_vline(xintercept=0) +
            scale_color_manual('Effect Size', values=pal_thr) + scale_fill_manual('Effect Size', values=pal_thr) +
            ylab('Size of Activation') + xlab('Total Disability') + ggtitle(title_hem) +
            guides(fill=FALSE) + theme_few() + theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)))

  } #end loop over robustness checks
  
  #HERE

  dat_pred_slowfast <- filter(dat_pred_slowfast, group=='ALS', !(threshold2 %in% c('FWER','FDR')))
  dat_pred_slowfast$progressors <- factor(dat_pred_slowfast$model_group, levels=c('fast','slow'), labels=c('fast', 'moderate'))

  pdf(paste0('Plots_2021/lmer_activation_vs_diability_',hem,'_fastslow.pdf'), width=5, height=4)

  # BAYESIAN GLM
  # xmax=0.5 is approximately equal to 90th quantile of Burden Hand Motor
  df_tmp <- filter(dat_pred_slowfast, BurdenOther==0)
  print(ggplot(df_tmp,  aes(x=BurdenHandMotor)) +
          geom_line(aes(y=Area_pred, linetype=progressors), size=1) +
          geom_ribbon(aes(ymin=(SE_lo), ymax=(SE_hi), group=progressors), alpha=0.2) +
          facet_grid(. ~ threshold2) +
          ylab('Size of Activation') + xlab('Hand Motor Disability') +
          theme_few() + theme(legend.position='bottom'))
  # xmax=0.35 is approximately equal to 90th quantile of Burden Other
  df_tmp <- filter(dat_pred_slowfast, BurdenHandMotor==0, BurdenOther <= 0.35)
  print(ggplot(df_tmp,  aes(x=BurdenOther)) +
          geom_line(aes(y=Area_pred, linetype=progressors), size=1) +
          geom_ribbon(aes(ymin=(SE_lo), ymax=(SE_hi), group=progressors), alpha=0.2) +
          facet_grid(. ~ threshold2) +
          ylab('Size of Activation') + xlab('Other Disability') +
          theme_few() + theme(legend.position='bottom'))
  dev.off()

} #end loop over hemispheres


### PREDICT HAND MOTOR DISABILITY WITH ACTIVATIONS?

active_areas_df <- dplyr::filter(active_areas_df, method=='excur', threshold %in% c(0,1,2))
active_areas_ALSFRS <- left_join(filter(active_areas_df, group=='ALS'),
                                 select(ALSFRS, subject, visit, HandMotor))
active_areas_ALSFRS$BurdenHandMotor <- (3*4 - active_areas_ALSFRS$HandMotor)/(3*4)
active_areas_ALSFRS <- filter(active_areas_ALSFRS, !(subject=='A14' & visit == 1)) #filter out earlier visits to see greater change in ALSFRS
active_areas_ALSFRS <- dplyr::filter(active_areas_ALSFRS, subject!='A04')

ggplot(active_areas_ALSFRS, aes(x=area, y=BurdenHandMotor)) +
  geom_point() + geom_smooth() +
  facet_grid(hemisphere ~ threshold2, scales='free') + theme_few()
active_areas_ALSFRS <- filter(active_areas_ALSFRS, !(hemisphere=='rh' & threshold==1 & (area > 1500)))
active_areas_ALSFRS <- filter(active_areas_ALSFRS, !(hemisphere=='rh' & threshold==2 & (area > 750)))

dat_ALS_lh0 = select(filter(active_areas_ALSFRS, threshold==0, hemisphere=='lh'), subject, visit, BurdenHandMotor, area)
dat_ALS_lh1 = select(filter(active_areas_ALSFRS, threshold==1, hemisphere=='lh'), subject, visit, BurdenHandMotor, area)
dat_ALS_lh2 = select(filter(active_areas_ALSFRS, threshold==2, hemisphere=='lh'), subject, visit, BurdenHandMotor, area)
#dat_ALS_rh0 = select(filter(active_areas_ALSFRS, threshold==0, hemisphere=='rh'), subject, visit, BurdenHandMotor, area)
#dat_ALS_rh1 = select(filter(active_areas_ALSFRS, threshold==1, hemisphere=='rh'), subject, visit, BurdenHandMotor, area)
#dat_ALS_rh2 = select(filter(active_areas_ALSFRS, threshold==2, hemisphere=='rh'), subject, visit, BurdenHandMotor, area)
names(dat_ALS_lh0)[4] <- 'area_lh0'
names(dat_ALS_lh1)[4] <- 'area_lh1'
names(dat_ALS_lh2)[4] <- 'area_lh2'
#names(dat_ALS_rh0)[4] <- 'area_rh0'
#names(dat_ALS_rh1)[4] <- 'area_rh1'
dat_ALS_lh <- inner_join(dat_ALS_lh0, dat_ALS_lh1, by=c('subject','visit','BurdenHandMotor'))
dat_ALS_lh <- inner_join(dat_ALS_lh, dat_ALS_lh2, by=c('subject','visit','BurdenHandMotor'))
#dat_ALS_rh <- inner_join(dat_ALS_rh0, dat_ALS_rh1, by=c('subject','visit','BurdenHandMotor'))
#dat_ALS <- inner_join(dat_ALS_lh, dat_ALS_rh, by=c('subject','visit','BurdenHandMotor'))

#lmer(area ~ 1 + days_since_visit1 + (1 | subject), data = dat_thr_HC)
predALS <- lm(BurdenHandMotor ~ 1 + ns(area_lh0, df=3) + ns(area_lh1, df=3) + ns(area_lh2, df=3), data=dat_ALS_lh)
subjects_ALS <- unique(dat_ALS_lh$subject)
dat_pred <- NULL
for(ss in subjects_ALS){
  predALS_s <- lm(BurdenHandMotor ~ 1 + ns(area_lh0, df=3)*ns(area_lh1, df=3), data=filter(dat_ALS_lh, subject != ss))
  dat_ALS_lh_s <- filter(dat_ALS_lh, subject == ss)
  dat_ALS_lh_s$pred <- predict(predALS_s, newdata=dat_ALS_lh_s)
  dat_pred <- rbind(dat_pred, dat_ALS_lh_s)
}

cor_pred <- cor(dat_pred$BurdenHandMotor, dat_pred$pred)
cor_pred_text <- paste0('cor = ',round(cor_pred,3))
ggplot(dat_pred, aes(x=BurdenHandMotor, y=pred)) +
  geom_point(aes(color=subject)) + geom_smooth(method='lm', se=FALSE) +
  annotate("text", x = 0.1, y = 0.6, label = cor_pred_text) +
  scale_color_manual(values=pal_ALS) +
  theme_few() + xlab('True Hand Motor Disability') + ylab('Predicted Hand Motor Disability')
#
#
# coefs_df_all$CI_lwr <- coefs_df_all$est - 1.96*coefs_df_all$SE
# coefs_df_all$CI_upr <- coefs_df_all$est + 1.96*coefs_df_all$SE
# coefs_df_all$tval <- coefs_df_all$est/coefs_df_all$SE
# coefs_df_all$pval <- 2*pnorm(abs(coefs_df_all$tval), lower.tail=FALSE)
#
# filter(coefs_df_all, group=='HC')
# filter(coefs_df_all, group=='ALS')

######################################################################################################
### COMPUTATION TIME
######################################################################################################

### TIME TO IDENTIFY ACTIVATIONS

#sum over hemispheres
comptime_act_df2 <- comptime_act_df %>% group_by(subject, visit, threshold) %>% summarize(comptime = sum(comptime))
#average over visits
comptime_act_df3 <- comptime_act_df2 %>% group_by(subject, threshold) %>% summarize(comptime_avg = mean(comptime), num_visits=n())
comptime_act_df3$threshold <- factor(comptime_act_df2$threshold, levels=0:2, labels=paste0(0:2,'%'))


### TIME FOR MODEL ESTIMATION

load(file=file.path(main_dir,'Results_2021','comptime.RData')) #comptime_df -- model estimation time
#sum over hemispheres
comptime_df2 <- comptime_df %>% group_by(subject, num_visits) %>% summarize(comptime = sum(comptime))


pdf(file.path(main_dir,'Plots_2021/comptime.pdf'), width=4, height=4)
ggplot(comptime_df2, aes(x=num_visits, y=comptime/60)) +
  geom_point() + geom_smooth() +
  #geom_line(aes(color=subject, group=interaction(subject, hemisphere))) +
  #geom_boxplot(aes(group=num_visits), alpha=0.3) +
  theme_few() + ggtitle('Computation Time to Estimate Model') +
  ylab('Minutes') + xlab('Number of Sessions in Model')
ggplot(comptime_act_df3, aes(x=num_visits, y=comptime_avg/60)) +
  geom_point() + geom_smooth() + #geom_smooth(aes(group=threshold, color=threshold)) +
  #geom_line(aes(group=subject, color=num_visits)) +
  #geom_violin(aes(group=threshold), alpha=0.5) +
  #scale_color_viridis_c() +
  theme_few() + #theme(legend.position = 'bottom') +
  ggtitle('Computation Time to Identify Activations') +
  ylab('Minutes (per Session & Threshold)') + xlab('Number of Sessions in Model')
dev.off()





