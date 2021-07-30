######################################################
# INSTALL AND LOAD DEPENDENCIES ----

## fMRIscrub package
#library(devtools)
#install_github('mandymejia/fMRIscrub', ref='3.0') 
library(fMRIscrub) #clever()

## CRAN packages --- install each with install.packages()
library(gifti) #readgii()
library(ggplot2)
library(ggthemes)

## BayesfMRI package + dependencies (ciftiTools + INLA + INLA dependencies)
#library(devtools)
#install_github('mandymejia/ciftiTools', ref='2.1')
#install.packages('sp') #necessary INLA dependency
#install.packages('foreach') #necessary INLA dependency
#install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=FALSE)
library(INLA)
#inla.pardiso() #download PARDISO license for fast computation
inla.setOption(pardiso.license = '~/pardiso.lic')
inla.pardiso.check() #verify PARDISO license
#install_github('mandymejia/BayesfMRI', ref='1.8')
library(BayesfMRI)

######################################################
# SETUP ----

main_dir <- '/my/project/dir' #path to repo contents
data_dir <- file.path(main_dir,'data')
result_dir <- file.path(main_dir,'results') 
setwd(data_dir)
load(file='subjects.Rdata') #subjects
groups <- ifelse(grepl('A', subjects), 'ALS', 'HC')
table(groups)
task <- 'RRC' #analyze "real" right hand clenching task

#read in design and 1st derivative
X <- as.matrix(read.delim('motor.mat', header=FALSE))[,1:2] #third col is NA
design_long <- data.frame(column = rep(c('HRF','dHRF'),each=nrow(X)),
                          volume = rep(1:nrow(X), 2),
                          value = c(X[,1],X[,2]))
design_long$column <- factor(design_long$column, levels=c('HRF','dHRF'), labels=c('HRF','HRF Derivative'))

### SUPPLEMENTARY FIGURE

ggplot(design_long, aes(x=volume, y=value, group=column, linetype=column)) +
  geom_line() + theme_few() + theme(legend.position = 'bottom', legend.title = element_blank()) +
  scale_x_continuous(breaks=seq(0, 140, 40)) + ylab('Arbitrary Units') + xlab('Volume Index')

######################################################
### RUN BAYESIAN & CLASSICAL MODELS 

thresholds <- c(0,1,2) #activation thresholds (percent signal change)
U <- length(thresholds)
 
# comptime_df <- comptime_act_df <- NULL
# active_areas_df <- NULL 
# n_visits_all <- n_visits_all2 <- rep(NA, length(subjects)) #before and after QC

s <- 'C07' # example subject data in repo
#for(s in subjects){
  isubj <- which(subjects==s)
  setwd(file.path(data_dir,s))

  ### GET LIST OF VISITS
  
  files_s <- list.files(pattern = paste('*',task,'lh','fmri.csv',sep='_'))
  visits_s <- sort(substr(files_s, start=5, stop=6)) #get visit numbers
  #limit max number of visits to analyze
  if(length(visits_s) > 10){
    which.visits <- round(seq(1, length(visits_s), length.out=10)) #roughly equally spaced, including first and last visit
    visits_s <- visits_s[which.visits]
  }
  n_visits <- length(visits_s)
  #n_visits_all[isubj] <- n_visits

  exclude <- NULL #any visits to exclude due to excessive outliers
  out_s <- vector('list', length=n_visits)
  for(v in 1:n_visits){
    #read in data
    setwd(file.path(data_dir,s))
    fname_BOLD_sv_lh <- paste(s,visits_s[v],task,'lh','fmri.csv',sep='_')
    fname_BOLD_sv_rh <- paste(s,visits_s[v],task,'rh','fmri.csv',sep='_')
    if(!(fname_BOLD_sv_lh %in% list.files())) break
    if(!(fname_BOLD_sv_rh %in% list.files())) break
    BOLD_sv_lh <- read.csv(fname_BOLD_sv_lh, header=FALSE)[,1:136]
    BOLD_sv_rh <- read.csv(fname_BOLD_sv_rh, header=FALSE)[,1:136]
    BOLD_sv <- rbind(BOLD_sv_lh, BOLD_sv_rh)

    #compute PCA leverage for scrubbing and QC
    PClev_sv <- clever(t(BOLD_sv), measures='leverage')
    out_s[[v]] <- PClev_sv$outlier_flags$leverage__PCA_kurt
    pct_out <- mean(out_s[[v]])
    if(pct_out > 0.25){ exclude <- c(exclude, v); next } #if more than 1/4 outliers, exclude visit (after loop)
  }
  save(visits_s, out_s, exclude, file=paste0('visits_',s,'.RData'))
  #load(file=paste0('visits_',s,'.RData'))
  n_visits <- length(visits_s)
  #n_visits_all[isubj] <- n_visits

  #exclude bad visits
  if(!is.null(exclude)){
    visits_s <- visits_s[-exclude]
    out_s <- out_s[-exclude]
    n_visits <- length(visits_s)
  }
  #n_visits_all2[isubj] <- n_visits

  #loop over hemispheres
  setwd(file.path(data_dir,s))
  results_s <- vector('list', length=2) #save Bayesian GLM results across hemispheres
  results0_s <- vector('list', length=2) #save classical GLM results across hemispheres
  for(h in 1:2){

    hem <- c('lh','rh')[h]
    cat(paste0('~~~~~~~~~~~~~~~~~~~ hemisphere ',hem,' ~~~~~~~~~~~~~~~~~~~\n'))

    #READ IN PIAL SURFACE FOR SPATIAL MODELING
    surf_fname <- paste0(s,'_',hem,'.pial.10K.surf.gii')
    surf_sh <- readgii(surf_fname)$data
    vertices_sh <- surf_sh$pointset
    faces_sh <- (surf_sh$triangle + 1)
    rm(surf_sh)

    ### READ IN MOTOR MASK
    mask_fname <- paste0(s,'_',hem,'.aparc.ascii.10K.motor.label.csv')
    mask_sh <- (as.matrix(read.csv(mask_fname, header=FALSE))==1)

    ### ORGANIZE DATA IN "SESSIONS" FORMAT
    data_sh <- vector('list', n_visits)
    names(data_sh) <- visits_s
    for(v in 1:n_visits){

      print(paste0('visit ',v,' of ',n_visits))

      #remove outliers from design
      keep_sv <- !out_s[[v]]
      design_sv <- X[keep_sv,]

      #create nuisance regressors & remove outliers
      fname_motion_sv <- paste(s,visits_s[[v]],task,'motion.txt',sep='_')
      motion_sv <- as.matrix(read.table(fname_motion_sv, header=FALSE))
      dmotion_sv <- apply(motion_sv, 2, gradient)[1:136,] #remove extra volumes collected for one subject
      motion_sv <- motion_sv[1:136,] #remove extra volumes collected for one subject
      trend <- (1:136)/136 #for linear and quadratic time trends
      Z <- cbind(motion_sv, dmotion_sv, trend, trend^2)
      nuisance_sv <- Z[keep_sv,]

      #read BOLD data & remove outliers
      fname_BOLD_svh <- paste(s,visits_s[v],task,hem,'fmri.csv',sep='_')
      BOLD_svh <- read.csv(fname_BOLD_svh, header=FALSE)[,1:136]
      BOLD_svh <- BOLD_svh[,keep_sv] #remove outlier volumes

      ### ORGANIZE DATA INTO SESSIONS
      data_svh <- list(BOLD = t(BOLD_svh), design=design_sv, nuisance=nuisance_sv)
      if(!is.session(data_svh)) break
      data_sh[[v]] <- data_svh
    } #end loop over visits

    ### ESTIMATE CLASSICAL AND BAYESIAN GLM
    results0_s[[h]] <- classicalGLM(data=data_sh, mask=mask_sh, avg_sessions = FALSE, scale_design = FALSE)
    tmp <- system.time( #5-10 minutes per hemisphere
      results_s[[h]] <- BayesGLM(data=data_sh,
                                 beta_names = c('RRC','dRRC'),
                                 vertices=vertices_sh,
                                 faces=faces_sh,
                                 mask=mask_sh,
                                 num.threads = 4,
                                 avg_sessions = FALSE,
                                 scale_design=FALSE))

    #comptime_df_ii <- data.frame(subject=s, hemisphere=hem, num_visits=n_visits, comptime=tmp[3])  #elapsed time (sec)
    #comptime_df <- rbind(comptime_df, comptime_df_ii)
    #save(comptime_df, file=file.path(result_dir,'comptime.RData'))

    ### IDENTIFY AREAS OF ACTIVATION

    #compute size of areas of activation
    mask_sh_used <- results_s[[h]]$mask #for Bayesian GLM
    mask_sh_used0 <- results0_s[[h]][[1]]$mask #for classical GLM
    mesh_sh <- results_s[[h]]$mesh
    mesh_sh <- inla.mesh.create(loc = as.matrix(vertices_sh), tv = as.matrix(faces_sh))
    areas_all_sh <- diag(inla.fmesher.smorg(mesh_sh$loc,mesh_sh$graph$tv, fem = 0, output = list("c0"))$c0)
    areas_sh <- areas_all_sh[mask_sh_used]
    areas_sh0 <- areas_all_sh[mask_sh_used0]
    
    alpha <- 0.05

    #classical GLM activations
    act0_combined <- matrix(0, length(mask_sh_used0), n_visits)
    for(v in 1:n_visits){
      act0_FWER <- id_activations.classical(model_obj=results0_s[[h]], field_inds = 1, session_name = visits_s[v], threshold=0, alpha=alpha, correction='FWER')
      act0_FDR <- id_activations.classical(model_obj=results0_s[[h]], field_inds = 1, session_name = visits_s[v], threshold=0, alpha=alpha, correction='FDR')
      act0_nocorr <- id_activations.classical(model_obj=results0_s[[h]], field_inds = 1, session_name = visits_s[v], threshold=0, alpha=alpha, correction='none')
      act0_combined[mask_sh_used0,v] <- act0_nocorr$active + act0_FDR$active + act0_FWER$active
    }
    save(act0_combined, file=file.path(result_dir, paste0('act0_',s,'_',hem,'.RData')))

    #Bayesian GLM activations (~2 min per visit)
    for(v in 1:n_visits){
      print(paste0('visit ',v,' of ',n_visits))
      excur_sv <- vector('list', length=U)
      for(u in 1:U){
        thr_u <- thresholds[u]
        time_u <- system.time(
          excur_sv[[u]] <- id_activations.posterior(model_obj=results_s[[h]],
                                                    field_names = task,
                                                    session_name = visits_s[v],
                                                    threshold=thr_u,
                                                    alpha=alpha))
        area_svu <- sum(areas_sh[excur_sv[[u]]$active==1]) #compute size of active area
        # active_areas_df_ii <- data.frame(subject=s,
        #                                  hemisphere=hem,
        #                                  visit=visits_s[v],
        #                                  method='excur',
        #                                  threshold=thr_u,
        #                                  area = area_svu)
        # active_areas_df <- rbind(active_areas_df, active_areas_df_ii)

        # comptime_act_df_ii <- data.frame(subject=s,
        #                                    hemisphere=hem,
        #                                    visit=visits_s[v],
        #                                    threshold=thr_u,
        #                                    comptime=time_u[3])  #elapsed time (sec)
        # comptime_act_df <- rbind(comptime_act_df, comptime_act_df_ii)
      } #loop over thresholds

      fname <- paste0('excur_',s,'_',visits_s[v],'_',hem,'.RData')
      save(excur_sv, file=file.path(result_dir,fname))
    } #loop over visits
    #save(active_areas_df, comptime_act_df, file=file.path(result_dir,'activeareas_df.Rdata'))
 
  } #loop over hemispheres
  
  save(results0_s, file=file.path(result_dir, paste0('GLM_classical_',s,'.RData')))
  save(results_s, file=file.path(result_dir, paste0('GLM_Bayesian_',s,'.RData')))
#}

# sum(n_visits_all[groups=='ALS'])
# sum(n_visits_all[groups=='HC'])
# sum(n_visits_all2[groups=='ALS'])
# sum(n_visits_all2[groups=='HC'])

