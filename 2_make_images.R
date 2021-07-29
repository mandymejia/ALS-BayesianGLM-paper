######################################################
# INSTALL AND LOAD DEPENDENCIES ----

#CRAN packages --- install each with install.packages()
library(gifti) #readgii

#install_github('mandymejia/ciftiTools', ref='2.1') #1.2
library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/')


######################################################
# SETUP ----

main_dir <- '/my/project/dir' #path to repo contents
data_dir <- file.path(main_dir,'data')
result_dir <- file.path(main_dir,'results') 
setwd(data_dir)
load(file='subjects.Rdata') #subjects
task <- 'RRC'

s <- 'C07' # example subject data in repo
#for(s in subjects){
  isubj <- which(subjects==s)
  setwd(file.path(data_dir,s))

	#GET NUMBER OF VISITS
  load(file=file.path(data_dir,s,paste0('visits_',s,'.RData'))) #visits_s, out_s, exclude
  n_visits <- length(visits_s)
  #exclude bad visits
  if(!is.null(exclude)){
    visits_s <- visits_s[-exclude]
    out_s <- out_s[-exclude]
    n_visits <- length(visits_s)
  }

  #READ IN INFLATED SURFACE FOR VISUALIZATION
  setwd(file.path(data_dir,s))
  surfs_s <- vector('list', length=2)
  for(h in 1:2){
    hem <- c('lh','rh')[h]
    surf_fname <- paste0(s,'_',hem,'.inflated.10K.surf.gii')
    surf_sh <- readgii(surf_fname)$data
    vertices_sh <- surf_sh$pointset
    faces_sh <- (surf_sh$triangle + 1)
    surfs_s[[h]] <- make_surf(list(pointset=vertices_sh, triangle=faces_sh))
  }

  ######################################################
  # VISUALIZE ESTIMATES OF ACTIVATION ----

  load(file.path(result_dir,paste0('GLM_classical_',s,'.RData')))
  load(file.path(result_dir,paste0('GLM_Bayesian_',s,'.RData')))
  
  #make a matrix with beta estimates for RRC task, with visits in columns of matrix
  estimates0L <- sapply(results0_s[[1]], function(x) return(x$estimates[,1])) #classical GLM
  estimates0R <- sapply(results0_s[[2]], function(x) return(x$estimates[,1])) #classical GLM
  estimatesL <- sapply(results_s[[1]]$beta_estimates, function(x) return(x[,1])) #Bayesian GLM
  estimatesR <- sapply(results_s[[2]]$beta_estimates, function(x) return(x[,1])) #Bayesian GLM

  #classical GLM
  xifti_betas0_s <- as.xifti(cortexL = estimates0L,
                             surfL = surfs_s[[1]],
                             cortexR = estimates0R,
                             surfR = surfs_s[[2]],
                             HCP_32k_auto_mwall = FALSE,
                             validate = FALSE)

  #dir.create(file.path(main_dir,'images',s))
  png_fname = file.path(main_dir,'images',s,paste0(s,'_',visits_s,'_',task,'_classical.png'))
  plot(xifti_betas0_s, idx=1:n_visits, zlim=c(-2,2), fname=png_fname, title=paste0('visit ',visits_s,' (Classical GLM)'))

  #Bayesian GLM
  xifti_betas_s <- as.xifti(cortexL = estimatesL,
                            surfL = surfs_s[[1]],
                            cortexR = estimatesR,
                            surfR = surfs_s[[2]],
                            HCP_32k_auto_mwall = FALSE,
                            validate = FALSE)

  png_fname = file.path(main_dir,'images',s,paste0(s,'_',visits_s,'_',task,'_Bayes.png'))
  plot(xifti_betas_s, idx=1:n_visits, zlim=c(-2,2), fname=png_fname, title=paste0('visit ',visits_s,' (Bayesian GLM)'))


  ######################################################
  ### VISUALIZE AREAS OF ACTIVATION ----

  pal_thr <- c(rgb(128,0,38,maxColorValue=255), rgb(227,26,28,maxColorValue=255), rgb(254,178,76,maxColorValue=255))
  pal_corr <- c('gray5','darkgray')
  #pal_act <- c('red3','orange','yellow')
  #pal_act <- brewer.pal(9, 'YlOrRd')

  act_Bayes <- vector('list', length=2)
  for(h in 1:2){
    hem <- c('lh','rh')[h]

    #classical
    load(file=file.path(result_dir,paste0('act0_',s,'_',hem,'.RData'))) #act0_combined #0 = not significant, 1 = significant with no correction, 2 = survives FDR correction, 3 = survives FWER correction
    act0_combined[act0_combined <= 1] <- 0 #remove voxels that do not survive FDR or FWER correction
    act0_combined[act0_combined == 2] <- 1 #FDR
    act0_combined[act0_combined == 3] <- 2 #FWER
    act_classical[[h]] <- act0_combined

    #Bayesian
    mask_sh_used <- results_s[[h]]$mask
    act_Bayes_h <- matrix(0, length(mask_sh_used), n_visits)
    for(v in 1:n_visits){
      fname <- paste0('excur_',s,'_',visits_s[v],'_',hem,'.RData')
      load(file=file.path(result_dir,fname)) #excur_sv
      act_Bayes_h[mask_sh_used,v] <- excur_sv[[1]]$active + excur_sv[[2]]$active + excur_sv[[3]]$active #combine activations over 0, 1 and 2% signal change
    }
    act_Bayes[[h]] <- act_Bayes_h
  }

  xifti_act_s <- as.xifti(cortexL = act_classical[[1]],
                             surfL = surfs_s[[1]],
                             cortexR = act_classical[[2]],
                             surfR = surfs_s[[2]],
                             HCP_32k_auto_mwall = FALSE,
                             validate = FALSE)
  png_fname = file.path(main_dir,'images',s,paste0(s,'_',visits_s,'_',task,'_act_classical.png'))
  plot(xifti_act_s, idx=1:n_visits, colors=c('white',rev(pal_corr)), color_mode='qual', fname=png_fname, title=paste0('visit ',visits_s,' (classical GLM)'))

  xifti_act_s <- as.xifti(cortexL = act_Bayes[[1]],
                           surfL = surfs_s[[1]],
                           cortexR = act_Bayes[[2]],
                           surfR = surfs_s[[2]],
                           HCP_32k_auto_mwall = FALSE,
                           validate = FALSE)
  png_fname = file.path(main_dir,'images',s,paste0(s,'_',visits_s,'_',task,'_act_Bayes.png'))
  plot(xifti_act_s, idx=1:n_visits, colors=c('white',pal_thr), color_mode='qual', fname=png_fname, title=paste0('visit ',visits_s,' (Bayesian GLM)'))

#}

