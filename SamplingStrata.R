#------------------------------------------------------
# Performance comparison between current SamplingStrata
# and SamplingStrata with Rcpp modules:
# - aggrStrata_RcppOpen
# - Bethel_Rcpp
# - stdev_Rcpp
# N.B.: does not work with Bethel_RcppOpen
#------------------------------------------------------


#------------------------------------------------------
# Current SamplingStrata
#------------------------------------------------------
install.packages("SamplingStrata")
library(SamplingStrata)


#-------------------------
# Test data (all regions)
#-------------------------
data(swissmunicipalities)
swissmun <- swissmunicipalities[,c("REG","COM","Nom","HApoly",
                                  "Surfacesbois","Surfacescult",
                                  "Airbat","POPTOT")]
swissmun$HApoly.cat <- var.bin(swissmun$HApoly,15)
table(swissmun$HApoly.cat)
swissmun$POPTOT.cat <- var.bin(swissmun$POPTOT,15)
table(swissmun$POPTOT.cat)

#------------------------------
# Preparation for atomic method
#------------------------------
frame1 <- buildFrameDF(df = swissmun,
                           id = "COM",
                           X = c("POPTOT.cat","HApoly.cat"),
                           Y = c("Airbat","Surfacesbois"),
                           domainvalue = "REG")

strata1 <- buildStrataDF(frame1, progress=F)
head(strata1)

#--------------------------------------------------------------------
# Preparation for continuous method (to be used with bethel_RcppOpen)
#--------------------------------------------------------------------
# swissmun$progr <- c(1:nrow(swissmun))
# frame2 <- buildFrameDF(df = swissmun,
#                        id = "COM",
#                        X = c("POPTOT","HApoly"),
#                        Y = c("Airbat","Surfacesbois"),
#                        domainvalue = "REG")
# set.seed(1234)
# init_sol2 <- KmeansSolution2(frame=frame2,
#                              errors=cv,
#                              maxclusters = 10)  
# nstrata2 <- tapply(init_sol2$suggestions,
#                    init_sol2$domainvalue,
#                    FUN=function(x) length(unique(x)))
# nstrata2


#----------------------------------
# Precision constraints
#----------------------------------

ndom <- length(unique(swissmun$REG))
cv <- as.data.frame(list(DOM=rep("DOM1",ndom),
                         CV1=rep(0.1,ndom),
                         CV2=rep(0.1,ndom),
                         domainvalue=c(1:ndom) ))
cv

checkInput(errors = checkInput(errors = cv, 
                               strata = strata1, 
                               sampframe = frame1))


#----------------------------------------
# Performance with current SamplingStrata 
#----------------------------------------

#----------------------------------------
# Performance atomic method no parallel 
#----------------------------------------
set.seed(1234)
system.time(
solution1 <- optimStrata(method = "atomic",
                        errors = cv, 
                        nStrata = rep(10,ndom),
                        framesamp = frame1,
                        iter = 50,
                        pops = 10,
                        parallel=F)
)
# utente   sistema trascorso 
#  63.22      5.93     72.83

#----------------------------------------
# Performance atomic method parallel 
#----------------------------------------
system.time(
  solution1 <- optimStrata(method = "atomic",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame1,
                           iter = 50,
                           pops = 10,
                           parallel=T)
)
# *** Starting parallel optimization for  7  domains using  5  cores
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=17s  
# 
# *** Sample size :  498
# *** Number of strata :  54
# ---------------------------   
# utente   sistema trascorso 
#   0.50      0.03     18.49 

sum(solution1$aggr_strata$SOLUZ)
# [1] 502.3473
expected_CV(solution1$aggr_strata)
# cv(Y1)    cv(Y2)
# DOM1 0.0985073 0.1000000
# DOM2 0.1000000 0.1000000
# DOM3 0.0974722 0.0977469
# DOM4 0.0999351 0.0967996
# DOM5 0.1000000 0.1000000
# DOM6 0.0992327 0.0990466
# DOM7 0.0974012 0.0998922

#--------------------------------------
# Performance with SamplingStrata Rcpp
#--------------------------------------
detach("package:SamplingStrata", unload = TRUE)
install.packages("D:/Google Drive/Sampling/SamplingStrata C++/SamplingStrata_1.5-4.tar.gz", repos = NULL, type = "source")
library(SamplingStrata)


#----------------------------------------
# Performance atomic method no parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution2 <- optimStrata(method = "atomic",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame1,
                           iter = 50,
                           pops = 10,
                           parallel=F)
)
# utente   sistema trascorso 
# 25.97      2.55     30.51 

#----------------------------------------
# Performance atomic method parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution2 <- optimStrata(method = "atomic",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame1,
                           iter = 50,
                           pops = 10,
                           parallel=T)
)
# *** Starting parallel optimization for  7  domains using  5  cores
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13s  
# 
# *** Sample size :  508
# *** Number of strata :  57
# ---------------------------   
# utente   sistema trascorso 
#   0.38      0.02     13.04

sum(solution2$aggr_strata$SOLUZ)
# [1] 507.6798
expected_CV(solution2$aggr_strata)
#         cv(Y1)    cv(Y2)
# DOM1 0.1000000 0.1000000
# DOM2 0.0998728 0.0998161
# DOM3 0.1000000 0.1000000
# DOM4 0.0974066 0.0996974
# DOM5 0.0993837 0.0997432
# DOM6 0.1000000 0.1000000
# DOM7 0.0996361 0.0988654

