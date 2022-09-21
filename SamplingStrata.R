#------------------------------------------------------
# Performance comparison between current SamplingStrata
# and SamplingStrata with Rcpp modules:
# - aggrStrata_RcppOpen
# - Bethel_Rcpp
# - stdev_Rcpp
# N.B.: bethelRcppOpen gives wrong results
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
swissmun <- swissmunicipalities[swissmunicipalities$REG < 3,c("REG","COM","Nom","HApoly",
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

#------------------------------
# Preparation for continuous method
#------------------------------
frame2 <- buildFrameDF(df = swissmun,
                       id = "COM",
                       X = c("Airbat","Surfacesbois"),
                       Y = c("Airbat","Surfacesbois"),
                       domainvalue = "REG")


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

all1 <- bethel(strata1,cv[1,])
sum(all1)
# [1] 482

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
#  16.18      1.58     18.92 

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
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=17s  
# 
# *** Sample size :  498
# *** Number of strata :  54
# ---------------------------   
# utente   sistema trascorso 
#   0.18      0.03      8.67 

sum(solution2$aggr_strata$SOLUZ)
# [1] 265.4858
expected_CV(solution2$aggr_strata)
# cv(Y1)    cv(Y2)
# DOM1 0.0966098 0.1000000
# DOM2 0.0967295 0.0997367


#----------------------------------------
# Performance continuous method no parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution3 <- optimStrata(method = "continuous",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame2,
                           iter = 50,
                           pops = 10,
                           parallel=F)
)
# utente   sistema trascorso 
# 25.45      1.72     28.14 

outstrata <- solution3$aggr_strata
sum(outstrata$SOLUZ)
# [1] 52.8884
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0989425 0.096506
# DOM2 0.0935584 0.098520

#----------------------------------------
# Performance continuous method parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution4 <- optimStrata(method = "continuous",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame2,
                           iter = 50,
                           pops = 10,
                           parallel=T)
)
# *** Starting parallel optimization for  7  domains using  5  cores
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13s  
# 
# *** Sample size :  51
# *** Number of strata :  16
# ---------------------------   
#   utente   sistema trascorso 
# 0.00      0.03     11.00 

outstrata <- solution4$aggr_strata
sum(outstrata$SOLUZ)
# [1] 51.26343
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0996368 0.0996492
# DOM2 0.0979906 0.0986974


#--------------------------------------
# Performance with SamplingStrata Rcpp
#--------------------------------------
detach("package:SamplingStrata", unload = TRUE)
# Better to restart R
install.packages("D://Google Drive//Sampling//SamplingStrata-C-//SamplingStrata_1.5-4.tar.gz", repos = NULL, type = "source")

library(SamplingStrata)

all2 <- bethelRcppOpen(strata1,cv[1,],realAllocation = TRUE)
sum(all2)
# [1] 482

#----------------------------------------
# Performance atomic method no parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution5 <- optimStrata(method = "atomic",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame1,
                           iter = 50,
                           pops = 10,
                           parallel=F)
)
# utente   sistema trascorso 
# 3.03      0.16      3.84  

outstrata <- solution5$aggr_strata
colnames(outstrata)[10] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1]  143.4412
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.2925264 0.3065943
# DOM2 0.0999229 0.0999656

#----------------------------------------
# Performance atomic method parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution6 <- optimStrata(method = "atomic",
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
# *** Sample size :  0
# *** Number of strata :  42
# ---------------------------   
#   utente   sistema trascorso 
#  0.21      0.03      1.88   

outstrata <- solution6$aggr_strata
colnames(outstrata)[10] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1] 173.5117
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.2422364 0.2703787
# DOM2 0.0980734 0.0997199

#----------------------------------------
# Performance continuous method no parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution7 <- optimStrata(method = "continuous",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame2,
                           iter = 50,
                           pops = 10,
                           parallel=F)
)
# utente   sistema trascorso 
# 19.44      1.15     21.50  

outstrata <- solution7$aggr_strata
colnames(outstrata)[ncol(outstrata)] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1] 45.16294
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.3635650 0.3872813
# DOM2 0.0957606 0.0990580

#----------------------------------------
# Performance continuous method parallel 
#----------------------------------------
set.seed(1234)
system.time(
  solution8 <- optimStrata(method = "continuous",
                           errors = cv, 
                           nStrata = rep(10,ndom),
                           framesamp = frame2,
                           iter = 50,
                           pops = 10,
                           parallel=T)
)
# *** Starting parallel optimization for  7  domains using  5  cores
# |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=13s  
# 
# *** Sample size :  0
# *** Number of strata :  42
# ---------------------------   
#   utente   sistema trascorso 
# 0.01      0.03      8.41 

outstrata <- solution8$aggr_strata
colnames(outstrata)[10] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1] 91
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.1343402 0.1235310
# DOM2 0.0792397 0.0835106
