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
# [1] 947

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
#  57.23      5.45     66.16 
sum(solution1$aggr_strata$SOLUZ)
# [1] 502.3473
expected_CV(solution1$aggr_strata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0985073 0.1000000
# DOM2 0.1000000 0.1000000
# DOM3 0.0974722 0.0977469
# DOM4 0.0999351 0.0967996
# DOM5 0.1000000 0.1000000
# DOM6 0.0992327 0.0990466
# DOM7 0.0974012 0.0998922

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
#   0.40      0.11     16.30 

sum(solution2$aggr_strata$SOLUZ)
# [1] 501.4645
expected_CV(solution2$aggr_strata)
# cv(Y1)    cv(Y2)
# DOM1 0.0992139 0.0997360
# DOM2 0.0995165 0.0999033
# DOM3 0.0988785 0.0972292
# DOM4 0.0986402 0.0828498
# DOM5 0.0999727 0.0994666
# DOM6 0.0995148 0.0980787
# DOM7 0.0999553 0.0993688


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
# 82.92      5.13     91.78  

outstrata <- solution3$aggr_strata
sum(outstrata$SOLUZ)
# [1] 145.5547
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0975718 0.0984605
# DOM2 0.0923348 0.0963319
# DOM3 0.0934600 0.0974126
# DOM4 0.0958892 0.0786081
# DOM5 0.0975011 0.0998418
# DOM6 0.0953433 0.0887643
# DOM7 0.0907852 0.0917804

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

#   utente   sistema trascorso 
# 0.12      0.14     22.68 

outstrata <- solution4$aggr_strata
sum(outstrata$SOLUZ)
# [1] 150.6363
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0992080 0.0993171
# DOM2 0.0981495 0.0988630
# DOM3 0.0981118 0.0934174
# DOM4 0.0979480 0.0974263
# DOM5 0.0973675 0.0970349
# DOM6 0.0946572 0.0899318
# DOM7 0.0961231 0.0984867

save.image(file="run.Rdata")
#--------------------------------------
# Performance with SamplingStrataC 
#--------------------------------------
detach("package:SamplingStrata", unload = TRUE)
# Better to restart R
# install.packages("D://Google Drive//Sampling//SamplingStrata-C-//SamplingStrataC_1.5-4.tar.gz", repos = NULL, type = "source")

library(SamplingStrataC)

all2 <- bethelRcppOpen(strata1,cv[1,],realAllocation = TRUE)
sum(all2)
# [1] 947  (ok!!!)

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
# 8.31      1.51     11.22 

outstrata <- solution5$aggr_strata
colnames(outstrata)[10] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1]  520.1812
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0985074 0.1000000
# DOM2 0.1000000 0.1000000
# DOM3 0.1000000 0.1000000
# DOM4 0.0917426 0.0990130
# DOM5 0.0999852 0.0999813
# DOM6 0.0990527 0.0994976
# DOM7 0.1000000 0.1000000

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
#   utente   sistema trascorso 
#  0.44      0.05      3.52    

outstrata <- solution6$aggr_strata
colnames(outstrata)[10] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1] 476.3474
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.1000000 0.1000000
# DOM2 0.0993866 0.0999889
# DOM3 0.1004838 0.0987443
# DOM4 0.0686157 0.0967519
# DOM5 0.1000000 0.1000000
# DOM6 0.0994908 0.0991029
# DOM7 0.0999291 0.0989241

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
# 65.29      4.56     74.97  

outstrata <- solution7$aggr_strata
colnames(outstrata)[ncol(outstrata)] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1] 144.2907
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.0989425 0.0965060
# DOM2 0.0935584 0.0985200
# DOM3 0.0918125 0.0945989
# DOM4 0.1103029 0.1000580
# DOM5 0.0967879 0.0953595
# DOM6 0.0938451 0.0952151
# DOM7 0.0893192 0.0948410

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
# 0.13      0.06     18.13 

outstrata <- solution8$aggr_strata
colnames(outstrata)[10] <- "SOLUZ"
sum(outstrata$SOLUZ)
# [1] 275
expected_CV(outstrata)
#         cv(Y1)    cv(Y2)
# DOM1 0.1189196 0.1066284
# DOM2 0.0996083 0.1268673
# DOM3 0.0458874 0.0720247
# DOM4 0.1018769 0.1124923
# DOM5 0.0732544 0.0761021
# DOM6 0.0874267 0.0980321
# DOM7 0.0798213 0.1314007
