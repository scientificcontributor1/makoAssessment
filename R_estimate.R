library(pacman)
p_load(demogR, devtools, matlib, grid, egg)

###################
### Parameters  ###
###################

### Age at maturity
Tmat_low    <- 18
Tmat_high   <- 21
T_mat_best  <- 18

### Life span
TermAge_low <- 29
TermAge_high <- 32
TermAge_best <- 32

### Growth
Linf_best <- 350.3
Linf_err  <- 0.1
k_best    <- 0.064
k_err     <- 0.1
t0_best   <- -3.09
t0_err    <- 0.1

### Pup term size
pupSize_best <- 73
pupSize_err  <- 0.1

### W-L relationship for youngs
a_WLYoung_best <- 8.198
a_WLYoung_err  <- 0.1
b_WLYoung_best <- 3.117
b_WLYoung_err  <- 0.1

### W-L relationship for adults
a_WLAdult_best <- 0.0000052432
a_WLAdult_err  <- 0.1
b_WLAdult_best <- 3.1407
b_WLAdult_err  <- 0.1


rEstimator <- function(Tmat, TermAge, Linf, k, t0, pupSize, a_WLYoung, b_WLYoung, a_WLAdult, b_WLAdult ){
  #####################
  ### Pup mortality ###
  #####################
  ### Calculate Length at age for pups
  Lyoung <- pupSize*exp(-4*exp(-0.01*0:(30*18)))
  ### Calculate weight-at-age for pups
  Wyoung <- (a_WLYoung*(Lyoung/100)^b_WLYoung)*1000
  ### Calculate mortality-at-age using McGurk 1986 relationship
  Myoung <- 0.00022*(Wyoung*0.15)^(-0.85)
  M0     <- sum(Myoung[1:365])
  #######################
  ### Adult mortality ###
  #######################
  ### Natural mortality Then et al. (independant for age)
  M_T = 4.899*TermAge^(-0.916)
  ### Calculate Length at age for adults
  LSMA   <- Linf*(1-exp(-k*((1:TermAge)-t0)))
  ### Natural mortality using Chen and Watanabe (1989)
  M_CW                   <- NULL
  M_CW[1:Tmat]           <-  k/(1-exp(-k*((0:(Tmat-1)) - t0)))
  a0                     <- (1-exp(-k*(Tmat - t0)))
  a1                     <- k*exp(-k*(Tmat-t0))
  a2                     <- -0.5*k^2*exp(-k*(Tmat-t0))
  M_CW[(Tmat+1):TermAge] <- k/(a0 + a1*(((Tmat+1):TermAge)-Tmat) + a2*(((Tmat+1):TermAge)-Tmat)^2)
  M_CW <- c(M0,M_CW)
  ### Natural mortality using McGurk / Peterson-Wroblewski
  WSMA <- a_WLAdult*LSMA^b_WLAdult
  M_PW <- (0.00526*(WSMA*1000*0.15)^(-0.25))*365
  M_PW <- c(M0, M_PW)
  #################
  ### Fecundity ###
  #################
  LS = 0.810*(LSMA/100)^2.346
  LS[1:(Tmat-1)] <- 0
  LS <- c(0,LS)
  #LS <- LS[2:length(LS)]
  #####################################
  ### Generate the Leslie matrices  ###
  #####################################
  # Survival at age for Then et al.
  mlt_T <- c(M0, rep(exp(-M_T), TermAge))
  # Survival at age for Chen and Watanabe (1989)
  mlt_CW <- exp(-M_CW)
  # Survival at age for McGurk / Peterson-Wroblewski
  mlt_PW <- exp(-M_PW)
  # Fecundity at age
  mx <- LS/2
  ### r estimate for Then et al.
  mad_T <- leslie.matrix(lx=mlt_T, mx=mx, SRB=1, infant.class = T)
  mad_T[(1:dim(mad_T)[1]),(1:dim(mad_T)[1])] <- 0
  mad_T[1,] <- LS[2:(dim(mad_T)[1]+1)]/2
  for (i in 1:(dim(mad_T)[1]-1)){mad_T[(i+1),(i)] <- mlt_T[i]}
  rpow_T <- log(powerMethod(mad_T)$value)
  r_T    <- log(eigen.analysis(mad_T)$lambda1)
  ### r estimate for Chen and Watanabe (1989)
  mad_CW <- leslie.matrix(lx=mlt_CW, mx=mx, SRB=1, infant.class = T)
  mad_CW[(1:dim(mad_CW)[1]),(1:dim(mad_CW)[1])] <- 0
  mad_CW[1,] <- LS[2:(dim(mad_CW)[1]+1)]/2
  for (i in 1:(dim(mad_CW)[1]-1)){mad_CW[(i+1),(i)] <- mlt_CW[i]}
  rpow_CW <- log(powerMethod(mad_CW)$value)
  r_CW    <- log(eigen.analysis(mad_CW)$lambda1)
  ### r estimate for McGurk / Peterson-Wroblewski
  mad_PW <- leslie.matrix(lx=mlt_PW, mx=mx, SRB=1, infant.class = T)
  mad_PW[(1:dim(mad_PW)[1]),(1:dim(mad_PW)[1])] <- 0
  mad_PW[1,] <- LS[2:(dim(mad_PW)[1]+1)]/2
  for (i in 1:(dim(mad_PW)[1]-1)){mad_PW[(i+1),(i)] <- mlt_PW[i]}
  rpow_PW <- log(powerMethod(mad_PW)$value)
  r_PW    <- log(eigen.analysis(mad_PW)$lambda1)
  ### Generate an output dataframe
  output <- data.frame(Tmat      = Tmat,
                       TermAge   = TermAge,
                       Linf      = Linf,
                       k         = k,
                       t0        = t0,
                       pupSize   = pupSize,
                       a_WLYoung = a_WLYoung,
                       b_WLYoung = b_WLYoung,
                       a_WLAdult = a_WLAdult,
                       b_WLAdult = b_WLAdult,
                       rpow_T    = rpow_T,
                       r_T       = r_T,
                       rpow_CW   = rpow_CW,
                       r_CW      = r_CW,
                       rpow_PW   = rpow_PW,
                       r_PW      = r_PW
                       )
}

### Bootstrap 1000 times with uncertainties on parameters
databaseR <- NULL

for (boot in 1:10000){
  print(boot)
  outputs <- rEstimator(Tmat      = floor(runif(1, Tmat_low, Tmat_high)),
                        TermAge   = floor(runif(1, TermAge_low, TermAge_high)),
                        Linf      = Linf_best*rnorm(1, mean=1, sd= Linf_err),
                        k         = k_best*rnorm(1,mean=1, sd=k_err),
                        t0        = t0_best*rnorm(1, mean=1, sd=t0_err),
                        pupSize   = pupSize_best*rnorm(1, mean=1, sd= pupSize_err),
                        a_WLYoung = a_WLYoung_best*rnorm(1, mean=1, sd= a_WLYoung_err),
                        b_WLYoung = b_WLYoung_best*rnorm(1, mean=1, sd= b_WLYoung_err),
                        a_WLAdult = a_WLAdult_best*rnorm(1, mean=1, sd= a_WLAdult_err),
                        b_WLAdult = b_WLAdult_best*rnorm(1, mean=1, sd= b_WLAdult_err))
  databaseR <- rbind(databaseR, outputs)
}

r_pow_PW_med <- median(databaseR$rpow_PW[which(databaseR$rpow_PW>0)])
r_pow_PW_sd  <- sd(databaseR$rpow_PW[which(databaseR$rpow_PW>0)])
r_pow_CW_med <- median(databaseR$rpow_CW[which(databaseR$rpow_CW>0)])
r_pow_CW_sd  <- sd(databaseR$rpow_CW[which(databaseR$rpow_CW>0)])
r_pow_T_med <- median(databaseR$rpow_T[which(databaseR$rpow_T>0)])
r_pow_T_sd  <- sd(databaseR$rpow_T[which(databaseR$rpow_T>0)])
cv_pow_VW <-  r_pow_CW_sd / r_pow_CW_med 

meltedData <- melt(databaseR[, colnames(databaseR)[seq(11,16, by=2)]])
p <- ggplot(meltedData, aes(x=value))+
      geom_histogram(binwidth=0.005, fill="grey85")+
     geom_density(aes(y=0.005 * ..count..), colour="dark red")+
      facet_wrap(~variable, scales = "free")+
      theme_light()
p

meltedDataPos <- meltedData[which(meltedData$value>0),]
medianValues <- ddply(meltedDataPos, .(variable), function(x) median(x$value))
method_names <- c(
  `rpow_T` = "r Then et al. 2015",
  `rpow_CW` = "r Chen and Watanabe (1989)",
  `rpow_PW` = "r McGurk 1986"
)

p <- ggplot(meltedDataPos, aes(x=value))+
  geom_histogram(binwidth=0.005, fill="grey85")+
  geom_density(aes(y=0.005 * ..count..), colour="dark red")+
  geom_vline(data=medianValues, mapping = aes(xintercept = V1), colour="dark blue")+
  facet_wrap(~variable, scales = "free_y", labeller=as_labeller(method_names))+
  labs(x=" r value", y ="Counts")+
  theme_light()
p



  
  
Tmat_low    <- 18
Tmat_high   <- 21
T_mat_best  <- 18

### Life span
TermAge_low <- 29
TermAge_high <- 32
TermAge_best <- 32

### Growth
Linf_best <- 350.3
Linf_err  <- 0.1
k_best    <- 0.064
k_err     <- 0.1
t0_best   <- -3.09
t0_err    <- 0.1

### Pup term size
pupSize_best <- 73
pupSize_err  <- 0.1

### W-L relationship for youngs
a_WLYoung_best <- 8.198
a_WLYoung_err  <- 0.1
b_WLYoung_best <- 3.117
b_WLYoung_err  <- 0.1

### W-L relationship for adults
a_WLAdult_best <- 0.0000052432
a_WLAdult_err  <- 0.1
b_WLAdult_best <- 3.1407
b_WLAdult_err  <- 0.1






Tmat=18
TermAge=32
Linf= 350.3
k= 0.064
t0=-3.09
pupSize=73
a_WLYoung=8.198
b_WLYoung=3.117
a_WLAdult= 0.0000052432
b_WLAdult=3.1407


###Natural mortality estimate

### Van Bertalanffy growth (https://www.researchgate.net/publication/329634097_Age_and_Growth_of_the_Shortfin_Mako_Shark_in_the_Southern_Indian_Ocean)
Linf <- 267.6
k <- 0.123
t0 <- -2.487
### Rosa et al. 2017
Linf <- 350.3
k <- 0.064
t0 <- -3.09
LSMA <- Linf*(1-exp(-k*((0:TermAge)-t0)))


## M for the first 12 month
a_W = 8.198
b_W = 3.117
Lyoung <- 73*exp(-4*exp(-0.01*0:(30*18)))
#Lyoung[1:23] <- 3
Wyoung <- (a_W*(Lyoung/100)^b_W)*1000
### M from McGurk 1986
Myoung <- 0.00022*(Wyoung*0.15)^(-0.85)
M0 <- sum(Myoung[1:365])

### Natural mortality Then et al.
M_T = 4.899*TermAge^(-0.916)

### Natural mortality Chen and Watanabe (1989)
M <- NULL
M[1:Tmat] = k/(1-exp(-k*((0:(Tmat-1)) - t0)))
a0 <- (1-exp(-k*(Tmat - t0)))
a1 <- k*exp(-k*(Tmat-t0))
a2 <- -0.5*k^2*exp(-k*(Tmat-t0))
M[(Tmat+1):TermAge] <- k/(a0 + a1*(((Tmat+1):TermAge)-Tmat) + a2*(((Tmat+1):TermAge)-Tmat)^2)
M_CW <- M
M_CW <- c(M0,M_CW)

### McGurk / Peterson-Wroblewski
a_WL <- 0.0000052432
b_WL <- 3.1407
WSMA <- a_WL*LSMA^b_WL

M_PW <- (0.00526*(WSMA*1000*0.15)^(-0.25))*365
M_PW <- c(M0, M_PW)

plot(M_CW, ylim=c(0.05,M0+0.05), type="l", col="dark blue",
     ylab=expression(Natural ~ Mortality ~ (years^{-1} )), xlab="Age (in years)",
     cex.axis=1.5, cex.lab=1.5, lwd=2)
lines(M_PW, col="dark red", lwd=2)
abline(h=M_T, col="dark green", lwd=2)
legend("topright", 
       legend=c("Then et al (2015)", "Chen & Watanabe (1989)", "Peterson & Wroblewski (1984)"),
       col=c("dark blue", "dark red", "dark green"),
       lwd=2, bty="n", cex=2)
       

### Fecundity
LS = 0.810*(LSMA/100)^2.346
LS[1:Tmat] <- 0
LS <- LS[2:length(LS)]


# Survival at age
mlt <- exp(-M_PW)
# Fecundity at age
mx <- LS/2

mad <- leslie.matrix(lx=mlt, mx=mx, SRB=1, infant.class = T)
mad[(1:31),(1:31)] <- 0
mad[1,] <- LS[1:31]/2
for (i in 1:30){mad[(i+1),(i)] <- mlt[i]}

rpow <- log(powerMethod(mad)$value)
rpow
r <- log(eigen.analysis(mad)$lambda1)
r

### W-L relationship for elbryos /MolletPrattetalmakorepro2000.pdf
a_W = 8.198
b_W = 3.117

aexp= 73
kg=0.123
Ti=180
curve(aexp*(-exp(-kg*(x-Ti))), -60,360)
curve(73*exp(-4*exp(-0.01*x)), 0,18*30)

Wegg <- rep(0.0005, 4) ### 4 days in egg 


Wyoung <- (exp(-3)*(Lyoung[1:30]*10)^3.117)/1000
## When size > 36 cm (174 days)  MolletPrattetalmakorepro2000.pdf
Wyoung <- a_W*(Lyoung[31:length(Lyoung)]/100)^b_W


Md <- 0.00022*(Wyoung*0.15)^-0.85
iwhWeight <- (Wyoung*0.15)>0.00504
Md[iwhWeight] <- 0.00526*(Wyoung[iwhWeight]*0.15)^-0.25
sum(Md)

1-exp(-Md)



