#library(devtools)
#install_github("jabbamodel/JABBA")
library(JABBA)
library(plyr)
library(reshape2)
library(ggplot2)

setwd("/home/sbonhomm/Nextcloud/rstudioGithub/makoAssessment/data")
# tata <- read.csv("IOTC-2020-WPEB16-DATA04-CELongline.csv", header=T)
# tata <- tata[, -seq(12,57, by=2)]
# 
# toto <- melt(tata, measure.vars = colnames(tata)[seq(12,34, by=1)])
# toto$variable <- as.character(toto$variable)
# tot_SWO <- ddply(toto, .(Year), function(x) sum(x$value[which(x$variable=="SWO.MT" & x$Gear=="ELL")], na.rm=T))
# 
# ###
# ratioTuna <- 0.05
# ratioSWO  <- 0.56
# 
# tot_TUN <- ddply(toto, .(Year), function(x) data.frame(YFT=sum(x$value[which(x$variable== c("YFT.MT") 
#                                                               & x$Gear %in% c("LL", "FLL"))], na.rm=T),
#                                                        BET=sum(x$value[which(x$variable== c("BET.MT") 
#                                                                              & x$Gear %in% c("LL", "FLL"))], na.rm=T) ,
#                                                        SKJ=sum(x$value[which(x$variable== c("SKJ.MT") 
#                                                                              & x$Gear %in% c("LL", "FLL"))], na.rm=T),
#                                                        ALB=sum(x$value[which(x$variable== c("ALB.MT") 
#                                                                              & x$Gear %in% c("LL", "FLL"))], na.rm=T),
#                                                        TOT=sum(x$value[which(x$variable%in% c("YFT.MT", "BET.MT", "ALB.MT", "SKJ.MT") 
#                                                                              & x$Gear %in% c("LL", "FLL"))], na.rm=T)))
# 
# tot_SWOSWO <- ddply(toto, .(Year), function(x) data.frame(SWO=sum(x$value[which(x$variable== c("SWO.MT") 
#                                                                              & x$Gear %in% c("LL", "FLL"))], na.rm=T)))
# 
# 
# tot_TUN2 <- melt(tot_TUN, id.vars = "Year")
# p <- ggplot(tot_TUN2, aes(x=Year, y= value, colour=variable))+geom_line()
# p
# tot_SMA <- data.frame(Year=tot_SWO$Year, Catch= ratioTuna * tot_TUN$TOT+ratioSWO*tot_SWO$V1)
# #tot_SMA$Catch2 <- c(tot_SMA$Catch[1:64],rep(NA, 3))
# p <- ggplot(tot_SMA, aes(x=Year, y= Catch))+geom_line()
# p

iotc <- NULL
iotc$sma <- NULL
#iotc$sma$catch <- tot_SMA[tot_SMA$Year>1992,]

catchNC <- read.csv("/home/sbonhomm/Downloads/SMA_NC_1964_2018.csv", header=T)

#iotc$sma$catch <- catchNC[catchNC$Year>1992,]
iotc$sma$catch <- catchNC[catchNC$Year>1963,]

### CPUE

#cpueChosen <- c("Year", "Japan", "Taiwan", "Spain_wt", "Portugal") ## out of c("Year", "Japan", "Spain_wt", "Spain_no", "Taiwan", "Portugal")
#cpueChosen <- c("Year", "Spain_wt") ## out of c("Year", "Japan", "Spain_wt", "Spain_no", "Taiwan", "Portugal")
#cpueChosen <- c("Year", "Spain_wt") ## out of c("Year", "Japan", "Spain_wt", "Spain_no", "Taiwan", "Portugal")
cpueChosen <- c("Year", "Taiwan", "Spain_wt", "Portugal") ## out of c("Year", "Japan", "Spain_wt", "Spain_no", "Taiwan", "Portugal")

#iotc$sma$cpue <- as.data.frame(apply(cpueSeries[, cpueChosen], 2, scale, center=T, scale = T))+5
#iotc$sma$cpue[,1] <- cpueSeries[,"Year"]
cpueSelect <- function(cpueChosen){
  output <- NULL
  cpueSeries <- read.csv(file = "/home/sbonhomm/Nextcloud/rstudioGithub/makoAssessment/data/CPUE_SMA_IOTC - all_scaled.csv", header = T)
  output$cpue <- cpueSeries[, cpueChosen]
  if (min(output$cpue$Year) > min(iotc$sma$catch$Year)){
    diffYear <- min(output$cpue$Year)-min(iotc$sma$catch$Year)
    addedDataFrame <- NULL
    for (i in 1:diffYear){addedDataFrame <- rbind(addedDataFrame, output$cpue[1,])}
    addedDataFrame[2:dim(addedDataFrame)[2]] <- NA
    addedDataFrame$Year <- min(iotc$sma$catch$Year):(min(output$cpue$Year)-1)
    output$cpue <- rbind(addedDataFrame, output$cpue)
  }
  cpueSeSeries <- read.csv(file = "/home/sbonhomm/Nextcloud/rstudioGithub/makoAssessment/data/CPUE_SMA_IOTC - all_standardized_se.csv", header = T)
  output$se   <- cpueSeSeries[, cpueChosen]
  if (min(output$se$Year) > min(iotc$sma$catch$Year)){
    diffYear <- min(output$se$Year)-min(iotc$sma$catch$Year)
    addedDataFrame <- NULL
    for (i in 1:diffYear){addedDataFrame <- rbind(addedDataFrame, output$se[1,])}
    addedDataFrame[2:dim(addedDataFrame)[2]] <- NA
    addedDataFrame$Year <- min(iotc$sma$catch$Year):(min(output$se$Year)-1)
    output$se <- rbind(addedDataFrame, output$se)
  }
  return(output)
}

runJABBA <- function(cpueChosen1, modelType, number){
  cpueData      <- cpueSelect(cpueChosen1)
  iotc$sma$cpue <- cpueData$cpue
  iotc$sma$se   <- cpueData$se
  scenario = paste(length(cpueChosen1), modelType, sep="_")
  print(scenario, i, j)
  jbinput = build_jabba(catch = iotc$sma$catch,
                        cpue  = iotc$sma$cpue,
                        se    = iotc$sma$se,
                        assessment = assessment,
                        scenario = scenario,
                        model.type = modelType,
                        sigma.est = TRUE,
                        igamma = c(0.001,0.001),
                        fixed.obsE = 0.01,
                        r.dist = "lnorm",
                        r.prior = c(0.03098706, 0.2),
                        KOBE.type = "IOTC")
  JABBApar = paste0("smainput", number)
  assign(JABBApar, jbinput)

  modelFit = fit_jabba(jbinput)
  JABBAout = paste0("smaoutput", number)
  assign(JABBAout, modelFit)
}

#iotc$sma$se[,2:5] <- 1

File = "/home/sbonhomm/Desktop/SMA_SA" # LINK to your folder of choice here
assessment = "SMAiotc"
output.dir = file.path(File,assessment)
dir.create(output.dir,showWarnings = F)
setwd(output.dir)

### Select CPUEs

cpueChosen <- NULL
cpueChosen$a <- c("Year", "Japan", "Taiwan", "Spain_wt", "Portugal")
cpueChosen$b <- c("Year", "Japan")
cpueChosen$c <- c("Year", "Taiwan")
cpueChosen$d <- c("Year", "Spain_wt")
cpueChosen$e <- c("Year", "Portugal")
cpueChosen$f <- c("Year", "Taiwan", "Spain_wt", "Portugal")
cpueChosen$g <- c("Year", "Japan", "Taiwan", "Spain_wt")

modelType <- c("Schaefer", "Fox", "Pella", "Pella_m")
number <- 38

for (i in 1:7){
  for (j in 1:4){
    number <- number+1
    cpueChosenModel <- eval(parse(text=paste0("cpueChosen$", letters[i])))
    #runJABBA(cpueChosen=cpueChosenModel, modelType[j], number)
    cpueData      <- cpueSelect(cpueChosenModel)
    iotc$sma$cpue <- cpueData$cpue
    iotc$sma$se   <- cpueData$se
    scenario = number
    print(paste(scenario, i, j))
    jbinput = build_jabba(catch      = iotc$sma$catch,
                          cpue       = iotc$sma$cpue,
                          se         = iotc$sma$se,
                          assessment = assessment,
                          scenario   = scenario,
                          model.type = modelType[j],
                          sigma.est  = TRUE,
                          igamma     = c(0.001,0.001),
                          fixed.obsE = 0.01,
                          r.dist     = "lnorm",
                          r.prior    = c(0.113, 0.2),
                          KOBE.type  = "IOTC")
    JABBApar = paste0("smainput", number)
    assign(JABBApar, jbinput)
    save(jbinput, file=paste0(JABBApar, ".Rdata"))
    
    modelFit = fit_jabba(jbinput, save.jabba = TRUE, output.dir=output.dir)
    JABBAout = paste0("smaoutput", number)
    assign(JABBAout, modelFit)
    #save(modelFit, file=paste0(JABBAout,".Rdata"))
  }
}

K       <- NULL
bmsy    <- NULL
fmsy    <- NULL
msy     <- NULL
r_est   <- NULL
stock   <- NULL
harvest <- NULL

for (i in 1:number){
  modelFit <- eval(parse(text=paste0("smaoutput", i)))
  K        <- c(K, modelFit$refpts$k[1])
  bmsy     <- c(bmsy, modelFit$refpts$bmsy[1])
  fmsy     <- c(fmsy, modelFit$refpts$fmsy[1])
  msy      <- c(msy, modelFit$refpts$msy[1])
  r_est    <- c(r_est, modelFit$pars$Median[2])
  stock    <- c(stock, median(modelFit$kobe$stock))
  harvest  <- c(harvest, median(modelFit$kobe$harvest))
}

globalResults <- data.frame(Scenario   = 1:number,
                            CPUE_Type  = c(rep("JPN+TWN+EU.SPAIN+EU.POR", 4), rep("JPN", 4), rep("TWN", 4), rep("EU.SPAIN", 4), rep("EU.POR",4), rep("TWN+EU.Spain+EU.POR",4), rep("JPN+TWN+EU.SPAIN", 4)),
                            Model_Type = rep(c("Schaefer", "Fox", "Pella", "Pella_m"), 7),
                            K          = K, 
                            Bmsy       = bmsy, 
                            Fmsy       = fmsy,
                            MSY        = msy, 
                            r          = r_est, 
                            B_BMSY     = stock,
                            F_FMSY     = harvest)

meltResults <- melt(globalResults, id.vars=c("Scenario", "CPUE_Type", "Model_Type"))

p_model_comp <- ggplot(meltResults, aes(x=Scenario, y=value, colour=Model_Type, shape=CPUE_Type))+
                  geom_point()+
                  facet_wrap(~variable, scales = "free")+
                  scale_shape_manual(values=c(15:19, 7:10))+
                  labs(y="")+
                  theme_light()
p_model_comp

sma1 <- modelFit
# Make individual plots
jbplot_catch(sma1)
jbplot_catcherror(sma1)
jbplot_ppdist(sma1)
jbplot_mcmc(sma1)
jbplot_residuals(sma1)
jbplot_cpuefits(sma1)
jbplot_runstest(sma1)
jbplot_logfits(sma1)
jbplot_procdev(sma1)

# Plot Status Summary
par(mfrow=c(2,2),mar = c(3.5, 3.5, 0.5, 0.1))
jbplot_trj(sma1,type="B",add=T)
jbplot_trj(sma1,type="F",add=T)
jbplot_trj(sma1,type="BBmsy",add=T)
jbplot_trj(sma1,type="FFmsy",add=T)

jbplot_spphase(sma1,add=T)
jbplot_kobe(sma1)

sma1$refpts
quantile(sma1$pars_posterior$K, c(0.025, 0.5, 0.975))
quantile(sma1$pars_posterior$r, c(0.025, 0.5, 0.975))
quantile(sma1$refpts$fmsy, c(0.025, 0.5, 0.975))
quantile(sma1$refpts$bmsy, c(0.025, 0.5, 0.975))
quantile(sma1$refpts$msy, c(0.025, 0.5, 0.975))
quantile(sma1$kobe$stock, c(0.025, 0.5, 0.975))
quantile(sma1$kobe$harvest, c(0.025, 0.5, 0.975))


#----------------------------------------------------------------
# Conduct Retrospective Analysis and Hind-Cast Cross-Validation
#----------------------------------------------------------------
# Organize folders by creating a "retro" subfolder
retro.dir = file.path(output.dir,"retro")
dir.create(retro.dir,showWarnings = F)

# Note as reference input 
refinput = smainput3
# Run hindcasts for reference model (set plotall = TRUE if you want to save plots from all runs)
hc3 = jabba_hindcast(refinput, 
                    save.hc=T, 
                    plotall=T,
                    output.dir = retro.dir, 
                    peels = 0:5)
# Note as reference input 
refinput = smainput23
# Run hindcasts for reference model (set plotall = TRUE if you want to save plots from all runs)
hc23 = jabba_hindcast(refinput, 
                     save.hc=T, 
                     plotall=T,
                     output.dir = retro.dir, 
                     peels = 0:5)

# Note as reference input 
refinput = smainput27
# Run hindcasts for reference model (set plotall = TRUE if you want to save plots from all runs)
hc27 = jabba_hindcast(refinput, 
                     save.hc=T, 
                     plotall=T,
                     output.dir = retro.dir, 
                     peels = 0:5)

save(hc3, file="hc3.Rdata")
save(hc23, file="hc23.Rdata")
save(hc27, file="hc27.Rdata")

# Retro Analysis Summary plot
jbplot_retro(hc23, as.png = T,single.plots = T,output.dir = retro.dir)
# Save plot and note Mohn's rho statistic
mohnsrho23 = jbplot_retro(hc23,as.png = F,single.plots = F,output.dir = retro.dir, col=rainbow(5))
# Zoom-in
mohnsrho23 = jbplot_retro(hc23,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2010,2018))
# eval mohnsrho
mohnsrho23

# Do Hindcast Cross-Validation (hcxval) 
# show multiplot
jbplot_hcxval(hc23,single.plots = F,as.png = F,output.dir=retro.dir,col=rainbow(8))
# Zoom-in
jbplot_hcxval(hc23,single.plots = F,as.png = F,output.dir=retro.dir,minyr=2010)
# save as png and note summary stats 
mase23 = jbplot_hcxval(hc23,single.plots = F,as.png = TRUE,output.dir=retro.dir)
#check stats
mase23

# Retro Analysis Summary plot
jbplot_retro(hc3,as.png = F,single.plots = F,output.dir = retro.dir)
# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc3,as.png = T,single.plots = F,output.dir = retro.dir)
# Zoom-in
mohnsrho3 = jbplot_retro(hc3,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2010,2018))
# eval mohnsrho
mohnsrho3

# Do Hindcast Cross-Validation (hcxval) 
# show multiplot
jbplot_hcxval(hc3,single.plots = F,as.png = F,output.dir=retro.dir,col=rainbow(8))
# Zoom-in
jbplot_hcxval(hc3,single.plots = F,as.png = F,output.dir=retro.dir,minyr=2010)
# save as png and note summary stats 
mase3 = jbplot_hcxval(hc3,single.plots = F,as.png = TRUE,output.dir=retro.dir)
#check stats
mase3

# Retro Analysis Summary plot
jbplot_retro(hc27,as.png = F,single.plots = F,output.dir = retro.dir)
# Save plot and note Mohn's rho statistic
mohnsrho = jbplot_retro(hc27,as.png = T,single.plots = F,output.dir = retro.dir)
# Zoom-in
mohnsrho27 = jbplot_retro(hc27,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2010,2018))
# eval mohnsrho
mohnsrho27

# Do Hindcast Cross-Validation (hcxval) 
# show multiplot
jbplot_hcxval(hc27,single.plots = F,as.png = F,output.dir=retro.dir,col=rainbow(8))
# Zoom-in
jbplot_hcxval(hc27,single.plots = F,as.png = F,output.dir=retro.dir,minyr=2010)
# save as png and note summary stats 
mase27 = jbplot_hcxval(hc27,single.plots = F,as.png = TRUE,output.dir=retro.dir)
#check stats
mase27





# Retro Analysis Summary plot
jbplot_retro(hc38, as.png = T,single.plots = T,output.dir = retro.dir)
# Save plot and note Mohn's rho statistic
mohnsrho38 = jbplot_retro(hc38,as.png = F,single.plots = F,output.dir = retro.dir, col=rainbow(5))
# Zoom-in
mohnsrho38 = jbplot_retro(hc38,as.png = F,single.plots = F,output.dir = retro.dir,Xlim=c(2010,2018))
# eval mohnsrho
mohnsrho38

jbplot_hcxval(hc38,single.plots = F,as.png = F,output.dir=retro.dir,col=rainbow(8))
# Zoom-in
jbplot_hcxval(hc38,single.plots = F,as.png = F,output.dir=retro.dir,minyr=2010)
# save as png and note summary stats 
mase27 = jbplot_hcxval(hc38,single.plots = F,as.png = TRUE,output.dir=retro.dir)
#check stats
mase27



load(file="hc3.Rdata")

#---------------------------
# Run REF with projections
#---------------------------
Scenarios=c("3", "23")
jbplot_summary(assessment=assessment,scenarios = Scenarios, mod.path = output.dir, cols=terrain.colors(3))
jbplot_summary(assessment=assessment,scenarios = Scenarios, plotCIs=FALSE, save.summary = T, as.png = T)
jbplot_summary(assessment=assessment,scenarios = Scenarios, plotCIs=FALSE, save.summary = T, as.png = F)

projinput <- smainput23

TAC  <- round(seq(0.6,1.4, by=0.1)*catchNC$sma_nc[dim(catchNC)[1]])
projections = build_jabba(catch      = projinput$data$catch,
                          cpue       = projinput$data$cpue,
                          se         = projinput$data$se,
                          assessment = assessment,
                          scenario   = projinput$settings$scenario,
                          model.type = projinput$settings$model.type,
                          sigma.est  = TRUE,
                          igamma     = c(0.001,0.001),
                          fixed.obsE = 0.01,
                          r.dist     = "lnorm",
                          r.prior    = c(0.03098706, 0.2),
                          KOBE.type  = "IOTC",
                          projection = TRUE,
                          TACs       = TAC)


smaprj23 = fit_jabba(projections,output.dir=output.dir, save.csvs = T)
# plot with CIs (80% for projections)
jbplot_prj(smaprj23,type="BBmsy", CIs=FALSE)
jbplot_prj(smaprj23,type="BB0", CIs=FALSE)

# or without CIs (80% for projections) 
jbplot_prj(smaprj23,type="FFmsy", CIs=FALSE)

