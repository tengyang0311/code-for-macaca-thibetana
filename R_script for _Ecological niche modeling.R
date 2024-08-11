#Notes:as we have two species used in our analysis, for simplicity, this script will just use Parus monticolus as an example.
#This script used for Ecological niche modeling
#The first part used for maxent parameter evaluation using ENMwizard, and the second part used for ensemble modeling using biomod2
#Install the dependencies package
#install.packages("devtools")
#devtools::install_github("HemingNM/ENMwizard")
#install.packages("biomod2", dependencies = TRUE)
rm(list = ls())
library(ENMwizard)#version 0.3.9
library(biomod2)#version 3.5.1
############
#first part#
############

#Load worldclim dataset
#select uncorrelated bio-factor
c <- c(2,9,14,15,18)
wc_cur <- dir("D:\\GIS\\bio19\\wc_cur_2.5m_bio",full.names = T,pattern = ".tif$")
wc_cur <- stack(wc2.1[c])
handle_fut <- function(x,y,z){
  fut <- getData('CMIP5', var='bio', res=2.5, rcp=x, model=y, year=z,
                  path="D:\\Work\\cczm-final\\SDM/rasters")
  fut <- fut[[bio]]
  projection(fut)<- crs("+proj=longlat +datum=WGS84 +no_defs")
  names(fut) <- names(wc_cur)
  return(fut)
}
#4550
futcn4550 <- handle_fut(45,'CN',50)
futmc4550 <- handle_fut(45,'MC',50)
futcc4550 <- handle_fut(45,'CC',50)
futmp4550 <- handle_fut(45,'MP',50)
#8550
futcn8550 <- handle_fut(85,'CN',50)
futmc8550 <- handle_fut(85,'MC',50)
futcc8550 <- handle_fut(85,'CC',50)
futmp8550 <- handle_fut(85,'MP',50)
#4570
futcn4570 <- handle_fut(45,'CN',70)
futmc4570 <- handle_fut(45,'MC',70)
futcc4570 <- handle_fut(45,'CC',70)
futmp4570 <- handle_fut(45,'MP',70)
#8570
futcn8570 <- handle_fut(85,'CN',70)
futmc8570 <- handle_fut(85,'MC',70)
futcc8570 <- handle_fut(85,'CC',70)
futmp8570 <- handle_fut(85,'MP',70)
#
predictors.l <- list(current = wc_cur,
                     futcn4550 = futcn4550,
                     futmc4550 = futmc4550,
                     futcc4550 = futcc4550,
                     futmp4550 = futmp4550,
                     futcn8550 = futcn8550,
                     futmc8550 = futmc8550,
                     futcc8550 = futcc8550,
                     futmp8550 = futmp8550,
                     futcn4570 = futcn4570,
                     futmc4570 = futmc4570,
                     futcc4570 = futcc4570,
                     futmp4570 = futmp4570,
                     futcn8570 = futcn8570,
                     futmc8570 = futmc8570,
                     futcc8570 = futcc8570,
                     futmp8570 = futmp8570)
#Load locations
gbt_cold_dry <- read.csv("./cold_dry.csv",header = T)
gbt_cold_dry.list <- list(gbt.cold_dry = gbt_cold_dry)
gbt_cold_dry.poly <- set_calibarea_b(gbt_cold_dry.list,plot = T,save.pts = T)
gbt_cold_dry.b <- buffer_b(gbt_cold_dry.poly,width = 2)
gbt_cold_dry.cut <- cut_calibarea_b(gbt_cold_dry.b, wc_cur,numCores= 5)
#thinned.dataset.batch.cold_dry <- thin_b(loc.data.lst = gbt_cold_dry.list)$filter locations if needed
#gbt_cold_dry.locs <- load_thin_occ(thinned.dataset.batch.cold_dry)
#load projected area
gbt_cold_dry_int <- rgdal::readOGR("./gbt_cold_dry.shp")
gbt_cold_dry_int.b <- raster::mask(current, buffer(gbt_cold_dry_int, width=2))
gbt_cold_dry_all<- raster::merge(gbt_cold_dry_int.b, gbt_cold_dry.cut$gbt.cold_dry)
gbt_cold_dry.cut <- list(gbt.cold_dry = gbt_cold_dry_all)
#Model tuning for maxent
ENMeval.res.lst.cold_dry <- ENMevaluate_b(gbt_cold_dry.locs, gbt_cold_dry.cut, RMvalues = seq(0.5, 5, 0.5),
                                          fc = c("L", "P", "Q", "H", "LP", "LQ", "LH", "PQ", "PH", "QH", "LPQ", "LPH", "LQH",
                                                 "PQH", "LPQH"),
                                      method="checkerboard2", algorithm="maxent.jar",
                                      numCores = 12,clamp = T)
#Model fitting for maxent
mxnt.mdls.preds.lst.cold_dry <- calib_mdl_b(ENMeval.o.l = ENMeval.res.lst.cold_dry, 
                                        a.calib.l = gbt_cold_dry.cut,
                                        mSel = c("LowAIC"),arg1 = "jackknife=true")

a.proj.l.cold_dry <- cut_projarea_rst_mscn_b(gbt_cold_dry.cut$gbt.cold_dry, predictors.l, gbt_cold_dry.poly,
                                         mask = F)
#Model projections for maxent
mxnt.mdls.preds.cf.cold_dry <- proj_mdl_b(mxnt.mdls.preds.lst.cold_dry, a.proj.l= a.proj.l.cold_dry,
                                      numCores=12)
#extract the 10% points that have lowest probability of presence
lo5 <- raster::extract(mxnt.mdls.preds.cf.cold_dry$gbt.cold_dry$mxnt.preds$current, 
                       ENMeval.res.lst.west[["gbt.cold_dry"]]@occ.pts) 
num5 <- round(nrow(lo5)*0.1, digits = 0)#check number to remove
lo5 <- cbind(ENMeval.res.lst.west[["gbt.cold_dry"]]@occ.pts, lo5)
lo5 <- arrange(lo5, lo5$Mod.LowAIC)[1:num5,][,1:2]
pointk <- setdiff(ENMeval.res.lst.west[["gbt.cold_dry"]]@occ.pts, lo5)
write.csv(pointk,"./gbt_cd_filter.csv")

#############
#second part#
#############

DataSpecies <- read.csv("./gbt_cd_filter.csv", header = T)[,2:3]
myRespName <- 'GBT'
myResp <- as.numeric(rep(1,nrow(DataSpecies)))
myExpl <- gbt_cold_dry_all
myExpl_proj <- mask(myExpl, gbt_cold_dry_int) %>% stack
#data formating, selecte PA randomly
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl ,
                                     resp.xy = DataSpecies,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 10000,
                                     PA.strategy = 'random'
)
#set maxent parameter based on the result of model tuning
myBiomodOption_2 <- BIOMOD_ModelingOptions(MAXENT.Phillips = list( path_to_maxent.jar = paste(system.file(package="dismo"), "/java", sep=''),
                                                                   memory_allocated = 4096,
                                                                   background_data_dir = 'default',
                                                                   maximumbackground = 'default',
                                                                   maximumiterations = 200,
                                                                   visible = FALSE,
                                                                   linear = F,#based on the tuning result
                                                                   quadratic = TRUE,#based on the tuning result
                                                                   product = TRUE,#based on the tuning result
                                                                   threshold = F,
                                                                   hinge = TRUE,#based on the tuning result
                                                                   lq2lqptthreshold = 80,
                                                                   l2lqthreshold = 10,
                                                                   hingethreshold = 15,
                                                                   beta_threshold = -1,
                                                                   beta_categorical = -1,
                                                                   beta_lqp = -1,
                                                                   beta_hinge = -1,
                                                                   betamultiplier = 2.0,
                                                                   defaultprevalence = 0.5))
#modeling
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                     models = c('MAXENT.Phillips','GBM','GAM',
                                                'MARS'), 
                                     models.options = myBiomodOption_2, 
                                     NbRunEval=5, 
                                     DataSplit=70, 
                                     Yweights=NULL, 
                                     VarImport=3, 
                                     models.eval.meth = c('TSS','ROC'),
                                     SaveObj = TRUE,
                                     rescal.all.models = T,
                                     do.full.models = FALSE,
                                     modeling.id='gbt_cd')
model.eva<- get_evaluations(myBiomodModelOut)
#EnsembleModeling
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by='all',
                                      eval.metric = c('TSS','ROC'),
                                      eval.metric.quality.threshold = c(0.6,0.8),
                                      models.eval.meth = c('TSS','ROC'),
                                      prob.mean = T,
                                      prob.cv = F,
                                      prob.ci = F,
                                      prob.ci.alpha = 0.05,
                                      prob.median = T,
                                      committee.averaging = F,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional' )
model.eva.EM <- get_evaluations(myBiomodEM)
#Projection
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl_proj,
  proj.name = 'current', 
  selected.models = 'all', 
  binary.meth = 'TSS', 
  compress = 'gzip', 
  clamping.mask = F, 
  output.format = '.grd',
  do.stack=T)
mod_projPres <- get_predictions(myBiomodProj)
#EnsembleModeling projection
myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj,
                                         binary.meth ="TSS",
                                         proj.name = 'current_EM')
projPresEnsemble <- get_predictions(myBiomodEF)