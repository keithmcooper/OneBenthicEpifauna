################################################################################
######### RANDOM FOREST MODELS - EPIFAUNA  METRICS    ##########################
######### 25/11/2025                                  ##########################
################################################################################

#### Note - This code is a simplified version of a previous script made to loop
#### through multiple response variables and sets of predictor variables that
#### included factor variables - this means parts of the code are redundant
####  and/or repetitive but have not been edited out to save time


######### 1. Set up libraries, directories and functions #######################

## Load libraries
require(data.table)
require(sf)
require(ggplot2)
require(gridExtra)
require(raster)
require(dplyr)
require(tidyr)
require(kableExtra)
require(summarytools)
require(ggcorrplot)
require(randomForest)
require(pdp)
require(patchwork)
require(caret)
require(ggh4x)

## Set Working Directory
setwd("WORKING DIRECTORY PATH HERE")

## Input and output paths
in_path = "PATH TO INPUT DATA"
out_path = "PATH FOR OUTPUTS"

## Predictor raster folder
rs ='PATH TO PREDICTOR RASTERS'

## Response data file
## Column 1: 'sample' Unique sample identifier
## Column 2: 'x' Longitude (dd.00)
## Column 3: 'y' Latitude (dd.00)    
## Column 4: 'metric' variable: Cluster class / S / N   
## Column 5: 'type' metric format: character / numeric   
## Column 6: 'value' metric value   
rd = 'Trawl_Epifauna_25112025.csv'

## File with labels to use for environmental variables
## Column 1: 'Variable' with predictor column names
## Column 2: 'Label' with labels to use for plotting
## Column 3: 'Include' with Yes/No for inclusion in preliminary models
rn = 'RasterNamesAndLabels2024.csv'

## Colour palette
colpal <- c("#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00", "#9aff9a",
            "#b40202","#ff0000","#ff8c00","#ffff00","#b4b404") # Define specific colours

## Highlight colours for plotting
hlp <- c("#9aff9a","#ffff00","#b4b404")


## Null model vif function
corvif =  function(dataz) {
  dataz <- as.data.frame(dataz)
  
  # vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(data.frame(vif=car::vif(lm_mod)))
}


## Objects to collect standard outputs for combining later

# Input response data
RD <- NULL
# An  object for appending all variable importances to
VIs <- NULL
# An object for appending model performance  table data
PSs <- NULL
# Observed vs predicted values across test runs
OvsPs <- NULL
# Partial Response Curve data
PRCs <- NULL

# Table for collecting all model performance statistics
forest.resultsSQ <- data.frame(Name=character(0),
                               Run=character(0),
                               N=numeric(0),
                               ExtRMSE=numeric(0),
                               ExtRMSEpctR=numeric(0),
                               ExtR2=numeric(0),
                               stringsAsFactors =F)


######### 2. Read in Data ######################################################

### Response point data  ----

# Read in csv, pivot to wide format and rename coordinate columns
response <- as.data.table(read.csv(paste0(in_path,rd))) %>%
  filter(type=='numeric') %>%
  select(sample,x,y,metric,value) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  mutate_at(vars(-1), as.numeric) 
summary(response)

# convert to spatial
resp_sf <- st_as_sf(response,coords=c('x','y'),crs=CRS('+proj=longlat +datum=WGS84 +no_defs'))

### Environmental predictor rasters  ----

# Read in all rasters in folder
f <- list.files(path=rs, pattern='.tif$', full.names=T)
s1 <- stack(f,RAT=F)
s1
plot(s1)


### Environmental variable labels -----

# Read in label file
varlabs <- read.csv(paste0(in_path,rn))

# Convert to named vector
envlab <- varlabs$Label
names(envlab) <- varlabs$Variable
envlab

# List variables to keep
# This removes some of the range variables at this point - see varlabs
keeprast <- varlabs %>%
  filter(Include=='Yes') %>%
  pull(Variable)

# Labels for ggplot axes
label_map <- list(
  x = expression(Longitude),
  y = expression(Latitude),
  Bathymetry = expression(Depth~(m)),
  VD = expression(Valley~depth~(m)),
  CDP0 = expression(Closed~depressions),
  CNBL = expression(Ch.~network~baselevel~(m)),
  CND = expression(Ch.~network~distance~(m)),
  LSF = expression(LS - factor),
  RSP = expression(Rel.~slope~pos.),
  Wave_veloc = expression(Wave~velocity~(m~s^{-1})),
  Current_Sp = expression(Current~speed~(m~s^{-1})),
  Predicted_Gravel_Fraction = expression(Gravel~"%"),
  Predicted_Mud_Fraction = expression(Mud~"%"),
  SPM_MEAN = expression(Mean~SPM~(g~m^{-3})),
  SPM_SUMMER = expression(Summer~SPM~(g~m^{-3})),
  SPM_WINTER = expression(Winter~SPM~(g~m^{-3})),
  dfe_mean = expression(Diss.~Iron~(mmol~m^{-3})),
  o2_mean = expression(Diss.~Oxygen~(mmol~m^{-3})),
  no3_mean = expression(Nitrate~(mmol~m^{-3})),
  thetao_mean = expression(Bottom~temp.~degree*C),
  thetao_range = expression(Bottom~temp.~range~degree*C),
  ph_mean = expression(pH),
  po4_mean = expression(Phosphate~(mmol~m^{-3})),
  so_mean = expression(Salinity~mean~(ppt)),
  so_range = expression(Salinity~range~(ppt)),
  si_mean = expression(Silicate~(mmol~m^{-3})),
  chl_mean = expression(Chlorophyll~(mmol~m^{-3})),
  KDPAR_mean_mean = expression(KD~PAR~(m^{-1})),
  phyc_mean = expression(Phytoplankton~(mmol~m^{-3}))
)


### Combine species and environmental data ----

# Extract raster values to points
env <- raster::extract(s1,resp_sf) %>%
  as.data.table()

# Combine data and reorder columns
pshp <- response %>%
  bind_cols(env %>%
              select(all_of(keeprast))) %>% # selected columns
  drop_na()                                 # Keep complete cases

## Final data as spatial
pshp_sf <- st_as_sf(pshp,coords=c('x','y'),crs=CRS('+proj=longlat +datum=WGS84 +no_defs'))

### Variable name lists as objects ----

# Predictor variable names including latitude and longitude
prnames <- pshp %>%
  select(x,y,Bathymetry:last_col()) %>%
  names()
facvars <- '' # Insert names of any factor variables
numvars <- prnames[prnames != facvars]


### Save model data and variables ----
save(pshp,pshp_sf,s1,envlab,prnames,facvars,numvars, 
     file = paste0(out_path,'EpifaunaMetricsModelData.RData'))

########## FROM HERE REPEAT FOR EACH RESPONSE VARIABLE #########################

## Response variable 1: S_epifauna
## Response variable 2: N_trans_epifauna

######### 3. Data Exploration ##################################################

### Select variables ----

## Select response variable
tax = 'S_epifauna'

## Select columns with response variable and environmental variables 
## to include as potential predictors
cols <- c(tax,prnames)
nrow(pshp)
temp <- copy(pshp %>%
               select(all_of(cols)))
temp <- temp[complete.cases(temp),]

# Rename the response column
setnames(temp,1,'resp')

summary(temp)
nrow(temp)

# Remove unnecessary range variables
sdata <- temp %>%
  select(resp,x,y,all_of(keeprast))

summary(sdata)

# Reassign variable lists after removal of variables
prnames <- names(sdata)[-1] # All predictors
facvars <- '' # List predictors that are factor variables
numvars <- prnames[prnames != facvars] # List predictors that are numerical variables

### Summary statistics  ----

resp.long.num <- sdata %>%
  select(where(is.numeric),-resp) %>%
  pivot_longer(cols=everything(),names_to = 'variable')

sumstats <- resp.long.num %>%
  group_by(variable) %>%
  summarise(across(everything(), list(count = ~sum(!is.na(.)),
                                      min= ~round(min(., na.rm = TRUE),2),
                                      mean = ~round(mean(., na.rm = TRUE),2),
                                      median = ~round(median(., na.rm = TRUE),2),
                                      max= ~round(max(., na.rm = TRUE),2)))) %>%
  left_join(varlabs %>%
              select(Variable,Label),by=c('variable'='Variable')) %>%
  relocate(last_col(), .after = variable) %>%
  setnames(c('Variable','Label','N','Min','Mean','Median','Max')) %>%
  mutate(Label=case_when(is.na(Label) ~ Variable,
                         TRUE ~ Label))

sumstats %>%
  select(-Label) %>%
  kable("html", caption = "Table 1. Summary Statistics for environmental variables.") %>%
  kable_classic(full_width = F, html_font='Merriweather',position = "left",fixed_thead = T) %>%
  column_spec(1,  border_right = T) %>%
  row_spec(0, bold = T)


stview(dfSummary(sdata,
                 plain.ascii  = FALSE, 
                 style        = "grid", 
                 graph.magnif = 0.75, 
                 valid.col    = FALSE,
                 display.labels = FALSE,
                 varnumbers = FALSE,
                 na.col = FALSE,
                 tmp.img.dir  = "/tmp"))

## Build histogram
indat <- ggplot(sdata,(aes(x=resp))) +
  geom_histogram(fill=colpal[2], col=colpal[1]) +
  xlab(paste(tax)) +
  theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_text(vjust=-1, size=12,colour="black"),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="#f5f1f1"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"),
        plot.margin = ggplot2::margin(3,3,8,3))

# Get the locations of the panel grid lines for fitting table
hmax.y <- quantile(ggplot_build(indat)$layout$panel_params[[1]]$y.range,0.95)
hmax.x <- quantile(ggplot_build(indat)$layout$panel_params[[1]]$x.range,0.95)

indat +
     annotate("text",
               label = paste("N =", nrow(sdata)),
               x = hmax.x,
               y = hmax.y,
               size = 5,
               hjust = 1)

### Co-variance between environmental variables  ----

# Correlation
corr <- sdata %>%
  select(all_of(numvars),-resp) %>%
  cor()
# Correlation plot
corplot <- ggcorrplot(corr, 
                      method='square',
                      type = 'upper',
                      lab = TRUE,
                      digits = 1,
                      show.legend = FALSE,
                      hc.order = TRUE,
                      sig.level = 0.95) +
  scale_y_discrete(labels = function(x) sapply(x, function(i) label_map[[i]])) +
  scale_x_discrete(labels = function(x) sapply(x, function(i) label_map[[i]]))
corplot


######### 4. Preliminary model #################################################

## Preliminary single model using all environmental variables are run first 
## to select best variables.Variables are dropped from the models based on 
## redundancy (high correlation with a more important variable), or poorly 
## defined relationship with the response variable (based on response curves).


###  Build Full model ----
prelRF <- randomForest(resp~.,
                       data=sdata, 
                       ntrees=1000,
                       replace=FALSE,
                       importance=TRUE,
                       nPerm=5)
prelRF

### Check Variable Importance ----
full.importance <- data.table(Predictor=rownames(prelRF$importance),prelRF$importance)
full.importance <- full.importance[order(full.importance[,IncNodePurity],decreasing=T),]
full.importance

### Plot Partial Dependence for All Variables ----

plotdata <- NULL        ## Object to store data
predselnf <- numvars    ## names of numeric variables

## Loop through all numeric predictor variables and combine response plot
## data into one table
for (j in 1:length(predselnf)){
  pdata <- partial(prelRF,pred.var = predselnf[j],plot = FALSE,train=sdata,grid.resolution=50)
  predname <- predselnf[j]
  temp <- data.frame(predvar=rep(predname,length(pdata[[1]])),x=pdata[[1]],y=pdata[[2]])
  plotdata <- rbind(plotdata,temp)
  print(paste(predselnf[j],'done'))
  next
}


## Set up list object for plots
fullRP.list <- list()

## Loop through each variable and add partial response plot to list
for (i in numvars) {
  
  fullRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y)) +
    geom_smooth(linewidth=0.8,se=FALSE,span = 0.3,col=hlp[3]) +
    facet_wrap(~ predvar,scales = "free", ncol=3) +
    ylim(c(min(plotdata$y),max(plotdata$y))) +
    theme(axis.title.y = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
          axis.title.x  = element_blank(),
          plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
          strip.background = element_rect(fill="grey90"),
          strip.text.x = element_text(size=12, face="bold"),
          panel.grid.major = element_line(colour="grey80"),
          panel.grid.minor = element_line(colour="grey80"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=12,colour="black"))
}


fullRP <- wrap_plots(fullRP.list,ncol = 4,guides = 'collect') + plot_layout(axis_titles='collect')
fullRP


### Variable selection ----

## Select uncorrelated variables to keep
# Variable list
vl <- full.importance[[1]]
# Numeric variables only, leaving out coordinates
vl <- vl[!vl %in% c("x", "y")]
# Correlation matrix with variables in order of full model importance
cr <- cor(sdata %>%
            select(all_of(vl)))
# Remove variables correlated to a higher importance variable
for(j in 1:length(cr[1,])){
  if (j == 1){
    pl <- c(names(cr[j,][1]),names( cr[j,][sqrt((cr[j,])^2)<0.6]))
    pl1 <- pl
  } else if (names(cr[j,])[j] %in% pl1){
    rem <- names(cr[j,-c(1:j)][sqrt((cr[j,-c(1:j)])^2)>0.6])
    if (length(rem) != 0L){  
      pl <- pl[!pl %in% rem]
    }
  }
  next
}

## Calculate Variance Inflation Factors (vif)
crval <- as.data.frame(pl)
crval[2] <- corvif(sdata %>%
                     select(all_of(pl)))
crval

## Choosing final set of predictor variables
predsel <- crval[[1]]

## List of all variables including response
clms <- c(names(sdata)[1],predsel)

### Data to use in model ----
mdata <- sdata %>%
  select(all_of(clms))
summary(mdata)


### Save preliminary model, model data  ----
save(prelRF,full.importance,fullRP.list, 
     file = paste0(out_path,'Preliminary_Model_',tax,'.RData'))
save(sdata,mdata,clms,predsel,
     file = paste0(out_path,'Final_Data_',tax,'.RData'))


######### 5. Repeated subset crossvalidation models ############################

### The final random forest models with the subset of predictor variables 
### defined in the previous step. Variable importance is determined through 
### a multiple permutation procedure during the model run. 
### Cross-validation via repeated sub-sampling is done to evaluate the 
### robustness of the model estimate and predictions to data sub-setting. 
### This consists of 10 random split samples with 75% used to train and 25% 
### to test models. The final model outputs the majority vote of all 10 runs. 
### A confidence map layer is made using the product of frequency (runs out  
### of 10), and the mean probability of the majority class

###  Prepare repeated model runs ----

## Raster data as a dataframe
s1d <- raster::as.data.frame(s1)
s1d <- s1d[complete.cases(s1d),predsel]

## Insert number of crossvalidation runs required
nruns <- 10 

## Variable names
preds <-   predsel
facvars <- NULL
predselnf <- predsel[!predsel %in% facvars] 

## List objects for saving results from repeated runs
ffs <- list() # Random Forest models
imps <- list() # Variable importance
res <- list() # Results

## Objects for saving data for plotting
plotdata <- NULL
predvsobs <- NULL

## Table for collecting all model performance statistics
results <- data.frame(Name=character(0),Run=character(0),
                      N=numeric(0),ExtRMSE=numeric(0),
                      ExtRMSEpctR=numeric(0),ExtR2=numeric(0),
                      stringsAsFactors =F)


###  Prepare repeated subsample datasets ----

## Set up lists for training and test datasets
train.sets <- list()
test.sets <- list()

## Splits for 10 random subsets (list of row numbers)
trainIndex <- createDataPartition(mdata$resp, p = .75,
                                  times = nruns)
## Create training and tests sets
for (j in 1:10){
  
  train.sets[[j]] <- mdata[trainIndex[[j]],]
  test.sets[[j]] <- mdata[-trainIndex[[j]],]
  
  next}

## Save training and test data subsets
save(train.sets,test.sets,
     file = paste0(out_path,'TrainTest_',tax,'.RData'))


###  Loop through crossvalidation models for each subset  ----

for (j in 1:10){
  
  train <- train.sets[[j]]
  test <- test.sets[[j]]
  
  ffs[j] <- list(forest <- randomForest(resp ~.,data=train,
                                        ntree=500, replace=FALSE,importance=T, keep.forest= T))
  
  
  predicted <-  predict(ffs[[j]],test,'response')
  observed <- test[[1]]
  
  plot(predicted,observed)
  predvsobs <- rbind(predvsobs,cbind(predicted,observed))
  
  results[j,1] <- names(train)[1]
  results[j,2] <- j
  results[j,3] <- nrow(train)
  results[j,4] <- round(postResample(predicted,observed)[[1]],2)
  results[j,5] <- round(postResample(predicted,observed)[[1]]/(max(train[,1])-min(train[,1])),2)
  results[j,6] <- round(postResample(predicted,observed)[[2]],2)
  
  
  imps[[j]] <- list(round(randomForest::importance(ffs[[j]]), 2))
  
  require(pdp)
  
  for (i in 1:length(predselnf)){
    pdata <- partial(ffs[[j]],pred.var = predselnf[i],plot = FALSE,train=train,grid.resolution=50)
    predname <- predselnf[i]
    temp <- data.frame(run=j,predvar=rep(predname,length(pdata[[1]])),x=pdata[[1]],y=pdata[[2]])
    plotdata <- rbind(plotdata,temp)
    next
  }
  
  print(paste('Run', j, 'done!')) 

  next
}


results$Name <- tax
forest.resultsSQ <- rbind(forest.resultsSQ,results)

###  Save model list, plotting data, importances and validation results  ----
save(ffs,plotdata,forest.resultsSQ,imps,
     file = paste0(out_path,'ForestResult_',tax,'.RData'))

#### Compile model performance statistics tables ----

# Performance table
model.perf <- data.table(Response=tax,
                       N=mean(results$N),
                       RMSE= paste(round(mean(results$ExtRMSE),2), '\u00B1',round(sd(results$ExtRMSE),2)),
                       'Relative RMSE' = paste(round(mean(results$ExtRMSEpctR),2), '\u00B1', round(sd(results$ExtRMSEpctR),2)),
                       'R^2^'= paste(round(mean(results$ExtR2),2), '\u00B1',round(sd(results$ExtR2),2)))
pander::pandoc.table(model.perf)


#### Model performance plots ----

### Observed vs. predicted plot

## Prepare data
predvsobs <- as.data.table(predvsobs)
limvals <- max(plyr::round_any(predvsobs$observed,10),plyr::round_any(predvsobs$predicted,10))

## Create plot
obprplot <- ggplot(predvsobs,aes(x=predicted,y=observed)) +
  geom_point(shape = 16, colour = hlp[1]) +
  geom_smooth(inherit.aes=FALSE,aes(x=predicted,y=observed),
              method='lm',linewidth=0.9,se=TRUE,span = 0.5,col= hlp[3]) +
  coord_equal(xlim = c(min(predvsobs),signif(max(predvsobs),1)),
              ylim = c(min(predvsobs),signif(max(predvsobs),1))) +
  scale_x_continuous(breaks = scales::breaks_extended(n = 7)) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 7)) +
  ylab("Observed") +
  xlab("Predicted") +
  theme(axis.title.y = element_text(vjust=4,  size=12,colour="black"),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_text(vjust=-4, size=12,colour="black"),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"),
        plot.margin = ggplot2::margin(0.5, 0.5, 1, 0.5, "cm"))


ggsave(obprplot,
       file = paste0(out_path,'Observed_vs_Predicted_',tax,'.png'),
       device = 'png', width = 20, height=22, units='cm', dpi=300, scale=1)


### Variable importance plot

## Set up data
imppl <- data.table(Var=rownames(imps[[1]][[1]]))

for (i in 1:10){
  
  imppl <- cbind(imppl,as.data.table(imps[[i]][[1]])[,1])
  
}

setnames(imppl,c('Var','Imp1','Imp2','Imp3','Imp4','Imp5','Imp6','Imp7','Imp8','Imp9','Imp10'))

imppl[,
      c("Mean",'Sd','Se') := 
        .(rowMeans(.SD, na.rm = TRUE), 
          apply(.SD, 1, sd, na.rm = TRUE),
          apply(.SD, 1, plotrix::std.error, na.rm = TRUE)), 
      .SDcols = 2:11]

imppl[,Var:=factor(Var,levels=Var[order(Mean)])]

## Create plot
impplot <-  ggplot(imppl,(aes(x=Var,y=Mean))) +
  geom_bar(stat = 'identity',fill=colpal[2], col=colpal[1]) +
  scale_x_discrete(labels = function(x) sapply(x, function(i) label_map[[i]])) +
  geom_linerange(inherit.aes=FALSE,
                 aes(x=Var, ymin=Mean-Se, ymax=Mean+Se), 
                 colour='#172957', alpha=0.9, linewidth=1.3) +
  ylab(label = '% Increase in MSE') +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5,hjust = 1, size=12,colour="black"),
        axis.text.x  = element_text(vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_text(vjust=-4, size=12,colour="black"),
        plot.title =  element_text(size=12,colour="white", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey90"),
        strip.text.x = element_text(size=12, face="bold"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.text = element_text(size=12,colour="black"),
        plot.margin = ggplot2::margin(0.5, 0.5, 1, 0.5, "cm"),)

impplot

## Write plot
ggsave(impplot,
       file = paste0(out_path,'Variable_Importane_',tax,'.png'),
       device = 'png', width = 30, height=18, units='cm', dpi=300, scale=1)


#### Partial Response Plots ----

## Prepare data
plotdata <- plotdata %>%
  mutate(across(y, \(x) round(x, 1)))


## List object for storing plots
cvRP.list <- list()

## Labels for facet strips (used as x-axis titles)
envlabeller <- as_labeller(envlab,default = label_parsed)
ylims <- plotdata %>%
  summarise(
    min_y = min(y, na.rm = TRUE),
    max_y = max(y, na.rm = TRUE)
  )

## Loop through variables and plot
for (i in predsel) {
  
  cvRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y,group=run)) +
    geom_smooth(method='loess',linewidth=0.01,se=FALSE,span = 0.2,col=hlp[1],) +
    geom_smooth(inherit.aes=FALSE,aes(x=x,y=y),
                method='loess',linewidth=0.9,se=FALSE,span = 0.2,col=hlp[3]) +
    facet_wrap(~ predvar,scales = "free_x",ncol =3,
               labeller = envlabeller,
               strip.position = 'bottom') +
    scale_y_continuous(limits = c(min(plotdata$y),max(plotdata$y))) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y  = element_text(vjust=0.5, size=10,colour="black"),
          axis.text.x  = element_text(vjust=0.5, size=10,colour="black"),
          legend.title = element_blank(),
          legend.key = element_rect(fill = NA),
          legend.text = element_text(size=10,colour="black"),
          strip.placement = 'outside',
          strip.background = element_blank(),
          strip.text = element_text(size=10,face = 'bold'),
          plot.margin = unit(c(0.5, 0.5, 0,0), "cm"))

}


## Layout of top six predictor variable plots
cvRP <-  wrap_plots(cvRP.list) + plot_layout(ncol=3) 
cvRP

## Write plot
ggsave(cvRP,
       file = paste0(out_path,'Partial_Response_',tax,'.png'),
       device = 'png', width = 30, height=30, units='cm', dpi=300, scale=1)


#### Save main results ----

## Validation results and plots
save(results,model.perf,impplot,cvRP,cvRP.list,obprplot,
     file = paste0(out_path,'Model_Results_',tax,'.RData'))

## Validation results to csv
write.csv(model.perf, file = paste0(out_path,'Model_Performance_',tax,'.csv'))

## Set objects for collating model results for all response variables
# Model data 
RD <- rbind(RD,data.table(Tax=tax,Metric='diversity',value=sdata[,'resp']))
# Importance plot
VIs <- rbind(VIs,data.table(Tax=tax,imppl[,.SD,.SDcols=c('Var','Mean','Se')]))  
## Performance stats---
PSs <- rbind(PSs,data.table(Tax=tax,transpose(mod.perf,keep.names = 'Stat')))
# Pred vs. obs
OvsPs <- rbind(OvsPs,data.table(Tax=tax,predvsobs))
# Save objects
save(VIs,forest.resultsSQ,RD,PSs,OvsPs,PRCs,envlab,
     file = paste0(out_path,'Model_Result_Objects',tax,'.RData'))


######### 6. Make model predictions ############################################

#### Drop unnecessary predictor layers  ----
dr <- names(s1)
dr <- dr[!dr %in% names(ffs[[1]]$forest$xlevels)]
s2 <- dropLayer(s1, dr)

#### Set up objects for saving predictions ----
cvpred <- NULL


#### Loop through each model and predict to raster ----

for (i in 1:length(ffs)){
  
  # Predict
  
  if (is.null(cvpred)){
    cvpred  <- raster::predict(s2,ffs[[i]])
  } else {
    cvpred <- addLayer(cvpred,predict(s2,ffs[[i]]))
  }
  
  next
}


#### Final outputs ----

## Calculate mean and standard deviation for predictions from CV runs
predres <- stack(calc(cvpred,fun=mean),
                 calc(cvpred,fun=sd))
names(predres) <- c('meanpred','sdpred')
## Calculate coefficient of variance
predres$relvar <- predres$sdpred/predres$meanpred

## Summary statistics for the rasters
summary(predres$meanpred)


#### Export Rasters ----

## Mean of predicted values
writeRaster(predres$meanpred,
            file = paste0(out_path,'RF_Mean_',tax,'.tif'),
            format="GTiff",overwrite=T)
## Standard deviation of predicted values
writeRaster(predres$sdpred,
            file = paste0(out_path,'RF_SD_',tax,'.tif'),
            format="GTiff",overwrite=T)
## Coefficient of variance of predicted values
writeRaster(predres$relvar,
            file = paste0(out_path,'RF_CV_',tax,'.tif'),
            format="GTiff",overwrite=T)

################################################################################
