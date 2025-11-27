################################################################################
######### RANDOM FOREST MODELS - TRAWL EPIFAUNA CLUSTERS   #####################
#########               25/11/2025                         #####################
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
colpal  <- c("#0000ee","#00ffff","#05aac1","#9a32cd","#00cd00", "#9aff9a",
             "#b40202","#ff0000","#ff8c00","#ffff00","#b4b404") # Define specific colours

## Highlight colours for confusion matrix plot
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


# Tables for collecting all model performance statistics
# Class specific validation statistics
forest.class.res <- data.frame(Name=character(0),
                               Run=character(0),
                               Class=character(0),
                               ClassN=numeric(0),
                               Sens=numeric(0),
                               Spec=numeric(0),
                               BA=numeric(0),
                               stringsAsFactors =F)

# Validation statistics for whole model
forest.class.res.all <- data.frame(Name=character(0),
                                   Run=character(0),
                                   N=numeric(0),
                                   Acc=numeric(0),
                                   NIR=numeric(0),
                                   P=numeric(0),
                                   Kappa=numeric(0),
                                   Q=numeric(0),
                                   A=numeric(0),
                                   stringsAsFactors =F)

# Combined validations statistic
forest.res.mat.all <- data.frame(Name=character(0),
                                 ModRun = numeric(),
                                 Comb=character(),
                                 Pred=character(),
                                 Obs=character(),
                                 Vals=numeric(),
                                 stringsAsFactors = FALSE)


######### 2. Read in Data ######################################################

### Response point data  ----

# Read in csv, pivot to wide format and rename coordinate columns
response <- as.data.table(read.csv(paste0(in_path,rd))) %>%
  filter(metric=='assemblage_epifauna') %>%
  mutate_at(c('value'), as.factor) %>%
  select(value,x,y) %>%
  rename(c(`Cluster`='value'))
summary(response)

# Convert to spatial
resp_sf <- st_as_sf(response,coords=c('x','y'),crs=CRS('+proj=longlat +datum=WGS84 +no_defs'))


### Environmental predictor rasters  ----

# Read in all rasters in folder
f <- list.files(path=rs, pattern='.tif$', full.names=T)
s1 <- stack(f,RAT=F)
s1
plot(s1)


### Environmetal variable labels -----

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
            select(x:last_col()) %>%
            names()
facvars <- '' # Insert names of any factor variables
numvars <- prnames[prnames != facvars]

# Response variable names
rspnames <- c("Cluster")

# Number of classes in the response variable
numclass <- nlevels(pshp[[1]])


### Save model data and variables ----
save(pshp,pshp_sf,s1,envlab,prnames,facvars,numvars,rspnames, 
     file = paste0(out_path,'EpifaunaClusterModelData.RData'))


######### 3. Data Exploration ##################################################

### Select variables ----

## Select response variable
tax = 'Cluster'

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

resp.long.num <- pshp %>%
  select(where(is.numeric)) %>%
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


stview(dfSummary(pshp,
                  plain.ascii  = FALSE, 
                  style        = "grid", 
                  graph.magnif = 0.75, 
                  valid.col    = FALSE,
                  display.labels = FALSE,
                  varnumbers = FALSE,
                  na.col = FALSE,
                  tmp.img.dir  = "/tmp"))

## Plot class data distribution
mxy <- max(summary(sdata$resp))

## Build histogram
indat <- ggplot(sdata,(aes(x=resp,fill=resp))) +
  geom_bar() +
  annotate(geom='text',label = paste('N =',nrow(sdata)), x = 1,
                                 y = mxy,size=5,hjust=0) +
  scale_fill_manual(values = colpal) +
  xlab('Class') +
  ylab('Number of samples') +
  theme(axis.title.y = element_text(vjust=0.5,size=12,colour="black"),
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

indat

### Co-variance between environmental variables  ----

# Correlation
corr <- pshp %>%
  select(all_of(numvars)) %>%
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
full.importance <- full.importance[order(full.importance[,MeanDecreaseGini],decreasing=T),]
full.importance

### Plot Partial Dependence for All Variables ----

plotdata <- NULL        ## Object to store data
predselnf <- numvars    ## names of numeric variables

## Loop through all numeric predictor variables and combine response plot
## data into one table
for (j in 1:length(predselnf)) {
  
  for (c in levels(sdata$resp)) {
    
    pdata <- partial(prelRF,pred.var = predselnf[j],which.class = c,
                     plot = FALSE,train=sdata,grid.resolution=50,prob = TRUE)
    predname <- predselnf[j]
    temp <- data.frame(predvar=predselnf[j],class=c,x=pdata[[1]],y=pdata[[2]])
    plotdata <- rbind(plotdata,temp)
    
    
  }
  print(paste(predselnf[j],'Done!'))

}


## Set up list object for plots
fullRP.list <- list()

## Loo through each variable and add partial response plot to list
for (i in predselnf) {
  
  fullRP.list[[i]] <- ggplot(plotdata[plotdata$predvar==i,],aes(x=x,y=y,col=class)) +
    geom_smooth(linewidth=0.8,se=FALSE,span = 0.3) +
    facet_wrap(~ predvar,scales = "free_x", ncol=3) +
    scale_colour_manual(values = colpal) +
    ylim(c(min(c(0,plotdata$y)),max(plotdata$y))) +
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


fullRP <- wrap_plots(fullRP.list,ncol = 5,guides = 'collect')
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
    pl <- c(names(cr[j,][1]),names( cr[j,][sqrt((cr[j,])^2)<0.7]))
    pl1 <- pl
  } else if (names(cr[j,])[j] %in% pl1){
    rem <- names(cr[j,-c(1:j)][sqrt((cr[j,-c(1:j)])^2)>0.7])
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

## Object for saving data for plotting
plotdata <- NULL

## Tables for collecting all model performance statistics
# Class specific validation statistics
class.res <- data.frame(Name=character(0),Run=character(0),Class=character(0),ClassN=character(0),
                        Sens=numeric(0),Spec=numeric(0),BA=numeric(0),
                        stringsAsFactors =F)

# Validation statistics for whole model
class.res.all <- data.frame(Name=character(0),
                            Run=character(0),
                            N=character(0),
                            Acc=numeric(0),
                            NIR=numeric(0),
                            P=numeric(0),
                            Kappa=numeric(0),
                            Q=numeric(0),
                            A=numeric(0),
                            stringsAsFactors =F)

# Combined validations statistic
res.mat.all <- data.frame(Name=character(0),ModRun = numeric(),Comb=character(),
                          Pred=character(),Obs=character(),Vals=numeric(),
                          stringsAsFactors = FALSE)


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
  
  # Get data subset j
  train <- train.sets[[j]]
  test <- test.sets[[j]]
  
  # Random Forest model j
  ffs[[j]] <- randomForest(resp ~.,data=train,
                           ntree=500, replace=FALSE,importance=T, keep.forest= T)
  
  # Set data frame with up observed values for j
  results <- as.data.frame(1:nrow(test))
  results$actual <- test[[1]]
  # Predict class with model j
  results$predicted <- as.data.table(predict(ffs[j],test))[,V1]
  names(results) <- c("id", "actual", "predicted")
  
  # Calculate confusion matrix for predictions by model j
  results.matrix <- confusionMatrix(results$predicted, results$actual)
  results.matrix
  
  # Get the number of objects predicted into each class by model j
  temp.pred <- NULL
  s1dnrow <- nrow(s1d)
  rs <- data.frame(start=c(1,ceiling(s1dnrow/3),2*ceiling(s1dnrow/3)),
                   end=c(ceiling(s1dnrow/3)-1,2*ceiling(s1dnrow/3)-1,s1dnrow))
  for (e in 1:nrow(rs)){
    
    ts1d <- s1d[rs[e,'start']:rs[e,'end'],]
    temp.pred <- rbind(temp.pred,as.data.table(predict(ffs[j],ts1d,'response')))
    
  } 
  
  pctObj <- as.data.frame(summary(temp.pred[,V1])/sum(summary(temp.pred[,V1])))
  
  # Get values from confusion matrix for model j into a data frame for future plotting
  rnx <- 1
  res.mat.0 <- data.frame(Name=character(), ModRun = numeric(),Comb=character(),
                          Pred=character(),Obs=character(),Vals=numeric(),
                          ValsP=numeric(),
                          stringsAsFactors = FALSE)
  
  for (r in 1:dim(results.matrix$table)[1]){
    
    for (o in 1:dim(results.matrix$table)[2]){
      rn <- rnx+o-1
      res.mat.0[rn,2] <- j
      labt <- paste(r,o,sep="") 
      res.mat.0[rn,3] <- labt
      res.mat.0[rn,4] <- levels(train[[1]])[r]
      res.mat.0[rn,5] <- levels(train[[1]])[o]
      res.mat.0[rn,6] <- results.matrix$table[r,o]
      next
    }
    
    CurSel <- res.mat.0[(rn-o+1):rn,]
    PrVals <- res.mat.0[(rn-o+1):rn,6]
    totPrCL <- sum(res.mat.0[(rn-o+1):rn,6])
    
    res.mat.0[(rn-o+1):rn,7] <- round((PrVals/totPrCL) * pctObj [[1]][r],3)
    CurSel <- res.mat.0[(rn-o+1):rn,]
    
    labt <- paste(r,'P',sep="") 
    rn <- rn+1
    res.mat.0[rn,2] <- j
    res.mat.0[rn,3] <- labt
    res.mat.0[rn,4] <- levels(train[[1]])[r]
    res.mat.0[rn,5] <- 'P'
    res.mat.0[rn,6] <- round(100*(results.matrix$table[r,r]/rowSums(results.matrix$table)[r]),1)
    res.mat.0[rn,7] <- round(100*(CurSel[CurSel[4]==CurSel[5],7]/sum(CurSel[7])),3)
    
    rnx <- rnx+dim(results.matrix$table)[2]+1
    next
  }
  
  
  
  for (o in 1:dim(results.matrix$table)[2]){
    ObsClassVal <- res.mat.0[res.mat.0[2]==j & res.mat.0[5]==levels(train[[1]])[o] ,]
    rn <- rnx+o-1
    res.mat.0[rn,2] <- j
    labt <- paste('U',o,sep="") 
    res.mat.0[rn,3] <- labt
    res.mat.0[rn,4] <- 'U'
    res.mat.0[rn,5] <- levels(train[[1]])[o]
    res.mat.0[rn,6] <- round(100*(results.matrix$table[o,o]/colSums(results.matrix$table)[o]),1)
    res.mat.0[rn,7] <- round(100*(ObsClassVal[ObsClassVal[4]==ObsClassVal[5],7]/sum(ObsClassVal[7])),1)
    next
  }
  
  RunValAll <- res.mat.0[res.mat.0[2]==j,]
  RunValAll$Pred <- factor(RunValAll$Pred,levels = c(levels(train[[1]]),"U","P"))
  RunValAll$Obs <- factor(RunValAll$Obs,levels = c(levels(train[[1]]),"U","P"))
  RunValAll$ValsP[is.nan(RunValAll$ValsP) & RunValAll$Pred !="U" & RunValAll$Pred !="P"] <- 0
  RunValNum <- RunValAll[RunValAll[4]!="U" & RunValAll[5]!="P",]
  RunValCor <- RunValNum[RunValNum[4]==RunValNum[5],]
  Psum <- aggregate(ValsP~Pred,data=RunValNum,sum,na.action=na.pass)
  Osum <- aggregate(ValsP~Obs,data=RunValNum,sum,na.action=na.pass)
  
  res.mat.0[rn+1,2] <- j
  labt <- paste('U','P',sep="") 
  res.mat.0[rn+1,3] <- labt
  res.mat.0[rn+1,4] <- 'U'
  res.mat.0[rn+1,5] <- 'P'
  res.mat.0[rn+1,6] <- round(100*results.matrix[[3]][[1]],1)
  res.mat.0[rn+1,7] <- round(100*(sum(RunValCor[7])/sum(RunValNum[7])),3)
  
  res.mat.all <- rbind(res.mat.all,res.mat.0)
  
  require(matrixStats)
  sum(2*rowMins(as.matrix(cbind(Osum[2]-RunValCor[7],Psum[2]-RunValCor[7]))))/2
  
  # Get overall accuracy measures for model validation run j
  class.res.all[j,2] <- j
  class.res.all[j,3] <- nrow(test)
  class.res.all[j,4] <- results.matrix[[3]][[1]]
  class.res.all[j,5] <- results.matrix[[3]][[5]]
  class.res.all[j,6] <- results.matrix[[3]][[6]]
  class.res.all[j,7] <- results.matrix[[3]][[2]]
  class.res.all[j,8] <- sum(abs(Osum[2]-Psum[2]))/2
  class.res.all[j,9] <- sum(2*rowMins(as.matrix(cbind(Osum[2]-RunValCor[7],
                                                      Psum[2]-RunValCor[7]))))/2
  
  class.res.0 <- data.frame(Name=character(0),Run=character(0),Class=character(0),ClassN=character(0),
                            Sens=numeric(0),Spec=numeric(0),BA=numeric(0),
                            stringsAsFactors =F)
  
  # Get class-specific accuracy measures for model validation run j
  for (i in 1:numclass){
    
    class.res.0[i,2] <- j
    class.res.0[i,3] <- row.names(as.data.frame(results.matrix[[4]]))[i]
    class.res.0[i,4] <- sum(results.matrix$table[,i])
    class.res.0[i,5] <- as.data.frame(results.matrix[[4]])[i,1]
    class.res.0[i,6] <- as.data.frame(results.matrix[[4]])[i,2]
    class.res.0[i,7] <- as.data.frame(results.matrix[[4]])[i,11]
    
    next
  }
  
  class.res <- rbind(class.res,class.res.0)
  
  
  ## Store Validation Results in main tables
  class.res$Name <- tax
  class.res
  
  
  class.res.all$Name <- tax
  class.res.all
  
  
  res.mat.all$Name <- tax
  res.mat.all
  
  
  imps[[j]] <- list(round(randomForest::importance(ffs[[j]]), 2))
  
  require(pdp)
  
  for (p in 1:length(predsel)) {
    
    for (c in levels(mdata$resp)) {
      
      pdata <- partial(ffs[[j]],pred.var = predsel[p],which.class = c,
                       plot = FALSE,train=mdata,grid.resolution=50,prob = TRUE)
      predname <- predsel[p]
      temp <- data.frame(Name='Cluster',predvar=predsel[p],class=c,x=pdata[[1]],y=pdata[[2]])
      plotdata <- rbind(plotdata,temp)
      
    }
    
  }
  
  if (!is.null(facvars)) {
    
    for (f in 1:length(facvars)){
      
      for (c in levels(mdata$resp)) {
        
        pdata <- partial(ffs[[j]],pred.var = facvars[f],which.class = c,
                         plot = FALSE,train=mdata,grid.resolution=50,prob = TRUE)
        predname <- facvars[j]
        temp <- data.frame(Name='Cluster',predvar=facvars[f],class=c,x=pdata[[1]],y=pdata[[2]])
        plotdata2 <- rbind(plotdata2,temp)
        
      }
      
    }
    
    
  }
  
  print(paste('Run',j,'done!'))
  
  next
}

forest.class.res <- rbind(forest.class.res,class.res)
forest.class.res.all <- rbind(forest.class.res.all,class.res.all)  
forest.res.mat.all <- rbind(forest.res.mat.all,res.mat.all)

###  Save model list, plotting data, importances and validation results  ----
save(ffs,plotdata,class.res,class.res.all,res.mat.all,imps,
     file = paste0(out_path,'ForestResult_',tax,'.RData'))

#### Compile model performance statistics tables ----

## Separate the main matrix and producers, users and overall accuracy into
## different data frames
res.mat.num <- res.mat.all[res.mat.all[4]!="U" & res.mat.all[5]!="P",]
res.mat.num <- transform(res.mat.num, 
                         Pred = factor(Pred, levels=levels(mdata[[1]])), 
                         
                         Obs = factor(Obs, levels=levels(mdata[[1]])))
res.mat.up <- res.mat.all[res.mat.all[4]=="U" | res.mat.all[5]=="P",]
res.mat.up <- transform(res.mat.up, 
                        Pred = factor(Pred, levels = c(levels(mdata[[1]]),"U","P")),
                        Obs = factor(Obs, levels = c(levels(mdata[[1]]),"U","P")))
res.mat.up[res.mat.up$Comb=="UP","Comb"] <- "OA"
highlights <- res.mat.num[res.mat.num[4]==res.mat.num[5],][1:numclass,4:6]


## Compare between classes

# Rename classes for class specific results for consistency
str(class.res)
class.res$Class <- as.factor(class.res$Class)
levels(class.res[[3]])
levels(class.res[[3]]) <- levels(mdata[[1]])
class.res$Run <- as.numeric(class.res$Run)
class.res$ClassN <- as.numeric(class.res$ClassN)

# Calculate overall true positives and sensitivity for each model run
TP <- aggregate(Vals~ModRun,data=res.mat.num[res.mat.num$Pred == res.mat.num$Obs,],sum)
SEN <- TP[2]/aggregate(Vals~ModRun,data=res.mat.num,sum)[2]

# Calculate overall true negatives and specififity for each model run
negs <- rep(0,length(unique(res.mat.num[1])[[1]]))
tnegs <- rep(0,length(unique(res.mat.num[1])[[1]]))

for (i in unique(res.mat.num[4])[[1]]){
  
  inegs <- aggregate(Vals~ModRun,data=res.mat.num[res.mat.num$Obs != i,],sum)
  negs <- negs+inegs[[2]]
  itnegs <- aggregate(Vals~ModRun,data=res.mat.num[res.mat.num$Obs != i & res.mat.num$Pred != i ,],sum)
  tnegs <- tnegs+itnegs[[2]]
  
  
  next
}

SPE <- tnegs/negs

# Calculate overall balanced accuracy for each model run   
BA <- (SEN+SPE)/2

# Combine values to a matrix
classrestot <- data.frame('Cluster',1:length(SPE),"Overall",nrow(train),SEN,SPE,BA)
names(classrestot) <- names(class.res)

# Add overall accuracy values to class specific table
classvalB <- rbind(class.res,classrestot)

classvalB
str(classvalB)

# Calculate averages and standard deviations for validation statistics across model runs
cavevalsB <- aggregate(x = classvalB[4:7], by = list(classvalB$Class), FUN = "mean")
csdvalsB <- aggregate(x = classvalB[4:7], by = list(classvalB$Class), FUN = "sd")

# Combine values in a table
BCvalsT <- data.frame(Name=cavevalsB[1],
                      N=cavevalsB[2],
                      SENSmean=round(cavevalsB[3],2),
                      SENSsd=round(csdvalsB[3],2),
                      SPECmean=round(cavevalsB[4],2),
                      SPECsd=round(csdvalsB[4],2),
                      BAmean=round(cavevalsB[5],2),
                      BAsd=round(csdvalsB[5],2))
# Rename columns
names(BCvalsT) <- c("Name","N","SENSmean","SENSsd","SPECmean","SPECsd",
                    "BAmean","BAsd")
# Print table
BCvalsT

asg.cl.perf <- data.table(Cluster=BCvalsT$Name,
                          N =BCvalsT$N,
                          'Sensitivity'= paste(BCvalsT$SENSmean, '\u00B1',BCvalsT$SENSsd),
                          'Specificity' =  paste(BCvalsT$SPECmean, '\u00B1',BCvalsT$SPECsd),
                          'Balanced Accuracy'= paste(BCvalsT$BAmean, '\u00B1',BCvalsT$BAsd))

asg.cl.perf[-9, ] %>%
  kbl('html',digits = 2,escape = FALSE, align=c('l',rep('r', 4)),
      caption='Class-specific performance') %>%
  kable_classic(full_width = F, position = "left",fixed_thead = T) %>%
  row_spec(0, bold = T)  %>%
  column_spec(1:2, width = "1.5cm") %>%
  column_spec(3:5, width = "3cm")

### Performance for whole classification
## Prepare data
classvalB.all <- data.frame(class.res.all,Sens=SEN[[1]],Spec=SPE,BA=BA[[1]])
classvalB.all$N <- as.numeric(classvalB.all$N)
# Calculate averages and standard deviations for validation statistics
callavevalsB <- colMeans(classvalB.all[,4:12])
callsdvalsB <- colSds(as.matrix(classvalB.all[,4:12]))
# Combine values in a table
BCallvalsT <- data.frame(Accmean=round(callavevalsB[1],2),
                         Accsd=round(callsdvalsB[1],2),
                         Pmean=round(callavevalsB[3],2),
                         Psd=round(callsdvalsB[3],2),
                         Kmean=round(callavevalsB[4],2),
                         Ksd=round(callsdvalsB[4],2),
                         Qmean=round(callavevalsB[5],2),
                         Qsd=round(callsdvalsB[5],2),
                         Amean=round(callavevalsB[6],2),
                         Asd=round(callsdvalsB[6],2),
                         Sensmean=round(callavevalsB[7],2),
                         Senssd=round(callsdvalsB[7],2),
                         Specmean=round(callavevalsB[8],2),
                         Specsd=round(callsdvalsB[8],2),
                         BAmean=round(callavevalsB[9],2),
                         BAsd=round(callsdvalsB[9],2))

# Rename columns
names(BCallvalsT) <- c("Accmean","Accsd","Pmean","Psd","Kmean","Ksd",
                       "Qmean","Qsd","Amean","Asd","Sensmean","Senssd",
                       "Specmean","Specsd","BAmean","BAsd")


# Print table
BCallvalsT

asg.perf <- data.table(N = nrow(train),
                       'Sensitivity'= paste(BCallvalsT$Sensmean, '\u00B1',BCallvalsT$Senssd),
                       'Specificity' =  paste(BCallvalsT$Specmean, '\u00B1',BCallvalsT$Specsd),
                       'Kappa' = paste(BCallvalsT$Kmean, '\u00B1',BCallvalsT$Ksd) ,
                       'Balanced Accuracy'= paste(BCallvalsT$BAmean, '\u00B1',BCallvalsT$BAsd),
                       'Quantity Disagreement'=paste(BCallvalsT$Qmean, '\u00B1',BCallvalsT$Qsd),
                       'Allocation Disagreement'=paste(BCallvalsT$Amean, '\u00B1',BCallvalsT$Asd))

asg.perf[, data.table(t(.SD), keep.rownames=TRUE),] %>%
  kbl('html',digits = 2,escape = FALSE, col.names = c('Statistic','Mean \u00B1 SD'),
      caption='Class-specific performance') %>%
  kable_classic(full_width = F, position = "left",fixed_thead = T) %>%
  row_spec(0, bold = T)  %>%
  column_spec(1:2, width = "3cm") 


#### Model performance plots ----

### Performance by class
## Prepare data
plotterBC <- classvalB %>%
                  select(Class, Sens, Spec, BA) %>%
                  pivot_longer(
                  cols = c(Sens, Spec, BA),
                  names_to = "variable",
                  values_to = "value"
                ) %>%
              mutate(variable = factor(variable, 
                                       labels = c("Sensitivity", "Specificity",
                                                  "Balanced Accuracy")))


## Create plot
ssbp <- ggplot(plotterBC,(aes(x=Class,y=value))) +
  geom_boxplot(width=0.7,position=position_dodge(width=0.71)) +
  facet_wrap(~ variable) +
  ylim(c(0,1)) +
  ggtitle("Accuracy Statistics Per Class") +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(vjust=0.5, size=12,colour="black"),
        axis.text.x  = element_text(angle=90,vjust=0.5, size=12,colour="black"),
        axis.title.x  = element_blank(),
        plot.title =  element_text(size=14,colour="black", face = "bold",vjust=2),
        strip.background = element_rect(fill="grey40"),
        strip.text.x = element_text(size=14, face="bold",colour="white"),
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_line(colour="grey80"),
        panel.background = element_rect(fill="white"))
ssbp

ggsave(ssbp,
       file = paste0(out_path,'Class_Performance_',tax,'.png'),
       device = 'png', width = 30, height=14, units='cm', dpi=300, scale=1)

### Plot of confusion matrix

## Prepare data
bgvals <- res.mat.num %>%
  group_by(ModRun,Pred) %>%
  mutate(sum= sum(Vals),
         percent = Vals / sum(Vals) * 100)

## Colours for facet strips
strip <- strip_themed(background_x = elem_list_rect(fill = colpal),
                      background_y = elem_list_rect(fill = colpal))

## Theme for swapping axes in plot
theme_swap_y <- function() {
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x.top = element_blank(), # remove ticks/text on labels
        axis.ticks.x.top = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.title.x.bottom = element_blank(), # remove titles
        axis.title.y.left = element_blank())
}

## Confusion matrix data table
res.mat.num2 <- res.mat.all %>%
  filter(Pred!="U" & Obs!="P")


## Percent of values predicted in class
bgvals2 <- res.mat.num2 %>%
  group_by(ModRun,Pred) %>%
  mutate(sum= sum(Vals),
         percent = Vals / sum(Vals) * 100)
## Highlight classes
highlights <- bgvals2 %>%
  group_by(Pred,Obs) %>%
  summarise(PctMean=mean(percent)) %>%
  mutate(fillcol= case_when(
    PctMean>=10 & PctMean<15  ~ '10-15%',
    PctMean>=15 & PctMean<30 ~ '15-30%',
    PctMean>=30  ~ '> 30%',
    TRUE~NA),
    fillcol=factor(fillcol,
                   levels=c('10-15%','15-30%','> 30%'))) %>%
  tidyr::drop_na()


## Define the plot of main confusion matrix
confmtx <-  ggplot(bgvals2) +
  geom_rect(data = highlights,aes(fill = fillcol),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf,alpha = 0.2) +
  geom_boxplot(aes(x=Obs,y=percent),width=0.7,fill='grey30') +
  scale_fill_manual(values = hlp,name='Proportion of predicted') +
  scale_y_continuous(position = 'right', sec.axis = dup_axis(),limits = c(0,100)) + 
  scale_x_discrete(position = "top") +
  facet_grid2(Pred~Obs, scales = 'free',strip = strip) +
  xlab("OBSERVED") +
  ylab("PREDICTED") +
  theme_swap_y() +
  theme(axis.title.y.right =element_text(size=12,colour="black",face = "bold",hjust=0,
                                         margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 10)),
        axis.text.y  = element_text(vjust=0.5, size=10,colour="black"),
        axis.text.x  = element_blank(),
        axis.title.x.top  = element_text(size=12,colour="black", face = "bold",hjust = 1,
                                         margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0)),
        strip.text.x = element_text(size=12, face="bold",colour="black"),
        strip.text.y = element_text(size=12, face="bold",colour="black"),
        strip.background = element_rect(colour="black"),
        strip.placement = 'outside',
        panel.grid.major = element_line(colour="grey80"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white",color="grey30", 
                                        linewidth=0.5, linetype="solid"),
        legend.position = 'top',
        legend.title.position = 'top',
        legend.title = element_text(hjust=0.5,face='bold',size=10),
        plot.margin = unit(c(0.5, 0.5, 0.5,0.7), "cm"))

confmtx

## Write plot
ggsave(confmtx,
       file = paste0(out_path,'Confidence_Matrix_',tax,'.png'),
       device = 'png', width = 25, height=27, units='cm', dpi=300, scale=1)

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


impplot <-  ggplot(imppl,(aes(x=Var,y=Mean))) +
  geom_bar(stat = 'identity',fill=colpal[3], col=colpal[3],width = 0.8) +
  scale_x_discrete(labels = function(x) sapply(x, function(i) label_map[[i]])) +
  geom_linerange(inherit.aes=FALSE,
                 aes(x=Var, ymin=Mean-Se, ymax=Mean+Se), 
                 colour='#172957', alpha=0.9, linewidth=1.3) +
  ylab(label = 'Mean decrease in Gini coefficient') +
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
       device = 'png', width = 20, height=15, units='cm', dpi=300, scale=1)


#### Partial Response Plots ----

## Prepare data
pplotdata <- plotdata %>%
  mutate(across(y, round, 2))

## Get top six variablers
top6 <- imppl %>%
  slice_max(order_by = Mean, n = 6) %>%
  pull(Var)

## List object for storing plots
cvRP.list <- list()

## Lables for facet strips (used as x-axis titles)
envlabeller <- as_labeller(envlab,default = label_parsed)
top6ylims <- pplotdata %>%
  filter(predvar %in% top6) %>%
  summarise(
    min_y = min(y, na.rm = TRUE),
    max_y = max(y, na.rm = TRUE)
  )



## Loop through variables and plot
for (i in top6) {
  
  pdata <- pplotdata[pplotdata$predvar==i,]
  cvRP.list[[i]] <- ggplot(pdata,aes(x=x,y=y,col=class)) +
    geom_smooth(method='loess',linewidth=0.5,se=FALSE,span = 0.2) +
    facet_wrap(~ predvar,scales = "free_x",ncol =3,
               labeller = envlabeller,
               strip.position = 'bottom') +
    scale_y_continuous(limits = c(top6ylims$min_y,top6ylims$max_y),
                       labels = scales::number_format(accuracy = 0.01),
                       n.breaks = 7,minor_breaks = NULL) +
    scale_x_continuous(n.breaks = 5,minor_breaks = NULL) +
    scale_colour_manual(values = colpal) +
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

length(cvRP.list)


## Layout of top six predictor variable plots
cvRP <-  wrap_plots(cvRP.list) + plot_layout(ncol=3,guides = 'collect') 
cvRP

## Write plot
ggsave(cvRP,
       file = paste0(out_path,'Partial_Response_',tax,'.png'),
       device = 'png', width = 30, height=15, units='cm', dpi=300, scale=1)

#### Save main results ----

save(results,asg.perf,asg.cl.perf,impplot,cvRP,cvRP.list, 
     file = paste0(out_path,'Model_Results_',tax,'.RData'))

write.csv(asg.cl.perf, file = paste0(out_path,'Model_Performance_',tax,'.csv'),
          fileEncoding = "UTF-8")


######### 6. Make model predictions ############################################

#### Drop unnecessary predictor layers  ----
dr <- names(s1)
dr <- dr[!dr %in% names(ffs[[1]]$forest$xlevels)]
s2 <- dropLayer(s1, dr)

#### Set up objects for saving predictions ----
cvpred <- NULL
cvpred.cps <- list()


#### Loop through each model and predict to raster ----

for (i in 1:length(ffs)){
  
  ## Run number
  rnn <-  paste0('Run',i)
  
  ## Class prediction
  if (is.null(cvpred)){
    cvpred  <- stack(raster::predict(s2,ffs[[i]]))
    names(cvpred) <- rnn
  } else {
    tmpl <- predict(s2,ffs[[i]])
    names(tmpl) <- rnn
    cvpred <- addLayer(cvpred,tmpl)
  }
  
  ## Probabilities for each class
  cvpred.cps[[rnn]] <- predict(s2,ffs[[i]],type='prob',index=1:numclass)
  
}

#### Final outputs ----

## Create a raster stack for confidence results
ROutput <- stack()

## Calculate most frequent class and its frequency
# Most frequent class
MaxClass <- modal(cvpred,freq=FALSE)
ROutput <- addLayer(ROutput,MaxClass)

# Frequency of most frequent class (fraction of runs)
MaxClassF <- modal(cvpred,freq=TRUE)/nruns
ROutput <- addLayer(ROutput,MaxClassF)

### Calculate average probabilities for classes
classsums <- Reduce("+", cvpred.cps)
AvePclass <- classsums / nruns

## Find average probability of maximum frequency class
MaxClassAveProb <- stackSelect(AvePclass, MaxClass)
ROutput <- addLayer(ROutput,MaxClassAveProb)

## Calculate new layer for frequency x probability
CombConf <- MaxClassF * MaxClassAveProb
ROutput <- addLayer(ROutput,CombConf)

## Rename layers
names(ROutput) <- c("MaxClass","MaxClassF","MaxClassAveProb","CombConf")

# Factor levels for the most frequent class
ROutput$MaxClass <- as.factor(ROutput$MaxClass)
rat <- levels(ROutput$MaxClass)[[1]]
rat$Cluster <- levels(train$resp)
levels(ROutput$MaxClass) <- rat

## Plot layers 
plot(ROutput)


#### Export Rasters ----

## Most frequent class over model runs
writeRaster(ROutput$MaxClass,
            file = paste0(out_path,'Most_Frequent_Class_',tax,'.tif'),
            format="GTiff",overwrite=T)
## Frequency of most frequent class over model runs
writeRaster(ROutput$MaxClassF,
            file = paste0(out_path,'Most_Frequent_Class_Frequency_',tax,'.tif'),
            format="GTiff",overwrite=T)
## Predicted probability of most frequent class over model runs
writeRaster(ROutput$MaxClassAveProb,
            file = paste0(out_path,'Most_Frequent_Class_Probability_',tax,'.tif'),
            format="GTiff",overwrite=T)
## Model prediction confidence (Frequency x probability of most frequent class)
writeRaster(ROutput$CombConf,
            file = paste0(out_path,'Most_Frequent_Class_Confidence_',tax,'.tif'),
            format="GTiff",overwrite=T)

