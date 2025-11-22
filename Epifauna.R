################################################################################
####            BIG DATA APPROACHES REVEAL LARGE-SCALE PATTERNS IN          ####
####            MARINE EPIFAUNA (ICES JOURNAL OF MARINE SCIENCE)            ####
################################################################################

## This script relates to the work in Cooper, K.M., Curtis, M., Downie, A-L.,
# Bolam, S.G. (2025). Big data approaches reveal large-scale spatial patterns in
# marine epifauna. ICES Journal of Marine Science.

# Data used in the script is sourced from the OneBenthic
# (https://rconnect.cefas.co.uk/onebenthic_portal/) trawl database using a sql 
# query. For users without direct access to this database, data
# can be sourced using the OneBenthic Data Extraction tool (Trawl): 
# https://rconnect.cefas.co.uk/onebenthic_dataextractiontrawl/

# This script includes code for running Random Forest models (a quick look-see),
# but the accompanying files: 
#  ClusterModel_2025.R and
#  ContinuousVariablesModel_2025.R 
# should be used for final modelling.

# The input data file for modelling is produced in this script and is named 
# 'poseidon_trawl_metrics_4_modelling.csv'

#_______________________________________________________________________________
#### GET DATA ####

## Set working directory
setwd("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl")

## Load packages
library(pool)
library(DBI)
library (RPostgres)
library(dplyr)
library(sf)
library(geojsonio)

## Connect to OneBenthic Trawl database
Sys.setenv(R_CONFIG_ACTIVE = "one_benthic")
dw <- config::get()
pool <- dbPool(drv = dbDriver(dw$driver),
               dbname = dw$database,
               host = dw$server,
               port =  dw$port,
               user = dw$uid,
               password = dw$pwd)

## Query OneBenthic trawl database
data = dbGetQuery(pool,"
SELECT su.surveyname,
s.samplecode,
s.samplelatstart,
s.samplelongstart,
s.samplelatend,
s.samplelongend,
s.date,
s.month,
s.year,
s.starttime,
s.endtime,
s.waterdepthstart,
s.waterdepthend,
g.geartype_geartype,
g.gearname,
s.gear_gearcode,
s.vesselspeed,
s.distance,
s.trawlvolume,
s.macrosieve,
ts.worrmstaxa_taxonname,
ts.taxaqual_qualifier,
w.validname,
w.rank,
tq.qualifiername,
w.validaphiaid,
w.family,
w.phylum,
ts.abund,
su.datapubliclyavailable,
s.samplecode2,
o.ownername,
su.programme,
o.source,
o.sector,
geom_line,
geom_line_length,
su.metadata

FROM 
associations.survey as su
INNER JOIN associations.surveysample as ss ON ss.survey_surveyname = su.surveyname 
INNER JOIN samples.sample as s ON ss.sample_samplecode = s.samplecode
INNER JOIN gear.gear as g ON s.gear_gearcode = g.gearcode
INNER JOIN faunal_data.taxasample as ts ON s.samplecode= ts.sample_samplecode 
LEFT JOIN faunal_data.taxaqual as tq ON ts.taxaqual_qualifier = tq.qualifier 
INNER JOIN faunal_data.worrms as w ON w.aphiaid = ts.worrms_aphiaid
INNER JOIN associations.sampleowner as so ON so.sample_samplecode = s.samplecode
INNER JOIN associations.owner as o ON so.owner_ownername = o.ownername
WHERE g.geartype_geartype = 'Trawl' AND
w.validaphiaid NOT IN (3,51,101,105,558,830,852,882,883,913,924,948,1066,1071,1128,1248,1267,1410,2081,2824,10194,11676,11707,11734,14712,21263,105695,105711,105729,105766,105814,105815,105821,105822,105869,105876,105883,105885,105887,105891,105923,106331,106358,106386,110671,110673,110690,110698,111795,111796,111807,111808,119822,119950,120017,120020,120087,120180,120203,120206,120208,125451,125464,125469,125470,125471,125508,125546,125591,125716,125739,125741,125742,125743,125885,126281,126285,126375,126417,126421,126425,126426,126435,126436,126437,126438,126439,126440,126441,126444,126445,126446,126448,126449,126450,126456,126457,126458,126461,126484,126501,126715,126716,126736,126822,126975,127023,127066,127312,127419,127427,129291,129352,129413,129414,129646,129868,129914,130103,130491,130494,130500,130508,130537,130544,130977,130980,131100,131117,131127,131141,131436,131458,135220,135301,135302,135304,135306,135994,138139,140270,140271,140625,143723,143834,144129,144199,145541,145546,145551,145984,146142,146420,148823,148824,151265,152203,152292,152352,153131,154385,154527,154747,182739,271509,272030,272392,300719,367297,368687,382226,399563,400076,400587,400596,1021266,1528231,999999001) AND
su.surveyname NOT IN ('North Sea Benthos_2m Beam Trawl_Cirolana00/05') AND
s.samplecode NOT IN ('BT1_2007', 'Tyne transect_6', 'Tyne transect_5', 'CEND0809_2MB_C47', 'CEND0707_N3/2MB/1','CEND0910_2MB_C1', 'LINCS_OFFSHORE_2014_BT01') AND s.id < 3926
ORDER by su.surveyname, s.samplecode,ts.abund desc;") # Check Plaice record with Matt # samples excluded where potential positional errors

## Inspect data
unique(data$surveyname)# Surveys x110
length(unique(data$samplecode))# Samples x3799
range(data$year)# Years: 1987 to 2023
dim(data)# 95565    38
#View(data)
#_______________________________________________________________________________
#### PREPARE DATA: EXCLUDE POINTS OUTSIDE RASTER ####

## Drop stations outside the extent of the raster predictor variables used for modelling

## Load raster extent polygon (predictor variables)
poly <-  st_read(pool, query ="SELECT * from spatial.raster_extent_polygon;")

## Plot raster extent
plot(poly)

## Change points to a sf dataframe
data_sf_df <- st_as_sf(x = data,                         
               coords = c("samplelongstart", "samplelatstart"),
               remove=F,
               crs = 'WGS 84')

## Plot  all location records
plot(data_sf_df)

## Check polygon and points data have same CRS
st_crs(poly)# EPSG:4326
st_crs(data_sf_df)# none

## Give points_sf_df the same crs as poly
data_sf_df <- data_sf_df %>% st_set_crs(st_crs(poly))

## Keep only the points inside the polygon
sf_use_s2(FALSE)
kept_data <- st_intersection(poly, data_sf_df)
dim(kept_data)# 90157    42
head(kept_data)

## Plot points kept
plot(kept_data)

## Turn back into a df
names(kept_data)
library(dplyr)
kept_data <- kept_data %>%st_drop_geometry()# Drop geom column
kept_data2 <-kept_data[,c(4:40)]# Drop unwanted rows at start (i.e. those associated with poly)
data <- kept_data2
#_______________________________________________________________________________
#### DATA EXPLORATION: MAPS AND ASSOCIATED HISTOGRAMS ####

## Load packages
library(sf)
library(reshape2)
library(scales)
library(ggplot2)
require(cowplot)

## Bring in country polygons for use in maps
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))

## Create copy of col 'geom_line_length' (will be used to create binned values)
data$geom_line_length2 <- data$geom_line_length

## Create dataframe for use in mapping and histograms with one row for each sample
pos <- distinct(data, samplelatstart,samplelongstart,sector,source,programme,year,month,gear_gearcode,macrosieve,vesselspeed,geom_line_length2)#,trawlvolume
#View(pos)

## Change extreme values to NA (believed to be errors). Arbitrary decision for values >5 (there are 24)
#pos$geom_line_length2[pos$geom_line_length2 >5] <- NA

## Confirm values >5km have been removed
#max(pos$geom_line_length2,na.rm=T)

## Create bins for col 'geom_line_length2'
pos$geom_line_length2 <- cut(pos$geom_line_length2,breaks = c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,10,15),
    labels = c(0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,10.0, 15.0))

## Melt data into form for faceting
facdat2=melt(pos,c("samplelatstart","samplelongstart"))
#View(facdat3)

## Insert NAs into empty cells in col 'value'
facdat2$value[facdat2$value==""]<- NA

## Add new col for concatenated variable and value cols - this is necessary to ensure unique
# values in maps, otherwise a 1 for month and sieve size will have the same colour
facdat3=transform(facdat2,varval=paste0(variable,value))

## Make varval a factor
facdat3$varval <- as.factor(facdat3$varval)

## Show levels
#levels(facdat3$varval)

## Return a specified number of colours.
#scales::hue_pal()(12) # Enter required number of colours in brackets

## Use this code to allow facet labels to be modified (changing col 'Sieve' to 'Sieve (mm)')
labels=c(gear_gearcode="Gear",macrosieve="Mesh Size (mm)",year="Year",month="Month",programme = "Programme",source = "Source",sector = "Sector", vesselspeed = "Tow Speed (kt)",geom_line_length2="Tow Length (km)")#,trawlvolume="Sample Volume (litres)"
#View(facdat3)

## Create map plots
sam.fac=ggplot()+
 #geom_polygon(data=defra.dem.df, aes(x=long, y=lat, group=group,fill=Depth), size=0.15)+
 #scale_fill_manual(values=grey)+
 geom_sf(data=countries, fill ="black",col ="black")+ 
 geom_point(data=facdat3,aes(samplelongstart,samplelatstart,col=factor(varval)), size=0.6,
 show.legend = FALSE, alpha=0.4)+
 coord_sf(xlim = c(-12, 9),ylim = c(47, 62))+
 labs(x="Longitude",y="Latitude")+
 theme_bw(base_size=14)+#was 15,16
 facet_wrap(~variable,labeller=labeller(variable=labels),ncol=3)+
 theme(strip.text.x = element_text(size = 14, face = "bold"))+
 scale_colour_manual(values = c( 
"#F8766D",#gear_gearcode2BT
"#00BFC4",#gear_gearcode4BT

"#F8766D",#geom_line_length20.25
"#EB8335",#geom_line_length20.5
"#DA8F00",#geom_line_length20.75 
"#C49A00",#geom_line_length21
"#A9A400",#geom_line_length21.25
"#86AC00",#geom_line_length21.5 
"#53B400",#geom_line_length21.75 
"#FF63B9",#geom_line_length210
"#FF6B96",#geom_line_length215
"#00BA38",#geom_line_length22
"#00BE6D",#geom_line_length22.25
"#00C094",#geom_line_length22.5
"#00C0B5",#geom_line_length22.75
"#00BDD2",#geom_line_length23
"#00B6EB",#geom_line_length23.25 
"#00ABFD",#geom_line_length23.5
"#619CFF",#geom_line_length23.75
"#A58AFF",#geom_line_length24
"#D078FF",#geom_line_length24.25
"#EC69EF",#geom_line_length24.75
"#FB61D7",#geom_line_length25
"#808080",#geom_line_lengthNA
   
#"#F8766D",#geom_line_length21 
#"#A3A500",#geom_line_length22 
#"#00BF7D",#geom_line_length23 
#"#00B0F6",#geom_line_length24
#"#E76BF3",#geom_line_length25 
#"#808080",#geom_line_length2NA

"#00BFC4",#macrosieve18
"#F8766D",#macrosieve4
"#C77CFF",#macrosieve40
"#7CAE00",#macrosieve5
"#F8766D",#month1
"#C77CFF",#month10
"#F564E3",#month11
"#FF64B0",#month12
"#DE8C00",#month2
"#B79F00",#month3
"#7CAE00",#month4
"#00BA38",#month5
"#00C08B",#month6
"#00BFC4",#month7
"#00B4F0",#month8
"#619CFF",#month9
#"#808080",#monthNA

"#F8766D",#programmeCSEMP
"#C49A00",#programmeEIA
"#53B400",#programmeLCM
"#00C094",#programmeICES
"#00B6EB",#programmeMPA
"#A58AFF",#programmeR&D
"#FB61D7",#programmeREC

"#F8766D", #sectorGovernment  
"#00BA38",#sectorIndustry
"#619CFF",#sectorIndustry\n

"#F8766D",#sourceCentre for Environment, Fisheries and Aquaculture Science (CEFAS)
"#B79F00",#sourceDepartment for Environment, Food and Rural Affairs (DEFRA)
"#00BA38",#sourceMarine Aggregates (MA)
"#00BFC4",#sourceNuclear (N)
"#619CFF",#sourceOffshore Wind (OW)
"#F564E3",#sourceOil & Gas (OG)

"#F8766D",#vesselspeed0.5
"#F27C56",#vesselspeed0.7 
"#EC823A",#vesselspeed0.8 
"#E58700",#vesselspeed0.9 
"#DC8D00",#vesselspeed1 
"#D39200",#vesselspeed1.05 
"#C99800",#vesselspeed1.09 
"#BD9C00",#vesselspeed1.1 
"#B1A100",#vesselspeed1.11
"#A3A500",#vesselspeed1.12 
"#93AA00",#vesselspeed1.13 
"#81AD00",#vesselspeed1.15
"#6BB100",#vesselspeed1.2 
"#4EB400",#vesselspeed1.25 
"#15B700",#vesselspeed1.27
"#00BA38",#vesselspeed1.28
"#00BC54",#vesselspeed1.3
"#00BE6A",#vesselspeed1.4 
"#00BF7D",#vesselspeed1.49 
"#00C08F",#vesselspeed1.5 
"#00C19F",#vesselspeed1.53 
"#00C0AF",#vesselspeed1.6 
"#00C0BD",#vesselspeed1.7 
"#00BECB",#vesselspeed1.75 
"#00BCD8",#vesselspeed1.8 
"#00B9E3",#vesselspeed1.9
"#00B5EE",#vesselspeed2 
"#00B0F6",#vesselspeed2.1 
"#00AAFE",#vesselspeed2.2
"#18A3FF",#vesselspeed2.3 
"#619CFF",#vesselspeed2.4 
"#8694FF",#vesselspeed2.5 
"#A28BFF",#vesselspeed2.6
"#B983FF",#vesselspeed2.7
"#CB7AFF",#vesselspeed2.8 
"#DB72FB",#vesselspeed3 
"#E76BF3",#vesselspeed3.3
"#F166E9",#vesselspeed3.4 
"#F862DE",#vesselspeed3.5
"#FD61D1",#vesselspeed3.7 
"#FF61C3",#vesselspeed3.9 
"#FF63B4",#vesselspeed4 
"#FF67A4",#vesselspeed4.3 
"#FF6B93",#vesselspeed4.6
"#FC7181",#vesselspeed5.9
"#808080",#vesselspeedNA

"#F8766D",#year1987
"#EF7F49",#year1988
"#E58700",#year1989
"#D89000",#year1992
"#C99800",#year1993
"#B79F00",#year1994
"#A3A500",#year1995
"#8AAB00",#year1996
"#6BB100",#year1997
"#39B600",#year2000
"#00BA38",#year2001
"#00BD5F",#year2002
"#00BF7D",#year2003
"#00C097",#year2004
"#00C0AF",#year2005
"#00BFC4",#year2006
"#00BCD8",#year2007
"#00B7E9",#year2008
"#00B0F6",#year2009
"#00A7FF",#year2010
"#619CFF",#year2011
"#9590FF",#year2012
"#B983FF",#year2013
"#D376FF",#year2014
"#E76BF3",#year2015
"#F564E3",#year2016
"#FD61D1",#year2017
"#FF62BC",#year2018
"#FF67A4",#year2019
"#FE6E8A"#year2023
))

## Plot histograms to accompany above maps (samples by factor)
p1= ggplot(pos, aes(x=reorder(sector, -table(sector)[sector]),fill=sector))+
 geom_bar(fill=c("#F8766D", "#00BA38","#619CFF"))+
 theme_classic(base_size=12)+
 ggtitle("Sector")+
 guides(fill=FALSE)+
 scale_x_discrete(labels=c("Gov.","Ind."))+
 theme(plot.title = element_text(size=14))+
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.text.x = element_text(hjust=0.5,vjust=0.5, size=12))+
 theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))# centres the plot title

p2= ggplot(pos, aes(x=reorder(source, -table(source)[source]),fill=source))+
 geom_bar()+
 theme_classic(base_size=12)+
 ggtitle("Source")+
 guides(fill=FALSE)+
 scale_x_discrete(labels=c("DEFRA","OW","MA","CEFAS","N","OG"))+
 xlab(label="Source")+
 theme(plot.title = element_text(size=14))+
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=12))+
  theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))

p3= ggplot(pos, aes(x=reorder(programme, -table(programme)[programme]),fill=programme))+
 geom_bar()+
 theme_classic(base_size=12)+
 scale_fill_hue(na.value = "lightsteelblue1")+
 ggtitle("Programme")+
 guides(fill=FALSE)+
 xlab(label="Programme")+
 theme(plot.title = element_text(size=14))+
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=12))+
  theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))

p4= ggplot(pos, aes(year, fill=as.factor(year)))+
 geom_bar()+
 theme_classic(base_size=12)+
 ggtitle("Year")+
 guides(fill=FALSE)+
 xlab(label="year")+
 theme(plot.title = element_text(size=14))+
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=12))+
  theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))

p5= ggplot(pos, aes(month, fill=as.factor(month)))+
 geom_bar()+
 theme_classic(base_size=12)+
 ggtitle("Month")+
 guides(fill=FALSE)+
 xlab(label="month")+
 theme(plot.title = element_text(size=14))+
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5, size=12))+
  theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))+ 
scale_x_continuous(breaks=breaks_pretty())

p6= ggplot(pos, aes(x=reorder(gear_gearcode, -table(gear_gearcode)[gear_gearcode]),fill=gear_gearcode))+
 geom_bar()+
 ggtitle("Gear")+
 guides(fill=FALSE)+
 scale_x_discrete(labels=c("2BT","4BT"))+
 theme_classic(base_size=12)+#was 15,16
 theme(plot.title = element_text(size=14))+
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.text.x = element_text(hjust=0.5,vjust=0.5, size=12))+
 theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))# centres the plot title

p7= ggplot(pos, aes(x=reorder(macrosieve, -table(macrosieve)[macrosieve]),fill=macrosieve))+
 geom_bar(fill = c("#F8766D" ,"#7CAE00", "#00BFC4", "#C77CFF"))+
 ggtitle("Mesh Size (mm)")+
 guides(fill=FALSE)+
 scale_x_discrete(labels=c("4","5","18","40"))+
 theme_classic(base_size=12)+#was 15,16
 theme(plot.title = element_text(size=14))+
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.text.x = element_text(hjust=0.5,vjust=0.5, size=12))+
 theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))# centres the plot title

p8= ggplot(pos, aes(vesselspeed,fill=as.factor(vesselspeed)))+
 geom_histogram(na.rm = FALSE)+
 theme_classic(base_size=12)+
 ggtitle("Tow Speed (kt)")+
 guides(fill=FALSE)+
 #xlab(label="Speed")+
 theme_classic(base_size=12)+#was 15,16
 theme(plot.title = element_text(size=14))+
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.text.x = element_text(hjust=0.5,vjust=0.5, size=12))+
 theme(plot.title = element_text(hjust = 0.5))+
 # scale_fill_manual(p8col,na.value="red")+
 theme(axis.text.y = element_text(hjust=1, size=12))# centres the plot title

p9= ggplot(pos, aes(geom_line_length2, fill=as.numeric(geom_line_length2)))+
 geom_bar(fill=c(
"#F8766D",#geom_line_length20.25
"#EB8335",#geom_line_length20.5
"#DA8F00",#geom_line_length20.75 
"#C49A00",#geom_line_length21
"#A9A400",#geom_line_length21.25
"#86AC00",#geom_line_length21.5 
"#53B400",#geom_line_length21.75 

"#00BA38",#geom_line_length22
"#00BE6D",#geom_line_length22.25
"#00C094",#geom_line_length22.5
"#00C0B5",#geom_line_length22.75
"#00BDD2",#geom_line_length23
"#00B6EB",#geom_line_length23.25 
"#00ABFD",#geom_line_length23.5
"#619CFF",#geom_line_length23.75
"#A58AFF",#geom_line_length24
"#D078FF",#geom_line_length24.25
"#EC69EF",#geom_line_length24.75
"#FB61D7",#geom_line_length25
"#FF63B9",#geom_line_length210
"#FF6B96",#geom_line_length215
"#808080"#geom_line_lengthNA  
   
 ))+
 theme_classic(base_size=12)+
 ggtitle("Tow Length (km)")+
 guides(fill=FALSE)+
 #xlab(label="year")+
 theme(plot.title = element_text(size=14))+
 theme(axis.title.x = element_blank())+ # remove x-axis title
 theme(axis.title.y = element_blank())+ # remove y-axis title
 theme(axis.text.x = element_text(angle=90,  size=12))+#angle=90,hjust=1,vjust=0.5,
  theme(plot.title = element_text(hjust = 0.5))+
 theme(axis.text.y = element_text(hjust=1, size=12))+
  scale_x_discrete(breaks = seq(0, 5, by = 0.5))

## Combine individual plots into one plot
histo=plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,align = "hv",nrow=3)

## Add labels
plott <- plot_grid(sam.fac,histo, labels = c("a)","b)"),nrow = 1,rel_widths = c(1.2, 1),scale = 0.98,
label_size = 20,hjust = c(-2,1), vjust = c(1,1),align = "hv") #SCALE 0.95, RW 1:1.

## Save plot
ggsave(plot = plott,
       filename = paste0("C:\\Users\\kmc00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_S1.tif"),
       width = 50,height = 30,units = "cm",pointsize = 36,device = "tiff",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### PREPARE DATA: ESTIMATE MISSING TOW DISTANCES BASED ON TOW TIME ####

## Load package
library(tidyverse)
library(ggplot2)
library(plotly)

## Col names
names(data)
str(data)

## Subset for 2m BT and where sample processed using 4/5mm sieve (samples likely to have been deployed in similar way)

data <- data %>% filter(data$gear_gearcode=='2BT', data$macrosieve == 4 | macrosieve ==5)

## Change geom_line_length values of 0 to NA
data$geom_line_length[data$geom_line_length == 0] <- NA

## Create copy of col 'geom_line_length' (will be used to create binned values)
#data$geom_line_length2 <- data$geom_line_length

## Create date/time field for trawl start and end. Necessary to calculate total trawl time
data$datetime_start <- as.POSIXct(paste(data$date, data$starttime), format="%Y-%m-%d %H:%M:%S")
data$datetime_end <- as.POSIXct(paste(data$date, data$endtime), format="%Y-%m-%d %H:%M:%S")

## Calculate tow time
data$tow_time <-difftime(data$datetime_end,data$datetime_start);
data$tow_time2 <- abs(as.numeric(data$tow_time))# convert to numeric
#View(data)

## Number of samples with tow time but not distance
#dim(distinct(subset(data[,c("samplecode","geom_line_length","tow_time2")],!is.na(tow_time2) & is.na(geom_line_length))))#75  3

## These samples have a trawl time of 0 - these are errors which need correcting
time_errors <- data[which(data$tow_time2==0),]
write.csv(time_errors, "C:\\Users\\kmc00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\trawls_zero_time_error.csv", row.names=FALSE)
#View(time_errors)

## Create a subset of the data for regressing time and distance
distime <- distinct(data, samplecode, date, tow_time,tow_time2,geom_line_length)
head(distime)

## Drop rows where time is 0 (this step shouldn't be necessary once errors corrected - see above)
distime2<-subset(distime, tow_time2!=0)

## Drop rows with missing data
distime3 <- distime2[complete.cases(distime2), ]
#View(distime3)

## Find samples with negative trawl time
distime3[distime3$tow_time < 0, ]

## Update trawl times for above samples collected over midnight. Note the date/time field doesn't deal with these samples as the nominal date is the same.
distime3$tow_time[distime3$samplecode == 'DGRB413_A2_BT' ] <- 6
distime3$tow_time2[distime3$samplecode == 'DGRB413_A2_BT' ] <- 6
distime3$tow_time[distime3$samplecode == '0205B_3' ] <- 15
distime3$tow_time2[distime3$samplecode == '0205B_3' ] <- 15
distime3$tow_time[distime3$samplecode == 'STN216_T1_35_CEND9/09' ] <- 16
distime3$tow_time2[distime3$samplecode == 'STN216_T1_35_CEND9/09' ] <- 16
distime3$tow_time[distime3$samplecode == 'STN366_T2_13_CEND9/09' ] <- 17
distime3$tow_time2[distime3$samplecode == 'STN366_T2_13_CEND9/09' ] <- 17
distime3$tow_time[distime3$samplecode == 'STN59_T1_51_CEND9/09' ] <- 20
distime3$tow_time2[distime3$samplecode == 'STN59_T1_51_CEND9/09' ] <- 20

## Remove columns where time value is negative (not necessary if above corrections made)
distime4 <- distime3[distime3$tow_time2 >= 0, ]
#View(distime4)

## Check initial model fit

plot <- ggplot(data = distime4, aes(x = tow_time2, y = geom_line_length)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +# Regression line
  ggtitle("Scatter Plot with Regression Line") +
  xlab("Tow time (mins)") +
  ylab("Tow distance (km)")
ggplotly(plot)

## Add an 'outlier' to allow for subsetting of good' data (decide on appropriate cut-off values after inspecting initial model fit)
distime4$outlier_len<- with(distime4, ifelse(geom_line_length > 1.75, TRUE, FALSE))# Remove length outliers
distime4$outlier_time<- with(distime4, ifelse(tow_time2 > 30, TRUE, FALSE))# Remove time outliers
#View(distime4)

## Remove outliers
distime5 <- distime4[which(distime4$outlier_len==FALSE & distime4$outlier_time==FALSE),]
#View(distime5)

## Plot tow time versus tow distance
plot <- ggplot(data = distime5, aes(x = tow_time2, y = geom_line_length)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +# Regression line
  ggtitle("Scatter Plot with Regression Line") +
  xlab("Tow time (mins)") +
  ylab("Tow distance (km)")
ggplotly(plot)

## Model 
model <- lm(geom_line_length ~ tow_time2, data=distime5)

## Use predict function to calculate 'tow distance' from 'tow time'
pred <- predict(model,data.frame(tow_time2 = c(600,800)))#0.6835246km 0.9375890km

## Find records with tow time but no geom_tow_length
data_time_only  <- subset(data, !is.na(tow_time2) & is.na(geom_line_length))
#View(data_time_only)

## Create a dataframe with the required cols
data_time_only2 <- data_time_only[,c(2,36,42)]# subset for required cols
data_time_only2 <- distinct(data_time_only, samplecode, tow_time2)# one row for each sample

## Calc predicted distance
data_time_only2$pred_dist <- predict(model,data_time_only2)
data_time_only2 <- data_time_only2[,c(1,3)]# keep only cols 'samplecode' and 'predicted 'geom-line_length'
colnames(data_time_only2) <- c("samplecode","geom_line_length")
#View(data_time_only2)

## Update value in col 'geom_line_length' (df 'data') where samplecode present in above df (data_time_only2)
# see https://stackoverflow.com/questions/63463265/update-few-values-in-a-data-frame-column-with-another-data-frame-column-in-r
#View(data[which(data$samplecode=='IDFBT09_64_B1'),])# sample with only tow time
data2 <- data %>% rows_update(data_time_only2, by = "samplecode")
#View(data[which(data$samplecode=='IDFBT09_64_B1'),]) check same sample now has a value for 'geom-line_length'
names(data2)
#_______________________________________________________________________________
#### PREPARE DATA: SAMPLE LOCATIONS (FIGURE 1) ####

## Load packages
library(sf)
library(rasterVis)# in order to use raster in ggplot
library(raster)
library(ggnewscale)
library(scales)
library(ggplot2)
library(dplyr)
library(shadowtext)

## Load countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))
ni_border <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\ni_border.shp"))

## Load DEM
dem <- raster("C:\\Users\\kmc00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\GEBCO_15_Nov_2023_be52dc9c9b2d\\gebco_2023_n61.0_s48.0_w-11.0_e11.0.tif")

## Remove elevations above sea level
dem[dem>0] <- 0

## Create slope and hillshade
slope = terrain(dem, opt='slope')
aspect = terrain(dem, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)

dem_spdf <- as(dem, "SpatialPixelsDataFrame")
dem_spdf <- as.data.frame(dem_spdf)
colnames(dem_spdf) <- c("value", "x", "y")

hill_spdf <- as(hill, "SpatialPixelsDataFrame")
hill_spdf <- as.data.frame(hill_spdf)
colnames(hill_spdf) <- c("value", "x", "y")

## Get unique position coordinates
unique_pos <- data2 %>% dplyr::select(samplelongstart, samplelatstart)%>%unique

## Convert latitude and longitude into geometries using st_as_sf(). 
points <- unique_pos %>%
  st_as_sf(coords = c("samplelongstart", "samplelatstart"), crs = 4326)

## st_coordinates() extracts the lon/lat values as a data frame with X and Y columns so you can use geom_point (ability to resize points etc)
points2 <- st_coordinates(points)
points2 <- as.data.frame(st_coordinates(points))
head(unique_pos)

## Create plot (map of sample locations, place names and background bathymetry)
# For bathy colours and breakpoints: https://stackoverflow.com/questions/70739780/is-there-a-scale-function-with-which-i-can-use-4-breaks-points
PSam2=ggplot()+ 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7)+
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0)))+
  geom_point(data = points2,aes(x = X, y = Y), fill="yellow",alpha = 0.6,colour="yellow",size=0.75)+#size = 0.1,color='blue'
  geom_sf(data=countries, fill ="black",col ="black")+ 
  geom_sf(data=ni_border, ,col ="grey",linewidth = 0.2)+ 
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw()+
  xlab("Longitude") +
  ylab("Latitude")+
  annotate("text",x=c(-1.4),y=c(52.5),label=c("UNITED \nKINGDOM"),color="white", size=3)+
  annotate("text",x=c(-6.7),y=c(51.2),label=c("Celtic \nSea"),color="white", size=5)+
  annotate("text",x=c(-1.1),y=c(50.19),label=c("English Channel"),color="white", size=5)+
  annotate("text",x=c(-4.7),y=c(53.8),label=c("Irish Sea"),color="white", size=4)+
  annotate("text",x=c(3),y=c(56.5),label=c("North Sea"),color="white", size=5)+
  annotate("text",x=c(6.3),y=c(54.5),label=c("German Bight"),color="white", size=3)+
  annotate("text",x=c(0),y=c(58.5),label=c("Fladden Ground"),color="white", size=3)+
  annotate("text",x=c(2.33),y=c(54.9),label=c("Dogger Bank"),color="white", size=3)+
  annotate("text",x=c(5.6),y=c(58.3),label=c("Norwegian Trench"),color="white", size=3, angle = -40)+
  annotate("text",x=c(-1),y=c(56.5),label=c("Scalp \nBank"),color="white", size=3)+
  annotate("text",x=c(-4.8),y=c(51.4),label=c("Bristol Channel"),color="white", size=3)+
  annotate("text",x=c(3),y=c(58),label=c("Ling Bank"),color="white", size=3)+
  annotate("text",x=c(0.95),y=c(53.5),label=c("Inner \nSilver \nPit"),color="white", size=3)+
  annotate("text",x=c(5.1),y=c(57),label=c("Fisher \nBanks"),color="white", size=3)+
  annotate("text",x=c(1.5),y=c(51),label=c("Strait of Dover"),color="white", size=3, angle = 45)+
  annotate("text",x=c(-5.8),y=c(52),label=c("St George's \nChannel"),color="white", size=2.5)+
  annotate("text",x=c(-7.6),y=c(53.2),label=c("IRELAND"),color="white", size=3)+
  annotate("text",x=c(1),y=c(49),label=c("FRANCE"),color="white", size=3)+
  annotate("text",x=c(6.5),y=c(52.5),label=c("NETHERLANDS"),color="white", size=3)+
  annotate("text",x=c(9.05),y=c(56),label=c("DENMARK"),color="white", size=3)+
  annotate("text",x=c(3.9),y=c(51),label=c("BELGIUM"),color="white", size=3)+
  annotate("text",x=c(7.6),y=c(59),label=c("NORWAY"),color="white", size=3)+
  annotate("text",x=c(9),y=c(53),label=c("GERMANY"),color="white", size=3)+
  annotate("text",x=c(-4.5),y=c(54.3),label=c("Isle \nof \nMan"),color="white", size=2)+
  annotate("text",x=c(-6.7),y=c(56.9),label=c("Hebrides"),color="white", size=3)+
  annotate("text",x=c(0.32),y=c(56.5),label=c("Devil's \nHole"),color="white", size=3)+
  annotate("text",x=c(-9),y=c(50),label=c("SW \nApproaches"),color="white", size=4)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 12),legend.title = element_text(color = "white", size = 12))+#make legend background transparent and text white
  theme(legend.position=c(0.9,0.18))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.key.size = unit(1.0, "cm"))+
  labs(fill = "Bathymetry (m)")

## Save Figure 1
ggsave(plot = PSam2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_1.tif"),
       width = 25,height = 26,units = "cm", pointsize = 60,
       device = "tiff",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### PREPARE DATA: COLLAPSE TO FAMILY LEVEL ####

## Load package
library(tidyr)
library(dplyr)

## Subset data for required columns
names(data)
data2test <- data2[,c(2,3,4,16,20,22,23,26,27,28,29,17,37)]
names(data2test)
dim(data2test)# 84493    13

## Collapse by validaphiaid
data2test2 <- data2test %>%
  group_by(samplecode, validaphiaid)%>%
  mutate(abund = sum(abund))%>%
  distinct()
head(data2test2)

## Remove rows with no validname
data2test2 <- data2test2[!(is.na(data2test2$validname) | data2test2$validname==""), ]
str(data2test2)

## Note missing values in col 'family'. These missing values indicate taxon name is at a higher level. Infill with appropriate name from col 'validname'
data2test2$family <- ifelse(data2test2$family == "", data2test2$validname, data2test2$family)
data2test2$family <- ifelse( is.na(data2test2$family), data2test2$validname, data2test2$family)
data2test2$family <- ifelse( data2test2$family=='[NULL]', data2test2$validname, data2test2$family)
head(data2test2)

## Collapse data by family name
data2test3 <- data2test2 %>% 
 group_by(samplecode, samplelatstart,samplelongstart,gear_gearcode,macrosieve, vesselspeed,geom_line_length,family) %>% #validname,validaphiaid,taxaqual_qualifier,
  summarise(abund = sum(abund))
head(data2test3)

## Change data2test3 from long to wide format: data, key (headers), values
#https://www.bing.com/videos/search?q=tidyr+long+to+wide&&view=detail&mid=9397360FBC77F8AE08099397360FBC77F8AE0809&&FORM=VRDGAR
data_wide <-spread(data2test3,family,abund)
head(data_wide)
names(data_wide)

## change NA to zero in faunal columns
data_wide[,8:ncol(data_wide)][is.na(data_wide[,8:ncol(data_wide)])] = 0
dim(data_wide)# 2800  464
class(data_wide)
#_______________________________________________________________________________
#### PREPARE DATA: CREATE A SUBSET OF COMPARABLE DATA ####

## Load package
library(dplyr)

## Convert from tibble to df
data_wide <- as.data.frame(data_wide)
head(data_wide)
dim(data_wide)# 2800  464

## Select 2m beam trawls with 4/5mm mesh where tow speed <2 knots and tow length less than 2km
data2 =data_wide %>%
  dplyr::filter(gear_gearcode=="2BT")%>%
  dplyr::filter(macrosieve %in% c(4,5))%>%
  dplyr::filter(is.na(vesselspeed) | vesselspeed < 2)%>%
  dplyr::filter(is.na(geom_line_length) | geom_line_length < 2)

## Check col names and dimensions of df 'data2'
names(data2)
dim(data2)# 2383  464

## Remove samples from the faunal data (df data3) where no fauna present
data3 <-  data2[ rowSums(data2[,8:ncol(data2)])!=0, ]
head(data3)
dim(data3)# 2383  464
class(data3)
#View(data3)

## Dataframe of sample positions after subsets (except spatial autocorrelation step)
data3_before_spat_AC_step <- data3
#_______________________________________________________________________________
#### PREPARE DATA: ASSESS SPATIAL AUTOCORRELATION ####

## Load packages
library(emon)
library(vegan)

## Prepare data for SAC check
sac_data <- data3[,c(1:3,9:ncol(data3))]# take cols samplecode, corrds, and faunal data
names(sac_data)

## Add col for richness
sac_data$richness <- specnumber(sac_data[,4:ncol(sac_data)])
names(sac_data)

## Drop faunal data cols
sac_data2 <- sac_data[,c(1:3,ncol(sac_data))]
head(sac_data2)

## Rename columns
colnames(sac_data2)[3] <- "x"#rename long
colnames(sac_data2)[2] <- "y"#rename lat
dim(sac_data2)# 2383    4 # check num of rows

## Get all longitudes onto a 0 to ? scale by adding max westerly degrees to all x values
min(sac_data2$x)#-11.1445
sac_data2$x <- sac_data2$x+11.1445
#data3 <- data2
head(sac_data2)

## Create seivariogram
x = sac_data2$x
y = sac_data2$y
rich = sac_data2$richness
u = seq(0, 5, 0.2)

semiv <- svariog(x=x, y=y, z=rich, u=u)

## Plot
par(mfrow=c(1,1))
plot(semiv$mid, semiv$cla, xlab='Distance', ylab='Classical')

gauss = function(C0, C1, a, sep) {
  sv = C0 + C1*(1 - exp(-a*sep^2))
  sv
}


f = function(par, z, h) {
  fitted = par[1] + par[2]*(1 - exp(-par[3]*h^2))
  error = sum((fitted - z)^2)
  error
}

## these are initial estimates for the optimising routine
C0 = 40
C1 = 100
a = 0.5
mid = semiv$mid
z = semiv$classical
inits = c(C0, C1, a)

gauss.opt = optim(par=inits, fn=f, z=z, h=mid)
gauss.par = gauss.opt$par                           # LS estimates for model
sep = seq(0,5,0.1)
gauss.fit = gauss(gauss.par[1], gauss.par[2], gauss.par[3], sep=sep)

tiff("C:\\Users\\kmc00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_S2.tiff", height=10, width=12, res=600, pointsize=9, units="cm")
plot(semiv$mid, semiv$cla, xlab='Distance (km)', ylab='Semi-Variogram')
lines(sep, gauss.fit, col=1, lwd=2)
legend("bottomright", legend=c("Gaussian"), col=c(2), lty=1, lwd=2)
dev.off()

## Parameter estimates
gauss.par
#_______________________________________________________________________________
#### PREPARE DATA: REMOVE SPATIALLY AUTOCORRELATED SAMPLES ####
#https://search.r-project.org/CRAN/refmans/sp/html/zerodist.html#:~:text=zerodist%20and%20zerodist2%20return%20a,to%20row%20pairs%20in%20obj%20.

## Load package
library(sp)

## Number of samples before removing those that are spatially autocorrelated
dim(data3)# 2383  464

## Set coordinates
coordinates(data3) <- c("samplelongstart", "samplelatstart")

## Set CRS as lat long (now distance can be defined in km using remove.duplicates)
proj4string(data3) <- CRS("+init=epsg:4326")
st_crs(data3)#none

## Drop 'replicates' zero = distance in km
data3norep  <- remove.duplicates(data3, zero = 2, remove.second = TRUE)
dim(data3norep)# 2km: 1400  462

## Number of samples dropped
length(data3)-length(data3norep) # 983

## Change class to df
data3norep2=data.frame(data3norep)
class(data3norep2)
names(data3norep2)

## Drop col 'optional'
data3norep2=data3norep2[,1:(ncol(data3norep2)-1)]
names(data3norep2)
dim(data3norep2)# 1400  464

## Find out what samples have been excluded
exc_samples <- as.data.frame(setdiff(data3$samplecode, data3norep2$samplecode))

## Update column name
colnames(exc_samples) <- 'samplecode'
exc_samples

## Check dimensions
dim(exc_samples)# 983   1

## Convert df back into 
data3 <- as.data.frame(data3norep2) 
#_______________________________________________________________________________
#### PREPARE DATA: PLOTS SHOWING RETAINED AND EXCLUDED SAMPLES ####

## Load packages
library(grid)
library(gridExtra)

## Get all sample locations
all_samples <-data3_before_spat_AC_step
dim(all_samples )# 2383  464

## Get excluded sample locations
exc_samples_local <- all_samples %>%  filter(samplecode %in% exc_samples$samplecode)
exc_samples_local
dim(exc_samples_local)# 983 464

#dim(points2)#2290    2

## Map showing all samples. Red points are those identified as autocorrelated
plot1 <- ggplot()+ 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7)+
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0)))+
    geom_point(data = data3_before_spat_AC_step,aes(x = samplelongstart, y = samplelatstart), fill="yellow",alpha = 0.4,colour="yellow",size=0.3)+#size = 0.1,color='blue' 
    geom_point(data = exc_samples_local,aes(x = samplelongstart, y = samplelatstart), fill="red",alpha = 0.4,colour="red",size=0.75)+#size = 0.1,color='blue'
  geom_sf(data=countries, fill ="black",col ="black")+ 
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  theme_bw()+
  ylab("Latitude")+
  theme(legend.position = "none",axis.title.x=element_blank())

## Samples remaining
dim(data3)# 1400  464

## Map following removal of autocorrelated samples
plot2 <- ggplot()+ 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7)+
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0)))+
 geom_point(data = data3,aes(x = samplelongstart, y = samplelatstart), fill="yellow",alpha = 0.4,colour="yellow",size=0.3)+#size = 0.1,color='blue' 
 geom_sf(data=countries, fill ="black",col ="black")+ 
 coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
 theme_bw()+
 theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  theme(legend.position = "none")

## Produce combined plot
excplot <- plot_grid(plot1, plot2, labels = c("a)","b)"),nrow = 1,label_size = 12,rel_widths = c(1.0,0.92))

## Add label to plot
x.grob <- textGrob("Longitude", 
                   gp=gpar( col="black", fontsize=12))# fontface="bold",

excplot2 <- grid.arrange(arrangeGrob(excplot, bottom = x.grob))

##Save plot
ggsave(plot = excplot2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_S3.tif"),
       height = 165, width =300, units = "mm", dpi = 500,
       device = "tiff",limitsize = FALSE,bg="white")#width =285
#_______________________________________________________________________________
#### PREPARE DATA: CREATE FAUNAL MATRIX ####

## Load package
library(dplyr)

## Faunal subset (ie remove Sample,Latitude_WGS84, Longitude_WGS84, month and year)
names(data3)
#View(data3)
data4=data3[,9:ncol(data3)]

## Remove taxa from the faunal data (df data4) where no abundance found
data5 <- data4 %>%
 dplyr::select(where(~ any(. != 0)))

## Check class of df data5
class(data5)# df
head(data5)

## Check dimensions of df 'data5'
dim(data5) #1400  429

## Check df 'data5' is just the faunal data
names(data5)# it is

## Change class of df data5 to a matrix
data6=data.matrix(data5)
str(data6)
#View(data6)
#_______________________________________________________________________________
#### PREPARE DATA: REMOVE EXTREME ABUNDANCES ####

## Turn matrix into a vector
vec <- c(data6)

## Order vector
vec_sort <- sort(vec, decreasing=T)
vec_sort

## Decide what max value should be and change higher values to this figure
data6[data6 >1000] <- 1000
head(data6)
dim(data6)# 1400  429
#_______________________________________________________________________________
#### PREPARE DATA: CREATE A DF OF POSITIONS FOR SAMPLES TAKEN FORWARD FOR FAUNAL ANALYSIS ####

## Create a df 'pos' for Sample, Latitude_WGS84 and Longitude_WGS84
pos2=data3[,c(1:3)]
names(pos2)
dim(pos2)# 1400    3

## Plot positions
plot(pos2$samplelongstart,pos2$samplelatstart)
#_______________________________________________________________________________
#### UNIVARIATE METRICS: CALCULATE BY SAMPLE ####

## Calculate univariate summary measures based on faunal abundance data in df 'data6'

## Call package
library(vegan)
library(reshape2)

## Number of species (S)
S = specnumber(data6) # Species Richness(S)
range(S)#1 63

## Number of individuals (N) (inc colonial taxa)
N =rowSums(data6) # Abundance
range(N)#1 8263

## Transformed number of individuals (N) (inc colonial taxa)
N_sqrt <- sqrt(sqrt(rowSums(data6)))

## Combine coordinates and metrics into one dataframe
univmeascoord=as.data.frame(cbind(pos2,S,N,N_sqrt))
class(univmeascoord)

## Update column names
colnames(univmeascoord) <- c('samplecode','samplelatstart','samplelongstart','S','N','N_sqrt')
str(univmeascoord)
head(univmeascoord)
dim(univmeascoord)

## Change Richness to numeric value
univmeascoord$S=as.numeric(univmeascoord$S)

## Drop the Sample column
#univmeascoord$samplecode <- NULL

## Change col names 
colnames(univmeascoord)[2] <- "Latitude_WGS84"
colnames(univmeascoord)[3] <- "Longitude_WGS84"

## Reshape data into long format for faceting
univmeascoord2=melt(univmeascoord,c("samplecode","Latitude_WGS84","Longitude_WGS84"))
head(univmeascoord2)
#_______________________________________________________________________________
#### UNIVARITE METRICS: PLOT SAMPLES WHERE S>1000 ####

# Load packages
library(ggplot2)
library(scales)

## Subset the dataframe for rows where value in column 'abund' is greater than 1000
subset_data <- subset(data, abund > 1000, select = c('samplelatstart', 'samplelongstart'))

## Display the subsetted dataframe
print(subset_data)

# Assuming 'hill_spdf' and 'dem_spdf' are already loaded in your R environment
# Create ggplot
ggplot() + 
  geom_raster(data = hill_spdf, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "white") +
  geom_raster(data = dem_spdf, aes(x = x, y = y, fill = value), alpha=0.7) +
  scale_fill_gradientn(
    colours = c("#051650","#02367b","#006ca5","#0496c7", "#04bade","#55e2e9"),
    limits  = c(-3800,0),
    values  = scales::rescale(c(-3800,-1000,-100,-50,-25,0), from = c(-3800,0))
  ) +
  geom_point(data = subset_data, aes(x = samplelongstart, y = samplelatstart), color = 'yellow') +
  labs(title = 'Locations where Variable S is 1000 or higher', x = 'Longitude', y = 'Latitude') +
  theme_minimal()
#_______________________________________________________________________________
#### SPATIAL MODELLING OF UNIVARIATE METRICS: PRE-TREATMENT OF DATA (IDENTIFY AND REMOVE OUTLIERS)  ####

## Load package
library(dplyr)

## Calculate medians for each biodiversity metric
univmeascoord2$variable <- as.character(univmeascoord2$variable)
univmeascoord2 %>%
  na.omit() %>%
group_by(variable)%>% 
summarise(median=median(value))

## Add medians
univmeascoord2 <- univmeascoord2 %>%
    mutate(median = case_when(
      variable == 'N' ~ 213,
      variable == 'N_sqrt' ~ 3.82,
      variable == 'S' ~ 21))

## Find abs difference between each value and the median
univmeascoord2$abs_diff_med <-abs( univmeascoord2$median -univmeascoord2$value)

## Find the Median Absolute Deviation (MAD)
univmeascoord2 %>%
  na.omit() %>%
group_by(variable)%>% 
summarise(mad=median(abs_diff_med))

## Add MAD column and insert values
univmeascoord2 <- univmeascoord2 %>%
    mutate(mad = case_when(
      variable == 'N' ~ 129,
      variable == 'N_sqrt' ~ 0.74,
      variable == 'S' ~ 7))
head(univmeascoord2)

## Find Modified Z-Score for each data value
univmeascoord2$mod_z_score = 0.6745*(univmeascoord2$abs_diff_med) / univmeascoord2$mad
head(univmeascoord2)
#View(univmeascoord2)

## Get samplecodes where modified z-score is <-3.5 or >3.5
pot_outliers <- univmeascoord2 %>%
  dplyr::filter(mod_z_score < -3.5 | mod_z_score > 3.5) %>%
  dplyr::filter(variable == 'N_sqrt') %>%
  #group_by(variable)%>%
  dplyr::select(samplecode, variable, mod_z_score)
head(pot_outliers)

## Save df 'pot_outliers'
#write.csv(pot_outliers, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\potential_outliers.csv", row.names=FALSE)

## Remove values with a modified Z-score less than -3.5 or greater than 3.5
univmeascoord2 = univmeascoord2 %>%
  na.omit() %>%
  group_by(variable) %>%
  filter(mod_z_score > -3.5 & mod_z_score < 3.5)# %>%
head(univmeascoord2)

## Remove unwanted cols
univmeascoord2 <-univmeascoord2[,1:5]

## Create df of biodiv metrics (long format) for modelling (outliers excluded)
univmeascoord2_mod <-as.data.frame(univmeascoord2)
head(univmeascoord2_mod)

## Update col names
colnames(univmeascoord2_mod)[4] <-'metric' 
colnames(univmeascoord2_mod)[5] <-'measurement'
#_______________________________________________________________________________
#### SPATIAL MODELLING OF UNIVARIATE METRICS: LOAD ENVIRONMENTAL PREDICTOR RASTERS ####

## Load package
library(raster)

## Load rasters
bathy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/bathy3.tif")
cur <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/cur3.tif")
gravel <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Gravel_Fraction.tif")
light <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/light3.tif")
mud <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Predicted_Mud_Fraction.tif")
oxy <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/oxy3.tif")
phyto <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/phyto3.tif")
sal <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sal3.tif")
sil <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/sil3.tif")
spm <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/SPM_MEAN.tif")
temp <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Final/temp3.tif")
wov <- raster("C:/Users/KMC00/OneDrive - CEFAS/R_PROJECTS/TraitsMapping/DATA/Rasters/Mitchell/Wave_veloc.tif")

## Update crs
crs(bathy) <- "+proj=longlat +datum=WGS84 +no_defs" 

## Create raster stack
predictors <- stack(bathy,cur,gravel,light,mud,oxy,phyto,sal,sil,spm,temp,wov)
#predictors

# Simple names for predictor variables
names(predictors)=c("Bathymetry","Current","Gravel","Light","Mud","Oxygen","Phyto","Salinity","Silicate","SPM","Temperature","WOV")#"Sand",
#names(predictors)

## Plot raster stack
#plot(predictors)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF UNIVARIATE METRICS: SELECT VARIABLE FOR MODELLING ####

## Select data by metric (do metrics one at a time)
metric <- univmeascoord2_mod %>% group_by(metric) %>% filter(metric=='S', na.rm = TRUE)
metric <- univmeascoord2_mod %>% group_by(metric) %>% filter(metric=='N_sqrt', na.rm = TRUE)
head(metric)

## Unload extract function from tidyr package (otherwise it won't work)
.rs.unloadPackage("tidyr")

## Extract predictor variables from raster stack
sdata <- raster::extract(predictors, metric[,3:2])

## Change from matrix to df
class(sdata)
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(metric$measurement,sdata)#metric$Sample

## Update name of col 1
colnames(sdata2)[1] <- "value"
head(sdata2)

## Remove rows with NA
sdata2 <- na.omit(sdata2)

## First check cols of correct type
str(sdata2)
dim(sdata2)#  1396   13
#____________________________________________________________________________________________________________________
#### SPATIAL MODELLING OF UNIVARIATE METRICS: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Load package
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 11 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## The training set
train = sdata2[ msk,] 
dim(train)# 1256   13
head(train)

## The test set
test  = sdata2[ !msk,]
dim(test) # 140  13
head(test)

## Check number of observations for train (TRUE) and test (FALSE) sets
print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
dim(sdata) # 1396   12
dim(train)+dim(test)# 1396   26
#_______________________________________________________________________________
#### SPATIAL MODELLING OF UNIVARIATE METRICS: MODELLING ####

# This code is intended for a quick look-see. For final models and associated RF outputs (inc associated model confidence maps), use R script in file xx

## Load package
library(randomForest)

## Model
model <- value ~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV

## Prepare training data
train2 <- train[complete.cases(train), ]
head(train2)
str(train2)

## Run model
rf2 <- randomForest(model, data=train2)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF UNIVARIATE METRICS: PRODUCE FULL COVERAGE RASTER FOR METRIC ####

## Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)
plot(pr)

## Save as .tiff
#writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\S.tif',overwrite=TRUE,format = "GTiff")
#writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\N2.tif',overwrite=TRUE,format = "GTiff")
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: SETUP ####

## Load packages
library(sf)
library(dplyr)
library(RColorBrewer)
library(terra)
library(raster)
library(stars)
library(ggplot2)
library(colorRamps)
library(RColorBrewer)
library(terra)
library(tidyterra)
library(raster)
library(dplyr)
library(stars)
library(ggplot2)
library(colorRamps)
library(climateStability)
library(ggpubr)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))

## Data required for modelling biodiversity metrics (long format, outliers removed, all metrics)
head(univmeascoord2_mod)
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: RICHNESS MAP (POINT SAMPLES) ####

## Select richness data only
S_data <- univmeascoord2_mod[which(univmeascoord2_mod$metric=='S'),]

## Produce point plot
S_point_plot <- ggplot()+
  geom_sf(data=countries,fill ="black",col ="black")+ 
  geom_point(data = S_data,aes(x = Longitude_WGS84, y = Latitude_WGS84,col=measurement), size = 1)+
  scale_color_gradientn(colors= matlab.like(50))+#, breaks=c(1,2000,5000,10000,18000)
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10))+#make legend background transparent and text white
  theme(legend.position=c(0.87,0.2))+#legend.position = "bottom","none"#legend.position = "bottom"
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y = element_text(colour = "white"),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = "S",
           parse = TRUE,
           size=5,
           colour = "#707070",
           fontface="bold",
           x = -9,
           y = 59.5)+
theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: RICHNESS MODEL MAP ####

## Load model raster
S_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Taxon Richness/Model/S_Mean.tif')

## Reduce size of raster
S_model_agg <- S_model
#S_model_agg <- aggregate(S_model, fact=2,fun = modal)

## Remove raster to save space
rm(S_model)

## Label
label1 =  "~italic(S)"

## Model plot
S_model_plot <- ggplot() +
  geom_spatraster(data = S_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10))+
  theme(legend.position=c(0.87,0.2))+
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y = element_text(colour = "white"),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=5,
           colour = "#707070",
           fontface="bold",
           x = -9,
           y = 59.5)+
theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: RICHNESS CONFIDENCE MAP ####

## Load confidence raster
S_conf <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Taxon Richness/Confidence/S_CV.tif')

## Reduce size of confidence raster
S_conf_agg <- S_conf
#S_conf_agg <- aggregate(S_conf, fact=7,fun = modal)

## Remove confidence raster to save space
rm(S_conf)

## Label
label1 =  "~italic(S)"

## Confidence plot
S_conf_plot <- ggplot() +
  geom_spatraster(data = S_conf_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10))+
  theme(legend.position=c(0.87,0.2))+
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = label1,
           parse = TRUE,
           size=5,
           colour = "#707070",
           fontface="bold",
           x = -9,
           y = 59.5)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: SQRT ABUNDANCE (POINT SAMPLES) ####

## Select abundance data only
N_data <- univmeascoord2_mod[which(univmeascoord2_mod$metric=='N_sqrt'),]

## Produce point plot
N_point_plot <- ggplot()+
  geom_sf(data=countries,fill ="black",col ="black")+ 
  geom_point(data = N_data,aes(x = Longitude_WGS84, y = Latitude_WGS84,col=measurement), size = 1)+
  scale_color_gradientn(colors= matlab.like(50))+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10))+
  theme(legend.position=c(0.87,0.2))+
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0.3), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.y = element_text(colour = "black"),
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"))+
  annotate(geom = "text",
           label = "N",
           parse = TRUE,
           size=5,
           colour = "#707070",
           fontface="bold",
           x = -9,
           y = 59.5)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: SQRT ABUNDANCE MODEL MAP ####

## Load model raster
N_model <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Square-root Abundance/Model/N_trans_Mean.tif')

## Reduce size of raster
N_model_agg <-N_model
#N_model_agg <- aggregate(N_model, fact=7,fun = modal)

## Remove raster to save space
rm(N_model)

## Model plot
N_model_plot <- ggplot() +
  geom_spatraster(data = N_model_agg) +
  scale_fill_gradientn(colors= matlab.like(50),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10))+
  theme(legend.position=c(0.87,0.2))+
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    panel.background = element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank())+
  annotate(geom = "text",
           label = as.character(expression(~italic(N)~(x^{0.25}))),
           parse = T,
           size=5,
           colour = "#707070",
           fontface="bold",
           x = -8.6,
           y = 59.5)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: SQRT ABUNDANCE CONFIDENCE MAP ####

## Load model raster
N_conf <- rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Square-root Abundance/Confidence/N_trans_CV.tif')

## Reduce size of confidence raster
N_conf_agg <-N_conf
#N_conf_agg <- aggregate(N_conf, fact=7,fun = modal)

## Remove confidence raster to save space
rm(N_conf)

## Confidence plot
N_conf_plot <- ggplot() +
  geom_spatraster(data = N_conf_agg) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent")+
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 12)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 10))+
  theme(legend.position=c(0.87,0.2))+
  theme(legend.title= element_blank())+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.5, "cm"))+#t, r, b, l
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
  annotate(geom = "text",
           label = as.character(expression(~italic(N)~(x^{0.25}))),
           parse = TRUE,
           size=5,
           colour = "#707070",
           fontface="bold",
           x = -8.6,
           y = 59.5)+
  theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())
#_______________________________________________________________________________
#### MAP UNIVARIATE METRICS: COMBINED PLOT (FIGURE 2) ####

## Stitch plots together
S_stitch <- egg::ggarrange(S_model_plot,S_conf_plot, labels = c("",""),nrow=1)#ggpubr
N_stitch <- egg::ggarrange(N_model_plot,N_conf_plot, labels = c("",""),nrow=1)#ggpubr

## Combine columns
figure2 <- ggpubr::ggarrange(S_stitch,N_stitch,labels = c("a)", "b)"),nrow=2,font.label=list(color="black",size=12,face='plain'),align="v",widths = c(0.5,0.5,1))

## Add annotations
fig2 <- annotate_figure(
  figure2, 
  top = text_grob("         Model                                                                  Confidence", 
  color = "black", face = "bold", size = 12),
  bottom = text_grob("          Longitude", 
  color = "black", face = "plain", size = 12))

## Save plot
ggsave(plot = fig2,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_2.tif"),
       height = 200, width =210, units = "mm", dpi = 500,
       device = "tiff",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### ASSEMBLAGES: PREPARE DATA ####

## Transform the data (fourth-root transformation)
datat=data6^(0.25)
#_______________________________________________________________________________
#### ASSEMBLAGES: ELBOW PLOT (FIGURE 4A) ####

## Load package
library(factoextra)

## Generate plot
elbow <- fviz_nbclust(datat,kmeans, method = "wss",k.max = 30,linecolor = "black")+
  geom_vline(xintercept = c(11), linetype = 2)+theme_classic(base_size = 14)
plot(elbow)

#_______________________________________________________________________________
#### ASSEMBLAGES: K-MEANS CLUSTERING ####

## Perform Kmeans clusterinig of data. Results (cluster group) to the object 'results'
set.seed(1234)
results=kmeans(datat,11,algorithm="MacQueen",iter.max=100,nstart=25)

## Number of samples belonging to each cluster group
results$size #  179 110 100  42 144  39 111 346 131  68 130
#_______________________________________________________________________________
#### ASSEMBLAGES: DENDROGRAM (FIGURE 4B) ####

## Load packages
library(ggplot2)
library(ggdendro)
library(ggplot2)
library(dplyr)
library(dendextend)

## Function to calculate absolute differences between cluster centres over all variables.
nclusters = 11
absdiff = matrix(0, nrow=nclusters, ncol=nclusters)
centers = results$centers
for (j in 1:nclusters) {
  for (k in 1:nclusters) {
    absdiff[j,k] = sum(abs(centers[j,] - centers[k,]))
  }
}
d=round(absdiff, 1)

## Find distance matrix
d1 <- dist(as.matrix(d))

## Produce dendrogram
d2 <- d1%>% hclust %>% as.dendrogram %>%set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 9) %>%  # node point size
  set("leaves_col", c("blue2","cyan1","#05aae1","purple","#b40202","red1","darkorange","yellow","#b4b404","green3","palegreen1")) %>%#"darkorchid3",,"plum2"
  set("labels_cex",0.8)%>% 
  set("labels", c('     Epi_A1', '     Epi_B1a', '     Epi_B1b', '     Epi_C1', '     Epi_D1', '     Epi_D2a', '     Epi_D2b',  '     Epi_D2c','     Epi_D2d','     Epi_E1a','     Epi_E1b'))%>% 
  set("branches_lwd", 0.7)


## Change dendrogram into a ggplot
ggd1 <- as.ggdend(d2)
dendrogram <-  ggplot(ggd1, horiz = T)+theme_classic(base_size = 14)+theme(axis.title.y=element_blank(),
                                                               axis.text.y=element_blank(),
                                                               axis.ticks.y=element_blank(),
                                                              axis.line.y=element_blank())+labs(y='Height')
dendrogram
#_______________________________________________________________________________
#### ASSEMBLAGES: COMBINED ELBOW & DENDROGRAM (FIGURE 4) ####

## Combined elbow plot and dendrogram (hased out otherwise '## png' message displayed in mrkdown)
tiff("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/Figure_4.tif", width = 33, height = 20, units = "cm", res = 800,pointsize = 12)
ggpubr::ggarrange(elbow,NULL,dendrogram,labels = c("a)","", "b)"),nrow=1,widths = c(1, 0.05, 1))
dev.off()
#_______________________________________________________________________________
####ASSEMBLAGES: CREATE DATAFRAME FOR PRODUCING FAUNAL CLUSTER MAPS ####

## Add cluster group from kmeans results file to df 'pos' which includes 'Sample',
# 'Latitude_WGS84' and 'Longitude_WGS84'
faunal.cluster=cbind(pos2,results$cluster)

## Change name of col 'results$cluster' to 'ClusterNum'
names(faunal.cluster)[4]<-paste("ClusterNum")

## Add a new empty col 'FaunalCluster' to df 'faunal.cluster
faunal.cluster["FaunalCluster"]=NA

## Populate FaunalCluster col with new names (see dendrogram from Step 21)
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 5] <- "Epi_E1b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 1] <- "Epi_E1a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 9] <- "Epi_D2d"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 3] <- "Epi_D2c"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 8] <- "Epi_D2b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 2] <- "Epi_D2a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 7] <- "Epi_D1"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 4] <- "Epi_C1"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 11] <- "Epi_B1b"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 10] <- "Epi_B1a"
faunal.cluster$FaunalCluster[faunal.cluster$ClusterNum == 6] <- "Epi_A1"
names(faunal.cluster)
head(faunal.cluster)
#_______________________________________________________________________________
#### ASSEMBLAGES: OUTPUT SAMPLES WITH THEIR CLUSTER GROUP ID (TO FIND REPRESESENTATIVE SAMPLE IMAGES ####

## Load package
library(dplyr)

dim(faunal.cluster)#1400 5

## Join cluster results to data to access surveynames
result <- faunal.cluster %>%
  left_join(data, by = "samplecode")%>% 
  dplyr::select(surveyname,samplecode,FaunalCluster)%>%
distinct()
#View(result)

## Save files
write.csv(result, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\result.csv", row.names=FALSE)
#_______________________________________________________________________________
#### DATA FOR USE WITH MODELLING SCRIPTS ####
#S, N sqrt - from object univmeascoord2.mod' and cluster from obect 'faunal.cluster'

## Subset univariate metric results
unidat <- univmeascoord2_mod[which(univmeascoord2_mod$metric=='S'|univmeascoord2_mod$metric=='N_sqrt'),]

## Faunal cluster results
clusdat <- faunal.cluster
colnames(clusdat) <- c('samplecode', 'Latitude_WGS84', 'Longitude_WGS84', 'metric', 'measurement')# update column names
clusdat$metric <- 'cluster'# updat metric name

## Stitch together univariate and cluster results
data4mod <- rbind(unidat,clusdat)# 

## Update metric name
data4mod$metric[data4mod$metric =="N_sqrt"] <- "N_trans" # update metric name

## Inspect data
View(data4mod)

## Save file for use with modelling scripts
write.csv(data4mod, "C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\poseidon_trawl_metrics_4_modelling.csv", row.names=FALSE)
#_______________________________________________________________________________
#### ASSEMBLAGES: OUTPUT DATA FOR SIMPER ANALYSIS USING PRIMER ####

## Load package
library(dplyr)

## Add Sample code to transformed faunal data matrix
data4simper=cbind(faunal.cluster$samplecode,datat)

## Change name of col 1 from 'faunal.cluster$Sample' to 'samplecode'
colnames(data4simper)[1] <- "samplecode"

## Check df are same length
dim(data4simper)#1400 430
dim(faunal.cluster)#1400 5

## Check both dfs have a col named 'Samplecode'
colnames(data4simper)
colnames(faunal.cluster)

## Merge dataframes by col 'Sample
simper_merged_100 <- merge(faunal.cluster[,c(1,5)],data4simper, by='samplecode')
dim(simper_merged_100)#1400  431
#View(simper_merged_100)

## Remove 'Sample' column
simper_merged_100 <- simper_merged_100%>%
  select(-samplecode)

## Create df for treatment
simper_merged_100_cluster <- simper_merged_100[,1]

## Create df for fauna
simper_merged_100_fauna <- simper_merged_100%>%
  select(-FaunalCluster)
dim(simper_merged_100_fauna)#1400  429

## Export both objects as .csv files for use with PRIMER6
write.csv(simper_merged_100_fauna,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/DATAFORSIMPER_100.csv",row.names=FALSE)
write.csv(simper_merged_100_cluster,file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/FACTORFORSIMPER_100.csv",row.names=FALSE)
#_______________________________________________________________________________
#### ASSEMBLAGES: FAUNAL CLUSTER DISTRIBUTION MAP (POINTS) ####

## Load package
library(ggplot2)
library(cowplot)

## Produce map
p2= ggplot()+
geom_sf(data=countries, fill ="black",col ="black")+ 
 geom_point(data=faunal.cluster,aes(samplelongstart,samplelatstart,col=FaunalCluster),
 size=4,show.legend = TRUE)+
 scale_colour_manual(values = c("blue2","cyan1","#05aae1","purple","#b40202","red1","darkorange","yellow","#b4b404","green3","palegreen1"),name="Cluster")+
 guides(colour = guide_legend(override.aes = list(size=3)))+ # Change size of legend dots
 coord_sf(xlim = c(-12, 9),ylim = c(47, 62))+
 theme_bw(base_size = 18)+
 labs(x="Longitude",y="Latitude")+
 theme(legend.position=c(0.9,0.21))+
 theme(legend.background=element_blank(),
      legend.text = element_text(color="white",size= 12),
    legend.box.background = element_blank(),
    legend.key = element_blank())

fig4a <- p2 + guides(colour = guide_legend(override.aes = list(size=6)))

## Save legend
legendfclus <- get_legend(p2 + theme_bw(base_size=24)+ guides(colour = guide_legend(override.aes =
list(size=8))))
plot(legendfclus)
#_______________________________________________________________________________
#### ASSEMBLAGES: FAUNAL CLUSTER DISTRIBUTION (POINT - FACET BY FAUNAL CLUSTER GROUP) ####

## Produce map
p6= ggplot()+
geom_sf(data=countries, fill ="black",col ="black")+ 
 geom_point(data=faunal.cluster,aes(samplelongstart,samplelatstart,col=FaunalCluster),
 size=2,show.legend = FALSE)+
 coord_sf(xlim = c(-12, 9),ylim = c(47, 62))+
 theme_bw(base_size=18)+
 scale_colour_manual(values = c("blue2","cyan1","#05aae1","purple","#b40202","red1","darkorange","yellow","#b4b404","green3","palegreen1"),name="Cluster")+
  theme(
    axis.title.x=element_blank(),
    axis.text.x=element_text(colour = "white"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank())+
 guides(colour = guide_legend(override.aes = list(size=5)))+ # Change size of legend dots
 labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0,0,0.2), "cm"))+
 facet_wrap(~FaunalCluster)

fig4b=p6+guides(fill=FALSE)# remove bathy legend
#_______________________________________________________________________________
#### ASSEMBLAGES: CLUSTER GROUP CHARACTERISTICS ####

# Create taxon (family) names with phylum code

## Load package
library(dplyr)

## Create a df for cluster centers
cluster_centers <- as.data.frame(results$centers)
#View(cluster_centers)

## See taxon names from df 'cluster_centers'
names(cluster_centers)

## Remove unwanted column
cluster_centers <- cluster_centers %>% dplyr::select(-"EPIBENTHIC.MIX.UNIDENTIFIED.")

## Identify 'family' and 'phylum' columns from raw data
names(data2test2)

## Take just the 'family' and 'phylum' cols
fam_phy <-unique(data2test2[,c(9,10)])
head(fam_phy)
#View(fam_phy)

## Remove row where value in Family is 'EPIBENTHIC MIX UNIDENTIFIED'
fam_phy <- subset(fam_phy, family != 'EPIBENTHIC MIX UNIDENTIFIED')

## Check for presence of individual families
#fam_phy[fam_phy$family == "Aspidosiphonidae", ]

## Find Phyla
unique(fam_phy$phylum)

## Create new col for phylum code
fam_phy$phy_code <- fam_phy$phylum

## Create Phylum codes
fam_phy$phy_code[fam_phy$phy_code == 'Mollusca'] <- 'Mol'
fam_phy$phy_code[fam_phy$phy_code == 'Echinodermata'] <- 'Ech'
fam_phy$phy_code[fam_phy$phy_code == 'Annelida'] <- 'Ann'
fam_phy$phy_code[fam_phy$phy_code == 'Arthropoda'] <- 'Art'
fam_phy$phy_code[fam_phy$phy_code == 'Cnidaria'] <- 'Cni'
fam_phy$phy_code[fam_phy$phy_code == 'Sipuncula'] <- 'Si'
fam_phy$phy_code[fam_phy$phy_code == 'Nemertea'] <- 'N'
#fam_phy$phy_code[fam_phy$phy_code == 'Phoronida'] <- 'Phor'
#fam_phy$phy_code[fam_phy$phy_code == 'Hemichordata'] <- 'Hem'
fam_phy$phy_code[fam_phy$phy_code == 'Bryozoa'] <- 'Bry'
fam_phy$phy_code[fam_phy$phy_code == 'Platyhelminthes'] <- 'Pla'
fam_phy$phy_code[fam_phy$phy_code == 'Chordata'] <- 'Cho'
fam_phy$phy_code[fam_phy$phy_code == 'Porifera'] <- 'Por'
fam_phy$phy_code[fam_phy$phy_code == 'Entoprocta'] <- 'Ent'
fam_phy$phy_code[fam_phy$phy_code == 'Brachiopoda'] <- 'Bra'
fam_phy$phy_code[fam_phy$phy_code == 'Ctenophora'] <- 'Cte'
#fam_phy$phy_code[fam_phy$phy_code == 'Tracheophyta'] <- 'Tra'
#fam_phy$phy_code[fam_phy$phy_code == 'Ochrophyta'] <- 'Och'
#fam_phy$phy_code[fam_phy$phy_code == 'Myzozoa'] <- 'Myz'
#fam_phy$phy_code[fam_phy$phy_code == 'Priapulida'] <- 'Pri'
#fam_phy$phy_code[fam_phy$phy_code == 'Foraminifera'] <- 'For'
#fam_phy$phy_code[fam_phy$phy_code == 'Ciliophora'] <- 'Cil'
#fam_phy$phy_code[fam_phy$phy_code == 'Nemertina'] <- 'Ne'
#fam_phy$phy_code[fam_phy$phy_code == 'Sipunculida'] <- 'Sip'
#fam_phy$phy_code[fam_phy$phy_code == 'Pogonophora'] <- 'Pog'
#fam_phy$phy_code[fam_phy$phy_code == 'Nematomorpha'] <- 'Nem'
#fam_phy$phy_code[fam_phy$phy_code == 'Gastrotricha'] <- 'Gas'
#fam_phy$phy_code[fam_phy$phy_code == 'Coelenterata'] <- 'Coe'

## Create mew 'taxon' col which is the family name and bracketed phylum code
fam_phy$taxon <- paste(fam_phy$family," (", fam_phy$phy_code, ")", sep = "")

# Remove duplicate rows
fam_phy2 <- fam_phy %>% distinct()

dim(fam_phy2)# 458   4
head(fam_phy2)
#View(fam_phy2)

## Remove these records (where phyla has been revised and thus there are two values)
remove_values1 <- c('Aspidosiphonidae (Si)','Cnidaria (Coe)','Golfingiidae (Si)','Siboglinidae (Pog)','NA (Mol)','EPIBENTHIC MIX UNIDENTIFIED  (NA)')#
fam_phy3 <- fam_phy2 %>% dplyr::filter(!taxon %in% remove_values1)
#View(fam_phy3)

## Remove these records (phylum not present in the data)
remove_values2 <- c('Nemertina','Sipunculida')
fam_phy3 <- fam_phy3 %>% filter(!phylum %in% remove_values2)
head(fam_phy3)

## Remove rows with NA
fam_phy3 <- na.omit(fam_phy2)

## check above values deleted
fam_phy3[fam_phy3$taxon == "Cnidaria (Coe)", ]

# Create a vector of column names from df 'cluster_centers'
column_names <- as.data.frame(colnames(cluster_centers))
dim(column_names)# 428 1
head(column_names)

## Rename column
colnames(column_names)[1] <- 'family'

# Replace dots with spaces in the 'family' column
#column_names$family <- gsub("\\.", " ", column_names$family)

# Merge data frames by 'family', keeping only rows from df 'column_names'
merged_df <- left_join(column_names, fam_phy3, by = "family")
head(merged_df)
dim(merged_df)# 429   4
head(fam_phy3)
head(merged_df)

## Remove extra row from merged_df. Note the extra row stems from two records for Golfinidae (Ann) and Golfinidae (Si)
merged_df <- merged_df[merged_df$taxon != "Golfingiidae (Si)", ]

## Update column names (family and phy code)in df 'cluster_centers'
colnames(cluster_centers) <- merged_df$taxon

dim(cluster_centers) #11 428
#View(cluster_centers)
length(merged_df$taxon)#428
#_______________________________________________________________________________
#### ASSEMBLAGES: DATAFRAME FOR USE WITH SIMPER (VEGAN) ####
## Save df ' cluster_centres' as 'cluster_centres_4_simper' (see later in script). This is for simper based on vegan
cluster_centers_4_simper <- cluster_centers
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - PREPARE DATA FOR PRODUCING PHYLA PIE CHARTS ####

## Load package
library(tidyr)

## Inspect df 'cluster_centers'
head(cluster_centers)
#View(cluster_centers)
#class(cluster_centers)

## Transpose df 'cluster_centers' so taxa are rows
cluster_centers_t <- as.data.frame(t(cluster_centers))
head(cluster_centers_t)

# Redefine row and column names
rownames(cluster_centers_t) <- colnames(cluster_centers)
colnames(cluster_centers_t) <- rownames(cluster_centers)

## Make column of taxon names based on row names
cluster_centers_t$taxon <- rownames(cluster_centers_t)

## Remove rownames
rownames(cluster_centers_t) <- NULL
head(cluster_centers_t)

## Update column (cluster group) names
colnames(cluster_centers_t) <- c('Epi_E1a','Epi_D2a','Epi_D2c','Epi_C1','Epi_E1b','Epi_A1','Epi_D1','Epi_D2b','Epi_D2d','Epi_B1a','Epi_B1b','taxon')
head(cluster_centers_t)
dim(cluster_centers_t)# 428 12

## Add family and phyla info
cluster_centers_t_wt_fam <- left_join(cluster_centers_t,fam_phy2, by='taxon')
dim(cluster_centers_t_wt_fam)# 428 15
head(cluster_centers_t_wt_fam)

## Drop cols you don't need and reorder
cluster_centers_t_wt_fam2 <- cluster_centers_t_wt_fam[,c(12,14,1:11)]
cluster_centers_t_wt_fam2

## Turn data from wide to long format
data_long <- gather(cluster_centers_t_wt_fam2, cluster, count, Epi_E1a:Epi_B1b, factor_key=TRUE)
data_long

## Reorder columns
data_long2 <- data_long[,c(3,1,2,4)]
data_long2

## Create a vector specifying the number of rows to select from each group
n_vector <- c('Epi_E1a' = 28,'Epi_D2a' = 15,'Epi_D2c' = 20,'Epi_C1' = 23,'Epi_E1b' = 28,'Epi_A1' = 30,'Epi_D1' = 20 ,'Epi_D2b' = 13,'Epi_D2d' = 19,'Epi_B1a' = 39,'Epi_B1b' = 29 )

## Order data by group with count from high to low 
data_long3 <- data_long2%>%
  group_by(cluster,taxon)%>%
  arrange(desc(count))
head(data_long3)

## Select number of rows by group according to above vector
library(purrr) #load package

selected_data <- data_long3 %>%
  group_by(cluster) %>%
  group_split() %>%
  map2_dfr(names(n_vector), ~ slice(.x, 1:n_vector[.y]))
head(selected_data)
#View(selected_data)
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - FIND MOST COMMON PHYLA ####

## Number of groups by phyla
selected_data%>%count(phylum)%>% arrange(desc(n))

#  phylum            n
#  <chr>         <int>
#1 Arthropoda       85
#2 Chordata         56
#3 Echinodermata    43
#4 Mollusca         37
#5 Cnidaria         23
#6 Bryozoa          10
#7 Annelida          9
#8 Porifera          1

## Update column names
colnames(selected_data)[1] <- "Cluster"

# Create a directory for the images
image_dir <- "pie_images"
dir.create(image_dir, showWarnings = FALSE)
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - ID COLOURS FOR PIE CHARTS ####

## Load package
library(scales)

#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(11)
hex
# "#F8766D" "#DB8E00" "#AEA200" "#64B200" "#00BD5C" "#00C1A7" "#00BADE" "#00A6FF" "#B385FF" "#EF67EB" "#FF63B6"
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - PRODUCE PIE CHARTS FOR THE 11 DIFFERENT CLUSTER GROUPS ####

## PIE Epi_A1
selected_data_Epi_A1 <- selected_data[which(selected_data$Cluster=='Epi_A1'),]
selected_data_Epi_A1

## Identify phyla for cluster
unique(selected_data_Epi_A1$phylum)

p_Epi_A1 <- ggplot(selected_data_Epi_A1, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
  theme(legend.position="none")
p_Epi_A1

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_A1.png", plot = p_Epi_A1, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_B1a
selected_data_Epi_B1a <- selected_data[which(selected_data$Cluster=='Epi_B1a'),]
selected_data_Epi_B1a

## Identify phyla for cluster
unique(selected_data_Epi_B1a$phylum)

p_Epi_B1a <- ggplot(selected_data_Epi_B1a, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_B1a

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_B1a.png", plot = p_Epi_B1a, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_B1b
selected_data_Epi_B1b <- selected_data[which(selected_data$Cluster=='Epi_B1b'),]
selected_data_Epi_B1b

## Identify phyla for cluster
unique(selected_data_Epi_B1b$phylum)

p_Epi_B1b <- ggplot(selected_data_Epi_B1b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_B1b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_B1b.png", plot = p_Epi_B1b, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_C1
selected_data_Epi_C1 <- selected_data[which(selected_data$Cluster=='Epi_C1'),]
selected_data_Epi_C1

## Identify phyla for cluster
unique(selected_data_Epi_C1$phylum)

p_Epi_C1 <- ggplot(selected_data_Epi_C1, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_C1

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_C1.png", plot = p_Epi_C1, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_D1
selected_data_Epi_D1 <- selected_data[which(selected_data$Cluster=='Epi_D1'),]
selected_data_Epi_D1

## Identify phyla for cluster
unique(selected_data_Epi_D1$phylum)

p_Epi_D1 <- ggplot(selected_data_Epi_D1, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_D1

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_D1.png", plot = p_Epi_D1, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_D2a
selected_data_Epi_D2a <- selected_data[which(selected_data$Cluster=='Epi_D2a'),]
selected_data_Epi_D2a

## Identify phyla for cluster
unique(selected_data_Epi_D2a$phylum)

p_Epi_D2a <- ggplot(selected_data_Epi_D2a, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_D2a

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_D2a.png", plot = p_Epi_D2a, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_D2b
selected_data_Epi_D2b <- selected_data[which(selected_data$Cluster=='Epi_D2b'),]
selected_data_Epi_D2b

## Identify phyla for cluster
unique(selected_data_Epi_D2b$phylum)

p_Epi_D2b <- ggplot(selected_data_Epi_D2b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_D2b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_D2b.png", plot = p_Epi_D2b, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_D2c
selected_data_Epi_D2c <- selected_data[which(selected_data$Cluster=='Epi_D2c'),]
selected_data_Epi_D2c

## Identify phyla for cluster
unique(selected_data_Epi_D2c$phylum)

p_Epi_D2c <- ggplot(selected_data_Epi_D2c, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_D2c

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_D2c.png", plot = p_Epi_D2c, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_D2d
selected_data_Epi_D2d <- selected_data[which(selected_data$Cluster=='Epi_D2d'),]
selected_data_Epi_D2d

## Identify phyla for cluster
unique(selected_data_Epi_D2d$phylum)

p_Epi_D2d <- ggplot(selected_data_Epi_D2d, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_D2d

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_D2d.png", plot = p_Epi_D2d, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_E1a
selected_data_Epi_E1a <- selected_data[which(selected_data$Cluster=='Epi_E1a'),]
selected_data_Epi_E1a

## Identify phyla for cluster
unique(selected_data_Epi_E1a$phylum)

p_Epi_E1a <- ggplot(selected_data_Epi_E1a, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_E1a

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_E1a.png", plot = p_Epi_E1a, width = 6, height = 4)
#_______________________________________________________________________________
## PIE Epi_E1b
selected_data_Epi_E1b <- selected_data[which(selected_data$Cluster=='Epi_E1b'),]
selected_data_Epi_E1b

## Identify phyla for cluster
unique(selected_data_Epi_E1b$phylum)

p_Epi_E1b <- ggplot(selected_data_Epi_E1b, aes(x = "", y = count, fill = phylum)) +#
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c('Annelida'='#F8766D','Arthropoda'='#DB8E00','Bryozoa'='#AEA200','Chordata'='#64B200','Cnidaria'='#00BD5C','Echinodermata'='#00C1A7','Entoprocta'='#00BADE','Nemertea'='#00A6FF','Mollusca'='#B385FF','Phoronida'='#EF67EB','Porifera'='#FF63B6'))+#
  coord_polar(theta = "y") +
  theme_void()+
theme(legend.position="none")

p_Epi_E1b

## Save plot
ggsave("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/pie_Epi_E1b.png", plot = p_Epi_E1b, width = 6, height = 4)
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - CONVERT PIE CHARTS TO IMAGES FOR USE IN GT TABLE ####

## Load packages
library(base64enc)
library(dplyr)

# Set the directory containing the .png files
directory <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/"

# List all .png files in the directory
png_files <- list.files(directory, pattern = "\\.png$", full.names = TRUE)

# Create a dataframe from the list of .png files
pie_images<- data.frame(files = png_files)

## Add column for Cluster
pie_images$Cluster <- c('Epi_A1','Epi_B1a','Epi_B1b','Epi_C1','Epi_D1','Epi_D2a','Epi_D2b','Epi_D2c','Epi_D2d','Epi_E1a','Epi_E1b')

##Update column order
pie_images <- pie_images[,2:1]

# Print the dataframe
print(pie_images)

# Set the directory containing the .png files
directory <- "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/R/pie_images/"

# List all .png files in the directory
png_files <- list.files(directory, pattern = "\\.png$", full.names = TRUE)

# Function to convert a file to base64
convert_to_base64 <- function(file) {
  base64encode(file)
}

# Apply the function to all .png files and create a dataframe
base64_df <- data.frame(
  file = png_files,
  base64 = sapply(png_files, convert_to_base64)
)

# Print the dataframe
View(base64_df)

## Create a df for rownames from df 'base64_df'. Names contain the cluster group name
Cluster <- rownames(base64_df)

## Remove rownames
rownames(base64_df) <- NULL

## Add cluster col to df 'base64_df'
pie_images <- cbind(Cluster,base64_df)

## Update column names
colnames(pie_images)[2] <- 'Image'
colnames(pie_images)[3] <- 'Image_base64'
#View(pie_images)

## Update cluster groups
pie_images$Cluster <- c('Epi_A1','Epi_B1a','Epi_B1b','Epi_C1','Epi_D1','Epi_D2a','Epi_D2b','Epi_D2c','Epi_D2d','Epi_E1a','Epi_E1b')

# Convert images to base64
pie_images2 <- pie_images %>%
  rowwise() %>%
  mutate(Image_base64 = base64enc::dataURI(file = Image, mime = "image/png"))
head(pie_images2)
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - IDENTIFY TAXA TO INCLUDE IN CHARACTERISTICS TABLE ####

## Analyse cluster centers
head(cluster_centers)
#View(cluster_centers)

# Function to return column names in order of values for each row
get_ordered_colnames <- function(row) {
  colnames(cluster_centers)[order(-row)]
}

# Apply the function to each row
result <- apply(cluster_centers, 1, get_ordered_colnames)

# Convert the result to a data frame for better readability
result_df <- data.frame(t(result))
colnames(result_df) <- paste0("Rank_", 1:ncol(result_df))
#View(result_df)

# Vector specifying the number of columns to keep for each row (this is the mean number of taxa)
N <- c(28,15,20,23,28,30,20,13,19, 39,29)

# Function to return top N column names in order of values for each row, highest first
get_top_n_colnames <- function(row, n) {
  colnames(cluster_centers)[order(-row)][1:n]
}
## Convert df 'cluster_centers' to a matrix
cluster_centers <- as.matrix(cluster_centers)

# Apply the function to each row with different N values
result <- mapply(get_top_n_colnames, split(cluster_centers, row(cluster_centers)), N, SIMPLIFY = FALSE)

# Convert the list to a data frame for better readability
result_df <-as.data.frame( do.call(rbind, lapply(result, function(x) {
  length(x) <- max(N)
  return(x)
})))

# Label the rows with row numbers
rownames(result_df) <- c('Epi_E1a','Epi_D2a','Epi_D2c','Epi_C1','Epi_E1b','Epi_A1','Epi_D1','Epi_D2b','Epi_D2d','Epi_B1a','Epi_B1b')

# Function to convert row values to a comma-separated string
row_to_string <- function(row) {
  paste(na.omit(row), collapse = ", ")
}

# Apply the function to each row
result2 <- apply(result_df, 1, row_to_string)

# Convert the result to a data frame for better readability
result_df2 <- data.frame(Row = paste0("Row_", 1:nrow(result_df)), Values = result2)
result_df2 <- data.frame(result2)
colnames(result_df2)[1] <- 'Taxa'

## Replace cluster labels
result_df2$Cluster <- c('Epi_E1a','Epi_D2a','Epi_D2c','Epi_C1','Epi_E1b','Epi_A1','Epi_D1','Epi_D2b','Epi_D2d','Epi_B1a','Epi_B1b')

## Remove row names
rownames(result_df2) <- NULL

# Move the last column to the first position
result_df2 <- result_df2[, c(ncol(result_df2), 1:(ncol(result_df2)-1))]
#View(result_df2)

## Get rows in correct order
chr_taxa <- result_df2[order(result_df2$Cluster, decreasing = F),]
#_______________________________________________________________________________
## ASSEMBLAGES: CHARACTERISTICS - UNIVARIATE SUMMARY MEASURES ####

## Load packages
library(vegan)
library(dplyr)

## Calculate univariate summary measures based on faunal abundance data in df 'data5'
Richness = specnumber(data6) # Species Richness(S)
Abundance=rowSums(data6) # Abundance

## Need cluster results before this works
uni <- cbind(pos2,Richness,Abundance,results$cluster)
colnames(uni)[6] <- "Cluster"

## Calculate metrics
uni2 <- uni %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(n = n(),r_mean = mean(Richness),r_sd=sd(Richness),a_mean=mean(Abundance),a_sd=sd(Abundance) )

## Replace cluster numbers for names
uni2$Cluster[uni2$Cluster == "1"] <- 'Epi_E1a'
uni2$Cluster[uni2$Cluster == "2"] <- 'Epi_D2a'
uni2$Cluster[uni2$Cluster == "3"] <- 'Epi_D2c'
uni2$Cluster[uni2$Cluster == "4"] <- 'Epi_C1'
uni2$Cluster[uni2$Cluster == "5"] <- 'Epi_E1b'
uni2$Cluster[uni2$Cluster == "6"] <- 'Epi_A1'
uni2$Cluster[uni2$Cluster == "7"] <- 'Epi_D1'
uni2$Cluster[uni2$Cluster == "8"] <- 'Epi_D2b'
uni2$Cluster[uni2$Cluster == "9"] <- 'Epi_D2d'
uni2$Cluster[uni2$Cluster == "10"] <- 'Epi_B1a'
uni2$Cluster[uni2$Cluster == "11"] <- 'Epi_B1b'

## Update column names
colnames(uni2)[1] <- "Cluster"
head(uni2)

## Order by r_mean
uni2 <- uni2[order(-uni2$r_mean),]
head(uni2)
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - ADD UNIVARIATE RESULTS TO TABLE OF CHARACTERISING TAXA ####

## Inspect existing table of characterising taxa
#View(chr_taxa)

## Merge dataframes
uni3 <- merge(uni2,chr_taxa,by="Cluster")
#View(uni3)

# Merge the base64 images data frame with the original data
merged_data <- left_join(uni3, pie_images2, by = "Cluster")
#View(merged_data)
#View(pie_images2)
#class(pie_images2)
#str(pie_images2)
#class(uni3)
#str(uni3)

## Rename
uni3 <- merged_data
names(uni3)

## Change colname from "Image_base64" to "Phylum" as this will be the col header in gt output
colnames(uni3)[9] <- 'Phylum'
#print(uni3)
#_______________________________________________________________________________
#### ASSEMBLAGES: CHARACTERISTICS - CREATE NICE TABLE FOR CLUSTER SUMMARY INFO ####

## Load packages
library(gt)
library(gtExtras)
library(dplyr)
library(stringr)
library(gtsummary)
library(webshot2)
library(chromote)

# Function to make selected text bold
make_bold <- function(text, bold_words) {
  for (word in bold_words) {
    text <- stringr::str_replace_all(text, word, paste0("<b>", word, "</b>"))
  }
  return(text)
}



# Apply the function to the 'Taxa' column (characterisinig taxa from SIMPER)
uni3$Taxa <- mapply(make_bold, uni3$Taxa, list(c("Pandalidae",#T_A1
                                                 "Polybiidae",
                                                 "Asteriidae",
                                                 "Crangonidae",
                                                 "Paguridae",
                                                 "Flustridae"),
                                               c("Parechinidae",#T_B1a
                                                 "Pectinidae",
                                                 "Paguridae",
                                                 "Inachidae",
                                                 "Asteriidae",
                                                 "Ophiuridae",
                                                 "Ophiotrichidae",
                                                 "Polybiidae",
                                                 "Galatheidae",
                                                 "Ascidiidae"),
                                               c("Paguridae",#T_B1b
                                                 "Pectinidae",
                                                 "Parechinidae",
                                                 "Inachidae",
                                                 "Pandalidae",
                                                 "Asteriidae",
                                                 "Polybiidae",
                                                 "Ophiotrichidae"),
                                               c("Crangonidae",#T_C1
                                                "Nuculidae",
                                                 "Nephropidae",
                                                 "Pandalidae",
                                                 "Processidae"),
                                               c("Asteriidae",#T_D1
                                                 "Soleidae",
                                                 "Polybiidae",
                                                 "Pleuronectidae",
                                                 "Astropectinidae"),
                                               c("Ophiuridae",#T_D2a
                                                 "Asteriidae",
                                                 "Paguridae"),
                                               c("Crangonidae",#T_D2b
                                                 "Paguridae",
                                                 "Asteriidae",
                                                 "Polybiidae"),
                                               c("Crangonidae",#T_D2c
                                                 "Ophiuridae",
                                                 "Gobiidae",
                                                 "Soleidae",
                                                 "Paguridae",
                                                 "Polybiidae"),
                                               c("Crangonidae",#T_D2d
                                                 "Ophiuridae",
                                                 "Paguridae",
                                                 "Asteriidae"),
                                               c("Paguridae",#T_E1a
                                                 "Crangonidae",
                                                 "Ophiuridae",
                                                 "Polybiidae",
                                                 "Inachidae",
                                                 "Gobiidae",
                                                 "Callionymidae"),
                                               c("Paguridae",#T_E1b
                                                 "Astropectinidae",
                                                 "Buccinidae",
                                                 "Asteriidae",
                                                 "Colidae",
                                                 "Crangonidae",
                                                 "Hormathiidae",
                                                 "Sertulariidae")
                                               ))
## Check text has bold markers added
head(uni3)
names(uni3)
#View(uni3)

## Drop std dev cols
  uni3 <- uni3[,c(1,2,3,5,9,7,8)]
  
## Update column names
colnames(uni3) <- c("Cluster","n","Richness","Abundance","Phylum","Taxa","Image")
 
#______________________________________________________________________________
### ASSEMBLAGES: CHARACTERISTICS - CALCULATE AREA OF EACH CLUSTER GROUP FOR INCLUSION IN CHARACTERISTICS TABLE ####

## Load packages
library(raster)
library(dplyr)

## Load raster data
#r <-raster('C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\DATA\\ANNA MODEL OUTPUTS\\TrawlAssemblageMaxClass_Sept24.tif')
r <-raster('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Assemblage/Model/TrawlAssemblageMaxClass_Sept24.tif')
#plot(r)

## Calculate the area of each cell in square kilometers
cell_areas <- area(r)

## Convert raster to data frame
raster_df <- as.data.frame(r, xy = TRUE)
raster_df$area_km2 <- cell_areas[]

##Filter out NA values, then summarize the total area by class
area_by_class <- raster_df %>%
  filter(!is.na(TrawlAssemblageMaxClass_Sept24)) %>%
  group_by(TrawlAssemblageMaxClass_Sept24) %>%
  summarize(total_area_km2 = sum(area_km2, na.rm = TRUE))

## Calculate the total area, ignoring NAs
total_area <- sum(area_by_class$total_area_km2, na.rm = TRUE)

## Add a percentage column
area_by_class <- area_by_class %>%
  mutate(percentage = (total_area_km2 / total_area) * 100)

## Add in cluster groups
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 1] <- "Epi_A1"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 2]<- "Epi_B1a"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 3] <- "Epi_B1b"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 4]<- "Epi_C1"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 5] <- "Epi_D1"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 6] <- "Epi_D2a"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 7] <- "Epi_D2b"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 8] <- "Epi_D2c"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 9] <- "Epi_D2d"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 10] <- "Epi_E1a"
area_by_class$TrawlAssemblageMaxClass_Sept24[area_by_class$TrawlAssemblageMaxClass_Sept24 == 11]<- "Epi_E1b"

## Update column names
colnames(area_by_class)[1] <- 'Cluster'

## View results
print(area_by_class)

## Join to other table data
uni4<- left_join(uni3,area_by_class, by='Cluster')
names(uni4)
uni5 <- uni4[,c(1,9,3,4,5,6,7)]
names(uni5)
View(uni5)
#______________________________________________________________________________
### ASSEMBLAGES: CHARACTERISTICS - CREATE TABLE 2 ####

## Load package
library(magick)

## Produce gt table
 uni_tab <-  uni5%>%gt()%>%
  fmt_number(
    columns = vars(percentage),
    decimals = 1
  )%>%
   # Reduce number of decimal places
   #fmt_number(columns = c(Richness,	Abundance), decimals = 0)%>%

   gt_plt_bar(column = Richness, scale_type = "number",  color = "grey",text_color = "black",width=60)%>% #keep_column = TRUE,,text_color = "black",,width = 80
   gt_plt_bar(column = Abundance, scale_type = "number", color = "grey",text_color = "black",width=60) %>%#keep_column = TRUE, scale_type = "number",

   #cols_hide(columns = c(n))%>%
   ## Add pie charts for phyla
   text_transform(
     locations = cells_body(columns = c(Phylum)),
     fn = function(x) {
       web_image(url = x, height = 100)
     }
   ) %>%
   cols_hide(columns = c(Image))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#9aff9a")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 11)
     ))  %>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#00cd00")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 10)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#b4b404")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 9)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#ffff00")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 8)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#ff8c00")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 7)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#ff0000"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 6)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#b40202"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 5)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#9a32cd"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 4)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#05aac1"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 3)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#00ffff")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 2)
     ))%>%
   # colour col header Clus2
   tab_style(
     style = list(
       cell_fill(color = "#0000ee"),
       cell_text(color = "white")
     ),
     location = list(
       cells_body(column = vars(Cluster),rows = 1)
     ))%>%
 # Apply HTML rendering (bold taxon names)
 fmt_markdown(columns = vars(Taxa))%>%
   ## Add a table caption
   #tab_header(
   #  title = md("<div style='text-align: left;'><b>Table 1</b>. Cluster group characteristics including mean richness, mean abundance, phylytic composition (pie charts) and characterising taxa (family level or above). Listed taxa are those with the highest mean centroid values, with the number of taxa reported based on the group mean richness value. Highlighted taxa are those identified by a SIMPER analysis as contributing to ~50% of the similarity between samples. Phyla codes are given in parenthesis (see Legend at foot of table).<br><br> </div>")
   #)%>%#: 'Ann' - Annelida, 'Art' - Arthropoda, 'Bra' - Brachiopoda, 'Bry' - Bryozoa, 'Cil' - Ciliophora, 'Cho' - Chordata, 'Coe' - Coelenterata, 'Cni' - Cnidaria, 'Ech' - Echinodermata, 'Ent' - Entoprocta, 'For' - Foraminifera, 'Gas' - Gastrotricha, 'Hem' - Hemichordata, 'Mol' - Mollusca, 'Myz' - Myzozoa, 'N' - Nemertea, 'Ne' - Nemertina, 'Nem' - Nematomorpha, 'Och' - Ochrophyta, 'Phor' - Phoronida, 'Pla' - Platyhelminthes, 'Pog' - Pogonophora, 'Por' - Porifera, 'Pri' - Priapulida, 'Si' - Sipuncula, 'Sip' - Sipunculida, 'Tra' - Tracheophyta)
   tab_options(
     table.border.top.style = "hidden"
     #heading.border.bottom.style = "hidden"
   )%>%tab_options(
     heading.title.font.size = px(18)  # Adjust the font size as needed
   )%>% 
   cols_label(
         Cluster = md("**Cluster**"),
         percentage = md("**% Area**"),
         Richness = md("**Richness (mean)**"),
         Abundance = md("**Abundance (mean)**"),
          Phylum = md("**Phyla**"),
         Taxa = md("**Taxa**"))%>%
   cols_align(
          align = "center",
     columns = everything())%>%
     ## Right align Taxa column
     cols_align(
      # uni_tab2,
       align = "left",
       columns = 'Taxa'
     )%>%
 ## Add pie legend
 
#tab_caption(
    #caption = md("**Legend:** 
    tab_footnote(
    footnote = md("**Pie Chart Legend:**
<span style='color:#DB8E00;font-size: 25px;'>■</span> Arthropoda (Art), 
<span style='color:#64B200;font-size: 25px;'>■</span> Chordata (Cho), 
<span style='color:#00C1A7;font-size: 25px;'>■</span> Echinodermata (Ech), 
<span style='color:#B385FF;font-size: 25px;'>■</span> Mollusca (Mol),
<span style='color:#00BD5C;font-size: 25px;'>■</span> Cnidaria (Cni),
<span style='color:#AEA200;font-size: 25px;'>■</span> Bryozoa (Bry), 
<span style='color:#F8766D;font-size: 25px;'>■</span> Annelida (Ann), 
<span style='color:#FF63B6;font-size: 25px;'>■</span> Porifera (Por)"
))%>%
  tab_options(
  table.font.size = "18px",
  column_labels.font.size = "18px"
)
 
  uni_tab 
 
## Save cluster characteristics table as .png
 uni_tab %>%
   gt::gtsave(
     "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/Table_2.png",vwidth = 1400, vheight = 1600)#,  zoom = 2, expand = 10

## Convert PNG to TIFF
 img <- image_read("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/Table_2.png")
 image_write(img, "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/Table_2.tif", format = "tiff")
#_______________________________________________________________________________
#### SPATIAL MODELLING OF ASSEMBLAGES: PREPARE DATA ####
 
# This code is intended for a quick look-see. For final model and associated RF outputs (inc associated model confidence map), use R script in file xx

#FaunalCluster <- read.csv("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\faunal.cluster6.csv",header=TRUE, stringsAsFactors=FALSE)
## Data used for modelling (from step 2.3.2.5)
BiodivCluster <- faunal.cluster[,c(1,3,2,5)]
head(BiodivCluster)
unique(BiodivCluster$FaunalCluster)

## Change names of cols
colnames(BiodivCluster)=c("Sample","lon","lat","cluster")
head(BiodivCluster)
dim(BiodivCluster)# 1400     4

## Unload extract function from tidyr package (otherwise it won't work)
#.rs.unloadPackage("tidyr")

## Extract variables from raster stack
sdata <- raster::extract(predictors, BiodivCluster[,2:3])
head(sdata)
#class(sdata)

## Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(BiodivCluster$Sample,BiodivCluster$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"
head(sdata2)
#View(sdata2)

## Remove rows with NA
sdata2 <- sdata2[complete.cases(sdata2), ]
unique(sdata2$Cluster)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF ASSEMBLAGES: MAKE TRAINING AND TESTING SET ####

## Equal splitting code, needs the caTools library;  sdata is the dataframe with cols for response and predictor variables

## Load package
library(caTools)

## Vector for response variable
Y <- sdata2[,2]

## Take a random 90% of samples in proportion across the 11 groups
set.seed(2)
msk= sample.split( Y, SplitRatio = 9/10, group = NULL )

## Check it was 90% 10% split

## The training set
train = sdata2[ msk,] 
#dim(train)#1260 14
#View(train)

## Remove station labels for train
train2 =train[,2:14]
#View(train2)

## The test set
test  = sdata2[ !msk,]
#dim(test)#140  14
#View(test)

## Remove station labels for test
test2 =test[,2:14]
#class(test2)
#View(test2)
#str(test2)

## Check number of observations for train (TRUE) and test (FALSE) sets
#print(table(Y, msk)) 

## Check number of samples in train and test sets is equal to total
#dim(sdata) # 1400 12
#dim(train)+dim(test)# 1400   28
#_______________________________________________________________________________
#### SPATIAL MODELLING OF ASSEMBLAGES: DO MODELLING ####

# This code is intended for a quick look-see. For final model and associated RF outputs (inc associated model confidence map), use R script in file xx

## Load package
library(randomForest)

## Model
model <- factor(Cluster)~Bathymetry+Current+Gravel+Light+Mud+Oxygen+Phyto+Salinity+Silicate+SPM+Temperature+WOV# cluster (multiple biodiv metrics)

## Run model
train3 <- train2[complete.cases(train2), ]
#str(train3)

## Make Cluster a factor
train3$Cluster <- as.factor(train3$Cluster)
#str(train3)
#head(train3)

## Run model
rf2 <- randomForest(model, data=train3)
#_______________________________________________________________________________
#### SPATIAL MODELLING OF ASSEMBLAGES: PRODUCE FULL COVERAGE RASTER FOR METRIC ####

##Use model to predict cluster group for each raster cell
pr <- predict(predictors, rf2)

## Basic plot
plot(pr)
#plot(pr, main='Random Forest, regression')#
#plot(wrld_simpl, add=TRUE, border='dark grey')

## Plot sample positions
#plot(FaunalCluster[,2:3])

## Save raster
writeRaster(pr,'C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\trawl_cluster.tif',overwrite=TRUE,format = "GTiff")
#_______________________________________________________________________________
#### MAP ASSEMBLAGES: SETUP ####

## Load packages
library(sf)
library(RColorBrewer)
library(terra)
library(tidyterra)
library(raster)
library(stars)
library(ggplot2)
library(ggpubr)

## Bring in countries polygon
countries <- st_read(file.path("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicHotSpots\\DATA\\EuropeLiteScoWal.shp"))

## Set CRS
st_crs(countries) = 4326
#_______________________________________________________________________________
#### MAP ASSEMBLAGES: MODEL ####

## Assemblage model layer
biodiv <-rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Assemblage/Model/TrawlAssemblageMaxClass_Sept24.tif')

## Reduce raster resolution
biodiv.agg  <-biodiv
#biodiv  <- aggregate(biodiv, fact=4,fun = modal)

## Remove NaN values
biodiv <- app(biodiv, fun=function(x) { ifelse(is.finite(x), x, NA) })

## Make cluster a factor
values(biodiv) <- as.factor(values(biodiv))

## Produce map
pbio <-ggplot() +
  geom_spatraster(data = biodiv) +#.agg
  scale_fill_manual(values = c("blue2","cyan1","#05aae1","purple","#b40202","red1","darkorange","yellow","#b4b404","green3","palegreen1"),labels = c("Epi_A1","Epi_B1a","Epi_B1b","Epi_C1","Epi_D1","Epi_D2a","Epi_D2b","Epi_D2c","Epi_D2d","Epi_E1a","Epi_E1b"),name="Cluster", na.translate = F)+#na.value="transparent" "red1","#b4b404"
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 24)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20))+
  theme(legend.position=c(0.9,0.21))+
  theme(axis.title.x = element_text(colour = "black"))+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.6, "cm"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
  
unique(faunal.cluster$FaunalCluster)
#_______________________________________________________________________________
#### MAP ASSEMBLAGES: CONFIDENCE ####

## Assemblage model confidence layer
biodiv_conf <-rast('C:/Users/kmc00/OneDrive - CEFAS/POSIDEN/POSEIDON RASTERS FOR AGE UPLOAD/Epifauna paper/Assemblage/Confidence/TrawlAssemblageConfidence_Sept24.tif')

## Reduce resolution of confidence layer
biodiv_conf.agg  <-biodiv_conf
#biodiv_conf.agg  <- aggregate(biodiv_conf, fact=4,fun = modal)

# stack
r_conf_list <- list(biodiv_conf.agg)

r_conf <- rast(r_conf_list)

## Rename
names(r_conf) <- c('CV')

## Convert to stars object
biodiv_conf_stars <- r_conf %>%
  st_as_stars()

## Produce map
pbio_conf <- ggplot() +
  geom_stars(data = biodiv_conf_stars) +
  scale_fill_gradientn(colors= brewer.pal(n = 5, name = "Greys"),na.value="transparent",name = "Condidence")+#Oranges
  geom_sf(data=countries, fill ="black",col ="black")+
  coord_sf(xlim = c(-10, 9),ylim = c(49, 60))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw(base_size = 24)+
  theme(legend.background=element_blank(),legend.text = element_text(color="white",size= 20))+
  theme(legend.position=c(0.9,0.2))+
  theme(
   axis.title.y =element_blank(),
   axis.text.y=element_blank())+
  labs(x="Longitude",y="Latitude")+
  theme(plot.margin = unit(c(0,0.2,0,0), "cm"),legend.key.size = unit(0.6, "cm"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
#_______________________________________________________________________________
#### MAP ASSEMBLAGES COMBINED MODEL & CONFIDENCE MAPS (FIGURE 5) ####

## Stitch plots together
model_stitch <- egg::ggarrange(pbio,pbio_conf, labels = c("",""),nrow=1)#ggpubr

## Annotate combined plot
fig4mod <- annotate_figure(
  model_stitch, 
  top = text_grob("               Model                                                                                 Confidence", 
  color = "black", face = "bold", size = 26),
  bottom = text_grob("          Longitude", 
  color = "black", face = "plain", size = 24),
  left = text_grob("Latitude", 
  color = "black", face = "plain", size = 24,rot = 90)
)

## Save plot
ggsave(plot = fig4mod,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_5.tif"),
       height = 260, width =510, units = "mm", dpi = 500,
       device = "tiff",limitsize = FALSE,bg="white")
#_______________________________________________________________________________
#### 2.5 EXPLAINING PATTERNS: PREPARE DATA

# Load packages
library(raster)

## Load biodiv data used for clustering (i.e. outliers removed, transformed, standardized, covariates removed - see step 2.3.1.5)
bio <- cbind(pos2,data6)

## Load factor data (faunal.cluster7 from step 2.3.2.5)
fac <- faunal.cluster[,c(1,3,2,5)]

## Change names of cols
colnames(fac)=c("Sample","lon","lat","cluster")

## Make cluster a factor
fac$cluster=as.factor(fac$cluster)


#############
## Get PHY data from rasters using variables used in RF modelling. Note sediment data comes from actual samples
so_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/so_mean.tif")# SALINITY MEAN
thetao_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/thetao_mean.tif")# BOTTOM TEMP
SPM_WINTER<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/SPM_WINTER.tif")#WINTER SPM
Current_Sp<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Current_Sp.tif")# CURRENT SPEED
gravel<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Gravel_Fraction.tif")#GRAVEL
Wave_veloc<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/wave_veloc.tif")# WAVE VELOCITY
vd<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/VD.tif")# VALLEY DEPTH
dfe_mean<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/dfe_mean.tif")# DISSOLVED IRON
RSP<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/RSP.tif")#RELATIVE SLOPE POSITION
LSF<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/LSF.tif")#LS-FACTOR
mud<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Mud_Fraction.tif")#MUD
CDP0<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/CDP0.tif")#CLOSED DEPRESSIONS

## Create raster stack
predictors <- stack(
  so_mean,
  thetao_mean,
  SPM_WINTER,
  Current_Sp,
  gravel,
  Wave_veloc,
  vd,
  dfe_mean,
  RSP,
  LSF,
  mud,
  CDP0
)

## Plot raster stack
plot(predictors)

## Update names for predictor variables
names(predictors)=c(
  'Salinity mean',
  'Bottom temp.',
  'Winter SPM',
  'Current speed',
  'Gravel',
  'Wave velocity',
  'Valley depth',
  'Diss. Iron',
  'Rel. slope pos.',
  'LS-factor',
  'Mud',
  'Closed depressions'
  
)

## Unload extract function from tidyr package (otherwise it won't work)
#.rs.unloadPackage("tidyr")

## Extract predictor variables from raster stack and store in df
sdata <- raster::extract(predictors, fac[,2:3])

##Change from matrix to df
sdata=as.data.frame(sdata)

## Combine response and predictor variables into one dataframe
sdata2=cbind(fac$Sample,fac$cluster,sdata)
colnames(sdata2)[1] <- "Sample"
colnames(sdata2)[2] <- "Cluster"
#head(sdata2)

## Merge response and predictor variables (df sdata2) with sediment data
phy <- sdata2
#phy <- merge(sdata2,sed_data4,by="Sample")
head(phy)
dim(phy)#1400 14
dim(bio)#1400 432
head(phy)
head(bio)

## Update col name
colnames(bio)[1] <- "Sample"

## Drop the cobbles column
#phy <- within(phy, rm(Cobbles))

## Merge the bio and phy data (already includes factor)
phy_bio <- merge(phy, bio, by = "Sample")  
head(phy_bio)
#View(phy_bio)
names(phy_bio)

## Drop rows with missing data
phy_bio <-  phy_bio[complete.cases(phy_bio), ]

## Change name of df to work with code below
best.data.ss.final <- phy_bio
names(best.data.ss.final)

## Create a df 'bestBIO' for biodiversity data (metrics used for clustering) and save
bestBIO=best.data.ss.final[,c(17:ncol(best.data.ss.final))]
head(bestBIO)
#names(best.data.ss.final)
#write.csv(bestBIO,file = "OUTPUTS/bestBIO.csv",row.names=TRUE)

## Create a df 'bestPHY' for physical variables and save.
bestPHY=best.data.ss.final[,c(3:14)]
head(bestPHY)
dim(bestPHY)#1400 12
#write.csv(bestPHY,file = "OUTPUTS/bestPHY.csv",row.names=TRUE)

## Create a df 'bestFAC' for factor cluster and save
bestFAC=as.data.frame(best.data.ss.final[,2])
head(bestFAC)
colnames(bestFAC)[1] <- 'Cluster'
#write.csv(bestFAC,file = "OUTPUTS/bestFAC.csv",row.names=TRUE)
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: CREATE A DF OF UNI METRICS FOR OVERLAYING IN MDS ORDINATION ####

## Unimetic values stored in this df
head(univmeascoord2_mod)
View(univmeascoord2_mod)

## This df contains data used for ordination
View(best.data.ss.final)

## Take only the samplecodes and cluster
bestFAC2 <- best.data.ss.final[,1:2]
head(bestFAC2)
#View(bestFAC2)
dim(bestFAC2) #1400    2

## Now create df for just the richness data
just_richness_values <- univmeascoord2_mod[which(univmeascoord2_mod$variable=='S'),]
just_richness_values2 <- just_richness_values[,c(1,5)]
colnames(just_richness_values2)[1] <- 'Sample'
head(just_richness_values2)

## Now add in the richness and abundance values
bestFAC3 <- merge(bestFAC2, just_richness_values2, by.x = 'Sample',by.y = 'Sample',all.x = TRUE, all.y = FALSE)
colnames(bestFAC3)[3] <- 'Richness'
head(bestFAC3 )
dim(bestFAC3)#1400    2
#View(bestFAC3)
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: DETECT COLINEARITY IN ENVIRONMENTAL VARIABLES ####

#Using Variance Inflation Factors (VIFs) from usdm package

## Load package
library(usdm)

## Get correlation coefficients for environmental variables
cor(bestPHY, use="all.obs", method="spearman")

## Get VIF scores for variables in df 'bestPHY'and remove highest value above 2.5, one variable at a time
vif(subset(bestPHY))
vif(subset(bestPHY, select = -Wave.velocity))
vif(subset(bestPHY, select = -c(Wave.velocity,Current.speed)))
#vif(subset(bestPHY, select = -c(Phyto, Oxygen,Bathymetry)))

## Create new df for variable with VIF SCORES <2.5
bestPHY2=subset(bestPHY, select = -c(Wave.velocity,Current.speed))
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: CHECK FOR SKEWNESS AND TRANSFORM AS NECESSARY ####

# Check for any skewness in the data. Negative skewness indicates that the mean of the data
# values is less than the median,and the data distribution is left-skewed. Positive skewness
# would indicate that the mean of the data values is larger than the median, and the data
# distribution is right-skewed.

summary(bestPHY2)# Positive skewness (mean>median) for Winter.SPM (3), Gravel (4), Rel..slope.pos. (7), LS.factor (8), Mud (9), Closed.depressions (10) 


## Transform relevant columns  Winter.SPM (3), Gravel (4), Rel..slope.pos. (7), LS.factor (8), Mud (9), Closed.depressions (10) 
bestPHY2[,c(3,4,7,8,9,10)]=log(bestPHY2[,c(3,4,7,8,9,10)]+0.1)
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: BIOENV ####

## Run bioenv with the transformed faunal (df 'bestBIO') and the selected env variables
#(df 'bestPHY2'). Note env variables will be normalised during best.

## Load package
library(vegan) # Load library

## Check dfs same length
dim(bestBIO)# 1400 429
dim(bestPHY2)# 1400  10

## Convert to matrix
bestBIO <- as.matrix(bestBIO)

## Perform BIOENV analysis
res<-bioenv(bestBIO, bestPHY2) 
res

## See all best results
summary(res)

## Output results as a dataframe
res2 <-data.frame(unclass(summary(res)), check.names = FALSE, stringsAsFactors = FALSE)
res3 <- res2[,c(1,3,2)]
class(res3)

## Reduce number of decimal places in correlation column
op = function(x, d=2) sprintf(paste0("%1.",d,"f"), x) 
res3$correlation <- op(res3$correlation, 4)
colnames(res3) <- c('Size','Variables','Correlation (ρ)')

## Create results table
library(knitr)
library(kableExtra)
library(magrittr)
library(webshot)
library(magick)
webshot::install_phantomjs()
kable(res3[1:5,], caption = "Table 3. Results of a 'best' analysis identifying the subset of environmental variables which are most correlated with the biodiversity data.")%>%
  kable_styling()%>%
  row_spec(5,bold=T,hline_after = T)%>%
  save_kable(file = "C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/Table_3.png")
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: ADONIS ####

# Use 'adonis' to quantify the variation in biodiversity data explained by the environmental predictors identified in 'best'

## Normalise env data prior to using adonis (subtract by pop mean and sd) of phy variable data 
#varephy3 <- scale(phy_dat2)# scale the transformed data

## varephy3 is a matrix so need to convert back to a df for use with adonis
#varephy3=as.data.frame(varephy3)

## Use Adonis to see how much of the variation is explained by the different variables. Enter
# phy variable which are important (see BEST results)
# http://www.talkstats.com/showthread.php/15912-Permanova-Adonis
adonis.res=adonis2(formula = bestBIO ~ Salinity.mean + Winter.SPM + Gravel + Rel..slope.pos. + Mud, data = bestPHY2, permutations = 999, method = "bray", by="terms")
adonis.res

#_______________________________________________________________________________
#### EXPLAINING PATTERNS: dbRDA ORDINATION WITH CLUSTER GROUP OVERLAY ####

# Use distance based redundancy analysis (dbRDA) ordination to visualise the  
# relationship between macrofaunal data and predictor variables 
#https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/capscale

# Ordination
vare.cap <- capscale(bestBIO ~ Salinity.mean + Winter.SPM + Gravel + Rel..slope.pos. + Mud, bestPHY2,dist="bray")#Gravel + WOV
vare.cap

## Set colours
palette(c("blue2","cyan1","#05aae1","purple","#b40202","red1","darkorange","yellow","#b4b404","green3","palegreen1"))

## Update column names for ordintion
colnames(bestPHY2)[colnames(bestPHY2) == "Salinity.mean"] <- "Salinity mean"
colnames(bestPHY2)[colnames(bestPHY2) == "Winter.SPM"] <- "Winter SPM"
colnames(bestPHY2)[colnames(bestPHY2) == "Rel..slope.pos."] <- "Rel. slope pos."

## Save plot to an image file (png or tiff)
tiff("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/Figure_6.tif",width = 18,height = 15,units = "cm", res = 800, pointsize = 12)
#tiff("OUTPUTS/FIGURE 6.tiff",width = 14,height = 12.5,units = "cm",res = 800,pointsize = 12)

## Plot empty ordination
plot(vare.cap,type="none",scaling=2,xlim = c(-2, 2),ylim = c(-4, 4),
     cex.lab=0.6,cex.axis=0.6)

## Fit environmental variables
#fit <- envfit(vare.cap, bestPHY2[,c("Gravel","Salinity","SPM","WOV")])
fit <- envfit(vare.cap, bestPHY2[,c('Salinity mean','Winter SPM','Gravel','Rel. slope pos.','Mud')])

# 50% transparency
my_colors <- c("blue2","cyan1","#05aae1","purple","#b40202","red1",
               "darkorange","yellow","#b4b404","green3","palegreen1")
transparent_colors <- adjustcolor(my_colors, alpha.f = 0.6)

## Add sites as points coloured by cluster group. Samples shown as points representing their mean species composition.
points(vare.cap, col=transparent_colors[bestFAC$Cluster],cex = 0.65,pch=16)

# Add environmental vectors
# Add explantory variable vectors pointing in the direction of max increase.  Long vectors equals strong trends
plot(fit, col="black", p.max=0.1, cex=0.8)   # p.max keeps only significant variables

## Add legend
legend(4.5,4, legend=c('Epi_A1a', 'Epi_B1a', 'Epi_B1b', 'Epi_C1', 'Epi_D1', 'Epi_D2a', 'Epi_D2b',  'Epi_D2c', 'Epi_D2d', 'Epi_E1a', 'Epi_E1b'),pch = 20,pt.cex=0.85,
       col=c("blue2","cyan1","#05aae1","purple","#b40202","red1","darkorange","yellow","#b4b404","green3","palegreen1") ,       box.lty=0,cex=0.6)

dev.off()

#anova(vare.cap)
#_______________________________________________________________________________
#### EXPLAINING PATTERNS: dbRDA ORDINATION WITH RICHNESS OVERLAY ####

# Optionally set colours using RColorBrewer
library(RColorBrewer)

#Create a function to generate a continuous color palette
#rbPal <- colorRampPalette(c('red','blue'))
rbPal <-colorRampPalette( matlab.like(50))

#This adds a column of color values
# based on the y values
bestFAC3$Col <- rbPal(10)[as.numeric(cut(bestFAC3$Richness,breaks = 10))]

legend_image <- as.raster(matrix(rev(matlab.like(50)), ncol=1))

## Save plot to an image file (png or tiff)
png("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/OneBenthicTrawl/OUTPUTS/ordination_richness.png",width = 18,height = 15,units = "cm", res = 800, pointsize = 12)
#tiff("OUTPUTS/FIGURE 6.tiff",width = 14,height = 12.5,units = "cm",res = 800,pointsize = 12)


plot(vare.cap,type="none",scaling=2,xlim = c(-2, 2),ylim = c(-4, 4),
     cex.lab=0.6,cex.axis=0.6)

## Add sites as points coloured by cluster group. Samples shown as points representing their
# mean species composition.
points(vare.cap, col=bestFAC3$Col,cex = 0.65,pch=16)


# Add explantory variable vectors pointing in the direction of max increase.  Long vectors 
#equals strong trends
text(vare.cap, dis="cn",col="black",cex=0.7)


## Legend position bottom right and top left
rasterImage(legend_image, 5.5, -4, 5,-2)

## Legend text x=dist along x-axis, 
text(x=5.8, y = seq(-4,-2,l=5), labels = seq(min(bestFAC3$Richness),max(bestFAC3$Richness),l=5),,col="black",cex=0.7)
dev.off()

#_______________________________________________________________________________
#### SUPPLEMENTARY INFORMATION: IMPORTANT VARIABLES PLOTS (QUANTILES) ####

#https://stackoverflow.com/questions/49910270/r-plot-raster-colorscheme-not-full-range

## Load packages
library(raster)
library(rasterVis)
library(classInt)
library(ggpubr)
library(RColorBrewer)

## Salinity mean
plotVar<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/so_mean.tif")# SALINITY MEAN
nColor <- 50
break1 <- classIntervals(plotVar[!is.na(plotVar)], n = nColor, style = "quantile")
lvp1 <- levelplot(plotVar, 
                 col.regions = colorRampPalette(rev(brewer.pal(9, 'RdYlGn'))), 
                 at = break1$brks, margin = FALSE, main="Salinity mean",colorkey = FALSE)

## Winter SMP
plotVar<-raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/SPM_WINTER.tif")#WINTER SPM
nColor <- 50
break1 <- classIntervals(plotVar[!is.na(plotVar)], n = nColor, style = "quantile")
lvp2 <- levelplot(plotVar, 
                 col.regions = colorRampPalette(rev(brewer.pal(9, 'RdYlGn'))), 
                 at = break1$brks, margin = FALSE, main="Winter SPM",colorkey = FALSE)

## Gravel plot
plotVar<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Gravel_Fraction.tif")#GRAVEL
nColor <- 50
break1 <- classIntervals(plotVar[!is.na(plotVar)], n = nColor, style = "quantile")
lvp3 <- levelplot(plotVar, 
                 col.regions = colorRampPalette(rev(brewer.pal(9, 'RdYlGn'))), 
                 at = break1$brks, margin = FALSE,main="Gravel",colorkey = FALSE)

## Rel. slope pos.
plotVar<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/RSP.tif")#RELATIVE SLOPE POSITION
nColor <- 50
break1 <- classIntervals(plotVar[!is.na(plotVar)], n = nColor, style = "quantile")
lvp4 <- levelplot(plotVar, 
                  col.regions = colorRampPalette(rev(brewer.pal(9, 'RdYlGn'))), 
                  at = break1$brks, margin = FALSE,main="Rel. slope pos.",colorkey = FALSE)

## Mud
plotVar<- raster("C:/Users/kmc00/OneDrive - CEFAS/R_PROJECTS/ModellingData/RasterPredictorsDec2024/Predicted_Mud_Fraction.tif")#MUD
nColor <- 50
break1 <- classIntervals(plotVar[!is.na(plotVar)], n = nColor, style = "quantile")
lvp5 <- levelplot(plotVar, 
                 col.regions = colorRampPalette(rev(brewer.pal(9, 'RdYlGn'))), 
                 at = break1$brks, margin = FALSE, main="Mud",colorkey = FALSE)

## Arrange plots together
import_var <- ggpubr::ggarrange(lvp1,lvp2,lvp3,lvp4,lvp5)
import_var

## Save combined plot
ggsave(plot = import_var,
       filename = paste0("C:\\Users\\KMC00\\OneDrive - CEFAS\\R_PROJECTS\\OneBenthicTrawl\\OUTPUTS\\Figure_S4.tif"),
       height = 250, width =350,units = "mm", dpi = 500,
       device = "tiff",limitsize = FALSE,bg="white")
