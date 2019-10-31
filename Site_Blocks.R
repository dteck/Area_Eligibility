library(httr)
library(stringr)

#ask users to set working directory
setwd(choose.dir())
workdir<-getwd()
#------

#--Custom function to call census geocoder
GPS_to_GEOID<- function(x,y){
  library(httr)
  library(stringr)
  x_san<- paste(str_trim(as.character(x), side = "both"))
  y_san<- paste(str_trim(as.character(y), side = "both"))
  URL<-paste("https://geocoding.geo.census.gov/geocoder/geographies/coordinates?x=", x_san, "&y=", y_san, "&benchmark=4&vintage=4&format=JSON", sep="") #build URL
  call<-try(GET(URL),silent = TRUE) #try to make call
  cont<-try(content(call), silent = TRUE) #parse returned data
  if(class(cont)=="try-error"){
    return(as.character(0))
  } else{
    return(as.character(cont$result$geographies$`2010 Census Blocks`[[1]]$GEOID))
  }
}
#------

#Load site lists and drop extra data
SFSP_2016 <-read_csv("Summer_Meal_Sites_2016.csv")
SFSP_2016 <-SFSP_2016[,c('X','Y','FNSID')]
  
SFSP_2017 <-read_csv("Summer_Meal_Sites_2017.csv")
SFSP_2017 <-SFSP_2017[,c('X','Y','FNSID')]

SFSP_2018 <-read_csv("Summer_Meal_Sites_2018.csv")
SFSP_2018 <-SFSP_2018[,c('X','Y','FNSID')]

SFSP_2019 <-read_csv("Summer_Meal_Sites_2019.csv")
SFSP_2019 <-SFSP_2019[,c('X','Y','FNSID')]
#------

#Merge all site lists, keep unique
Sites<-rbind.data.frame(SFSP_2016,SFSP_2017)
Sites<-rbind.data.frame(Sites,SFSP_2018)
Sites<-rbind.data.frame(Sites,SFSP_2019)

Sites<-Sites[!duplicated(Sites$FNSID), ]
#------

#Resolve the list of sites to 2010 census blocks
for(i in seq(1,length(Sites$X),1)){
  Sites$Block[i]<-GPS_to_GEOID(Sites$X[i],Sites$Y[i])
  print(paste(i,"of", length(Sites$X), "Result:",Sites$Block[i]))
}
#------

write.csv(Sites, "Sites.csv", row.names = FALSE)

















