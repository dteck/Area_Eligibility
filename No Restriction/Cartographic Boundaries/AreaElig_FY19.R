############################
# Author: Mark Richards    #
# Time to run: 29 hours    #
############################

#load libraries--
library(tigris)
library(leaflet)
library(leaflet.extras)
library(rgeos)
library(readxl)
library(RColorBrewer)
library(tidycensus)
library(rgdal)
library(rgeos)
library(htmltools)
library(spatstat)
library(spdep)
library(arrangements)
#------
start_time <- Sys.time() #start a timer
setwd('C:\\Users\\darkteck04\\Desktop\\AreaElig\\FY19')
#Download the Fiscal Year 2019 spreadsheet from the USDA to the working dir, then load it in.
download.file('http://data-cacfp-sfsp.opendata.arcgis.com/datasets/7a39c0ae571149c5a09c5790f0ff99f9_0.csv', method = 'auto',
              destfile = "./FY19.csv", quiet = TRUE)
FY19<-read.csv('./FY19.csv')
#------

#Separate CSV data into states
#Use sequence up to 56 to get states but not teritories and sotre as FIPS
StateData<-list() #list to hold data frames of data
for (i in seq(1,56,1)){ #subset CSV into state specific data
  StateData[[i]]<-subset.data.frame(FY19,FY19$STATEFP == i)
}
rm(FY19)
#------

#Create list of state outlines
StatesSHP<-states(cb=TRUE) #State shapefiles
StateOutline<-list() #list to hold state oulines
for(i in seq(1,56,1)){ #subset state shapes to indivigual states
  print(i)
  StateOutline[[i]]<-subset.data.frame(StatesSHP, as.numeric(StatesSHP$STATEFP)==i)
}
rm(StatesSHP) #drop unified state shapefile
#------

#create list of counties in each state
CountyOutline<-list() #create list to hold county shapes
for (i in seq(1,56,1)){ #pull county shapes for each state by FIPS code
  if (i %in% c(03,07,14,43,52)){ # list of non states if pulled it grabs *ALL* data, so exclude (159mb down to 29mb)
    print(paste('Do not pull County:', i))
    CountyOutline[[i]]<-0
  } else {
    print(i)
    CountyOutline[[i]]<-counties(i,cb=TRUE)
  }
}
#------

#Create list of block groups in each state
BGOutline<-list() #create list to hold Block Group shapes
for (i in seq(1,56,1)){ #pull county shapes for each state by FIPS code
  if (i %in% c(03,07,14,43,52)){ # list of non states, exclude
    print(paste('Do not pull BG:', i))
    CountyOutline[[i]]<-0
  } else {
    print(i)
    BGOutline[[i]]<-block_groups(i,cb=TRUE)
  }
}
#------

#merge state level data with block group shapes
BGData<-list()
for (i in seq(1,56,1)){
  a<-data.frame(StateData[[i]]) #pull out state data
  b<-BGOutline[[i]] #pull put state shape
  if (i %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", i))
  } else{
    if(nchar(b@data[["GEOID"]][1]) != nchar(a$GEOID[1])){ #check to see if FIPS lengths match
      print(paste("Adjust GEOID:", i))
      b@data[["GEOID"]]<-as.character(as.numeric(b@data[["GEOID"]])) #if not convert to numeric/character to drop zeros
    } else {
      print(i) #if yes continue
    }
    BGData[[i]]<-geo_join(b, a, by_sp='GEOID', by_df='GEOID', how = 'left') #join data and shape
  }
}
rm(a,b, StateData, BGOutline)
#------

#Build Lists of neighbors for each blockgroup in a state
Neighbors<-list()
for (i in seq(1,56,1)){
  if (i %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", i))
  } else{
    print(paste("Finding Neighbors for:", i))
    Neighbors[[i]]<-poly2nb(BGData[[i]], queen = TRUE)
  }
}
#------

#function to calc combos for under 18
Combo<- function(StateFIPS){ #function to calcuate combinatiosn for under 18
  dfcombo<-data.frame() #build df to hold results
  for(i in seq(1,length(BGData[[StateFIPS]]),1)){ #run through all group blocks
    if (BGData[[StateFIPS]]$ELIGFY19[i]=="No"){
      #print(paste(i, length(neighbor[[i]])))
      if (length(Neighbors[[StateFIPS]][[i]]) >1){ #if the # of neighboring BGs is more than 1
        comb<-arrangements::combinations(Neighbors[[StateFIPS]][[i]],2, replace = FALSE)
        for (n in seq(1,length(comb[,1]),1)){
          dfcomb2<-data.frame(
            GEOID=BGData[[StateFIPS]]$GEOID[i], #bg in question id
            ID18=BGData[[StateFIPS]]$Num18BG[i], #bg in question Identified under 18
            Tot18=BGData[[StateFIPS]]$TotPovUniv[i], #bg in question Total under 18
            ID12=BGData[[StateFIPS]]$Num12BG[i], #bg in question Identified under 12
            Tot12=BGData[[StateFIPS]]$TotPovUn_1[i], #bg in question Total under 12
            firstGEOID18=BGData[[StateFIPS]]$GEOID[comb[n,1]], #first calc block id
            B1_ID18=BGData[[StateFIPS]]$Num18BG[comb[n,1]], #first calc block Identified under 18
            B1_Tot18=BGData[[StateFIPS]]$TotPovUniv[comb[n,1]], #first calc block Total under 18
            B1_ID12=BGData[[StateFIPS]]$Num12BG[comb[n,1]], #first calc block Identified under 12
            B1_Tot12=BGData[[StateFIPS]]$TotPovUn_1[comb[n,1]],#first calc block Total under 12
            secondGEOID18=BGData[[StateFIPS]]$GEOID[comb[n,2]], #second calc block id
            B2_ID18=BGData[[StateFIPS]]$Num18BG[comb[n,2]], #second calc block Identified under 18
            B2_Tot18=BGData[[StateFIPS]]$TotPovUniv[comb[n,2]], #second calc block Total under 18
            B2_ID12=BGData[[StateFIPS]]$Num12BG[comb[n,2]], #second calc block Identified under 12
            B2_Tot12=BGData[[StateFIPS]]$TotPovUn_1[comb[n,2]])#second calc block Total under 12
          dfcomb2$weightnum18<-dfcomb2$ID18+dfcomb2$B1_ID18+dfcomb2$B2_ID18 #weight identified under 18
          dfcomb2$weightpov18<-dfcomb2$Tot18+dfcomb2$B1_Tot18+dfcomb2$B2_Tot18 #weight total under 18
          dfcomb2$weightPerc18<-(dfcomb2$weightnum18/dfcomb2$weightpov18)*100 #weight percent under 18
          dfcomb2$weightnum12<-dfcomb2$ID12+dfcomb2$B1_ID12+dfcomb2$B2_ID12 #weight identified under 12
          dfcomb2$weightpov12<-dfcomb2$Tot12+dfcomb2$B1_Tot12+dfcomb2$B2_Tot12 #weight total under 12
          dfcomb2$weightPerc12<-(dfcomb2$weightnum12/dfcomb2$weightpov12)*100 # weight percent under 12
          dfcombo<-rbind(dfcombo, dfcomb2)
          rm(dfcomb2)
        }
      } else { #if the # of neighboring BGs is not greater than 1
        if (Neighbors[[StateFIPS]][[i]][1]==0){
          dfcomb2<-data.frame(GEOID=BGData[[StateFIPS]]$GEOID[i],
                              ID18=BGData[[StateFIPS]]$Num18BG[i],
                              Tot18=BGData[[StateFIPS]]$TotPovUniv[i],
                              ID12=BGData[[StateFIPS]]$Num12BG[i],
                              Tot12=BGData[[StateFIPS]]$TotPovUn_1[i],
                              firstGEOID18=as.factor(0),
                              B1_ID18=0,
                              B1_Tot18=0,
                              B1_ID12=0,
                              B1_Tot12=0,
                              secondGEOID18=as.factor(0),
                              B2_ID18=0,
                              B2_Tot18=0,
                              B2_ID12=0,
                              B2_Tot12=0,
                              weightnum18=0,
                              weightpov18=0,
                              weightPerc18=0,
                              weightnum12=0,
                              weightpov12=0,
                              weightPerc12=0)
          dfcombo<-rbind(dfcombo, dfcomb2)
          rm(dfcomb2)
        } else {
          dfcomb2<-data.frame(
            GEOID=BGData[[StateFIPS]]$GEOID[i],
            ID18=BGData[[StateFIPS]]$Num18BG[i],
            Tot18=BGData[[StateFIPS]]$TotPovUniv[i],
            ID12=BGData[[StateFIPS]]$Num12BG[i],
            Tot12=BGData[[StateFIPS]]$TotPovUn_1[i],
            firstGEOID18=BGData[[StateFIPS]]$GEOID[Neighbors[[StateFIPS]][[i]][1]], #neighbor[[i]][1]
            B1_ID18=BGData[[StateFIPS]]$Num18BG[Neighbors[[StateFIPS]][[i]][1]], 
            B1_Tot18=BGData[[StateFIPS]]$TotPovUniv[Neighbors[[StateFIPS]][[i]][1]], 
            B1_ID12=BGData[[StateFIPS]]$Num12BG[Neighbors[[StateFIPS]][[i]][1]],
            B1_Tot12=BGData[[StateFIPS]]$TotPovUn_1[Neighbors[[StateFIPS]][[i]][1]],
            secondGEOID18=as.factor(0),
            B2_ID18=0,
            B2_Tot18=0,
            B2_ID12=0,
            B2_Tot12=0)
          dfcomb2$weightnum18<-dfcomb2$ID18+dfcomb2$B1_ID18 
          dfcomb2$weightpov18<-dfcomb2$Tot18+dfcomb2$B1_Tot18
          dfcomb2$weightPerc18<-(dfcomb2$weightnum18/dfcomb2$weightpov18)*100 
          dfcomb2$weightnum12<-dfcomb2$ID12+dfcomb2$B1_ID12
          dfcomb2$weightpov12<-dfcomb2$Tot12+dfcomb2$B1_Tot12
          dfcomb2$weightPerc12<-(dfcomb2$weightnum12/dfcomb2$weightpov12)*100
          dfcombo<-rbind(dfcombo, dfcomb2)
          rm(dfcomb2)
        }
      }
    }
    print(paste('FIPS:',StateFIPS,'Block Group:', i, "of", length(BGData[[StateFIPS]])))
  }
  return(dfcombo)
}


#calcuate all of the weighted averages, both under 18 and under 12
###############################################################################
############################WARNING COMPUTATION HEAVY##########################
###############################################################################
BGDataAll<<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste("Calculating Averages for:", r))
    BGDataAll[[r]]<-Combo(r)
  }
}
saveRDS(BGDataAll, "BGDataAll.rds") #save output to prevent having to rerun on error

#keep only one result above 50%
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste('Keeping one result per GEOID for:',r))
    BGDataAll[[r]]<-subset.data.frame(BGDataAll[[r]], BGDataAll[[r]]$weightPerc18 >= 50 | BGDataAll[[r]]$weightPerc12 >=50)
    BGDataAll[[r]]<-BGDataAll[[r]][!duplicated(BGDataAll[[r]]$GEOID),]
    }
}


#merge calualtions to shapefiles
BGCombined<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste("Joining data for:", r))
    BGCombined[[r]]<-geo_join(BGData[[r]],BGDataAll[[r]], by_sp='GEOID', by_df='GEOID', how='left')
  }
}
#------

#convert NA's to zeros for use in operator subsetting
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste('cleaning data for:', r))
    BGCombined[[r]]$weightnum18[is.na(BGCombined[[r]]$weightnum18)]<-0
    BGCombined[[r]]$weightpov18[is.na(BGCombined[[r]]$weightpov18)]<-0
    BGCombined[[r]]$weightPerc18[is.na(BGCombined[[r]]$weightPerc18)]<-0
    BGCombined[[r]]$weightnum12[is.na(BGCombined[[r]]$weightnum12)]<-0
    BGCombined[[r]]$weightpov12[is.na(BGCombined[[r]]$weightpov12)]<-0
    BGCombined[[r]]$weightPerc12[is.na(BGCombined[[r]]$weightPerc12)]<-0
  }
}
#------  

#separeate the results into yes/no/calc so they can be mapped separately
BGSubsets<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    YNC<-list()
    print(paste('Subsetting Yes/No/Calc for:', r))
    YNC[[1]]<-subset.data.frame(BGCombined[[r]],BGCombined[[r]]$ELIGFY19=="Yes")
    YNC[[2]]<-subset.data.frame(BGCombined[[r]],BGCombined[[r]]$ELIGFY19=="No" & BGCombined[[r]]$weightPerc18 < 50 & BGCombined[[r]]$weightPerc12 < 50)
    YNC[[3]]<-subset.data.frame(BGCombined[[r]], BGCombined[[r]]$ELIGFY19=="No" & BGCombined[[r]]$weightPerc18 >=50 | BGCombined[[r]]$weightPerc12 >=50)
    BGSubsets[[r]]<-YNC
    rm(YNC)
  }
}
#------

#Build leaflets for each state
BGLeaflets<-list()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste('Building Leaflet for:', r))
    BGLeaflets[[r]] <- leaflet() %>% enableTileCaching() %>%
      addProviderTiles(providers$OpenStreetMap.BlackAndWhite, group = "Grey")
    
    BGLeaflets[[r]]<- BGLeaflets[[r]] %>% addPolygons(data=BGSubsets[[r]][[1]],
                                                      weight = 1,
                                                      fill = TRUE,
                                                      fillOpacity = 0.3,
                                                      fillColor = 'red',
                                                      stroke = TRUE,
                                                      color = "black",
                                                      group = 'plotT',
                                                      label = BGSubsets[[r]][[1]]$GEOID, #paste("FY19Elig:", BGSubsets[[r]][[1]]$ELIGFY19),
                                                      labelOptions = labelOptions(textsize = "16px"),
                                                      highlight = highlightOptions(
                                                        weight = 5,
                                                        color = "red",
                                                        fillOpacity = 0.7,
                                                        bringToFront = TRUE),
                                                      popup = paste('<b>FY19Elig</b>: Eligable<br>',
                                                                    '<b>GEOID:</b>', BGSubsets[[r]][[1]]$GEOID, '<br>',
                                                                    '<b>Identified Under 18:</b>', BGSubsets[[r]][[1]]$Num18BG, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[1]]$TotPovUniv, '<br>',
                                                                    '<b>Identified Under 12:</b>', BGSubsets[[r]][[1]]$Num12BG, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[1]]$TotPovUn_1, '<br>'))
    BGLeaflets[[r]]<- BGLeaflets[[r]] %>% addPolygons(data=BGSubsets[[r]][[2]],
                                                      weight = 1,
                                                      fill = TRUE,
                                                      fillOpacity = 0.3,
                                                      fillColor = 'blue',
                                                      stroke = TRUE,
                                                      color = "black",
                                                      group = 'plotT',
                                                      label = BGSubsets[[r]][[2]]$GEOID, #paste("FY19Elig:", BGSubsets[[r]][[2]]$ELIGFY19),
                                                      labelOptions = labelOptions(textsize = "16px"),
                                                      highlight = highlightOptions(
                                                        weight = 5,
                                                        color = "red",
                                                        fillOpacity = 0.7,
                                                        bringToFront = TRUE),
                                                      popup = paste('<b>FY19Elig</b>: Not Eligable<br>',
                                                                    '<b>GEOID:</b>', BGSubsets[[r]][[2]]$GEOID, '<br>',
                                                                    '<b>Identified Under 18:</b>', BGSubsets[[r]][[2]]$Num18BG, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[2]]$TotPovUniv, '<br>',
                                                                    '<b>Identified Under 12:</b>', BGSubsets[[r]][[2]]$Num12BG, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[2]]$TotPovUn_1, '<br>'))
    BGLeaflets[[r]]<- BGLeaflets[[r]] %>% addPolygons(data=BGSubsets[[r]][[3]],
                                                      weight = 1,
                                                      fill = TRUE,
                                                      fillOpacity = 0.3,
                                                      fillColor = 'yellow',
                                                      stroke = TRUE,
                                                      color = "black",
                                                      group = 'plotT',
                                                      label = BGSubsets[[r]][[3]]$GEOID, #paste("FY19Elig: Calc Eligable"),
                                                      labelOptions = labelOptions(textsize = "16px"),
                                                      highlight = highlightOptions(
                                                        weight = 5,
                                                        color = "red",
                                                        fillOpacity = 0.7,
                                                        bringToFront = TRUE),
                                                      popup = paste('<b>FY19Elig</b>: Calculation Eligable<br>',
                                                                    '<b>GEOID:</b>', BGSubsets[[r]][[3]]$GEOID, '<br>',
                                                                    '<b>Identified Under 18:</b>', BGSubsets[[r]][[3]]$Num18BG, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[3]]$TotPovUniv, '<br>',
                                                                    '<b>Identified Under 12:</b>', BGSubsets[[r]][[3]]$Num12BG, '<br>',
                                                                    '<b>Total Under 18:</b>', BGSubsets[[r]][[3]]$TotPovUn_1, '<br><hr>',
                                                                    '<b><center>Calculation</center></b><hr>',
                                                                    '<table border="1">
                                                                      <tr>
                                                                      <th style="background-color:#c5d9d5;" align="center">GEOID</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">Under 18</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">18 pop</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">Under 12</th>
                                                                      <th style="background-color:#c5d9d5;" align="center">12 pop</th>
                                                                      </tr>
                                                                      <tbody>
                                                                      <tr>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$firstGEOID18,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_ID18,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_Tot18,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_ID12,'</td>
                                                                      <td align="center">',BGSubsets[[r]][[3]]$B1_Tot12,'</td>
                                                                      </tr>
                                                                      <tr>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$secondGEOID18, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_ID18, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_Tot18, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_ID12, '</td>
                                                                      <td align="center">', BGSubsets[[r]][[3]]$B2_Tot12, '</td>
                                                                      </tr>
                                                                      <tr>
                                                                      <td style="background-color:#f7daad;" align="center">Totals</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$Num18BG+BGSubsets[[r]][[3]]$B1_ID18+BGSubsets[[r]][[3]]$B2_ID18, '</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$TotPovUniv+BGSubsets[[r]][[3]]$B1_Tot18+BGSubsets[[r]][[3]]$B2_Tot18, '</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$Num12BG+BGSubsets[[r]][[3]]$B1_ID12+BGSubsets[[r]][[3]]$B2_ID12, '</td>
                                                                      <td style="background-color:#f7daad;" align="center">', BGSubsets[[r]][[3]]$TotPovUn_1+BGSubsets[[r]][[3]]$B1_Tot12+BGSubsets[[r]][[3]]$B2_Tot12, '</td>
                                                                      </tr>
                                                                      </tbody>
                                                                      </table>',
                                                                    '<b>Weighted % under 18:</b>', round(BGSubsets[[r]][[3]]$weightPerc18,2), '<br>',
                                                                    '<b>Weighted % under 12:</b>', round(BGSubsets[[r]][[3]]$weightPerc12,2), '<br>'))
    BGLeaflets[[r]] <- BGLeaflets[[r]] %>% addLegend(position = "bottomright", 
                                                     labels =  c("USDA Eligible","Calculated Eligible", "Ineligible"), 
                                                     colors =  c('red','yellow','blue'),
                                                     opacity = 0.5,
                                                     group = 'BlockG')
    
    BGLeaflets[[r]] <- BGLeaflets[[r]] %>% addSearchOSM(options = searchOptions(autoCollapse = TRUE, minLength = 2))
  }
}
#------
saveRDS(BGLeaflets, "BGLeaflets.rds") #save leaflet maps to prevent having to rerun

#export all data into folders
SumAll<-data.frame()
for (r in seq(1,56,1)){
  if (r %in% c(03,07,14,43,52)){
    print(paste("Skip non-State", r))
  } else{
    print(paste("Saving Map for:", r))
    mainDir <- getwd()
    subDir <- paste(mainDir,"/",gsub(" ","_",StateOutline[[r]]$NAME),sep="")
    if (file.exists(subDir)){
      setwd(file.path(subDir))
    } else {
      dir.create(file.path(subDir))
      setwd(file.path(subDir))
      
    }
    name<-paste(StateOutline[[r]]$NAME,"_FY19.html", sep="")
    name<-gsub(" ","_",name)
    htmlwidgets::saveWidget(BGLeaflets[[r]], name, selfcontained = TRUE)
    write.csv(BGCombined[[r]],paste(StateOutline[[r]]$NAME,".csv", sep = ""), row.names = FALSE)
    saveRDS(BGSubsets[[r]][[1]], "Elig_Yes.rds")
    saveRDS(BGSubsets[[r]][[2]], "Elig_No.rds")
    saveRDS(BGSubsets[[r]][[3]], "Elig_Calc.rds")
    Summary<-data.frame(
      Name=StateOutline[[r]]$NAME, #name of area
      AreasYes=length(BGSubsets[[r]][[1]]$GEOID), # number of USDA areas
      Under18Yes=sum(BGSubsets[[r]][[1]]$Num18BG), # Number of USDA identified under 18
      Under12Yes=sum(BGSubsets[[r]][[1]]$Num12BG), # Number of USDA identified under 12
      
      AreasCalc=length(BGSubsets[[r]][[3]]$GEOID), # number of calc eligable areas
      Under18Calc=sum(BGSubsets[[r]][[3]]$Num18BG), # Number of calc identified under 18
      Under12Calc=sum(BGSubsets[[r]][[3]]$Num12BG), # Number of calc identified under 12
      
      AreasDeltaPerc=(((length(BGSubsets[[r]][[1]]$GEOID)+length(BGSubsets[[r]][[3]]$GEOID))-length(BGSubsets[[r]][[1]]$GEOID))/length(BGSubsets[[r]][[1]]$GEOID))*100, # %change num of areas
      Under18DeltaPerc=(((sum(BGSubsets[[r]][[1]]$Num18BG)+sum(BGSubsets[[r]][[3]]$Num18BG))-sum(BGSubsets[[r]][[1]]$Num18BG))/sum(BGSubsets[[r]][[1]]$Num18BG))*100, # %change Ident under 18
      Under12DeltaPerc=(((sum(BGSubsets[[r]][[1]]$Num12BG)+sum(BGSubsets[[r]][[3]]$Num12BG))-sum(BGSubsets[[r]][[1]]$Num12BG))/sum(BGSubsets[[r]][[1]]$Num12BG))*100, # %change Ident under 12
      
      AreasNo=length(BGSubsets[[r]][[2]]$GEOID), #num of not eligable
      Under18No=sum(BGSubsets[[r]][[2]]$Num18BG) #num of under 18 still not served
    )
    SumAll<-rbind.data.frame(SumAll,Summary)
    write.csv(Summary,paste(StateOutline[[r]]$NAME,"_Summary.csv", sep = ""), row.names = FALSE)
    setwd(mainDir)
  }
}
#------
write.csv(SumAll, "Overall_Summary.csv", row.names = FALSE)
End_time <- Sys.time()
Run_Time<-End_time - start_time