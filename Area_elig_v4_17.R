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
library(tidyverse)
#------

#ask users to set working directory
setwd(choose.dir())
workdir<-getwd()
#------

#start Stopwatch
Start<-Sys.time()
#------

#Check to see if state data separation has already been done
if(file.exists(paste(workdir,'/Queue/Queue.rds',sep = ""))==TRUE){
  print("States already separated: loading Queue")
  Queue<-readRDS(paste(workdir,"/Queue/Queue.rds", sep = ""))
} else {
  print("Separating States")
  #read USDA area eligability file
  shapeData <- readOGR(paste(workdir,"/USDA/",list.files(paste(workdir,"/USDA", sep = ""), pattern = ".shp"),sep = "")) #read shapefile
  #------
  #Separate data into states
  StateData<-list() #create emplyt list to hold states
  for (i in seq(1,56,1)){ #separate by FIPS code
    if (i %in% c(3,7,14,43,52)){
      print(paste('Skip non-state:', i))
    } else {
      print(i)
      StateData[[i]]<-subset.data.frame(shapeData,shapeData$STATEFP == str_pad(paste(as.character(i)),2,side='left',pad = "0"))
      if(dir.exists(paste(workdir,'/Queue', sep = ""))){
        #print("dir exists")
      } else{
        print("creating queue")
        dir.create(paste(workdir,'/Queue', sep = ""))
      }
      saveRDS(StateData[[i]], paste(workdir,'/Queue/', StateData[[i]]$STATEFP[1],".rds",sep = "")) #save as indivigal states
    }
  }
  #saveRDS(StateData, paste(workdir,"/USDA/StateData.rds", sep = "")) #save RDS to USDA folder
  rm(shapeData)
  Queue<-list() #create a queue, with name, FIPS, and file location
  for (i in seq(1,56,1)){ #separate by FIPS code
    if (i %in% c(3,7,14,43,52)){
      print(paste('Skip non-state:', i))
    } else {
      print(paste('Queueing:', i))
      st<-data.frame(
        FIPS=str_pad(paste(as.character(i)),2,side='left',pad = "0"),
        Name=gsub(" ", x=as.character(unique(unlist(StateData[[i]]$State))), replacement = "_"),
        File=paste(workdir,'/Queue/', StateData[[i]]$STATEFP[1],".rds",sep = ""),
        Done=0
      )
      Queue[[i]]<-st
      rm(st)
    }
  }
  saveRDS(Queue, paste(workdir,'/Queue/Queue.rds',sep = ""))
  rm(StateData)
}
#------

calculations<-function(FIPS, state, Threshold){
  #generate a list of touching group blocks 
  neighbor<-poly2nb(state, queen = TRUE)
  #------
  
  #Drop polygons and keep data
  state<-data.frame(state)
  #------
  
  #Embed neighbor index values into Dataframe
  state$neighbors<-neighbor
  #------
  
  #set threshold value
  #Threshold=40
  #------
  
  #added geoid column
  state$GEOID<-state$GEOID_1
  #------
  
  #Create subsets of the dataframe to keep only items that are not eligable and meet the threshold
  state_Thresh<-subset.data.frame(state, state$Pct18BG>=Threshold | state$Pct12BG >=Threshold)
  state_Thresh_no<-subset.data.frame(state_Thresh, state_Thresh$ELIGFY17=="NO")
  #------
  
  #Drop neighbors that do not meet the threshold - keep as index values,  to change to GEOIDS use FY19_Thresh_no$neighborsGEOID[[i]]
  for(i in seq(1, length(state_Thresh_no$neighbors),1)){
    print(i)
    state_Thresh_no$neighborsThresh[i]<-list(state$GEOID[state_Thresh_no$neighbors[[i]]] %in% state_Thresh$GEOID_1)
    if(length(subset(state_Thresh_no$neighbors[[i]], state$GEOID[state_Thresh_no$neighbors[[i]]] %in% state_Thresh$GEOID_1))<1){
      state_Thresh_no$neighborsThresh[i]<-0
    } else{
      state_Thresh_no$neighborsThresh[i]<-list(subset(state_Thresh_no$neighbors[[i]], state$GEOID[state_Thresh_no$neighbors[[i]]] %in% state_Thresh$GEOID_1))
    }
  }
  #------
  
  #Create count of neighbors
  for(i in seq(1, length(state_Thresh_no$neighbors),1)){
    if(state_Thresh_no$neighborsThresh[[i]][1]==0){
      print(0)
      state_Thresh_no$Ncount[i]<-0
    } else{
      print(length(state_Thresh_no$neighborsThresh[[i]]))
      state_Thresh_no$Ncount[i]<-length(state_Thresh_no$neighborsThresh[[i]])
    }
  }
  #------
  
  #Subset threshold table to only those that have testable neighbors
  state_Thresh_Testable<-subset.data.frame(state_Thresh_no, state_Thresh_no$Ncount>0)
  #------
  
  #create list of combinations to test as well as a count for the number of combinations
  #added zero length support
  if(length(state_Thresh_Testable$Ncount)>0){
    for( i in seq(1, length(state_Thresh_Testable$Ncount),1)){
      #print(i)
      if (state_Thresh_Testable$Ncount[i] >=2){
        print(state_Thresh_Testable$Ncount[i])
        state_Thresh_Testable$test_Combo[i]<-list(arrangements::combinations(state_Thresh_Testable$neighborsThresh[[i]],2, replace = FALSE))
      } else{
        state_Thresh_Testable$test_Combo[i]<-list(0)
      }
      state_Thresh_Testable$test_NCombo[i]<-round(length(state_Thresh_Testable$test_Combo[[i]])/2,0)
    }
  } else{
    state_Thresh_Testable$test_Combo<-list()
  }

  #------
  
  #Preallocate rows to hold results
  dt1 <- data.frame(GEOID=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    ID18=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    Pop18=rep(0,sum(state_Thresh_Testable$Ncount)),
                    Perc18=rep(0,sum(state_Thresh_Testable$Ncount)),
                    ID12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    Pop12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    Perc12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    
                    firstGEOID=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    firstID18=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    firstPop18=rep(0,sum(state_Thresh_Testable$Ncount)),
                    firstPerc18=rep(0,sum(state_Thresh_Testable$Ncount)),
                    firstID12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    firstPop12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    firstPerc12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    
                    secondGEOID=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    secondID18=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    secondPop18=rep(0,sum(state_Thresh_Testable$Ncount)),
                    secondPerc18=rep(0,sum(state_Thresh_Testable$Ncount)),
                    secondID12=rep(0,sum(state_Thresh_Testable$Ncount)), 
                    secondPop12=rep(0,sum(state_Thresh_Testable$Ncount)),
                    secondPerc12=rep(0,sum(state_Thresh_Testable$Ncount))
  )
  #------
  
  #build list of single neighbor results
  #added zero support
  if(length(state_Thresh_Testable$Ncount)>0){
    for ( i in seq(1,length(state_Thresh_Testable$GEOID_1),1)){
      print(i)
      for(n in seq(1,state_Thresh_Testable$Ncount[i],1)){
        print(n)
        row<-which(dt1$GEOID==0)[1] #find first row that has 0 in GEOID
        dt1$GEOID[row]<-as.character(state_Thresh_Testable$GEOID[i])
        dt1$ID18[row]<-as.numeric(as.character(state_Thresh_Testable$Num18BG[i]))
        dt1$Pop18[row]<-as.numeric(as.character(state_Thresh_Testable$TotPovUniv[i]))
        dt1$Perc18[row]<-as.numeric(as.character(state_Thresh_Testable$Pct18BG[i]))
        dt1$ID12[row]<-as.numeric(as.character(state_Thresh_Testable$Num12BG[i]))
        dt1$Pop12[row]<-as.numeric(as.character(state_Thresh_Testable$TotPovUn_1[i]))
        dt1$Perc12[row]<-as.numeric(as.character(state_Thresh_Testable$Pct12BG[i]))
        
        dt1$firstGEOID[row]<-as.character(state$GEOID_1[state_Thresh_Testable$neighborsThresh[[i]][n]])
        dt1$firstID18[row]<-as.numeric(as.character(state$Num18BG[state_Thresh_Testable$neighborsThresh[[i]][n]]))
        dt1$firstPop18[row]<-as.numeric(as.character(state$TotPovUniv[state_Thresh_Testable$neighborsThresh[[i]][n]]))
        dt1$firstPerc18[row]<-as.numeric(as.character(state$Pct18BG[state_Thresh_Testable$neighborsThresh[[i]][n]]))
        dt1$firstID12[row]<-as.numeric(as.character(state$Num12BG[state_Thresh_Testable$neighborsThresh[[i]][n]]))
        dt1$firstPop12[row]<-as.numeric(as.character(state$TotPovUn_1[state_Thresh_Testable$neighborsThresh[[i]][n]]))
        dt1$firstPerc12[row]<-as.numeric(as.character(state$Pct12BG[state_Thresh_Testable$neighborsThresh[[i]][n]]))
        
        #dt1$secondGEOID<-0
        #dt1$secondID18<-0
        #dt1$secondPop18<-0
        #dt1$secondPerc18<-0
        #dt1$secondNum12<-0
        #dt1$secondPop12<-0
        #dt1$secondPerc12<-0
      }
      
    }
  } else{
    print("No neighbors")
  }
  
  
  #------
  
  #calculate weighted averages for single neighbors
  dt1$weight_ID18<-dt1$ID18+dt1$firstID18
  dt1$weight_Pop18<-dt1$Pop18+dt1$firstPop18
  dt1$weight_ID12<-dt1$ID12+dt1$firstID12
  dt1$weight_Pop12<-dt1$Pop12+dt1$firstPop12
  
  dt1$weight_Perc18<-dt1$weight_ID18/dt1$weight_Pop18
  dt1$weight_Perc12<-dt1$weight_ID12/dt1$weight_Pop12
  #------
  
  #Filter to only those items that are above 50% with percents above the threshold
  Single18<-subset.data.frame(dt1, dt1$weight_Perc18 >=0.5 & dt1$Perc18 >=Threshold & dt1$firstPerc18 >=Threshold)
  Single12<-subset.data.frame(dt1, dt1$weight_Perc12 >=0.5 & dt1$Perc12 >=Threshold & dt1$firstPerc12 >=Threshold)
  #------
  
  #Combine Single above thresh into one list and remove duplicates
  Single<-rbind.data.frame(Single18, Single12)
  Single<-Single[!duplicated(Single$GEOID),]
  #------
  
  #Create list of areas that were not found to be eligable 
  '%notin%'<- Negate(`%in%`) #create notin function to negate in operator
  secondCalc<- subset.data.frame(state_Thresh_Testable, state_Thresh_Testable$GEOID_1 %notin% Single$GEOID)
  #------
  
  #Preallocate rows to hold results
  dt2 <- data.frame(GEOID=rep(0,sum(secondCalc$test_NCombo)), 
                    ID18=rep(0,sum(secondCalc$test_NCombo)), 
                    Pop18=rep(0,sum(secondCalc$test_NCombo)),
                    Perc18=rep(0,sum(secondCalc$test_NCombo)),
                    ID12=rep(0,sum(secondCalc$test_NCombo)),
                    Pop12=rep(0,sum(secondCalc$test_NCombo)),
                    Perc12=rep(0,sum(secondCalc$test_NCombo)),
                    
                    firstGEOID=rep(0,sum(secondCalc$test_NCombo)), 
                    firstID18=rep(0,sum(secondCalc$test_NCombo)), 
                    firstPop18=rep(0,sum(secondCalc$test_NCombo)),
                    firstPerc18=rep(0,sum(secondCalc$test_NCombo)),
                    firstID12=rep(0,sum(secondCalc$test_NCombo)),
                    firstPop12=rep(0,sum(secondCalc$test_NCombo)),
                    firstPerc12=rep(0,sum(secondCalc$test_NCombo)),
                    
                    secondGEOID=rep(0,sum(secondCalc$test_NCombo)), 
                    secondID18=rep(0,sum(secondCalc$test_NCombo)), 
                    secondPop18=rep(0,sum(secondCalc$test_NCombo)),
                    secondPerc18=rep(0,sum(secondCalc$test_NCombo)),
                    secondID12=rep(0,sum(secondCalc$test_NCombo)), 
                    secondPop12=rep(0,sum(secondCalc$test_NCombo)),
                    secondPerc12=rep(0,sum(secondCalc$test_NCombo))
  )
  #------
  
  #build list of double neighbor results
  #added zero support
  if(length(state_Thresh_Testable$Ncount)>0){
    for ( i in seq(1,length(secondCalc$GEOID),1)){
      print(i)
      if(secondCalc$test_NCombo[i]==0){
        print('skip')
      } else{
        for(n in seq(1,secondCalc$test_NCombo[i],1)){
          print(n)
          row<-which(dt2$GEOID==0)[1] #find first row that has 0 in GEOID
          dt2$GEOID[row]<-as.character(secondCalc$GEOID[i])
          dt2$ID18[row]<-as.numeric(as.character(secondCalc$Num18BG[i]))
          dt2$Pop18[row]<-as.numeric(as.character(secondCalc$TotPovUniv[i]))
          dt2$Perc18[row]<-as.numeric(as.character(secondCalc$Pct18BG[i]))
          dt2$ID12[row]<-as.numeric(as.character(secondCalc$Num12BG[i]))
          dt2$Pop12[row]<-as.numeric(as.character(secondCalc$TotPovUn_1[i]))
          dt2$Perc12[row]<-as.numeric(as.character(secondCalc$Pct12BG[i]))
          
          dt2$firstGEOID[row]<-as.character(state$GEOID_1[secondCalc$test_Combo[[i]][n,1]])
          dt2$firstID18[row]<-as.numeric(as.character(state$Num18BG[secondCalc$test_Combo[[i]][n,1]]))
          dt2$firstPop18[row]<-as.numeric(as.character(state$TotPovUniv[secondCalc$test_Combo[[i]][n,1]]))
          dt2$firstPerc18[row]<-as.numeric(as.character(state$Pct18BG[secondCalc$test_Combo[[i]][n,1]]))
          dt2$firstID12[row]<-as.numeric(as.character(state$Num12BG[secondCalc$test_Combo[[i]][n,1]]))
          dt2$firstPop12[row]<-as.numeric(as.character(state$TotPovUn_1[secondCalc$test_Combo[[i]][n,1]]))
          dt2$firstPerc12[row]<-as.numeric(as.character(state$Pct12BG[secondCalc$test_Combo[[i]][n,1]]))
          
          dt2$secondGEOID[row]<-as.character(state$GEOID_1[secondCalc$test_Combo[[i]][n,2]])
          dt2$secondID18[row]<-as.numeric(as.character(state$Num18BG[secondCalc$test_Combo[[i]][n,2]]))
          dt2$secondPop18[row]<-as.numeric(as.character(state$TotPovUniv[secondCalc$test_Combo[[i]][n,2]]))
          dt2$secondPerc18[row]<-as.numeric(as.character(state$Pct18BG[secondCalc$test_Combo[[i]][n,2]]))
          dt2$secondID12[row]<-as.numeric(as.character(state$Num12BG[secondCalc$test_Combo[[i]][n,2]]))
          dt2$secondPop12[row]<-as.numeric(as.character(state$TotPovUn_1[secondCalc$test_Combo[[i]][n,2]]))
          dt2$secondPerc12[row]<-as.numeric(as.character(state$Pct12BG[secondCalc$test_Combo[[i]][n,2]]))
        }
      }
      
    }
  } else{
    print("No second calc Neighbors")
  }
  
  #------
  
  #calculate weighted averages for double neighbors
  dt2$weight_ID18<-dt2$ID18+dt2$firstID18+dt2$secondID18
  dt2$weight_Pop18<-dt2$Pop18+dt2$firstPop18+dt2$secondPop18
  dt2$weight_ID12<-dt2$ID12+dt2$firstID12+dt2$secondID12
  dt2$weight_Pop12<-dt2$Pop12+dt2$firstPop12+dt2$secondPop12
  
  dt2$weight_Perc18<-dt2$weight_ID18/dt2$weight_Pop18
  dt2$weight_Perc12<-dt2$weight_ID12/dt2$weight_Pop12
  #------
  
  #Filter to only those items that are above 50% with percents above the threshold
  Double18<-subset.data.frame(dt2, dt2$weight_Perc18 >=0.5 & dt2$Perc18 >=Threshold & dt2$firstPerc18 >=Threshold & dt2$secondPerc18 >=Threshold)
  Double12<-subset.data.frame(dt2, dt2$weight_Perc12 >=0.5 & dt2$Perc12 >=Threshold & dt2$firstPerc12 >=Threshold & dt2$secondPerc12 >=Threshold)
  #------
  
  #Combine double above thresh into one list and remove duplicates
  Double<-rbind.data.frame(Double18, Double12)
  Double<-Double[!duplicated(Double$GEOID),]
  #------
  
  
  #Combine all calcualted areas above the thresholds into one list
  All_calc<-rbind(Single, Double)
  #------
  
  #Create Combined data frame of all USDA data and newly calculated fields
  Merged<-merge.data.frame(state, All_calc, by.x = 'GEOID_1', by.y = 'GEOID', all.x = TRUE)
  Merged[is.na(Merged)] <- 0
  #------
  
  #get shapefiles for cartographic boundries
  STBlockGCB<-block_groups(FIPS,cb=TRUE)
  #------
  
  #Drop extra data
  STBlockGCB@data$STATEFP<-NULL
  STBlockGCB@data$COUNTYFP<-NULL
  STBlockGCB@data$TRACTCE<-NULL
  STBlockGCB@data$BLKGRPCE<-NULL
  STBlockGCB@data$AFFGEOID<-NULL
  STBlockGCB@data$NAME<-NULL
  STBlockGCB@data$LSAD<-NULL
  STBlockGCB@data$ALAND<-NULL
  STBlockGCB@data$AWATER<-NULL
  #------
  
  #Join calculated data and polygons
  StateCalc<-geo_join(STBlockGCB, Merged, by_sp= 'GEOID', by_df= 'GEOID_1', how='left')
  #------
  
  return(StateCalc)

}


BuildMap<-function(BGYes, BGNo, BGCalc){
  #build Map
  map <- leaflet() %>% 
    enableTileCaching() %>% addProviderTiles(providers$OpenStreetMap.BlackAndWhite, group = "Grey")
  map<- map %>% 
    addMapPane("plotT", zIndex = 410) %>% 
    addMapPane("QMi", zIndex = 420) %>% 
    addMapPane("BlockG", zIndex = 430) %>% 
    addMapPane("The", zIndex = 480) %>% 
    addMapPane("top", zIndex = 490)
  
  
  map<- map %>% addPolygons(data=BGYes,
                            weight = 1,
                            fill = TRUE,
                            fillOpacity = 0.2,
                            fillColor = "red",#~BlockGYes(Pct18BG),
                            stroke = TRUE,
                            color = "black",
                            label = paste("Est. Eligable: ", as.character(BGYes$Pct18BG),"%", sep=""),
                            options = leafletOptions(pane = "BlockG"),
                            popup = paste('<b>Elig</b>: Eligable<br>',
                                          '<b>GEOID:</b>', as.character(BGYes$GEOID), '<br>',
                                          '<b>Identified Under 18:</b>', as.character(BGYes$Num18BG), '<br>',
                                          '<b>Total Under 18:</b>', as.character(BGYes$TotPovUniv), '<br>',
                                          '<b>Under 18 Percent:</b>', as.character(BGYes$Pct18BG), '<br>',
                                          '<b>Identified Under 12:</b>', as.character(BGYes$Num12BG), '<br>',
                                          '<b>Total Under 12:</b>', as.character(BGYes$TotPovUn_1), '<br>',
                                          '<b>Under 12 Percent:</b>', as.character(BGYes$Pct12BG), '<br>'),
                            # popup = paste("<b>FIPS:</b> ",as.character(BGYes$GEOID),
                            #               "<br><b> Eligable:</b> ", as.character(BGYes$Pct18BG), "%",
                            #               "<br><b> Population:</b> ", BGYes$TotPovUniv,
                            #               "<br><b> Est. Under 18 Eligable:</b> ", BGYes$Num18BG,
                            #               #"<br><b> Coverage:</b> ", round(BGYes$CoverArea*100,2),"%",
                            #               '<br><b> Area elig:</b> ', as.character(BGYes$ELIGFY20), sep=""),
                            labelOptions = labelOptions(style = list("font-size" = "16px")),
                            group = 'BlockG',
                            highlight = highlightOptions(weight = 5,
                                                         color = "Red",
                                                         bringToFront = TRUE))
  map<- map %>% addPolygons(data=BGNo,
                            weight = 1,
                            fill = TRUE,
                            fillOpacity = 0.2,
                            fillColor = "blue",#~BlockGNo(Pct18BG),
                            stroke = TRUE,
                            color = "black",
                            label = paste("Est. Eligable: ", as.character(BGNo$Pct18BG),"%", sep=""),
                            options = leafletOptions(pane = "BlockG"),
                            popup = paste('<b>Elig</b>: Not Eligable<br>',
                                          '<b>GEOID:</b>', as.character(BGNo$GEOID), '<br>',
                                          '<b>Identified Under 18:</b>', as.character(BGNo$Num18BG), '<br>',
                                          '<b>Total Under 18:</b>', as.character(BGNo$TotPovUniv), '<br>',
                                          '<b>Under 18 Percent:</b>', as.character(BGNo$Pct18BG), '<br>',
                                          '<b>Identified Under 12:</b>', as.character(BGNo$Num12BG), '<br>',
                                          '<b>Total Under 12:</b>', as.character(BGNo$TotPovUn_1), '<br>',
                                          '<b>Under 12 Percent:</b>', as.character(BGNo$Pct12BG), '<br>'),
                            labelOptions = labelOptions(style = list("font-size" = "16px")),
                            group = 'BlockG',
                            highlight = highlightOptions(weight = 5,
                                                         color = "Red",
                                                         bringToFront = TRUE))
  map<- map %>% addPolygons(data=BGCalc,
                            weight = 1,
                            fill = TRUE,
                            fillOpacity = 0.3,
                            fillColor = "orange",
                            stroke = TRUE,
                            color = "black",
                            label = paste("Calc Eligable"),
                            options = leafletOptions(pane = "BlockG"),
                            # popup = paste("<b>FIPS:</b> ",as.character(BGCalc$GEOID),
                            #               "<br><b> Self Eligable:</b> ", as.character(BGCalc$Pct18BG), "%",
                            #               "<br><b> Population:</b> ", BGCalc$TotPovUniv,
                            #               "<br><b> Est. Under 18 Eligable:</b> ", BGCalc$Num18BG,
                            #               #"<br><b> Coverage:</b> ", round(BGCalc$CoverArea*100,2),"%",
                            #               '<br><b> Area elig:</b> ', "Calculated Yes", 
                            #               '<br><b> Calc Block 1:</b> ', BGCalc$firstGEOID,
                            #               '<br><b> Calc Block 2:</b> ', BGCalc$secondGEOID,
                            #               sep=""),
                            popup = paste('<b>Elig</b>: Calculation Eligable<br>',
                                          '<b>GEOID:</b>', as.character(BGCalc$GEOID), '<br>',
                                          '<b>Identified Under 18:</b>', as.character(BGCalc$Num18BG), '<br>',
                                          '<b>Total Under 18:</b>', as.character(BGCalc$TotPovUniv), '<br>',
                                          '<b>Under 18 Percent:</b>', as.character(BGCalc$Pct18BG), '<br>',
                                          '<b>Identified Under 12:</b>', as.character(BGCalc$Num12BG), '<br>',
                                          '<b>Total Under 12:</b>', as.character(BGCalc$TotPovUn_1), '<br>',
                                          '<b>Under 12 Percent:</b>', as.character(BGCalc$Pct12BG), '<br><hr>',
                                          '<b><center>Calculation</center></b><hr>',
                                          '<table border="1">
                                        <tr>
                                        <th style="background-color:#c5d9d5;" align="center">GEOID</th>
                                        <th style="background-color:#c5d9d5;" align="center">Under 18</th>
                                        <th style="background-color:#c5d9d5;" align="center">18 pop</th>
                                        <th style="background-color:#c5d9d5;" align="center">18 perc</th>
                                        <th style="background-color:#c5d9d5;" align="center">Under 12</th>
                                        <th style="background-color:#c5d9d5;" align="center">12 pop</th>
                                        <th style="background-color:#c5d9d5;" align="center">12 perc</th>
                                        </tr>
                                        <tbody>
                                        <tr>
                                        <td align="center">',BGCalc$firstGEOID,'</td>
                                        <td align="center">',BGCalc$firstID18,'</td>
                                        <td align="center">',BGCalc$firstPop18,'</td>
                                        <td align="center">',BGCalc$firstPerc18,'</td>
                                        <td align="center">',BGCalc$firstID12,'</td>
                                        <td align="center">',BGCalc$firstPop12,'</td>
                                        <td align="center">',BGCalc$firstPerc12,'</td>
                                        </tr>
                                        <tr>
                                        <td align="center">', BGCalc$secondGEOID, '</td>
                                        <td align="center">', BGCalc$secondID18, '</td>
                                        <td align="center">', BGCalc$secondPop18, '</td>
                                        <td align="center">', BGCalc$secondPerc18, '</td>
                                        <td align="center">', BGCalc$secondID12, '</td>
                                        <td align="center">', BGCalc$secondPop12, '</td>
                                        <td align="center">', BGCalc$secondPerc12, '</td>
                                        </tr>
                                        <tr>
                                        <td style="background-color:#f7daad;" align="center">Totals</td>
                                        <td style="background-color:#f7daad;" align="center">', BGCalc$weight_ID18, '</td>
                                        <td style="background-color:#f7daad;" align="center">', BGCalc$weight_Pop18, '</td>
                                        <td style="background-color:#f7daad;" align="center">', "", '</td>
                                        <td style="background-color:#f7daad;" align="center">', BGCalc$weight_ID12, '</td>
                                        <td style="background-color:#f7daad;" align="center">', BGCalc$weight_Pop12, '</td>
                                        <td style="background-color:#f7daad;" align="center">', "", '</td>
                                        </tr>
                                        </tbody>
                                        </table>',
                                          '<b>Weighted % under 18:</b>', round(BGCalc$weight_Perc18*100,2), '<br>',
                                          '<b>Weighted % under 12:</b>', round(BGCalc$weight_Perc12*100,2), '<br>'),
                            labelOptions = labelOptions(style = list("font-size" = "16px")),
                            group = 'BlockG',
                            highlight = highlightOptions(weight = 5,
                                                         color = "Red",
                                                         bringToFront = TRUE))
  
  
  map <- map %>% addLegend(position = "bottomright", 
                           labels =  c("USDA Eligible","Calculated Eligible", "Ineligible"), 
                           colors =  c('red','yellow','blue'),
                           opacity = 0.5,
                           title = 'Color Code',
                           group = 'BlockG')
  
  map <- map %>% addSearchOSM(options = searchOptions(autoCollapse = TRUE, minLength = 2))
  return(map)
}







#Read queue, find next item, load RDS, run calc function, generate output
Data<-data.frame()
for (i in seq(1,56,1)){ #separate by FIPS code
  if (i %in% c(3,7,14,43,52)){
    print(paste('Skip non-state:', i))
  } else {
    if(Queue[[i]]$Done[1]==1){
      print("Done")
    } else{
      print("To Be Done")
      Start_state<-Sys.time()
      state<-readRDS(as.character(Queue[[i]]$File[1]))
      FIPS<-as.character(Queue[[i]]$FIPS[1])
      Output<-calculations(FIPS, state, 40)
      Output$CalcElig<-ifelse(Output$ID18 >0, "1", "0")
      CalcAreas<-data.frame(Output[c('GEOID','CalcElig')])
      if(dir.exists(paste(workdir,'/Done', sep = ""))){
        #print("dir exists")
      } else{
        print("creating output directory")
        dir.create(paste(workdir,'/Done', sep = ""))
      }
      setwd(paste(workdir,'/done/',sep = ""))
      write.csv(CalcAreas, paste(as.character(Queue[[i]]$FIPS[1]),"_",as.character(Queue[[i]]$Name[1]),"_Calc.csv",sep = ""), row.names = FALSE)
      setwd(workdir)
      print('Saving calculation output')
      saveRDS(Output, paste(workdir,'/done/', as.character(Queue[[i]]$FIPS[1]),"_",as.character(Queue[[i]]$Name[1]),"_Calc.rds",sep = "")) #save calc output
      print("separating and building map")
      BGYes<-subset.data.frame(Output, Output$ELIGFY17=='YES')
      BGNo<-subset.data.frame(Output, Output$ELIGFY17=='NO')
      BGCalc<-subset.data.frame(BGNo, BGNo$weight_Perc18>=0.5 | BGNo$weight_Perc12>=0.5)
      BGNo<-subset.data.frame(BGNo, BGNo$weight_Perc18<0.5 & BGNo$weight_Perc12<0.5)
      
      maps<-BuildMap(BGYes,BGNo,BGCalc)
      print('Saving Map')
      setwd(paste(workdir,'/done/',sep = ""))
      htmlwidgets::saveWidget(maps,file = paste(as.character(Queue[[i]]$FIPS[1]),"_",as.character(Queue[[i]]$Name[1]),"_Calc.html",sep = ""), selfcontained = TRUE)
      setwd(workdir)
      End_state<-Sys.time()
      Queue[[i]]$BG_Yes_Areas[1]<-length(BGYes$Num18BG)
      Queue[[i]]$BG_Yes_ID18[1]<-sum(as.numeric(as.character(BGYes$Num18BG)))
      Queue[[i]]$BG_No_Areas[1]<-length(BGNo$Num18BG)
      Queue[[i]]$BG_No_ID18[1]<-sum(as.numeric(as.character(BGNo$Num18BG)))
      Queue[[i]]$BG_Calc_Areas[1]<-length(BGCalc$Num18BG)
      Queue[[i]]$BG_Calc_ID18[1]<-sum(as.numeric(as.character(BGCalc$Num18BG)))
      Queue[[i]]$ElapsedTime[1]<-as.numeric(difftime(End_state, Start_state, units = "secs"))
      rm(Output, BGYes, BGNo, BGCalc, maps)
      Queue[[i]]$Done<-1
      Data<-rbind(Data, Queue[[i]])
    }
  }
}

Data<-subset.data.frame(Data, is.na(Data$Name)==FALSE)
write.csv(Data, "Data.csv", row.names = FALSE)
#grep("B", colnames(df)) #ets column undex for match of text




#Stop stopwatch
End<-Sys.time()
print(End-Start)
#------

