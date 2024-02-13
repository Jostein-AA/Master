germany_map  <- read_sf('./Shapefiles/Germany')[,c('NAME_2', 'TYPE_2', 'ID_2', 'geometry')]
germany_map_2 <- read_sf('./Shapefiles/Germany_first_level_unpacked/')[,c('NAME_1', 'TYPE_1', 'ID_1', 'geometry')]
germany_border <- read_sf('./Shapefiles/Germany_border/')[, c('geometry')]
Population_data <- read_excel('./Data/germany_population.xlsx')

#Make sure no NA's in the map!
germany_map <- na.omit(germany_map)

#Remove water body from germany map
germany_map <- germany_map[germany_map$TYPE_2 != 'Water body', ]

#Omit population data that is not from Germany
Population_data <- Population_data[grepl('DE', Population_data$Code), ]

#Create TYPE_2 column for population data
NAME_2 = rep('', nrow(Population_data)); TYPE_2 = rep('', nrow(Population_data))
for(i in 1:nrow(Population_data)){
  if(grepl(',', Population_data$Name[i])){
    temp = strsplit(Population_data$Name[i], split = ', ')
    NAME_2[i] = temp[[1]][1]; TYPE_2[i] = temp[[1]][2]
  } else{NAME_2[i] =Population_data$Name[i]; TYPE_2[i] = ''}
}
Population_data$NAME_2 = NAME_2; Population_data$TYPE_2 = TYPE_2

#Merge population data to map of Germany
germany_map$Population <- rep(0, nrow(germany_map))
for(i in 1:nrow(Population_data)){
  #If name is in Germany map names
  if(Population_data$NAME_2[i] %in% germany_map$NAME_2){
    if(Population_data$TYPE_2[i] == ''){ #If landskreis, statdkreis, etc not specified merge directly
      germany_map$Population[
        germany_map$NAME_2 == Population_data$NAME_2[i]] <- Population_data$Pop[i]
    } else{ #If landskreis, statdkreis, etc specified, merge accordingly
      germany_map$Population[
        germany_map$NAME_2 == Population_data$NAME_2[i] &
          germany_map$TYPE_2 == Population_data$TYPE_2[i]] <- Population_data$Pop[i]
    }
  }
}
#Special cases
germany_map$Population[germany_map$NAME_2 == 'Dillingen an der Donau'] = Population_data$Pop[Population_data$NAME_2 == 'Dillingen a.d. Donau']
germany_map$Population[98] = 115436
germany_map$Population[germany_map$NAME_2 == 'Neumarkt in der Oberpfalz'] = Population_data$Pop[Population_data$NAME_2 == 'Neumarkt i. d. OPf.']
germany_map$Population[germany_map$NAME_2 == 'Neustadt an der Aisch-Bad Windsheim'] = Population_data$Pop[Population_data$NAME_2 == 'Neustadt a. d. Aisch-Bad Windsheim']
germany_map$Population[germany_map$NAME_2 == 'Neustadt an der Waldnaab'] = Population_data$Pop[Population_data$NAME_2 == 'Neustadt a. d. Waldnaab']
germany_map$Population[germany_map$NAME_2 == 'Pfaffenhofen an der Ilm'] = Population_data$Pop[Population_data$NAME_2 == 'Pfaffenhofen a. d. Ilm']
germany_map$Population[germany_map$NAME_2 == 'Weiden in der Oberpfalz'] = Population_data$Pop[Population_data$NAME_2 == 'Weiden i. d. Opf']
germany_map$Population[germany_map$NAME_2 == 'Wunsiedel im Fichtelgebirge'] = Population_data$Pop[Population_data$NAME_2 == 'Wunsiedel i. Fichtelgebirge']
germany_map$Population[193] = Population_data$Pop[Population_data$NAME_2 == 'Landkreis Rostock']
germany_map$Population[germany_map$NAME_2 == 'Friesland'] = Population_data$Pop[Population_data$NAME_2 == 'Friesland (DE)']
germany_map$Population[223] = Population_data$Pop[Population_data$NAME_2 == 'Oldenburg (Oldenburg)']
germany_map$Population[germany_map$NAME_2 == 'Osterode am Harz'] = 23993
germany_map$Population[382] = 42000
germany_map$Population[germany_map$NAME_2 == 'Ilm-Kreis'] = Population_data$Pop[Population_data$NAME_2 == 'Ilm-Kreis (NUTS 2021)']
germany_map$Population[germany_map$NAME_2 == 'Saalfeld-Rudolstadt'] = Population_data$Pop[Population_data$NAME_2 == 'Saalfeld-Rudolstadt (NUTS 2021)']
germany_map$Population[germany_map$NAME_2 == 'Schmalkalden-Meiningen'] = Population_data$Pop[Population_data$NAME_2 == 'Schmalkalden-Meiningen (NUTS 2021)']
germany_map$Population[germany_map$NAME_2 == 'Sonneberg'] = Population_data$Pop[Population_data$NAME_2 == 'Sonneberg (NUTS 2021)']
germany_map$Population[germany_map$NAME_2 == 'Suhl'] = 37000
germany_map$Population[germany_map$NAME_2 == 'Wartburgkreis'] = Population_data$Pop[Population_data$NAME_2 == 'Wartburgkreis (NUTS 2021)']

rm(list = c('temp', 'i', 'NAME_2', 'TYPE_2'))
