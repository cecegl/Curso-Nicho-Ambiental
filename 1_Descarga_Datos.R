##################################################

#CURSO ANÁLISIS DE NICHO AMBIENTAL

#Luis R. Pertierra y Celia González López

#SESIÓN 1: Confección del nicho y datos necesarios
##################################################

#Preparación del entorno de trabajo####

#Podemos limpiar previamente todos los elementos que tengamos en la sesión
rm(list = ls()) #Este paso es opcional, pero es apropiado para evitar solapamiento
                #con elementos previos que tuvieramos creados

#Instala los paquetes necesarios para la sesión
#install.packages("rgbif")
#install.packages("tidyverse")
#install.packages("geodata")
#install.packages("terra)
#install.packages("usdm")
#install.packages("dplyr")
#install.packages("CoordinateCleaner")

#Carga los paquetes
#(tendrás que hacer esto cada vez que inicies una nueva sesión en R)
library(rgbif)
library(tidyverse)
library(geodata)
library(terra)
library(raster)
library(usdm)
library(dplyr)
library(CoordinateCleaner)

#Elige un directorio de trabajo
#(esta será la carpeta donde guardaremos los datos)
results_folder <- "" #añadid vuestra ruta entre las comillas
setwd(results_folder)

#Añadimos nuestros credenciales de GBIF para la descarga de datos

GBIF_USER="" #vuestro nombre de usuario de GBIF entre comillas
GBIF_PWD="" #vuestra contraseña de GBIF (permitid que no haga públca la mía)
GBIF_EMAIL="" #vuestro correo asociado a GBIF

#Descarga de datos de especies####

#Especie (o especies) con la que vamos a trabajar
spp_names <- c("especie1", "especie2",
              "especie3", "especie4") #cambiadlo por el nombre de las especies que queráis

#Descargamos y limpiamos los datos

#Primero debemos obtener el identificador de cada una de las especies:
taxon_keys <- name_backbone_checklist(spp_names) %>%
    filter(!matchType == "NONE") %>%
    pull(usageKey)
  
#Y ahora procederemos a la descarga de los datos
#Podemos modificar la búsqueda en base a los datos que deseemos para nuestro estudio
gbif_download <- occ_download(
    pred_in("taxonKey", taxon_keys),
    pred_in("occurrenceStatus", "PRESENT"), #sólo presencias, obviar ausencias
    pred("hasCoordinate", TRUE), #sólo queremos datos que tengan coordenadas
    pred("hasGeospatialIssue", FALSE), #además, filtramos los que tengan algún problema
    #en sus datos geoespaciales
    #pred("continent","europe"), #si queremos datos más concretos podemos acotar
    #la búsqueda a un único país o continente
    pred_not(pred_in("BASIS_OF_RECORD", c("FOSSIL_SPECIMEN"))), #descartamos registros fósiles
    pred_gte("year", 1970), #podemos filtrar el año de inicio para eliminar registros muy antiguos
    user = GBIF_USER, pwd = GBIF_PWD, email = GBIF_EMAIL,
    format = "SIMPLE_CSV"
  ) #todos los argumentos de esta función podeis modificarlos en base a los datos que querais obtener
  
#Estas descargas llevan tiempo, puedes comprobar el estado de la tuya con el
#siguiente comando:
occ_download_wait(gbif_download)

  
#Cuando la descarga se completé, aparecerá un mensaje en la consola
#y entonces podrás guardarla en tu directorio con el siguiente comando:
presencias <- occ_download_get(gbif_download) |>
    occ_download_import()
  
#Es importante comprobar que no existan registros duplicados en nuestra descarga
flags <- clean_coordinates(x = presencias,
                                       lon = "decimalLongitude",
                                       lat = "decimalLatitude",
                                       species = "species",
                                       tests = c("validity", "equal",
                                                 "zeros", "duplicates",
                                                 "outliers"),
                                       outliers_td = 1)

summary(flags)

presencias <- clean_coordinates(x = presencias,
                                   lon = "decimalLongitude",
                                   lat = "decimalLatitude",
                                   species = "species",
                                   tests = c("validity", "equal",
                                             "zeros", "outliers",
                                             "duplicates"),
                                   outliers_td = 1,
                                   value = "clean")

presencias <- presencias %>%
  mutate(across(c("decimalLatitude","decimalLongitude"),round,3))

presencias <- clean_coordinates(x = presencias,
                            lon = "decimalLongitude",
                            lat = "decimalLatitude",
                            species = "species",
                            tests = c("validity", "equal",
                                      "zeros", "outliers",
                                      "duplicates"),
                            outliers_td = 1,
                            value = "clean")

#Descarga de variables ambientales####
  
clim <- worldclim_global(var = "bio", #todas las variables climáticas disponibles
                         res = 10, #resolución espacial de los datos
                         path = results_folder,
                         #country = "spain"
                         )

#Unir datos de presencia y variables ambientales####

#Vamos a cortar el espacio geográfico en base a nuestros datos de presencia
studyArea = extent(-180,180,10,80) #Modificad esta extensión en base a vuestros datos (lon,lon,lat,lat)

world_map <- world(resolution = 3,
                   path = results_folder)

my_map <- crop(x = world_map, y = studyArea)

#Dibujando el mapa podemos comprobar que todo está funcionando correctamente
plot(my_map, axes = TRUE,
     col = "grey95")

points(linces[linces$species=="Lynx pardinus",23:22],pch=20,cex=0.3,col="orange")
points(linces[linces$species=="Lynx lynx",23:22],pch=20,cex=0.3,col="lightblue")
points(linces[linces$species=="Lynx rufus",23:22],pch=20,cex=0.3,col="olivedrab3")
points(linces[linces$species=="Lynx canadensis",23:22],pch=20,cex=0.3,col="brown")


#Ahora vamos a cortar los datos de clima a la extensión deseada

bioclim_data <- crop(x = clim, y = studyArea)

#Comprobamos que ha funcionado dibujando en el mapa la primera de las variables climáticas
plot(bioclim_data[[1]])

#Hacemos un stack de todas las variables ambientales para tenerlas en un mismo raster
clima_completo <- stack(bioclim_data)

#plot(clima_completo)

#guardamos un raster con las variables deseadas para otras sesiones
writeRaster(clima_completo, "Climate_Layers.tif")

#Por último y como paso MÁS importante, tenemos que asignar valores de nuestras
#variables climáticas a cada uno de los puntos geográficos donde tenemos una presencia

linces_simple <- linces[,c("species","decimalLongitude","decimalLatitude")]
#nos quedamos sólo con los datos que nos interesan

columnas <- c("sp","lon","lat")

colnames(linces_simple) <- columnas


#extraemos las variables climáticas en las coordenadas de presencias
bioclim_extract <- raster::extract(clima_completo,
                           linces_simple[,c("lon","lat")])

#Las variables climáticas pueden estar correlacionadas entre sí, lo cuál podría
#afectar a nuestros análisis sobre o subestimando la influencia de algunas de ellas

#Para evitarlo vamos a descartar las variables correlacionadas mediante un análisis
#de inflación de varianza (VIF)
clima_df <- as.data.frame(bioclim_extract)

vif_select <- vifstep(clima_df, th = 2)

no_corr <- vif_select@results$Variables

variables_nocorr = subset(clima_df, select = no_corr)

#lo unimos a los datos de gbif (guardamos también un df completo sin descartar variables)
occ_todas_variables <- data.frame(cbind(linces_simple, clima_df))
occ_varnocorr <- data.frame(cbind(linces_simple, variables_nocorr))

summary(occ_todas_variables)
summary(occ_varnocorr) #Existen puntos donde nuestras variables seleccionadas no están
#medidas
# por tanto, debemos borrarlos (para posteriores análisis)
occ_todas_variables <- na.omit(occ_todas_variables)
occ_varnocorr <- na.omit(occ_varnocorr)

#y por ultimo, lo guardamos para usar más adelante
write.csv(occ_todas_variables, "occ_completo.csv")
write.csv(occ_varnocorr, "occ_variables_nocorr.csv")
