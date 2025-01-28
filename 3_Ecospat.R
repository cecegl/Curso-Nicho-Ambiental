#########################################

#CURSO ANÁLISIS DE NICHO AMBIENTAL

#Luis R. Pertierra y Celia González López

#SESIÓN 3: Solapamiento de nicho
#########################################


#Preparación del entorno de trabajo####

#Limpiamos los objetos de sesiones previas que puedan haber quedado activos
rm(list = ls())

#Instala los paquetes necesarios para la sesión
#install.packages("ecospat")
#install.packages("ade4")
#install.packages("terra")

#Carga los paquetes
library(ecospat)
library(ade4)
library(terra)


# Establece el directorio de trabajo
results_folder <- "C:/Users/Usuario/Desktop/cosas serias/MNCN/Curso Nicho Ambiental/Pruebas/2_Results"
setwd(results_folder)
#Carga de datos####

#Cargamos los datos de linces que descargamos en la primera sesión
linces <- read.csv("occ_variables_nocorr.csv")

summary(linces)

#Separamos por especie
lynx <- linces %>%
  filter(sp == "Lynx lynx")
rufus <- linces %>%
  filter(sp == "Lynx rufus")
pardinus <- linces %>%
  filter(sp == "Lynx pardinus")
canadensis <- linces %>%
  filter(sp == "Lynx canadensis")

lynx_samp <- lynx[sample(nrow(lynx),100),]
pardinus_samp <- pardinus[sample(nrow(pardinus),100),]

#Vamos a añadir datos aleatorios de falsas ausencias, ya que necesitamos comparar
#el espacio que ocupa la especie realmente frente al que no

#Recordad este código de la primera sesión
studyArea = extent(-10,60,25,75)

world_map <- world(resolution = 3,
                   path = results_folder)

my_map <- crop(x = world_map, y = studyArea)

#Vamos a crear 500 puntos aleatorios, para un análisis real no son muchos, pero
#en este caso no queremos que los tiempos de computación se alarguen mucho
bg_rand <- spatSample(my_map, 500, "random")

plot(my_map)

points(bg_rand)

#Extraemos las coordenadas
ceros <- data.frame(terra::geom(bg_rand)[,c("x","y")])
summary(ceros)
#Si recordais del primer día, hicimos una relación de las variables climáticas
#con los puntos de presencia, tenemos que hacer lo mismo con el background
clima_completo <- rast("Climate_Layers.tif") #las variables climáticas que guardamos
#el primer día

#Asociamos coordenadas de ausencias a variables climáticas (como el primer día)
bioclim_extract <- raster::extract(clima_completo,
                                   ceros[,c("x","y")])
print(names(linces)) #Comprobamos qué variables nos habíamos quedado previamente
#y las cogemos también de nuestra nueva tabla

print(names(bioclim_extract))

clima <- subset(bioclim_extract,select = c(Climate_Layers_2, Climate_Layers_4,
                                         Climate_Layers_8, Climate_Layers_13))

pseudoausencias <- data.frame(cbind(ceros, clima))

columnas <- c("lon","lat","bio2","bio4","bio8","bio13")

#Las próximas líneas hasta el principio de los análisis se encargan de encadenar
#los datos de las especies con las falsas ausencias que acabamos de crear
names(lynx_samp)
names(pardinus_samp)

lynx_samp <- subset(lynx_samp,select = -c(X,sp))
pardinus_samp <- subset(pardinus_samp,select = -c(X,sp))

colnames(pseudoausencias) <- columnas
colnames(lynx_samp) <- columnas
colnames(pardinus_samp) <- columnas

pseudoausencias$pa <- rep(0,length(pseudoausencias$lon))
lynx_samp$pa <- rep(1,length(lynx_samp$lon))
pardinus_samp$pa <- rep(1,length(pardinus_samp$lon))

lynxlynx <- rbind(lynx_samp,pseudoausencias)
lynxpardinus <- rbind(pardinus_samp,pseudoausencias)

#Análisis de solapamiento####

#Análisis de componentes principales (PCA) para las variables climáticas

#Vamos a comparar las dos especies de linces europeos entre ellas
#Esta tarde podréis probar vosotros con vuestras propias especies

#Necesitamos un sample aleatorio del mismo número de ocurrencias para ambas especies

pca.env <- dudi.pca(rbind(lynxlynx, lynxpardinus)[, 3:6], scannf = FALSE, nf = 2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)

# Puntajes de PCA para toda el área de estudio
scores.globclim <- pca.env$li

# Puntajes de PCA para la distribución de Lynx lynx
scores.sp.lynx <- suprow(pca.env, lynxlynx[which(lynxlynx[,7] == 1), 3:6])$li

# Puntajes de PCA para la distribución de Lynx pardinus
scores.sp.pardinus <- suprow(pca.env, lynxpardinus[which(lynxpardinus[,7] == 1), 3:6])$li

# Puntajes de PCA para todo el área de estudio para Lynx lynx
scores.clim.lynx <- suprow(pca.env, lynxlynx[, 3:6])$li

# Puntajes de PCA para todo el área de estudio para Lynx pardinus
scores.clim.pardinus <- suprow(pca.env, lynxpardinus[, 3:6])$li

# Calcular la Grilla de Densidades de Ocurrencia con ecospat.grid.clim.dyn()
# Para Lynx lynx
grid.clim.lynx <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.lynx, sp = scores.sp.lynx, R = 200, th.sp = 0)

# Para Lynx pardinus
grid.clim.pardinus <- ecospat.grid.clim.dyn(glob = scores.globclim, glob1 = scores.clim.pardinus, sp = scores.sp.pardinus, R = 200, th.sp = 0)

# Calcular el Solapamiento de Nicho con ecospat.niche.overlap()
D.overlap <- ecospat.niche.overlap(grid.clim.lynx, grid.clim.pardinus, cor = TRUE)$D
D.overlap

# Realizar el Test de Equivalencia de Nicho con ecospat.niche.equivalency.test() según Warren et al. (2008)
eq.test <- ecospat.niche.equivalency.test(grid.clim.lynx, grid.clim.pardinus, rep = 70)
# Test de Equivalencia de nicho H1: ¿Es el solapamiento entre el nicho de ambas especies mayor
#al solapamiento entre dos nichos aleatorios?

# Realizar el Test de Similaridad de Nicho con ecospat.niche.similarity.test()
sim.test <- ecospat.niche.similarity.test(grid.clim.lynx, grid.clim.pardinus, rep = 70, rand.type = 2)
# Test de Similaridad de nicho H1: ¿Es el solapamiento entre el nicho nativo y el invadido mayor a cuando el nicho invadido ha sido introducido aleatoriamente en el área invadida?

# Graficar los Test de Equivalencia/Similaridad de Nicho
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
ecospat.plot.overlap.test(sim.test, "I", "Similarity")

# Delimitación de categorías de nicho y cuantificación de la dinámica del nicho en climas análogos con ecospat.niche.dyn.index()
niche.dyn <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0)

####Cambio de colores por ecologia
ecospat.plot.niche.dyn( grid.clim.pardinus, grid.clim.lynx, 
  interest = 2, 
  title = "",
  name.axis1 = "PC1", 
  name.axis2 = "PC2",
)

#Desplazamiento de los centroides del nicho
ecospat.shift.centroids(scores.sp.lynx, scores.sp.pardinus, scores.clim.lynx, scores.clim.pardinus)

# Calcular los índices de dinámica del nicho con diferentes valores de intersección
# Nota: No se debe usar NA como valor para intersection
index_na <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0.1)
index_0 <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0)
index_0_05 <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0.05)
index_na
index_0
index_0_05

#Probad ahora a estudiar el solapamiento de nicho para una única variable ambiental
#(Pista: sólo tenéis que seleccionar una única columna de entre las variables
#en las líneas de código donde antes ponía [,3:6])
