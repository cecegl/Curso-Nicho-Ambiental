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

#Carga los paquetes
library(ecospat)
library(ade4)


# Establece el directorio de trabajo
setwd("C:/Users/Usuario/Desktop/cosas serias/MNCN/Curso Nicho Ambiental/Pruebas/2_Results")

#Carga de datos####

#Cargamos los datos de linces que descargamos en la primera sesión
linces <- read.csv("occ_completo.csv")

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

#Análisis de solapamiento####

#Análisis de componentes principales (PCA) para las variables climáticas

#Linces europeos
pca.env <- dudi.pca(rbind(lynx, pardinus)[, 5:8], scannf = FALSE, nf = 2)
ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)

# Puntajes de PCA para toda el área de estudio
scores.globclim <- pca.env$li

# Puntajes de PCA para la distribución de Lynx lynx
scores.sp.lynx <- suprow(pca.env, lynx[, 5:8])$li

# Puntajes de PCA para la distribución de Lynx pardinus
scores.sp.pardinus <- suprow(pca.env, pardinus[, 5:8])$li

# Puntajes de PCA para todo el área de estudio para Lynx lynx
scores.clim.lynx <- suprow(pca.env, lynx[, 5:8])$li

# Puntajes de PCA para todo el área de estudio para Lynx pardinus
scores.clim.pardinus <- suprow(pca.env, pardinus[, 5:8])$li

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


#Añadir el título con control de su posición
#title(main = "Total Niche Overlap", line = 3)  # Ajusta el valor de `line` para más separación

#Desplazamiento de los centroides del nicho
#ecospat.shift.centroids(scores.sp.C_alb, scores.sp.C_min, scores.clim.C_alb, scores.clim.C_min)

# Calcular los índices de dinámica del nicho con diferentes valores de intersección
# Nota: No se debe usar NA como valor para intersection
index_na <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0.1)
index_0 <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0)
index_0_05 <- ecospat.niche.dyn.index(grid.clim.lynx, grid.clim.pardinus, intersection = 0.05)
index_na
index_0
index_0_05
