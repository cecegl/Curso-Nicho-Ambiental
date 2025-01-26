##########################################

#CURSO ANÁLISIS DE NICHO AMBIENTAL

#Luis R. Pertierra y Celia González López

#SESIÓN 5: Correspondencias ambiente y genoma
##########################################

#Preparamos el entorno de trabajo####

#Limpiamos los objetos de sesiones previas que puedan haber quedado activos
rm(list = ls())

#Instalamos los paquetes necesarios
#install.packages("remotes")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("LEA")

#Cargamos los paquetes
library(remotes)
library(psych)
library(vegan)
library(LEA)
library(permute)
library(lattice)
library(dplyr)
library(corrplot)

setwd("C:/Users/Usuario/Desktop/cosas serias/MNCN/Curso Nicho Ambiental/5_Sesion")

#Cargamos los datos
Genotypes<-read.lfmm("./rock0-adapt.lfmm") #Leemos los datos genéticos
EnvVars<-read.csv("ClimValues_Rock3.csv") # Leemos las variables ambientales

#Comprobamos que no falten datos
dim(Genotypes)
sum(is.na(Genotypes))

EnvVars$MEM1<-as.character(EnvVars$MEM1)

Population<-EnvVars$MEM1
row.names(Genotypes)<- Population
row.names(Genotypes)
identical(rownames(Genotypes), EnvVars$Population)

EnvVars[, 6:31] <- lapply(EnvVars[, 6:31], function(x) as.numeric(as.character(x)))


M <- cor(EnvVars [,6:31])

corrplot(M, method = "circle")

Env2 <- EnvVars[,6:31]

corrplot(as.matrix(Env2))


SelectedVars<-EnvVars %>%
  dplyr::select(c(8,10,11,13,31))
colnames(SelectedVars) <- c("sst","chlo","velo","bio2","pH")


pairs.panels(SelectedVars, scale = TRUE)


HL_rda<- rda(Genotypes ~ ., data = SelectedVars, scale = TRUE)
HL_rda

Rsquared<-RsquareAdj(HL_rda)
Rsquared

summary(eigenvals(HL_rda, model = "constrained"))

Fullsig <- anova.cca(HL_rda, parallel=getOption("mc.cores = 7"))
Fullsig

Axissig <- anova.cca(HL_rda, by="axis", parallel=getOption("mc.cores = 4"))
Axissig


vif.cca(HL_rda)

EnvVars$Population = as.factor(EnvVars$Population)
Pop<- EnvVars$Population
bg <- c("#FE842B","#800E93","#16604E")

plot(HL_rda, type="n", scaling=3)
points(HL_rda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=bg[EnvVars$Population]) #displays individuals
text(HL_rda, scaling=2, display="bp", col="black", cex=1.0)  #displays predictors
legend("bottomright", legend=levels(Pop), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg) #adds legend


plot(HL_rda, type="n", scaling=3, choices = c(1,3))
points(HL_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3, choices = c(1,3)) #displays SNPs
points(HL_rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[Pop], choices = c(1,3)) #displays individuals
text(HL_rda, scaling=3, display="bp", col="#0868ac", cex=0.5, choices = c(1,3))  #displays predictors
legend("bottomright", legend=levels(Pop), bty="n", col="gray32", pch=21, cex=0.8, pt.bg=bg) #adds legend

