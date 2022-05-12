
# 1.Instalar paquetes estadisticos y cargar librerias

install.packages("TraMineR")
install.packages("RColorBrewer")
install.packages("TraMineRextras")
install.packages("WeightedCluster")
install.packages("cluster")

library(TraMineR)
library(RColorBrewer)
library(TraMineRextras)
library(WeightedCluster)
library(cluster)


# 2.Cargar y revision de base de datos 

# 2.1.Borrar todo en el "environment"
rm(list=ls())

# 2.2.Cargar datos ELSOC (base enviada titulada "data.elsoc")
load(".../data.elsoc.RData")

# 2.3.Revision de variables de base de datos
names(elsoc)

# 2.3.1.Estatus laboral (1:"Paid work", 2:"Retired", 3:"Out of the labor force")
prop.table(table(elsoc$lbfs.1))*100
prop.table(table(elsoc$lbfs.2))*100
prop.table(table(elsoc$lbfs.3))*100
prop.table(table(elsoc$lbfs.4))*100

# 2.3.2.Sintomas depresivos (1:"None", 2:"Mild", 3:"Moderate", 4:"Moderately Severe", 5:"Severe")
prop.table(table(elsoc$phq.1))*100
prop.table(table(elsoc$phq.2))*100
prop.table(table(elsoc$phq.3))*100
prop.table(table(elsoc$phq.4))*100


# 3.Analisis de secuencias laborales

# 3.1.Creacion de objetos con etiquetas de indicadores de variables
emp.lab <- c("Paid work","Retired","Out of the labor force")

# 3.2.Creacion de objetos con etiquetas abreviadas de indicadores de variables
emp.shortlab <- c("Work","Ret","OLF")

# 3.3.Creacion de objetos con colores para indicadores de variables
emp.cpal <- c("chartreuse4","red3","gray94")

# 3.4.Crear objeto con secuencias
names(elsoc)
emp.seq <- seqdef(elsoc, 2:5, labels=emp.lab, states=emp.shortlab, cpal=emp.cpal)

# 3.5.Medir distancias entre secuencias usando analisis de secuencias usando optimal matching analysis
dist.emp <- seqdist(emp.seq, method="OM",sm="CONSTANT")
View(dist.emp)

# 3.6.Construir tipos de trayectorias con analisis de cluster jerarquico Ward
ward.emp <- hclust(as.dist(dist.emp), method = "ward.D")

# 3.6.1.Comparar soluciones con diferentes numeros de trayectorias
ward.range.emp <- as.clustrange(ward.emp, diss = dist.emp, ncluster = 15)
round(summary(ward.range.emp, max.rank = 15), digits=4)
plot(ward.range.emp, stat=c("ASW","ASWw", "HG", "PBC", "HC"))
plot(ward.range.emp, stat=c("ASW","ASWw", "HG", "PBC", "HC"), norm="zscore")

# 3.6.2.Revision de soluciones con 7 tipos
ward.emp7 <- cutree(ward.emp, k = 7)
xtlab <-c("2016","2017","2018","2019")
seqdplot(emp.seq, group = ward.emp7, xtlab=xtlab, xlab="Tiempo")
seqIplot(emp.seq, group = ward.emp7, xtlab=xtlab, sortv="from.start", xlab="Tiempo")
seqmtplot(emp.seq, group = ward.emp7, xlab="Estatus")
seqHtplot(emp.seq, group = ward.emp7, xlab="Tiempo")

# 3.7. Juntar informacion de tipos con base inicial
ward.emp7 <- as.data.frame(ward.emp7)
elsoc <- cbind(elsoc, ward.emp7)
prop.table(table(elsoc$ward.emp7))*100


# 4.Analisis de secuencias de sintomas depresivos

# 4.1.Creacion de objetos con etiquetas de indicadores de variables
phq.lab <- c("None","Mild","Moderate","Moderately Severe","Severe")

# 4.2.Creacion de objetos con etiquetas abreviadas de indicadores de variables
phq.shortlab <- c("None","Mild","Mod","Mod Sev","Sev")

# 4.3.Creacion de objetos con colores para indicadores de variables
phq.cpal <- c("#005700","#639B64","#F7BF35","#EB8500","#8A1923")

# 4.4.Crear objeto con secuencias
names(elsoc)
phq.seq <- seqdef(elsoc, 6:9, labels=phq.lab, states=phq.shortlab, cpal=phq.cpal)

# 4.5.Medir distancias entre secuencias usando analisis de secuencias usando optimal matching analysis
dist.phq <-seqdist(phq.seq, method="OM",sm="CONSTANT")
View(dist.phq)

# 4.6.Construir tipos de trayectorias con analisis de cluster jerarquico Ward
ward.phq  <- hclust(as.dist(dist.phq), method = "ward.D")

# 4.6.1.Comparar soluciones con diferentes numeros de trayectorias
ward.range.phq <- as.clustrange(ward.phq, diss = dist.phq, ncluster = 15)
round(summary(ward.range.phq, max.rank = 15), digits=4)
plot(ward.range.phq, stat=c("ASW","ASWw", "HG", "PBC", "HC"))
plot(ward.range.phq, stat=c("ASW","ASWw", "HG", "PBC", "HC"), norm="zscore")

# 4.6.2.Revision de soluciones con 5 tipos
ward.phq7 <- cutree(ward.phq, k = 7)
seqdplot(phq.seq, group = ward.phq7, xtlab=xtlab, xlab="Tiempo")
seqIplot(phq.seq, group = ward.phq7, xtlab=xtlab, sortv="from.start", xlab="Tiempo")
seqmtplot(phq.seq, group = ward.phq7, xlab="Estatus")
seqHtplot(phq.seq, group = ward.phq7, xlab="Tiempo")

# 4.7. Juntar informacion de tipos con base inicial
ward.phq7 <- as.data.frame(ward.phq7)
elsoc <- cbind(elsoc, ward.phq7)
prop.table(table(elsoc$ward.phq7))*100


# 5.Analisis de secuencias laborales y de sintomas depresivos

# 5.1.Medir distancias entre secuencias en ambos dominios usando analisis de secuencias usando optimal matching analysis
mcdist.emp.phq <- seqdistmc(channels=list(emp.seq, phq.seq),method="OM",sm="CONSTANT",cweight=c(1,1))

# 5.2.Construir tipos de trayectorias con analisis de cluster jerarquico Ward
ward.emp.phq  <- hclust(as.dist(mcdist.emp.phq), method = "ward.D")

# 5.2.1.Comparar soluciones con diferentes numeros de trayectorias
ward.range.emp.phq <- as.clustrange(ward.emp.phq, diss = mcdist.emp.phq, ncluster = 15)
round(summary(ward.range.emp.phq, max.rank = 15), digits=4)
plot(ward.range.emp.phq, stat=c("ASW","ASWw", "HG", "PBC", "HC"))
plot(ward.range.emp.phq, stat=c("ASW","ASWw", "HG", "PBC", "HC"), norm="zscore")

# 5.3.Revision de soluciones con 6 tipos
ward.emp.phq6 <- cutree(ward.emp.phq, k = 6)
seqdplot(emp.seq, group = ward.emp.phq6, xtlab=xtlab, xlab="Tiempo")
seqdplot(phq.seq, group = ward.emp.phq6, xtlab=xtlab, xlab="Tiempo")
seqIplot(emp.seq, group = ward.emp.phq6, xtlab=xtlab, sortv="from.start", xlab="Tiempo")
seqIplot(phq.seq, group = ward.emp.phq6, xtlab=xtlab, sortv="from.start", xlab="Tiempo")
seqmtplot(emp.seq, group = ward.emp.phq6, xlab="Estatus")
seqmtplot(phq.seq, group = ward.emp.phq6, xlab="Estatus")
seqHtplot(emp.seq, group = ward.emp.phq6, xlab="Tiempo")
seqHtplot(phq.seq, group = ward.emp.phq6, xlab="Tiempo")

# 5.4.Juntar informacion de tipos con base inicial
ward.emp.phq6 <- as.data.frame(ward.emp.phq6)
elsoc <- cbind(elsoc, ward.emp.phq6)
prop.table(table(elsoc$ward.emp.phq6))*100





