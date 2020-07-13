##Análisis

#Librerías
library(readxl)
library(sp)
library(sf)
library(tidyverse)
library(rgdal)
library(RColorBrewer)
library(spdep)
library(tmap)
library(tmaptools)
library(spatialreg)
library(epitools)
library(DCluster)
library(plotrix)
library(MASS)
library(mgcv)
library(INLA)

indicadores <- read_excel("Datos/Datos Cantonales.xlsx")
cantones_sp <- read_sf(dsn="Datos",layer = "Cantones")

#Merge
datos_sf <- merge(cantones_sp,indicadores)
datos_sf <- datos_sf[,-c(1:4)]
rm(cantones_sp,indicadores)

#Quitar la Isla del Coco
new_bb = c(286803.0, 889158.2, 658864.2,1241118.1)
names(new_bb) = c("xmin", "ymin", "xmax", "ymax")
attr(new_bb, "class") = "bbox"
attr(st_geometry(datos_sf), "bbox") <- new_bb
rm(new_bb)

datos_sf2 <- datos_sf


# Primer ploteo -----------------------------------------------------------


#Estadísticas descriptivas
#pdf("Trabajo Escrito/F11.#pdf")

tm_shape(datos_sf) +
  tm_polygons("dengue", palette=c("lightgreen","tomato"), legend.hist=TRUE)

#dev.off()

#pdf("Trabajo Escrito/F12.#pdf")

tm_shape(datos_sf) +
  tm_polygons("dengue", palette=c("lightgreen","tomato"),style="quantile")

#dev.off()

#pdf("Trabajo Escrito/FA1.#pdf")

tm_shape(datos_sf) +
  tm_polygons("casos", palette=c("lightgreen","tomato"), legend.hist=TRUE)

#dev.off()

#pdf("Trabajo Escrito/FA2.#pdf")

tm_shape(datos_sf) +
  tm_fill("dengue",style="sd",palette=c("lightgreen","tomato")) +
  tm_borders()

#dev.off()

#pdf("Trabajo Escrito/FA3.#pdf")

tm_shape(datos_sf) +
  tm_polygons("tugurio",n=6, palette="-Spectral")

#dev.off()

#pdf("Trabajo Escrito/FA4.#pdf")

tm_shape(datos_sf) +
  tm_polygons("densidad",n=6, palette="-Spectral", style="quantile")

#dev.off()

#pdf("Trabajo Escrito/FA5.#pdf")

tm_shape(datos_sf) +
  tm_polygons("residuos",n=6, palette="Spectral")

#dev.off()

#pdf("Trabajo Escrito/FA6.#pdf")

tm_shape(datos_sf) +
  tm_polygons("acueducto",n=6, palette="Spectral")

#dev.off()


# Análisis Frecuentista ---------------------------------------------------

datos_sp <- as(datos_sf,"Spatial")
datos_sp@bbox <- matrix(c(286803.0, 889158.2, 658864.2,1241118.1),ncol = 2,byrow = F)

rm(datos_sf)

#Análisis Vecinos
coords <- coordinates(datos_sp)
id <-row.names(datos_sp) 

nb.1 <- poly2nb(datos_sp,queen = T)
nb.2 <- poly2nb(datos_sp,queen = F)
nb.3 <- knn2nb(knearneigh(coords, k=2), row.names=id)
nb.4 <- knn2nb(knearneigh(coords, k=4), row.names=id)

#pdf("Trabajo Escrito/F21.#pdf")

plot(datos_sp, axes=F, border="gray")
plot(nb.1,coords, pch = 20, cex = 0.6, add = T, col = "red")

#dev.off()

#pdf("Trabajo Escrito/F22.#pdf")

plot(datos_sp, axes=F, border="gray")
plot(nb.4,coords, pch = 20, cex = 0.6, add = T, col = "red")

#dev.off()

#Matrices de pesos
w.11 <- nb2listw(nb.1,style = "W")
w.12 <- nb2listw(nb.1,style = "B")
w.13 <- nb2listw(nb.1,style = "S")

w.21 <- nb2listw(nb.2,style = "W")
w.22 <- nb2listw(nb.2,style = "B")
w.23 <- nb2listw(nb.2,style = "S")

w.31 <- nb2listw(nb.3,style = "W")
w.32 <- nb2listw(nb.3,style = "B")
w.33 <- nb2listw(nb.3,style = "S")

w.41 <- nb2listw(nb.4,style = "W")
w.42 <- nb2listw(nb.4,style = "B")
w.43 <- nb2listw(nb.4,style = "S")

#Test de Moran
moran.test(datos_sp$dengue,listw=w.11)
moran.test(datos_sp$dengue,listw=w.12)
moran.test(datos_sp$dengue,listw=w.13)
moran.test(datos_sp$dengue,listw=w.21)
moran.test(datos_sp$dengue,listw=w.22)
moran.test(datos_sp$dengue,listw=w.23)
moran.test(datos_sp$dengue,listw=w.31)
moran.test(datos_sp$dengue,listw=w.32)
moran.test(datos_sp$dengue,listw=w.33)
moran.test(datos_sp$dengue,listw=w.41)
moran.test(datos_sp$dengue,listw=w.42)
moran.test(datos_sp$dengue,listw=w.43)

#Se elije Reina y matriz W
rm(nb.2,nb.3,nb.4,w.12,w.13,w.21,w.22,w.23,w.31,w.32,w.33,w.41,w.42,w.43,coords,id)

#Casos de influencia

#pdf("Trabajo Escrito/F31.#pdf")

#Supuestos

datos_sp@data <- datos_sp@data[,c(1:8)]

# Modelos

#Lineal
m1 <- lm(dengue~tugurio+densidad+residuos+acueducto,data = datos_sp)
summary(m1)
step(m1)
m1 <- lm(sqrt(dengue)~residuos+acueducto,data = datos_sp)
summary(m1)
#plot(m1)
lm.morantest(m1, listw = w.11)

#SAR
m2 <- spautolm(dengue~tugurio+densidad+residuos+acueducto,data = datos_sp,listw=w.11)
summary(m2)
m2 <- spautolm(dengue~residuos+acueducto,data = datos_sp,listw=w.11)
summary(m2)
moran.mc(residuals(m2),w.11, 999)

#CAR
m3 <- spautolm(dengue~tugurio+densidad+residuos+acueducto,data = datos_sp,listw=w.11,family = "CAR")
summary(m3)
m3 <- spautolm(dengue~residuos+acueducto,data = datos_sp,listw=w.11,family = "CAR")
summary(m3)
moran.mc(residuals(m3),w.11, 999)

#GLM
datos_sp$x<-coordinates(datos_sp)[,1]/1000
datos_sp$y<-coordinates(datos_sp)[,2]/1000
m4 <- gam(as.integer(casos)~+residuos+acueducto+offset(log(pob))+s(x,y), data=datos_sp, family= "quasipoisson")
summary(m4)
moran.mc(residuals(m4),w.11, 999)

#Residuales
datos_sp$Lineal <- residuals(m1)
datos_sp$SAR <- residuals(m2)
datos_sp$CAR <- residuals(m3)
datos_sp$GAM <- residuals(m4)

#pdf("Trabajo Escrito/F4.#pdf")

spplot(datos_sp,c("SAR", "CAR", "Lineal","GAM"), 
       at=c(seq(-500,2000,500)), 
       col.regions=colorRampPalette(gry)(7))

#dev.off()

#Epidem

datos_sp$observados <- datos_sp$casos
r <- sum(datos_sp$observados)/sum(datos_sp$pob)
datos_sp$esperados <- datos_sp$pob*r
datos_sp$SMR <- datos_sp$observados/datos_sp$esperados

#pdf("Trabajo Escrito/FA8.#pdf")

spplot(datos_sp,c("observados","esperados"), col.regions=rev(brewer.pal(7, "RdYlGn")), cuts=6)

#dev.off()

#pdf("Trabajo Escrito/F5.#pdf")

spplot(datos_sp,"SMR",col.regions=rev(brewer.pal(7, "RdYlGn")), cuts=6)

#dev.off()

int <- pois.exact(datos_sp$SMR)
int <- cbind(int,datos_sp$NOM_CANT_1)
col <- 1*(int$lower>1)
col <- ifelse(col==0,"grey","red")
linea <- ifelse(col=="grey",4,1) 

#pdf("Trabajo Escrito/F6.#pdf")

plotCI(x = 1:81, y = int$x, ui = int$upper,li = int$lower,pch=18,err="y",
       col=col,sfrac = 0,xlab="Cantones",ylab="Riesgo Relativo",xaxt="n")
abline(h=1,col="grey",lty=2,lwd=1.75)




# INLA --------------------------------------------------------------------

#Estructuras de vecinos y matriz W
datos_sp2 <- as(datos_sf2,"Spatial")
nb <- poly2nb(datos_sp2,queen = T)
rm(datos_sp2)

#Valores esperado
r <- sum(datos_sf2$casos)/sum(datos_sf2$pob)
datos_sf2$esperados <- datos_sf2$pob*r
rm(r)

datos_inla <- datos_sf2 %>%
  st_drop_geometry() %>% 
  mutate(region_1 = row_number(),
         region_2 = row_number()) %>% 
  dplyr::select(region_1,esperados,casos,tugurio,
                densidad,residuos,acueducto,region_2) 
rm(datos_sf2)

# Fórmulas INLA -----------------------------------------------------------

nb2INLA("dengue.graph", nb)


formula1 <- casos ~ f(region_2,model="besag",graph.file="dengue.graph",
                     param=c(1,0.00005))+ residuos + acueducto + f(region_1)

m_in1 <- inla(formula1 ,family="poisson",
                    data=datos_inla,E=esperados,
                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                    control.predictor = list(compute = TRUE))


formula2 <- casos ~ f(region_2,model="besagproper",graph.file="dengue.graph",
                      param=c(1,0.00005))+ residuos + acueducto + f(region_1)

m_in2 <- inla(formula2 ,family="poisson",
              data=datos_inla,E=esperados,
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
              control.predictor = list(compute = TRUE))

formula3 <- casos ~ f(region_2,model="bym",graph.file="dengue.graph",
                      param=c(1,0.00005))+ residuos + acueducto + f(region_1)

m_in3 <- inla(formula3 ,family="poisson",
              data=datos_inla,E=esperados,
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
              control.predictor = list(compute = TRUE))

formula4 <- casos ~ f(region_2,model="iid",graph.file="dengue.graph",
                      param=c(1,0.00005))+ residuos + acueducto + f(region_1)

m_in4 <- inla(formula4 ,family="poisson",
              data=datos_inla,E=esperados,
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
              control.predictor = list(compute = TRUE))

summary(m_in1)
summary(m_in2)
summary(m_in3)
summary(m_in4)

datos_sp$besag <- m_in1$summary.fitted.values[,"mean"]
datos_sp$besagproper <- m_in2$summary.fitted.values[,"mean"]
datos_sp$bym <- m_in3$summary.fitted.values[,"mean"]
datos_sp$iid <- m_in4$summary.fitted.values[,"mean"]

#Intervalos
int <- cbind(m_in1$summary.fitted.values[, "mean"],
             m_in1$summary.fitted.values[, "0.025quant"],
             m_in1$summary.fitted.values[, "0.975quant"])
int <- as.data.frame(int)
str(int)
int <- cbind(int,datos_sp$NOM_CANT_1)


colnames(int) <- c("x","lower","upper","nombre")

col <- 1*(int$lower>1)
col <- ifelse(col==0,"grey","red")
linea <- ifelse(col=="grey",4,1) 

pdf("Trabajo Escrito/FF1.pdf")
plotCI(x = 1:81, y = int$x, ui = int$upper,li = int$lower,pch=18,err="y",
       col=col,sfrac = 0,xlab="Cantones",ylab="Riesgo Relativo",xaxt="n")
abline(h=1,col="grey",lty=2,lwd=1.75)
dev.off()

pdf("Trabajo Escrito/FF2.pdf")
spplot(datos_sp, c("bym", "iid","besag","besagproper"))
dev.off()
#Ploteo de betas

marginal1.0 <- data.frame(inla.smarginal(m_in1$marginals.fixed$`(Intercept)`))

pdf("Trabajo Escrito/FF3.1.pdf")
ggplot(marginal1.0, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(beta[0]), y = "Density") +
  geom_vline(xintercept = 0, col = "blue") + theme_bw()
dev.off()

marginal1.1 <- data.frame(inla.smarginal(m_in1$marginals.fixed$residuos))
pdf("Trabajo Escrito/FF3.2.pdf")
ggplot(marginal1.1, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(beta[1]), y = "Density") +
  geom_vline(xintercept = 0, col = "blue") + theme_bw()
dev.off()

marginal1.2 <- data.frame(inla.smarginal(m_in1$marginals.fixed$acueducto))
pdf("Trabajo Escrito/FF3.3.pdf")
ggplot(marginal1.2, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(beta[2]), y = "Density") +
  geom_vline(xintercept = 0, col = "blue") + theme_bw()
dev.off()

#Fin del análisis