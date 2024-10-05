# Paquetes requeridos
library(picante)

# Ahora que entendemos filogenias, vamos a entender como podemos usar la información
# evolutiva para realizar inferencias sobre el ensamblaje de comunidades.

# Sin embargo, exploremos un poco la forma como se ven las filogenias en R
# Existen varias formas de cargar una filogenia. Por ejemplo, en el paquete APE
# "Analysis of Phylogenetics and Evolution". Por ejemplo, es posible tener 
# arboles evolutivos reconstruidos en otro software e importarlo a R para reproducirlo
# La forma más sencilla de representar un arbol evolutivo en R es a traves de 
# formato parentético. Este formato normalmente es denominado como "Newick". Sin embargo
# hay otros formatos como Nexus que tambien pueden ser leido por R a tavés de APE.

# Por ejemplo, si queremos recrear el arbol evolutivo de los grandes simios:

simios.new<-"((((Homo:0.21, Pongo:0.21):0.28, Macaca:0.49):0.13, 
                           Ateles:0.62):0.38, Galago:1.00);"
simios.tree<-read.tree(text=simios.new)

# Veamos la estructura del archivo del arbol

simios.tree

str(simios.tree)

plot(simios.tree)

# Ejercicio

# Siguiendo el ejemplo de simios.new, recrearun arbol filogenético de sus organismos 
# de interés usando la notación parentética y graficarlo.

#################################################################################
# Iniciemos entonces caracterizando la diversidad de una comunidad en términos filognéticos
# Para esto, vamos a utilizar comunidades hipotéticas creadas por Campbell Webb y 
# Steven Kembel en el año 2002 cuando crearon el software de Phylocom.

data(phylocom)

# Revisemos que datos tiene phylocom

str(phylocom)

# Vamos a nombrar cada componente separadamente

arbol<-phylocom$phylo
comms<-phylocom$sample
caracteres<-phylocom$traits

# Veamos en uno o varios arboles la presencia o ausencia de especies en comunidades

especies.muestreo<-colnames(comms)

offs<-seq(0,1,length=6)
my.cols<-brewer.pal(6,"Spectral")

par(mar=c(0,0,0,0))
plot(arbol,direction="upwards",show.tip.label=FALSE,y.lim=c(0,6)
     ,x.lim=c(-1,Ntip(arbol)))
for(i in 1:nrow(comms)){
  comm1<-especies.muestreo[comms[i,]>0]
  points(match(comm1,arbol$tip.label),rep(5+offs[i],length(comm1)),pch=15
         ,col=my.cols[i],cex=2)
  text(-0.5,5+offs[i],rownames(comms)[i])
}

# 1. Diversidad filogenética propuesta por Faith

faith.pd<-pd(comms,arbol)

# Como interpretamos estas diferencias en las comunidades?

barplot(faith.pd$PD,names.arg=rownames(faith.pd),space=0,border="white")
box(bty="l")

# Veamos como se calcula la diversidad filogenética
# En este caso de este árbol, todas las longitudes de ramas son iguales.

plot(arbol,direction="upwards",show.tip.label=FALSE,y.lim=c(0,6)
     ,x.lim=c(-1,Ntip(arbol)),edge.color=rep(c(my.cols[1],"black"),times=c(16,62-16))
     ,edge.width=rep(2:1,times=c(16,62-16)))
comm1<-especies.muestreo[comms[1,]>0]
points(match(comm1,arbol$tip.label),rep(5+offs[1],length(comm1)),pch=15
       ,col=my.cols[1],cex=2)


# Ejercicio, replicar este mismo código para colorear las ramas de las especies
# que componen la comunidad Clump2a


# Veamos en un ejemplo real como esta relacionada la diversidad filogenética con 
# la riqueza de especies.
jasper.arbol<-read.tree("jasper_tree.phy")
jasper.dat<-read.csv("jasper_data.csv",row.names=1)

jasper.pd<-pd(jasper.dat,jasper.arbol)

plot(jasper.pd$SR,jasper.pd$PD,pch=19,xlab="Riqueza de especies"
     ,ylab="Diversidad filogenética")

# Vemos entonces que hay una relación positiva muy estrecha. A medida que incrementa
# la riqueza de especies, incrementa la diversidad filogenética.

# Que información adicional entonces nos puede dar la diversidad filogenética que 
# no nos de la diversidad taxonómica?

# Otros indices de diversidad filogenética

# 1. Distancia promedio entre especies en la comunidad

# Volvamos al ejemplo hipotético

par(mar=c(0,0,0,0))
plot(arbol,direction="upwards",show.tip.label=FALSE,y.lim=c(0,5.2)
     ,x.lim=c(-1,Ntip(arbol)),edge.color=rep(c(my.cols[1],"black"),times=c(16,62-16))
     ,edge.width=rep(2:1,times=c(16,62-16)))
comm1<-especies.muestreo[comms[1,]>0]
points(match(comm1,arbol$tip.label),rep(5+offs[1],length(comm1)),pch=15
       ,col=my.cols[1],cex=2)

# Calculemos primero a mano para la comunidad 1
# Las filogenias son esencialmente representaciones de distancias entre especies
# Consecuentemente, una filogenia se puede transformar en una matriz de distancias.
# Para esto vamos a usar la función cophenetic.

arbol.dist<-cophenetic(arbol)

# Buscar entonces todos los posibles pares existentes en la communidad 1.
# Primero, creamos entonces una matriz de todos los posibles pares de especies
# que existen en la comunidad.

comm1.combn<-t(combn(comm1,2))
colnames(comm1.combn)<-c("Especie1","Especie2")

# Ahora busquemos la distancia por cada par de especies
comm1.dist<-rep(NA,nrow(comm1.combn))

for(i in 1:length(comm1.dist)){
  
  pos1<-which(rownames(arbol.dist)==comm1.combn[i,1])
  pos2<-which(colnames(arbol.dist)==comm1.combn[i,2])
  comm1.dist[i]<-arbol.dist[pos1,pos2]
  
}

mean(comm1.dist)

# Ejercicio.
# Calcular manualmente la distancia filogenética promedio para la comunidad Clump2


# El paquete picante tiene una función que permite calcular esta distancia para 
# muchas comunidades al mismo tiempo

comm.mpd<-mpd(comms,cophenetic(arbol),abundance.weighted=FALSE)

# 2. Distancia promedio entre parientes más cercanos.

par(mar=c(0,0,0,0))
plot(arbol,direction="upwards",show.tip.label=FALSE,y.lim=c(0,5.2)
     ,x.lim=c(-1,Ntip(arbol)),edge.color=rep(c(my.cols[1],"black"),times=c(16,62-16))
     ,edge.width=rep(2:1,times=c(16,62-16)))
comm1<-especies.muestreo[comms[1,]>0]
points(match(comm1,arbol$tip.label),rep(5+offs[1],length(comm1)),pch=15
       ,col=my.cols[1],cex=2)
tiplabels(frame = "none")

# Nuevamente, calculemoslo a mano
# Cortemos primero la matriz a las especies de la comunidad

col.pos<-match(comm1,colnames(arbol.dist))
row.pos<-match(comm1,rownames(arbol.dist))

comm1.arbol.dist<-arbol.dist[row.pos,col.pos]
diag(comm1.arbol.dist)<-NA

dtmc<-apply(comm1.arbol.dist,1,min,na.rm=TRUE)

mean(dtmc)

# Ejercicio
# Calcular manualmente la distancia filogenética promedio al taxón más cercano 
# para la comunidad Clump2

# Igual que para el mpd, para el mntd existe tambien una función precargada en el 
# paquete picante.

comm.mntd<-mntd(comms,arbol.dist)
comm.mntd

# Veamos para la comunidad real, la cual tiene más diversidad, como están relacionadas
# las diferentes medidas de diversidad. 

jasper.mpd<-mpd(jasper.dat,cophenetic(jasper.arbol))
jasper.mntd<-mntd(jasper.dat,cophenetic(jasper.arbol))

par(mfrow=c(3,2),mar=c(3.5,3.5,1,1),mgp=c(1.75,0.5,0))
plot(jasper.pd$SR,jasper.mpd,pch=19,xlab="Riqueza de Especies",ylab="mpd",bty="l")
plot(jasper.pd$SR,jasper.mntd,pch=19,xlab="Riqueza de Especies",ylab="mntd",bty="l")

plot(jasper.pd$PD,jasper.mpd,pch=19,xlab="PD",ylab="mpd",bty="l")
plot(jasper.pd$PD,jasper.mntd,pch=19,xlab="PD",ylab="mntd",bty="l")

plot(jasper.mntd,jasper.mpd,pch=19,xlab="mntd",ylab="mpd",bty="l")
plot(jasper.pd$SR,jasper.pd$PD,pch=19,xlab="Riqueza de Especies",ylab="PD",bty="l")

# Que podemos concluir de estos resultados?

# Que otras cosas nos podemos preguntar con estos datos?

# ¿Cual es el objetivo de tener todas estas medidas de diversidad?

# Podemos entonces empezar a hacer inferencias sobre los mecanismos que determinan
# el ensamblaje de las comunidades si hacemos algunos supuestos y si los patrones 
# observados difieren de lo que uno esperaría por azar.

# Que significa azar?

# Volvemos entonces al ejemplo con los datos de phylocom. Vamos a calcular la 
# diversidad filogenética que uno esperaría si simplemente las especies colonizaran
# las comunidades sin ningún mecanismo aparente. Solamente por puro azar.

# 

rnd.mpd<-rnd.mntd<-matrix(NA,ncol=6,nrow=1000,dimnames=list(1:1000,rownames(comms)))

for(i in 1:1000){
  
  ith.shuffle<-   sample(rownames(arbol.dist)) 
  ith.tree<-arbol.dist
  rownames(ith.tree)<-ith.shuffle
  colnames(ith.tree)<-ith.shuffle
  ith.mpd<-mpd(comms,ith.tree)
  ith.mntd<-mntd(comms,ith.tree)
  rnd.mpd[i,]<-ith.mpd
  rnd.mntd[i,]<-ith.mntd
  
}

par(mfrow=c(3,2))
for(i in 1:6){
  hist(rnd.mpd[,i],main=rownames(comms)[i],xlab="MPD",xlim=range(c(comm.mpd[i],rnd.mpd[,i])))
  abline(v=comm.mpd[i],col="red",lty=2)
}

# Interpretemos los resultados obtenidos en cada una de las comunidades.
# Que quieren decir estos resultados?

# De la misma forma, podemos calcular un tamaño de efecto estandarizado. El tamaño de efecto
# estandarizado simplemente nos dice que tan alejado esta nuestra observación de lo 
# esperado por azar. Para esto restamos el valor observado del promedio de la distribución 
# aleatoria. Luego estandarizamos por la desviación estandar de la distribución nula.

rnd.mpd.prom<-apply(rnd.mpd,2,mean)
rnd.mpd.sd<-apply(rnd.mpd,2,sd)
rnd.mntd.prom<-apply(rnd.mntd,2,mean)
rnd.mntd.sd<-apply(rnd.mntd,2,sd)

my.ses.mpd<-(comm.mpd-rnd.mpd.prom)/rnd.mpd.sd
my.ses.mntd<-(comm.mntd-rnd.mntd.prom)/rnd.mntd.sd

# Veamos los resultados graficamente

plot(1:6,my.ses.mpd,pch=19,xlab="Comunidad",ylab="Indice de diversidad filogenética",bty="l",xaxt="n",cex=2
     ,ylim=c(-11,3))
points(1:6,my.ses.mntd,pch=19,col="grey",cex=2)
abline(h=c(-1.96,0,1.96),lty=2,col="red")
axis(1,1:6,rownames(comms))

# Otra forma de calcular los resultados es por medio de las funciones pre establecidas
# en el paquete picante.
nri<-ses.mpd(samp=comms,arbol.dist,null.model="taxa.labels")
nti<-ses.mntd(samp=comms,arbol.dist,null.model="taxa.labels")

# Ejercicio.

# Calcular las métricas de diversidad filogenética para las comunidades de
# Jasper Ridge



# Sin embargo, esto es solamente una parte de la historia... Recordemos que el objetivo
# de utilizar este tipo de información es poder realizar inferencias sobre los mecanismos
# de ensamblaje de comunidades

# Recordemos que los patrones observados sobre la diversidad filogenética
#de las especies que ocupan un sitio dependen de dos procesos diferentes.

# 1.) Las características que permiten a las especies vivir en un ambiente particular
# 2.) Las tasas de evolución de los caracteres.

# En ese orden de ideas, tenemos que determinar la distribución de los caracteres
# en la comunidad y la forma como evolucionan cada uno de los caracteres.

# Para eso vamos a aprovechar todo el conocimiento que tenemos hasta el momento para
# poder calcular la estructura de los caracteres en la comunidad.

# Calculemos primero una matriz de distancias de los caracteres
# Podemos calcular una similitud por caracter y entender cuales son los caracteres que
# más influyen

caract.mpd<-list()

for(i in 1:ncol(caracteres)){
  
  ith.dist<-as.matrix(dist(caracteres[,i]))
  colnames(ith.dist)<-rownames(ith.dist)<-rownames(caracteres)
  ith.sesmpd<-ses.mpd(comms,ith.dist)
  caract.mpd[[i]]<-ith.sesmpd
}

# Miremos los resultados e interpretemos.

# Finalmente debemos calcular la señal filogenética que tiene una característica.
# La señal filogenética es simplemente la influencia que tienen las relaciones evolutivas
# sobre una característica. Según la teoría, se esepraría que especies más cercanamente
# emparentadas se parecieran más que especies mas lejanamente emparentadas. Esto
# quiere decir que la similitud en características entre dos especies es proporcional 
# tiempo de divergencia entre estas dos especies. Un proceso natural a través del cual 
# se puede estimar la evolución de los caracteres es denominado como movimiento browniano.
# Simon Blomberg, propuso una medida para determinar si las caracteristicas observadas
# evolucionan bajo un modelo de evolución browniana o mas o menos rapido que este modelo.
# Esta medida esta basada en la relación entre varianza observada de un caracter y la
# varianza esperada si los caracteres evolucionaran de forma neutral (browniana). K varia
# entre 0 e infinito. Un K menor que uno indica una evolución rápida de los caracteres, mientras
# que un K mayor que uno indica que los caracteres evolucionan más lentamente que 
# lo esperado por azar, algo que se ha denominado como conservatismo filogenético
# de nicho.

# Veamos entonces como se ven los caracteres en la filogenia.
traitA<-as.vector(caracteres[,1])
names(traitA)<-rownames(caracteres)

par(mfrow=c(2,1),mar=c(0,3,1,0))
plot(1:length(traitA),traitA,pch=19,axes=FALSE
     ,xlab="",ylab="")
axis(2)
par(mar=c(0,3,0,0))
plot(arbol,direction="upwards")

# Calculemos la señal filogenética
Kcalc(traitA,arbol)
phylosignal(traitA,arbol)
