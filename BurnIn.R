jpeg('acf1.jpg')
acf1<-acf(cadena[,1],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O1",xlab="Iteraciones")
dev.off()

jpeg('acf2.jpg')
acf2<-acf(cadena[,2],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O2",xlab="Iteraciones")
dev.off()

jpeg('acf3.jpg')
acf3<-acf(cadena[,3],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O3",xlab="Iteraciones")
dev.off()

jpeg('acf4.jpg')
acf4<-acf(cadena[,4],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O4",xlab="Iteraciones")
dev.off()

jpeg('acf5.jpg')
acf5<-acf(cadena[,5],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O5",xlab="Iteraciones")
dev.off()

jpeg('acf6.jpg')
acf6<-acf(cadena[,6],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O6",xlab="Iteraciones")
dev.off()

jpeg('acf7.jpg')
acf7<-acf(cadena[,7],lag.max =iteraciones,plot=TRUE,main="Autocorrelación cadena O7",xlab="Iteraciones")
dev.off()

plot(acf1[[4]],acf1[[1]])
abline(h=0,col="red")
