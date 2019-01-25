
U_exp<-matrix(data=scan(file="U_exp.txt"),nc=220,nr=4)
time<-matrix(data=scan(file="time.txt"),nc=220,nr=4)
age<-scan(file="age.txt")
vol<-scan(file="vol.txt")

pacientes=220
Aj<-vector(mode="logical",length=pacientes)
ka_j<-vector(mode="logical",length=pacientes)
ke_j<-vector(mode="logical",length=pacientes)
Uij<-matrix(nr=4,nc=pacientes)

#Armo la función de likelihood
likelihood <- function(param){
  o1 = param[1]
  o2 = param[2]
  o3 = param[3]
  o4 = param[4]
  o5 = param[5]
  o6 = param[6]
  o7 = param[7]
  
  Aj<- o1 + o2*vol
  for(j in 1:pacientes){
    ka_j[j]<- o3 + o4*vol[j] + o5*age[j]
  }
  for(j in 1:pacientes){
    ke_j[j]<- o6+ o7*age[j] 
  }
  for(j in 1:pacientes){
    Uij[,j]<-Aj[j]*(exp(-ke_j[j]*time[,j])-exp(-ka_j[j]*time[,j]))
  }
  pred = as.vector(Uij)
  cadalikelihood = dnorm(as.vector(U_exp), mean = pred, sd = 1, log = T)
  cadalikelihood<-replace(cadalikelihood,is.na(cadalikelihood),0)
  sumll = sum(cadalikelihood)
  return(sumll)   
}


#Armo la probabilidad a priori
priori <- function(param){
  o1 = param[1]
  o2 = param[2]
  o3 = param[3]
  o4 = param[4]
  o5 = param[5]
  o6 = param[6]
  o7 = param[7]

  o1priori = dunif(o1, min=0, max=100, log = T)
  o2priori = dunif(o2, min=0,max=1, log = T)
  o3priori = dunif(o3, min=0, max=1, log = T)
  o4priori = dunif(o4, min=0,max=1, log = T)
  o5priori = dunif(o5, min=-1, max=1, log = T)
  o6priori = dunif(o6, min=0,max=1, log = T)
  o7priori = dunif(o7, min=-1, max=1, log = T)

  return(o1priori+o2priori+o2priori+o4priori+o5priori+o6priori+o7priori)
}




posteriori <- function(param){
  return (likelihood(param) + priori(param))
}



proposicion <- function(param){
  return(rnorm(7,mean = param, sd= c(0.15,0.0039,0.0056,0.00014,6.4e-5,4.7e-5,7.9e-7)))
}

metropolis_MCMC <- function(valor_0, iteraciones){
  cadena = array(dim = c(iteraciones+1,7))
  cadena[1,] = valor_0
  for (i in 1:iteraciones){
    proposicion_i = proposicion(cadena[i,])
    
    probab = exp(posteriori(proposicion_i) - posteriori(cadena[i,]))
    if (runif(1) < probab){
      cadena[i+1,] = proposicion_i
    }else{
      cadena[i+1,] = cadena[i,]
    }
  }
  return(cadena)
}

valor_0 = c(60,0.17,0.5,0.003,-0.004,0.004,-3E-5)
iteraciones=10^6
cadena = metropolis_MCMC(valor_0, iteraciones)
write(cadena,file = "Cadena.txt",sep ="\n")

jpeg('cadena_o1.jpg')
plot(seq(0,iteraciones,1),cadena[,1],type="l",xlab="Iteraciones",ylab = "o1")
dev.off()
jpeg('cadena_o2.jpg')
plot(seq(0,iteraciones,1),cadena[,2],type="l",xlab="Iteraciones",ylab = "o2")
dev.off()
jpeg('cadena_o3.jpg')
plot(seq(0,iteraciones,1),cadena[,3],type="l",xlab="Iteraciones",ylab = "o3")
dev.off()
jpeg('cadena_o4.jpg')
plot(seq(0,iteraciones,1),cadena[,4],type="l",xlab="Iteraciones",ylab = "o4")
dev.off()
jpeg('cadena_o5.jpg')
plot(seq(0,iteraciones,1),cadena[,5],type="l",xlab="Iteraciones",ylab = "o5")
dev.off()
jpeg('cadena_o6.jpg')
plot(seq(0,iteraciones,1),cadena[,6],type="l",xlab="Iteraciones",ylab = "o6")
dev.off()
jpeg('cadena_o7.jpg')
plot(seq(0,iteraciones,1),cadena[,7],type="l",xlab="Iteraciones",ylab = "o7")
dev.off()

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

write(acf1[[1]],file = "acf1.txt",sep ="\n")
write(acf2[[1]],file = "acf2.txt",sep ="\n")
write(acf3[[1]],file = "acf3.txt",sep ="\n")
write(acf4[[1]],file = "acf4.txt",sep ="\n")
write(acf5[[1]],file = "acf5.txt",sep ="\n")
write(acf6[[1]],file = "acf6.txt",sep ="\n")
write(acf7[[1]],file = "acf7.txt",sep ="\n")

burnIn = 2.3e3

o1_final<-mean(cadena[-(1:burnIn),1])
o2_final<-mean(cadena[-(1:burnIn),2])
o3_final<-mean(cadena[-(1:burnIn),3])
o4_final<-mean(cadena[-(1:burnIn),4])
o5_final<-mean(cadena[-(1:burnIn),5])
o6_final<-mean(cadena[-(1:burnIn),6])
o7_final<-mean(cadena[-(1:burnIn),7])
O_final<-c(o1_final,o2_final,o3_final,o4_final,o5_final,o6_final,o7_final)
write(O_final,file = "o_final.txt",sep ="\n")

sd_o1<-sd(cadena[-(1:burnIn),1])
sd_o2<-sd(cadena[-(1:burnIn),2])
sd_o3<-sd(cadena[-(1:burnIn),3])
sd_o4<-sd(cadena[-(1:burnIn),4])
sd_o5<-sd(cadena[-(1:burnIn),5])
sd_o6<-sd(cadena[-(1:burnIn),6])
sd_o7<-sd(cadena[-(1:burnIn),7])


jpeg('dist_o1.jpg')
hist(cadena[-(1:burnIn),1],nclass=40, main="Distribución de o1", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o1_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o1_final, sd=sd_o1),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
jpeg('dist_o2.jpg')
hist(cadena[-(1:burnIn),2],nclass=40, main="Distribución de o2", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o2_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o2_final, sd=sd_o2),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
jpeg('dist_o3.jpg')
hist(cadena[-(1:burnIn),3],nclass=40, main="Distribución de o3", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o3_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o3_final, sd=sd_o3),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
jpeg('dist_o4.jpg')
hist(cadena[-(1:burnIn),4],nclass=40, main="Distribución de o4", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o4_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o4_final, sd=sd_o4),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
jpeg('dist_o5.jpg')
hist(cadena[-(1:burnIn),5],nclass=40, main="Distribución de o5", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o5_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o5_final, sd=sd_o5),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
jpeg('dist_o6.jpg')
hist(cadena[-(1:burnIn),6],nclass=40, main="Distribución de o6", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o6_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o6_final, sd=sd_o6),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
jpeg('dist_o7.jpg')
hist(cadena[-(1:burnIn),7],nclass=40, main="Distribución de o7", prob=TRUE, xlab="Valores",ylab = "Densidad")
abline(v = o7_final ,col="blue", lwd=2)
curve(dnorm(x, mean=o7_final, sd=sd_o7),col="red", lwd=2, add=TRUE, yaxt="n", lty="dotted")
dev.off()
