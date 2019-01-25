#Evaluación de la función
U_exp<-matrix(data=scan(file="U_exp.txt"),nc=220,nr=4)
time<-matrix(data=scan(file="time.txt"),nc=220,nr=4)
age<-scan(file="age.txt")
vol<-scan(file="vol.txt")
mage<-mean(age)
mvol<-mean(vol)

ka<- o3_final + o4_final*mvol + o5_final*mage
ke<- o6_final+ o7_final*mage
A<- o1_final + o2_final*mvol
tim<-seq(0,150,1)
U<-A*(exp(-ke*tim)-exp(-ka*tim))
jpeg('Promedio.jpg')
plot(tim, U, type="l", col="blue", lwd=3, ylim=c(0, 100),xlab="Tiempo [s]", ylab="% Absorción" )
points(time, U_exp, pch=8)
dev.off()

j=8
Aj<- o1_final + o2_final*vol[j]
ka_j<- o3_final + o4_final*vol[j] + o5_final*age[j]
ke_j<- o6_final+ o7_final*age[j]
Uj<-Aj*(exp(-ke_j*tim)-exp(-ka_j*tim))
jpeg('Caso8.jpg')
plot(tim, Uj, type="l", col="blue", lwd=3, ylim=c(0, 100),xlab="Tiempo [s]", ylab="% Absorción" )
points(time[,j], U_exp[,j], pch=8)
dev.off()


