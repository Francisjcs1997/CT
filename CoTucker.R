library(CA3variants)
#####
#TRATAMENTO DE DADOS

#Extração de dados
saliva=read.csv2("saliva.csv")
rownames(saliva)=saliva$ï..code
saliva=saliva[,-1]
colnames(saliva)=c("etanol","NAG","acetoina","ureia")

urina=read.csv2("urina.csv")
rownames(urina)=urina$ï..code
urina=urina[,-1]
colnames(urina)=c("Pn3G","P3G","X4.DEA","X3.HIBA","alanina","X4.DTA","colina","carnitina",
                  "glicina","creatina","creatinina","glicose","hipurico","GAA","X1.6.anidro",
                  "treonina","trigonelina","histidina","taurina","X2.KG","U3","U1","N5AC","U2")

#Normalização
saliva1=scale(saliva[1:7,],center=T,scale=F)
saliva2=scale(saliva[8:14,],center=T,scale=F)
saliva3=scale(saliva[15:21,],center=T,scale=F)
saliva_norm=rbind(saliva1,saliva2,saliva3)
saliva_norm=data.frame(saliva_norm)

urina1=scale(urina[1:7,],center=T,scale=T)
urina2=scale(urina[8:14,],center=T,scale=T)
urina3=scale(urina[15:21,],center=T,scale=T)
urina_norm=rbind(urina1,urina2,urina3)
urina_norm=data.frame(urina_norm)

#remove(saliva,urina)


#####
#CO-INERCIA

Dn=as.matrix(diag(1/7,nrow=7,ncol=7))
Z1=t(saliva1)%*%Dn%*%urina1
Z2=t(saliva2)%*%Dn%*%urina2
Z3=t(saliva3)%*%Dn%*%urina3

CO_I=rbind(Z1,Z2,Z3)
CO_I=data.frame(CO_I)


A=sum(diag(saliva1%*%t(saliva1)%*%Dn%*%urina1%*%t(urina1)%*%Dn))
B=sum(diag(saliva2%*%t(saliva2)%*%Dn%*%urina2%*%t(urina2)%*%Dn))
C=sum(diag(saliva3%*%t(saliva3)%*%Dn%*%urina3%*%t(urina3)%*%Dn))
barplot(c(A,B,C),names.arg = c("1ºT","2ºT","3ºT"),ylab="Interação")



#####
#ESCOLHA E CONSTRUÇÃO DO MODELO

#1º Modo
Xa=as.matrix(cbind(CO_I[1:4,],CO_I[5:8,],CO_I[9:12,]))

#2º Modo
Xb=t(CO_I)

#3º Modo
l1=as.numeric(t(CO_I[1:4,]));l2=as.numeric(t(CO_I[5:8,]));
l3=as.numeric(t(CO_I[9:12,]))
Xc=rbind(l1,l2,l3)


U=tucker(array(Xa,dim=c(4,24,3)),2,2,2)
A=U$a;B=U$b;C=U$cc
rownames(A)=colnames(saliva1)
rownames(B)=colnames(urina1)
rownames(C)=c("1T","2T","3T")
G=array(U$g,dim=c(2,2,2))
perc=sum(G^2)/sum(Xa^2)


#####
#VISUALIZAÇÃO DE COMPONENTES
par(mfrow=c(1,3))
#1º Modo
plot(NA,xlim=c(min(A[,1],0),max(A[,1],0)),ylim=c(min(A[,2]-0.1,0),max(A[,2]+0.1,0)),
     xlab="P1",ylab="P2",main="P")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
arrows(0,0,A[,1],A[,2],length = 0.1,lwd=2)
text(A[1,1]-0.05,A[1,2]-0.02,labels = rownames(A)[1],cex=2)#etanol
text(A[2,1]+0.08,A[2,2]+0.06,labels = rownames(A)[2],cex=2)#NAG
text(A[3,1]-0.18,A[3,2],labels = rownames(A)[3],cex=2)#acetoine
text(A[4,1]-0.05,A[4,2]+0.05,labels = rownames(A)[4],cex=2)#Urea

#2ª Modo
plot(NA,xlim=c(min(B[,1]-0.1,0),max(B[,1],0)),ylim=c(min(B[,2],0),max(B[,2],0)),
     xlab="Q1",ylab="Q2",
     main="Q",col="red")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(B[,1],B[,2],labels = rownames(B),cex=2, pos=3,col="red")

#3º Modo
plot(C[,1],C[,2],xlim=c(min(C[,1],0),max(C[,1],0)),ylim=c(min(C[,2],0),max(C[,2],0)),
     xlab="R1",ylab="R2",
     main="R",col="blue")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(C[1,1],C[1,2]+0.05,labels = rownames(C)[1],cex=2,col="blue")#1T
text(C[2,1],C[2,2]+0.05,labels = rownames(C)[2],cex=2,col="blue")#2T
text(C[3,1],C[3,2]-0.05,labels = rownames(C)[3],cex=2,col="blue")#3T
par(mfrow=c(1,1))
#####
#CONSTRUÇÃO DOS BIPLOTS

#R1
G1=svd(G[,,1])
Ax1=((4/24)^(1/4))*A%*%G1$u%*%diag(sqrt(G1$d))
Bx1=((24/4)^(1/4))*B%*%G1$v%*%diag(sqrt(G1$d))

xmin=min(Ax1[,1],Bx1[,1]);xmax=max(Ax1[,1],Bx1[,1])
ymin=min(Ax1[,2],Bx1[,2]);ymax=max(Ax1[,2],Bx1[,2])

plot(NA,col="red",xlim=c(xmin-0.02,xmax+0.01),ylim=c(ymin,ymax+0.005),
     main="R1",xlab="Primeira Componente",ylab="Segunda Componente")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(Bx1[,1], Bx1[,2], labels=rownames(Bx1), cex= 1.2, pos=3,col="red")
arrows(0,0,Ax1[,1],Ax1[,2],length = 0.1,lwd=2)
text(Ax1[,1], Ax1[,2], labels=rownames(Ax1), cex= 1.3, pos=3) 


#R2
G2=svd(G[,,2])
Ax2=((2/24)^(1/4))*A%*%G2$u%*%diag(sqrt(G2$d))
Bx2=((24/2)^(1/4))*B%*%G2$v%*%diag(sqrt(G2$d))

xmin=min(Ax2[,1],Bx2[,1]);xmax=max(Ax2[,1],Bx2[,1])
ymin=min(Ax2[,2],Bx2[,2]);ymax=max(Ax2[,2],Bx2[,2])

#2º componente do 3º modo
plot(NA,col="red",xlim=c(xmin-0.005,xmax),ylim=c(ymin,ymax+0.005),
     main="R2",xlab="Primeira Componente",ylab="Segunda Componente")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(Bx2[,1], Bx2[,2], labels=rownames(Bx2), cex= 1.2, pos=3,col="red")
arrows(0,0,Ax2[,1],Ax2[,2],length = 0.1,lwd=2)
text(Ax2[,1], Ax2[,2], labels=rownames(Ax2), cex= 1.2, pos=3)