#####
#TRATAMENTO DE DADOS
#(Dados não disponibilizados)

#Tratamento de dados
saliva1=as.matrix(saliva[8:14,]-saliva[1:7,]); #X2*-X1*
saliva1=scale(saliva1,center=T,scale=F);
rownames(saliva1)=nomes1
saliva2=as.matrix(saliva[15:21,]-saliva[1:7,]); #X3*-X1*
saliva2=scale(saliva2,center=T,scale=F);
rownames(saliva2)=nomes2
saliva3=as.matrix(saliva[15:21,]-saliva[8:14,]); #X3*-X2*
saliva3=scale(saliva3,center=T,scale=F);
rownames(saliva3)=nomes3

urina1=as.matrix(urina[8:14,]-urina[1:7,]); #Y2*-Y1*
urina1=scale(urina1,center=T,scale=T);
rownames(urina1)=nomes1;
urina2=as.matrix(urina[15:21,]-urina[1:7,]); #Y3*-Y1*
urina2=scale(urina2,center=T,scale=T);
rownames(urina2)=nomes2
urina3=as.matrix(urina[15:21,]-urina[8:14,]); #Y3*-Y2*
urina3=scale(urina3,center=T,scale=T);
rownames(urina3)=nomes3

saliva_norm=as.data.frame(rbind(saliva1,saliva2,saliva3))
urina_norm=as.data.frame(rbind(urina1,urina2,urina3))


#####
#CO-INERCIA

Dn=as.matrix(diag(1/7,nrow=7,ncol=7))
Z1=t(saliva1)%*%Dn%*%as.matrix(urina1)
Z2=t(saliva2)%*%Dn%*%as.matrix(urina2)
Z3=t(saliva3)%*%Dn%*%as.matrix(urina3)

CO_I=rbind(Z1,Z2,Z3)
CO_I=data.frame(CO_I)

#calculo co-inércia
A=sum(diag(as.matrix(saliva1)%*%t(saliva1)%*%Dn%*%as.matrix(urina1)%*%t(urina1)%*%Dn))
B=sum(diag(as.matrix(saliva2)%*%t(saliva2)%*%Dn%*%as.matrix(urina2)%*%t(urina2)%*%Dn))
C=sum(diag(as.matrix(saliva3)%*%t(saliva3)%*%Dn%*%as.matrix(urina3)%*%t(urina3)%*%Dn))
barplot(c(A,B,C),names.arg = c("1C","2C","3C"))
#
saliva1=data.frame(saliva1)
saliva2=data.frame(saliva2)
saliva3=data.frame(saliva3)
urina1=data.frame(urina1)
urina2=data.frame(urina2)
urina3=data.frame(urina3)

remove(Dn,Z1,Z2,Z3,A,B,C)

#####
#ALS em https://github.com/Francisjcs1997/CT/blob/main/ALS.R

#####
#VISUALIZAÇÃO DE COMPONENTES

#1º Modo
plot(NA,xlim=c(min(A[,1]-0.05,0),max(A[,1]+0.05,0)),ylim=c(min(A[,2]-0.1,0),max(A[,2]+0.1,0)),
     xlab="P1",ylab="P2",main="P")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
arrows(0,0,A[,1],A[,2],length = 0.1,lwd=2)
text(A[1,1],A[1,2]-0.05,labels = rownames(A)[1],cex=1.5)#etanol
text(A[2,1],A[2,2]-0.05,labels = rownames(A)[2],cex=1.5)#NAG
text(A[3,1]+0.1,A[3,2]+0.05,labels = rownames(A)[3],cex=1.5)#aceine
text(A[4,1],A[4,2]-0.05,labels = rownames(A)[4],cex=1.5)#urea

#2ª Modo
plot(NA,xlim=c(min(B[,1]-0.1,0),max(B[,1]+0.1,0)),ylim=c(min(B[,2],0),max(B[,2]+0.1,0)),
     xlab="Q1",ylab="Q2",
     main="Q",col="red")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(B[,1],B[,2],labels = rownames(B),cex=1.5, pos=3,col="red")


#####
#CONSTRUÇÃO DOS BIPLOTS
#R1
G1=svd(G[,,1])
Ax1=((4/24)^(1/4))*A%*%G1$u%*%diag(sqrt(G1$d))
Bx1=((24/4)^(1/4))*B%*%G1$v%*%diag(sqrt(G1$d))

xmin=min(Ax1[,1],Bx1[,1]);xmax=max(Ax1[,1],Bx1[,1])
ymin=min(Ax1[,2],Bx1[,2]);ymax=max(Ax1[,2],Bx1[,2])

plot(NA,col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
     main="R1",xlab="Primeira Componente",ylab="Segunda Componente")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(Bx1[,1], Bx1[,2], labels=rownames(Bx1), cex= 1.2, pos=3,col="red")
arrows(0,0,Ax1[,1],Ax1[,2],length = 0.1,lwd=2)
text(Ax1[,1], Ax1[,2], labels=rownames(Ax1), cex= 1.3, pos=3) 
