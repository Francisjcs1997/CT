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

nomes1=c("1_Cond1","2_Cond1","3_Cond1","4_Cond1","5_Cond1","6_Cond1","7_Cond1")
nomes2=c("1_Cond2","2_Cond2","3_Cond2","4_Cond2","5_Cond2","6_Cond2","7_Cond2")
nomes3=c("1_Cond3","2_Cond3","3_Cond3","4_Cond3","5_Cond3","6_Cond3","7_Cond3")

saliva1=as.matrix(saliva[8:14,]-saliva[1:7,]);
saliva1=scale(saliva1,center=T,scale=F);
rownames(saliva1)=nomes1
saliva2=as.matrix(saliva[15:21,]-saliva[1:7,]);
saliva2=scale(saliva2,center=T,scale=F);
rownames(saliva2)=nomes2
saliva3=as.matrix(saliva[15:21,]-saliva[8:14,]);
saliva3=scale(saliva3,center=T,scale=F);
rownames(saliva3)=nomes3

urina1=as.matrix(urina[8:14,]-urina[1:7,]);
urina1=scale(urina1,center=T,scale=T);
rownames(urina1)=nomes1;
urina2=as.matrix(urina[15:21,]-urina[1:7,]);
urina2=scale(urina2,center=T,scale=T);
rownames(urina2)=nomes2
urina3=as.matrix(urina[15:21,]-urina[8:14,]);
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

#
A=sum(diag(as.matrix(saliva1)%*%t(saliva1)%*%Dn%*%as.matrix(urina1)%*%t(urina1)%*%Dn))
B=sum(diag(as.matrix(saliva2)%*%t(saliva2)%*%Dn%*%as.matrix(urina2)%*%t(urina2)%*%Dn))
C=sum(diag(as.matrix(saliva3)%*%t(saliva3)%*%Dn%*%as.matrix(urina3)%*%t(urina3)%*%Dn))
#barplot(c(A,B,C),names.arg = c("1C","2C","3C"))
#
saliva1=data.frame(saliva1)
saliva2=data.frame(saliva2)
saliva3=data.frame(saliva3)
urina1=data.frame(urina1)
urina2=data.frame(urina2)
urina3=data.frame(urina3)

remove(Dn,Z1,Z2,Z3,A,B,C)

#####
#ESCOLHA E CONSTRUÇÃO DO MODELO

#1º Modo
Xa=as.matrix(cbind(CO_I[1:4,],CO_I[5:8,],CO_I[9:12,]))
decomp_Xa=svd(Xa%*%t(Xa))#diag(1/ncol(Xa),nrow=ncol(Xa),ncol=ncol(Xa))%*%

A=as.matrix(decomp_Xa$u[,1:2])
#rownames(A)=rownames(CO_I)[1:4];colnames(A)=c("P1","P2")

#2º Modo
Xb=t(as.matrix(CO_I))
decomp_Xb=svd(Xb%*%t(Xb))#diag(1/ncol(Xb),nrow=ncol(Xb),ncol=ncol(Xb))%*%


B=as.matrix(decomp_Xb$u[,1:2])
#rownames(B)=colnames(CO_I);colnames(B)=c("Q1","Q2")

#3º Modo
l1=as.numeric(t(CO_I[1:4,]));l2=as.numeric(t(CO_I[5:8,]));l3=as.numeric(t(CO_I[9:12,]))
Xc=rbind(l1,l2,l3)
decomp_Xc=svd(Xc%*%t(Xc))#diag(1/ncol(Xc),nrow=ncol(Xc),ncol=ncol(Xc))%*%

C=as.matrix(decomp_Xc$u[,1:2])
#rownames(C)=c("C1","C2","C3");colnames(C)=c("R1","R2")

#Cubo G
a=2;b=2;c=2
Ga=t(A)%*%Xa%*%(C%x%B)
G=array(c(Ga[,1:(b*c)]),dim=c(a,b,c))
G_prop=(G^2)/(sum(G^2))

perc=sum(G^2)/sum(CO_I^2)
G_prop=round(G_prop*perc*100,3)

#remove(decomp_Xa,decomp_Xb,decomp_Xc,Ga,Xa,Xb,Xc,a,b,c,l1,l2)

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

#plot(NA,col="red",xlim=c(xmin-0.005,xmax+0.005),ylim=c(ymin-0.05,ymax+0.02),
    # main="R1",xlab="Primeira Componente",ylab="Segunda Componente")
plot(NA,col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
     main="R1",xlab="Primeira Componente",ylab="Segunda Componente")
abline(v=0, col='black', lty=10, lwd=0.3)
abline(h=0, col='black', lty=10, lwd=0.3)
text(Bx1[,1], Bx1[,2], labels=rownames(Bx1), cex= 1.2, pos=3,col="red")
arrows(0,0,Ax1[,1],Ax1[,2],length = 0.1,lwd=2)
#text(Ax1[c(3,4),1], Ax1[c(3,4),2], labels=rownames(Ax1)[c(3,4)], cex= 1.3, pos=3) #ureia e acetonina
#text(Ax1[c(1,2),1], Ax1[c(1,2),2]-0.12, labels=rownames(Ax1)[c(1,2)], cex= 1.3, pos=3) 
text(Ax1[,1], Ax1[,2], labels=rownames(Ax1), cex= 1.3, pos=3) 



#####
sp=function(i,j){
  return(as.numeric(Ax1[i,1]*Bx1[j,1]+Ax1[i,2]*Bx1[j,2]))
}


for (i in 2:2){
  vec=numeric(24)
  for (j in 1:24){
    vec[j]=sp(i,j)
  }
  print(rownames(A)[i])
  print("Interacao Positiva")
  M=cbind(rownames(B)[order(vec,decreasing = T)],round(vec[order(vec,decreasing = T)],3))
  linha=1
  while(abs(as.numeric(M[linha,2]))>0.6){
    print(c(M[linha,1],M[linha,2]))
    print(CO_I[c(i,i+4,i+8),c(which(rownames(B)==as.character(M[linha,1])))])
    linha=linha+1
  }
  print("Interacao Negativa")
  M=cbind(rownames(B)[order(vec,decreasing = F)],round(vec[order(vec,decreasing = F)],3))
  linha=1
  while(abs(as.numeric(M[linha,2]))>0.6){
    print(c(M[linha,1],M[linha,2]))
    print(CO_I[c(i,i+4,i+8),c(which(rownames(B)==as.character(M[linha,1])))])
    linha=linha+1
  }
  print("--------------------------------------")
}


######
model=array(NA,dim=c(4,24,3))

d=0
for (i in 1:4){
  for (j in 1:24){
    for (k in 1:3){
      for (p in 1:2){
        for (q in 1:2){
          for (r in 1:2){
            d=d+A[i,p]*B[j,q]*C[k,r]*G[p,q,r]
          }
        }
      }
      model[i,j,k]=d
      d=0
    }
  }
}

model=A2%*%cbind(G2[,,1],G2[,,2],G2[,,3])%*%(t(C2)%x%t(B2))

(sum(model^2)/sum(W^2))

A=U$a;B=U$b;C=U$cc;G=array(U$g,dim=c(2,2,2))
A1=cbind(W[,,1],W[,,2],W[,,3])-model
A1=cbind(A1[,,1],A1[,,2],A1[,,3])

#eigen(t(A1)%*%A1)$values[1] #39.28193
norm(A1,"F")^2
 
#29.01   #76.53
#28.58   #77.37
#22.60
#19.02
G=t(A)%*%Xa%*%(C%x%B)

res=tucker(W,p=2,q=2,r=2)
sum(U$g^2)/sum(W^2)


G2=array(NA,dim=c(2,2,2))
d=0
for (p in 1:2){
  for (q in 1:2){
    for (r in 1:2){
      for (i in 1:4){
        for (j in 1:24){
          for (k in 1:3){
            d=d+A[i,p]*B[j,q]*C[k,r]*W[i,j,k]
          }
        }
      }
      G2[p,q,r]=d
      d=0
    }
  }
}
G=G2



O=as.numeric(which(Bx2[,1]<0 & Bx2[,2]<0))
CO_I[c(9,12),O]

######
#Gráficos

#R1
#acetoina 1 
par(mfrow=c(1,3))
a=saliva_norm$acetoina;b=urina_norm$X4.DTA
plot(1:7,a[1:7],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="1ºC")
lines(1:7,b[1:7],type="l",lty=10,col="blue")
legend("bottomleft", legend=c("acetoina", "X4.DTA"),
       col=c("black", "blue"), lty=1:2, cex=0.8)
plot(1:7,a[8:14],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="2ºC")
lines(1:7,b[8:14],type="l",lty=10,col="blue")
plot(1:7,a[15:21],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="3ºC")
lines(1:7,b[15:21],type="l",lty=10,col="blue")
par(mfrow=c(1,1))

#ureia
par(mfrow=c(1,3))
a=saliva_norm$ureia;b=urina_norm$X4.DTA
plot(1:7,a[1:7],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="1ºC")
lines(1:7,b[1:7],type="l",lty=10,col="blue")
legend("bottomleft", legend=c("ureia", "X4.DTA"),
       col=c("black", "blue"), lty=1:2, cex=0.8)
plot(1:7,a[8:14],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="2ºC")
lines(1:7,b[8:14],type="l",lty=10,col="blue")
plot(1:7,a[15:21],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="3ºC")
lines(1:7,b[15:21],type="l",lty=10,col="blue")
par(mfrow=c(1,1))

par(mfrow=c(1,3))
a=saliva_norm$NAG;b=urina_norm$N5AC
plot(1:7,a[1:7],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="1ºC")
lines(1:7,b[1:7],type="l",lty=10,col="blue")
legend("bottomleft", legend=c("NAG", "N5AC"),
       col=c("black", "blue"), lty=1:2, cex=0.8)
plot(1:7,a[8:14],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="2ºC")
lines(1:7,b[8:14],type="l",lty=10,col="blue")
plot(1:7,a[15:21],type="l",col="black",lwd=2,ylab="Metabolitos",xlab="Grávidas",
     ylim=c(min(a,b),max(a,b)),main="3ºC")
lines(1:7,b[15:21],type="l",lty=10,col="blue")
par(mfrow=c(1,1))
