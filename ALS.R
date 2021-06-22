#Depois do tratamento de dados realizado e calculo da CO_I:
c1=2;c2=2;c3=2 #número P,Q,R, respetivamente

Xa=as.matrix(cbind(CO_I[1:4,],CO_I[5:8,],CO_I[9:12,]))
Xb=t(as.matrix(CO_I))
l1=as.numeric(t(CO_I[1:4,]));l2=as.numeric(t(CO_I[5:8,]));l3=as.numeric(t(CO_I[9:12,]));
Xc=rbind(l1,l2,l3)

########################
#Inicialização 1 com coeficientes nulos
b=matrix(0,nrow=24,ncol=c2)
c=matrix(0,nrow=3,ncol=c3)

#Inicialização 2 com Wb Wb^T
decomp_Xb=svd(Xb%*%t(Xb))
b=as.matrix(decomp_Xb$u[,1:c2])
decomp_Xc=svd(Xc%*%t(Xc))
c=as.matrix(decomp_Xc$u[,1:c3])

erro=FALSE
control=1000000
k=0

while(erro==FALSE){
  k=k+1
  M1=Xa%*%(c%x%b)
  a=svd(M1)$u[,1:c1]
  
  M2=Xb%*%(c%x%a)
  b=svd(M2)$u[,1:c2]
  
  M3=Xc%*%(a%x%b)
  c=svd(M3)$u[,1:c3]
  
  A1=Xa-a%*%t(a)%*%Xa%*%(c%x%b)%*%(t(c)%x%t(b))
  G1=array(t(a)%*%Xa%*%(c%x%b),dim=c(c1,c2,c3))
  perc1=sum(G1^2)/sum(Xa^2)
  if(control-norm(A1,"F")^2<10^-6){
    A=a
    B=b
    C=c
    erro=TRUE
  }
  else{
    A=a
    B=b
    C=c
    control=norm(A1,"F")^2
  }
  print(k)
  print(norm(A1,"F")^2)
  print(round(perc1,3))
}

G=array(t(A)%*%Xa%*%(C%x%B),dim=c(c1,c2,c3))
perc=sum(G^2)/sum(Xa^2)



#ALS com função Tucker do R
library(CA3variants)

CO_I=as.matrix(CO_I)
W=array(c(CO_I[1:4,],CO_I[5:8,],CO_I[9:12,]),dim=c(4,24,3))
U=tucker(W,c1,c2,c3)

A=U$a;B=U$b;C=U$cc
G=array(U$g,dim=c(c1,c2,c3))
perc=sum(G^2)/sum(Xa^2)
