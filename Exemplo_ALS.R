library(ade4)
library(CA3variants)

#dados
data(meau)

dados1=meau$env
dados2=meau$spe

#centrar matrizes
env1=scale(dados1[1:6,],center=T,scale=F)
env2=scale(dados1[7:12,],center=T,scale=F)
env3=scale(dados1[13:18,],center=T,scale=F)
env4=scale(dados1[19:24,],center=T,scale=F)

spe1=scale(dados2[1:6,],center=T,scale=F)
spe2=scale(dados2[7:12,],center=T,scale=F)
spe3=scale(dados2[13:18,],center=T,scale=F)
spe4=scale(dados2[19:24,],center=T,scale=F)

#Co-in�rcia
T1=t(env1)%*%diag(1/6,nrow=6,ncol=6)%*%spe1
T2=t(env2)%*%diag(1/6,nrow=6,ncol=6)%*%spe2
T3=t(env3)%*%diag(1/6,nrow=6,ncol=6)%*%spe3
T4=t(env4)%*%diag(1/6,nrow=6,ncol=6)%*%spe4

CO_I=rbind(T1,T2,T3,T4)


Xa=as.matrix(cbind(CO_I[1:10,],CO_I[11:20,],CO_I[21:30,],CO_I[31:40,]))
Xb=t(as.matrix(CO_I))
l1=as.numeric(t(CO_I[1:10,]));l2=as.numeric(t(CO_I[11:20,]));l3=as.numeric(t(CO_I[21:30,]));
l4=as.numeric(t(CO_I[31:40,]))
Xc=rbind(l1,l2,l3,l4)


Final=matrix(NA,nrow=278,ncol=13)
linha=0

for (a1 in 1:10){
  for (b1 in 1:13){
    for (c1 in 1:4){
      if (a1<=b1*c1 & b1<=a1*c1 & c1<=b1*a1){
        linha=linha+1
        
        #1 inicializa��o com zeros
        b=matrix(0,nrow=13,ncol=b1)
        c=matrix(0,nrow=4,ncol=c1)
        
        #1
        erro=FALSE
        control=1000000
        k=0
        
        while(erro==FALSE){
          k=k+1
          M1=Xa%*%(c%x%b)
          a=svd(M1)$u[,1:a1]
          
          M2=Xb%*%(c%x%a)
          b=svd(M2)$u[,1:b1]
          
          M3=Xc%*%(a%x%b)
          c=svd(M3)$u[,1:c1]
          
          A1=Xa-a%*%t(a)%*%Xa%*%(c%x%b)%*%(t(c)%x%t(b))
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
        }
        G=array(t(A)%*%Xa%*%(C%x%B),dim=c(a1,b1,c1))
        Final[linha,1]=a1
        Final[linha,2]=b1
        Final[linha,3]=c1
        Final[linha,4]=a1+b1+c1
        Final[linha,5]=k
        Final[linha,6]=round(control,2)
        Final[linha,7]=round(sum(G^2)/sum(Xa^2),2)
        
        #2 inicializa��o com SVD
        decomp_Xb=svd(Xb%*%t(Xb))
        b=as.matrix(decomp_Xb$u[,1:b1])
        
        decomp_Xc=svd(Xc%*%t(Xc))
        c=as.matrix(decomp_Xc$u[,1:c1])
        
        erro=FALSE
        control=1000000
        k=0
        
        while(erro==FALSE){
          k=k+1
          M1=Xa%*%(c%x%b)
          a=svd(M1)$u[,1:a1]
          
          M2=Xb%*%(c%x%a)
          b=svd(M2)$u[,1:b1]
          
          M3=Xc%*%(a%x%b)
          c=svd(M3)$u[,1:c1]
          
          A1=Xa-a%*%t(a)%*%Xa%*%(c%x%b)%*%(t(c)%x%t(b))
          if(control-norm(A1,"F")^2<10^-6){
            erro=TRUE
          }
          else{
            A=a
            B=b
            C=c
            control=norm(A1,"F")^2
          }
        }
        
        G=array(t(A)%*%Xa%*%(C%x%B),dim=c(a1,b1,c1))
        Final[linha,8]=k
        Final[linha,9]=round(control,2)
        Final[linha,10]=round(sum(G^2)/sum(Xa^2),2)
        
        #Usando fun��o Tucker
        W=array(c(CO_I[1:10,],CO_I[11:20,],CO_I[21:30,],CO_I[31:40,]),dim=c(10,13,4))
        U=tucker(W,p=a1,q=b1,r=c1)
        Final[linha,11]=U$cont
        A1=Xa-U$a%*%t(U$a)%*%Xa%*%(U$cc%x%U$b)%*%(t(U$cc)%x%t(U$b))
        Final[linha,12]=round(norm(A1,"F")^2,2)
        Final[linha,13]=round(sum(U$g^2)/sum(W^2),2)
      }
    }
  }
}
#P,Q,R - n�mero de componentes do 1�,2� e 3� modo respetivamente. S=P+Q+R
#k1,k2,k3 - n�mero de itera��es
#Frob1, Frob2, Frob3 - Norma de Frobenius em cada inicializa��o
#%1,%2,%3 - propor��o de vari�ncia explicada
colnames(Final)=c("P","Q","R","S","K1","Frob1","%1","K2","Frob2","%2","K3","Frob3","%3") 
View(Final) #Tabela com todos os modelos

colMeans(Final[,c(5,8,11)])#itera��es medias em cada modelo

