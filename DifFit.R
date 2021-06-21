library(CA3variants)

#Depois do tratamento de dados e do cálculo das matrizes de covariância cruzada:
Xa=as.matrix(cbind(CO_I[1:4,],CO_I[5:8,],CO_I[9:12,]))
Xb=t(as.matrix(CO_I))
l1=as.numeric(t(CO_I[1:4,]));l2=as.numeric(t(CO_I[5:8,]));l3=as.numeric(t(CO_I[9:12,]))
Xc=rbind(l1,l2,l3)

M=matrix(0,nrow=19,ncol=5)
l=0

for (n in 1:19){
  l=l+1
  for (j in 1:24){
    for (i in 1:4){
      for (k in 1:3){
        if (i+j+k==n & i*j>=k & j*k>=i & k*i>=j){
          b=matrix(0,nrow=24,ncol=j)
          c=matrix(0,nrow=3,ncol=k)
          
          erro=FALSE           #faz ALS para cada modelo
          control=1000000
          
          while(erro==FALSE){
            M1=Xa%*%(c%x%b)
            a=svd(M1)$u[,1:i]
            
            M2=Xb%*%(c%x%a)
            b=svd(M2)$u[,1:j]
            
            M3=Xc%*%(a%x%b)
            c=svd(M3)$u[,1:k]
            
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
          
          G=array(t(A)%*%Xa%*%(C%x%B),dim=c(i,j,k))
          
          
          perc=(sum(G^2)/sum(Xa^2))
          if (perc>=M[l,5]){      #Escolhe o melhor modelo entre os modelos com S=s componentes
            M[l,1]=i
            M[l,2]=j
            M[l,3]=k
            M[l,4]=n
            M[l,5]=round(perc,4)
          }
        }
      }
    }
  }
}
M=M[c(3,5:19),] #retirar modelo com 4 componentes porque é impossível construir
M_diff=M

vec=numeric(16)
vec[1]=M_diff[1,5]

for (i in 2:16){ #calcular diferença com o modelo anterior
  vec[i]=round(M_diff[i,5]-M_diff[i-1,5],4)
}
M_diff=cbind(M_diff,vec)

vec=numeric(16)
for (i in 1:15){      #calcular rácio
  vec[i]=round(M_diff[i,6]/max(M_diff[(i+1):16,6]),3)
}
M_diff=cbind(M_diff,vec)


View(M_diff)
