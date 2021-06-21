Final=matrix(NA,nrow=278,ncol=13)
linha=0

for (a1 in 1:10){
  for (b1 in 1:13){
    for (c1 in 1:4){
      if (a1<=b1*c1 & b1<=a1*c1 & c1<=b1*a1){
        linha=linha+1
        
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
        
        decomp_Xb=svd(Xb%*%t(Xb))
        b=as.matrix(decomp_Xb$u[,1:b1])
        
        decomp_Xc=svd(Xc%*%t(Xc))
        c=as.matrix(decomp_Xc$u[,1:c1])
        
        #2
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
        
        U=tucker(W,p=a1,q=b1,r=c1)
        Final[linha,11]=U$cont
        A1=Xa-U$a%*%t(U$a)%*%Xa%*%(U$cc%x%U$b)%*%(t(U$cc)%x%t(U$b))
        Final[linha,12]=round(norm(A1,"F")^2,2)
        Final[linha,13]=round(sum(U$g^2)/sum(W^2),2)
      }
    }
  }
}
colnames(Final)=c("P","Q","R","S","K1","Frob1","%1","K2","Frob2","%2","K3","Frob3","%3")

length(which(Final[,7]==Final[,10] & Final[,10]==Final[,13]))#percentagem

######
tentar=function(x, p, q, r, test = 10^-6) 
{
  loss.old <- criter(x, 0)
  loss.new <- rep(0, 4)
  param <- init3(x, p, q, r)
  param <- step.g3(param)
  param <- newcomp3(param)
  loss.new[1] <- loss1.3(param, 0)
  paramp <- param
  cont <- 0
  while (abs(loss.old - loss.new[1]) > test) {
    cont <- cont + 1
    a.old <- paramp$a
    b.old <- paramp$b
    c.old <- paramp$cc
    loss.old <- loss.new[1]
    paramp <- stepi3(paramp)
    paramp <- step.g3(paramp)
    loss.new[3] <- loss2(paramp, b.old)
    paramp <- stepi3(paramp)
    paramp <- step.g3(paramp)
    loss.new[4] <- loss2(paramp, c.old)
    paramp <- stepi3(paramp)
    paramp <- step.g3(paramp)
    loss.new[2] <- loss2(paramp, a.old)
    loss.new[1] <- loss1.3(paramp, a.old)
  }
  paramp$cont <- cont
  paramp
}

criter(W, 0)
