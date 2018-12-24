##################@European Call Price@###################

BSM_c <- matrix(nrow = 500, ncol = 1)

for (t in 1:500) {
  
  
  T = 1 #Time to Maturity (Time to maturity is in years)
  S_o =100 #Initial Stock Price
  K=100  #Strike Price 
  sig=0.25 #Volatility 
  r=0.05 #Interest Rate
  q= 0.01 #Dividend Rate (q=Dividend rate because we are dealing with stocks.)
  
  
  N = t    
  Nj = N+1    
  dt <- T/N    
  dx <- sig*sqrt(3*dt)    
 
   
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0,S_t - K)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
  }
  
  BSM_c[t,1] <- Value_t[1,Nj+1]
  
}
write.csv(BSM_c, "Black Scholes Merton Finite Difference European Call Solution.csv")


############## @European Put Price@###############

BSM_p <- matrix(nrow = 500, ncol = 1)
for (t in 1:500) {
  
  N = t
  Nj = N+1
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0,K-S_t)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
  }
  
  BSM_p[t,1] <- Value_t[1,Nj+1]
  
}

write.csv(BSM_p, "Black Scholes Merton Finite Difference European Put Solution.csv")


###############################################################
BSM_c <- matrix(nrow = 500, ncol = 1)

for (t in 1:500) {
  
  N = t
  Nj = N+1
  S_o = 100
  K = 100
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_mat <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_mat[N+1,] <- pmax(0,S_t - K)
  
  for(i in N:1){
    Value_mat[i,2:(2*Nj)]<-pu*Value_mat[i+1,3:(2*Nj+1)] + pm*Value_mat[i+1,2:(2*Nj)] + pd*Value_mat[i+1,1:(2*Nj-1)]
    Value_mat[i,1]<- Value_mat[i,2]
    Value_mat[i,(2*Nj+1)]<-Value_mat[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
  }
  
  BSM_c[t,1] <-Value_mat[1,Nj+1]
} 
write.csv(BSM_c, "Finite Difference Method Call (To verify).csv")



##################@American Call Price@#####################


american_call <- matrix(nrow = 500, ncol = 3)

for (t in 1:500) {
  
  S_o = 100
  K=100  
  N = t
  T = 1 
  
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, S_t - K)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],S_t - K)
  }
  
  american_call[t,1] <- Value_t[1,Nj+1]
  
}

write.csv(american_call, "American Call Price.csv")

##################@American Call Price with increase in time to maturity@##################


american_call <- matrix(nrow = 500, ncol = 3)

for (t in 1:500) {
  
  S_o = 100
  K=100  
  N = 100
  T = t 
  
  dt <- T/100
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, S_t - K)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],S_t - K)
  }
  
  american_call[t,1] <- Value_t[1,Nj+1]
  
}
write.csv(american_call, "American Call Price with inrease in time to maturity.csv")

#######@Amercan Call Price with increasing Spot@########

american_call <- matrix(nrow = 500, ncol = 3)

for (t in 1:500) {
  
  S_o = 100+t  
  K=100 
  N = t
  T = 1 
  
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, S_t - K)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],S_t - K)
  }
  
  american_call[t,1] <- Value_t[1,Nj+1]
  american_call[t,2] <- S_o
  
}
write.csv(american_call, "American Call Price Finite Difference Method with Increasing Spot.csv")

#####@American Call with Increasing Strike@#####

american_call <- matrix(nrow = 500, ncol = 3)

for (t in 1:500) {
  
  S_o = 100 
  K=100+t 
  N = t
  T = 1 
  
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, S_t - K)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],S_t - K)
  }
  
  american_call[t,1] <- Value_t[1,Nj+1]
  american_call[t,3] <- K
  
}
write.csv(american_call, "American Call Price Finite Difference Method Increasing Strike.csv")

############################# @American Put Price@#############################################

T=1
N=500
Nj = N+1
K=100
american_put <- matrix(nrow = 500, ncol = 3)
for (t in 1:500) {
  
  S_o =100 
  N=t
  Nj=N+1
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, K - S_t)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],K-S_t)
  }
  
  american_put[t,1] <- Value_t[1,Nj+1]
  
}

write.csv(american_put, "American Put Price Finite Difference Method.csv")

####### @American Put with increase in time to maturity@  ###########################################################################################
T=1
N=500
Nj = N+1
K=100
american_put <- matrix(nrow = 500, ncol = 3)
for (t in 1:500) {
  
  S_o =100
  N=100
  T=t
  Nj = N+1
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, K - S_t)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],K-S_t)
  }
  
  american_put[t,1] <- Value_t[1,Nj+1]
  
}

write.csv(american_put, "American Put Price with increase in time to maturity.csv")

#######(@American Put Price with Increasing Spot@)########
T=1
N=500
Nj = N+1
K=100
american_put <- matrix(nrow = 500, ncol = 3)
for (t in 1:500) {
  
  S_o =100+ t 
  N=500
  Nj=N+1
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, K - S_t)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],K-S_t)
  }
  
  american_put[t,1] <- Value_t[1,Nj+1]
  american_put[t,2] <- S_o
}

write.csv(american_put, "American Put Price Finite Difference Method Increasing Spot.csv")


##########(@American Put with increasing Strike@)###########
T=1
N=500
Nj = N+1
K=100
american_put <- matrix(nrow = 500, ncol = 3)
for (t in 1:500) {
  
  S_o =100
  N=500
  Nj=N+1
  K=100+t
  dt <- T/N
  dx <- sig*sqrt(3*dt)
  nu <- r - q - 0.5*sig^2
  edx <- exp(dx)
  pu <- 0.5*dt*(sig^2/dx^2+nu/dx)
  pd <- 0.5*dt*(sig^2/dx^2-nu/dx)
  pm <- 1 - pu - pd - r*dt
  
  S_t<-S_o*exp((-Nj:Nj)*dx)
  Value_t <- matrix(ncol=2*Nj+1,nrow=N+1)
  Value_t[N+1,] <- pmax(0, K - S_t)
  
  for(i in N:1){
    Value_t[i,2:(2*Nj)]<-pu*Value_t[i+1,3:(2*Nj+1)] + pm*Value_t[i+1,2:(2*Nj)] + pd*Value_t[i+1,1:(2*Nj-1)]
    Value_t[i,1]<- Value_t[i,2]
    Value_t[i,(2*Nj+1)]<-Value_t[i,(2*Nj)] + (S_t[(2*Nj+1)]-S_t[(2*Nj)])
    Value_t[i,]<-pmax(Value_t[i,],K-S_t)
  }
  
  american_put[t,1] <- Value_t[1,Nj+1]
  american_put[t,3] <- K
}

write.csv(american_put, "American Put Price Finite Difference Method Increasing Strike.csv")
