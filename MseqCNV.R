Pos_Size <- read.csv("D:\\SampleData.csv")
chromosome_size =5000000
sample_size=40
segment_size=150
lambda_initial_estimate <-  13
T=1000
Insertion_Size=200
Read_Length=50
iteration=3
##############################################################################
segment_position=c()
current=1
num_segments=0
while (current<= (chromosome_size-segment_size)){
           segment_position=c(segment_position, c(current,current+segment_size-1))
           current=current+segment_size
           num_segments=num_segments+1
}
x <- list()
segment_reads_passing <- array(x,c(sample_size,num_segments))
segment_reads_flag    <- array(x,c(sample_size,num_segments))
segment_reads_count   <- matrix(rep(1,sample_size*num_segments),nrow=sample_size,ncol=num_segments)
########################################################################################################
d=dim(Pos_Size)
for  (i in 1:d[1]){
         Pos            =Pos_Size[i,1]
         Insert_Size    =Pos_Size[i,2]
         paired_flag    =Pos_Size[i,3]
         mapping_quality=Pos_Size[i,4]
         pair_num       =Pos_Size[i,5]
         sample_number  =Pos_Size[i,6]
         if  (paired_flag==2){
             paired_flag=18}   
         if ((paired_flag==18)&(Insert_Size<=12000)&(Insert_Size>0)&(mapping_quality>=25)){
            S=Pos+Read_Length
            E=Pos+Insert_Size+Read_Length 
            start_segment_row=ceiling(S/segment_size)
            end_segment_row  =ceiling(E/segment_size)
            for (j in start_segment_row:end_segment_row){
                L=length(segment_reads_passing[[sample_number,j]])
                segment_reads_passing[[sample_number,j]][L+1]=Insert_Size
                segment_reads_flag[[sample_number,j]][L+1]=paired_flag
            }
            segment_reads_count[sample_number,start_segment_row]=segment_reads_count[sample_number,start_segment_row]+1
            segment_reads_count[sample_number,end_segment_row] =segment_reads_count[sample_number,end_segment_row]  +1
         } 
}
###############################################################################
Mu <- function(Insertion_Size,k,s,t,observation, mu_sigma_beta_alpha_membership){ 
membership_pro = mu_sigma_beta_alpha_membership[[k,s,6]][t,]
estimated_mu1=0
estimated_mu2=0
sum1=0.0001
sum2=0.0001
sum3=0.0001
sum4=0.0001
if (length(observation)>0){
   for (i in 1:length(observation)){
        sum1 = sum1 +    membership_pro[i]*observation[i]
        sum2 = sum2 +    membership_pro[i]
        sum3 = sum3 + (1-membership_pro[i])*observation[i]
        sum4 = sum4 + (1-membership_pro[i])
   }
}
estimated_mu1 = sum1/sum2
estimated_mu2 = sum3/sum4
if (estimated_mu1 < 40){
   estimated_mu1  = 100000
}
if (estimated_mu2 < (Insertion_Size + 90)){
   estimated_mu2  = 100000
}
mu=c(Insertion_Size,estimated_mu2)
res <-list() 
res$mu <- mu
return(res)
}
###############################################################################
EmpericalDist <- function(mi,ma,length_of_interval,insertion_sizes){
   range=(ma-mi);
   INTERVALS=c()
   upper_bond=mi
   while (upper_bond <= (ma+4)){
         INTERVALS=rbind(INTERVALS,c(upper_bond,(upper_bond+length_of_interval)))
         upper_bond=upper_bond+length_of_interval
   }
   d=dim(INTERVALS)
   Pro_Intervals=rep(0,d[1])
   for (i in 1:length(insertion_sizes)){
          point=insertion_sizes[i]
          Found=0
          k=1
          d=dim(INTERVALS)
          while ((Found==0)&(k<=d[1])){
               if ((point>=INTERVALS[k,1])&(point<=INTERVALS[k,2])){
                  Found=1
                  Pro_Intervals[k]=Pro_Intervals[k]+1
               }
               k=k+1
          }
   }
Pro_Intervals=Pro_Intervals/sum(Pro_Intervals)
pro_intervals_limit=cbind(INTERVALS,Pro_Intervals)

res <-list() 
res$pro_intervals_limit <- pro_intervals_limit
return(res)
}
################################################################################
betaMuSigma <- function(Insertion_Size,C,X){
   X=matrix(X)
   Cluster_A=c()
   Cluster_B=c()
   for (i in 1:length(X)){
       if (X[i] <= (Insertion_Size+80)){
               Cluster_A=c(Cluster_A, X[i])
       }
       else if (X[i]> (Insertion_Size+80)){
               Cluster_B=c(Cluster_B, X[i])
       }
   }
   mu_A=-1
   mu_B=-1
   if (length(Cluster_A)>0){
       mu_A=mean(Cluster_A)
   }
   if (length(Cluster_B)>0){
       mu_B=mean(Cluster_B)
   }
   ##############################################################################
   beta=matrix(c(0.1,0.9),1,2)
   sigma=matrix(c(20,20),1,2)

   if (C==0){
      if (mu_B>0){
         mu   =c(Insertion_Size,mu_B)}
      else if (mu_A>0){
         mu   =c(Insertion_Size,100000)}
   }  else if (C==1){
         mu   =c(Insertion_Size,250)
   } else if (C==2){
         mu   =c(Insertion_Size,500)
   } else if (C==3){
         mu   =c(Insertion_Size,750)
   } else if (C==4){
         mu   =c(Insertion_Size,1000)
   } else if (C==5){
         mu   =c(Insertion_Size,1500)
   } else if (C==6){
         mu   =c(Insertion_Size,2200)
   } else if (C==7){
         mu   =c(Insertion_Size,2700)
   } else if (C==8){
         mu   =c(Insertion_Size,3200)
   } else if (C==9){
         mu   =c(Insertion_Size,3700)
   } else if (C==10){
         mu   =c(Insertion_Size,4200)
   } else if (C==11){
         mu   =c(Insertion_Size,4700)
   } else if (C==12){
         mu   =c(Insertion_Size,5200)
   } else if (C==13){
         mu   =c(Insertion_Size,5700)
   } else if (C==14){
         mu   =c(Insertion_Size,6200)
   } else if (C==15){
         mu   =c(Insertion_Size,6700)
   } else if (C==16){
         mu   =c(Insertion_Size,7200)
   } else if (C==17){
         mu   =c(Insertion_Size,7700)
   } else if (C==18){
         mu   =c(Insertion_Size,8200)
   } else if (C==19){
         mu   =c(Insertion_Size,8700)
   } else if (C==20){
         mu   =c(Insertion_Size,9200)
   } else if (C==21){
         mu   =c(Insertion_Size,9700)
   } else if (C==22){
         mu   =c(Insertion_Size,10200)}
   ##############################################################################
   mu=matrix(mu,1,2)
   it=2
   for (g in 1:it){
        p_yi_x=c(0,0)
        beta_new=c(0,0)
        mu_new=c(0,0)
        sigma_new=c(0,0)
        P_Y=c()
        for (k in 1:length(X)){
            p_x_k=beta[g,1]*dnorm(X[k],mu[g,1],sigma[g,1])+beta[g,2]*dnorm(X[k],mu[g,2],sigma[g,2])
            for (i in 1:2){
               p_yi_x[i]=beta[g,i]*dnorm(X[k],mu[g,i],sigma[g,i])/(p_x_k+0.0000000001)
            }
            p=rbind(p_yi_x[1],p_yi_x[2])
            P_Y=cbind(P_Y,p)
        }
        for (L in 1:2){
            beta_new[L] =sum(P_Y[L,])/length(X)
            mu_new[L]   =(P_Y[L,]%*%X)/(sum(P_Y[L,]+0.000001))
            sigma_new[L]=sqrt((P_Y[L,]%*%((X-mu_new[L])^2))/(sum(P_Y[L,])+0.000001))  

            if (beta_new[L]==0){
                beta_new[L]= 0.0001
            }
            if (sigma_new[L]<=1){
                sigma_new[L]=1
            }
            if (!(sigma_new[L]>0)){
                  sigma_new[L]=sigma[g,L]
            }
            if (!(mu_new[L]>0)){
                  mu_new[L]=mu[g,L]
            }
        }
        beta =rbind(beta, c(beta_new[1], beta_new[2]))
        mu   =rbind(mu  , c(Insertion_Size, mu_new[2]))
        sigma=rbind(sigma, c(sigma_new[1], sigma_new[2]))
    }
res <-list() 
res$beta <- beta
res$mu <- mu
res$sigma <- sigma
return(res)
}
##########################################################################
betaModification <- function(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g){
F=1
mu1=mu[g,1]
mu2=mu[g,2]
        if ((beta[g,1]!=0.0001)&&(beta[g,2]==0.0001)){
        mu_segment   =rbind(mu_segment,    c(mu1,100000))
        sigma_segment=rbind(sigma_segment, c(sigma[g,1],1))
        beta_segment =rbind(beta_segment,  c(1,0))
        }
         else if ((beta[g,1]==0.0001)&&(beta[g,2]!=0.0001)){
        mu_segment   =rbind(mu_segment,    c(100000,mu2))
        sigma_segment=rbind(sigma_segment, c(1,sigma[g,2]))
        beta_segment =rbind(beta_segment , c(0,1))
        } 
        else if (((beta[g,1]!=0.0001)&&(beta[g,2]!=0.0001))&&(abs(sum(beta[g,])-1)<0.00001)){       
        mu_segment   =rbind(mu_segment, c(mu1,mu2))
        sigma_segment=rbind(sigma_segment, c(sigma[g,1],sigma[g,2]))
        beta_segment =rbind(beta_segment,  c(beta[g,1], beta[g,2]))
        } 
        else { F=0 }
res <-list() 
res$F <- F
res$mu_segment <- mu_segment
res$sigma_segment <- sigma_segment
res$beta_segment <- beta_segment
return(res)
}
###########################################################################
insertion_sizes=matrix(rnorm(8000,mean = Insertion_Size, sd = 20),ncol=1)
mi=120
ma=280
length_of_interval=2
ProDensity=EmpericalDist(mi,ma,length_of_interval,insertion_sizes)$pro_intervals_limit
x <- list()
mu_sigma_beta_alpha_membership <- array(x,c(iteration,sample_size,7))
pro <- array(x,c(iteration,sample_size))

for (k in 1:iteration){
    for (s in 1:sample_size){
         mu_sigma_beta_alpha_membership[[k,s,1]] <- matrix(rep(0,T*2),T,2)
         mu_sigma_beta_alpha_membership[[k,s,2]] <- matrix(rep(0,T*2),T,2)
         mu_sigma_beta_alpha_membership[[k,s,3]] <- matrix(rep(0,T*2),T,2)
         mu_sigma_beta_alpha_membership[[k,s,4]] <- matrix(rep(0,T*6),T,6)
         mu_sigma_beta_alpha_membership[[k,s,5]] <- matrix(rep(0,T*6),T,6)
         MP_numbers <- 0
         for (t in 1:T){
              d=length(segment_reads_passing[[s,t]])
              if (d > MP_numbers){
                  MP_numbers <- d
              }
         }
         mu_sigma_beta_alpha_membership[[k,s,6]] <- matrix(rep(0,T*MP_numbers),T,MP_numbers)
         mu_sigma_beta_alpha_membership[[k,s,7]] <- matrix(rep(0,T*2),T,2)

         pro[[k,s]] <- matrix(rep(0,T*7),T,7)
    }
}

ini <- 0.02
for (s in 1:sample_size){
    for (t in 1:T){
         mu_sigma_beta_alpha_membership[[1,s,4]][t,] <- c(ini,ini,1-5*ini,ini,ini,ini)      
         mu_sigma_beta_alpha_membership[[1,s,5]][t,] <- c(0.001,0.75,0.999,0.999,0.999,0.999)
    }
}

lambda <- c()
lambda <- c(lambda, lambda_initial_estimate*rep(1,T))
gamma <-  c(2,2,10,2,2,2)
nu <- matrix(c(2,50,60,20,20,2,20,2,20,2,20,2), ncol=2, nrow=6,byrow=T)
teta <- c(0.05,0.5,1,1.5,2,2.5)

for (s in 1:sample_size){
     g <- 2
     mu_segment <- c()
     sigma_segment <- c()
     beta_segment <- c()
     for (t in 1:T){
                L_pre=length(beta_segment[ ,1])
                if (length(segment_reads_passing[[s,t]])> 0){
                        F=0;
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,0,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,1,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,2,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,3,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,4,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,5,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,6,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,7,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,8,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,9,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,10,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,11,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,12,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,13,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,14,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,15,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,16,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,17,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,18,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,19,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,20,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment
                            F=OUT$F  
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,21,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                        if (F==0){
                            OUT=betaMuSigma(Insertion_Size,22,segment_reads_passing[[s,t]])
                            sigma=OUT$sigma
                            beta =OUT$beta
                            mu   =OUT$mu
                            OUT=betaModification(beta,mu,sigma,mu_segment,sigma_segment,beta_segment,g)
                            mu_segment=OUT$mu_segment
                            sigma_segment=OUT$sigma_segment
                            beta_segment=OUT$beta_segment 
                            F=OUT$F 
                        }
                }
                else {
                            beta_segment =rbind(beta_segment, c(1,0));
                            mu_segment   =rbind(mu_segment  , c(Insertion_Size,100000)); 
                            sigma_segment=rbind(sigma_segment,c(20,1));
                } 
                L_post=length(beta_segment[ ,1])
                if (L_pre==L_post){
                     t
                }
     }
     beta_segment=round(beta_segment,2)
     mu_sigma_beta_alpha_membership[[1,s,1]]=mu_segment
     mu_sigma_beta_alpha_membership[[1,s,2]]=sigma_segment
     mu_sigma_beta_alpha_membership[[1,s,3]]=beta_segment
}
##################### initialize alpha, mu, and sigma based on the previous step ####################
for (s in 1:sample_size){
       for (t in 1:T){  
               mu_1          =mu_sigma_beta_alpha_membership[[1,s,1]][t,1]
               mu_2          =mu_sigma_beta_alpha_membership[[1,s,1]][t,2]
               min_prime     = mi+mu_2-mu_1
               if (length(segment_reads_passing[[s,t]])>0){
                     for (n in 1:length(segment_reads_passing[[s,t]])){
                           sample_insertion=segment_reads_passing[[s,t]][n]
                           p_z_jr=0
                           for (i in 1:6){          
                                 alpha_i =mu_sigma_beta_alpha_membership[[1,s,4]][t,i]
                                 beta_i  =mu_sigma_beta_alpha_membership[[1,s,5]][t,i]
                                 row=round(((sample_insertion-mi)+1)/length_of_interval)
                                 d=dim(ProDensity)
                                 if ((row>=1)&&(row<=d[1])){
                                    p_1=ProDensity[row,3]
                                 } else {
                                    p_1=0
                                 }  
                                 row=round(((sample_insertion-min_prime)+1)/length_of_interval)
                                 if ((row>=1)&&(row<=d[1])){
                                    p_2=ProDensity[row,3]
                                 } 
                                 else{
                                    p_2=0
                                 }
                                 if ((beta_i*p_1+(1-beta_i)*p_2)!=0){
                                    p_z_jr=p_z_jr + alpha_i*beta_i*p_1/(beta_i*p_1+(1-beta_i)*p_2)
                                
                                 }
                           }
                           mu_sigma_beta_alpha_membership[[1,s,6]][t,n]=p_z_jr
                    }
             }
      } 
}
#######################################################################       
for (k in 2:iteration){
           for (t in 1:T){
              for (s in 1:sample_size){
                    n    =length(segment_reads_passing[[s,t]]) 
                    n_s_1=floor(sum(mu_sigma_beta_alpha_membership[[k-1,s,6]][t,]))
                    n_s_2=ceiling(n-sum(mu_sigma_beta_alpha_membership[[k-1,s,6]][t,])) 
                    mu_sigma_beta_alpha_membership[[k,s,7]][t,1]=n_s_1
                    mu_sigma_beta_alpha_membership[[k,s,7]][t,2]=n_s_2 
                    p_state=c()
                    p_n_j=0;
                    for (i in 1:6){
                        alpha_i= mu_sigma_beta_alpha_membership[[k-1,s,4]][t,i] 
                        beta_i = mu_sigma_beta_alpha_membership[[k-1,s,5]][t,i]
                        if (k==2){
                            poisson_par=teta[i]*(lambda[t])
                        }
                        else{
                            poisson_par=teta[i]*(lambda[k-1,t])
                        }
                        read_count_pro= dpois(length(segment_reads_passing[[s,t]]),poisson_par)
                        mate_pair_pro = dbinom(n_s_1,n,beta_i)
                        p_state_i     = alpha_i*read_count_pro*mate_pair_pro
                        p_state=rbind(p_state , p_state_i)
                        p_n_j         = p_n_j + p_state_i;
                    }
                    p_state = rbind(p_state , p_n_j);
                    pro[[k,s]][t,]=p_state;   
              }  
           }
           #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% alpha update
           for (t in 1:T){
                  sum1=rep(0,7)
                  sum2=0
                  for (i in 1:6){
                      for (s in 1:sample_size){
                          sum1[i] = sum1[i] + pro[[k,s]][t,i]
                      }    
                  }
                  for (s in 1:sample_size){
                      sum2 = sum2 + pro[[k,s]][t,7]
                  }
                  new_alpha=c()
                  for (i in 1:6){
                      new_alpha_i = (sum1[i]+ (gamma[i]-1)*sum2)/(sum2*(1+sum(gamma)-6))
                      new_alpha=rbind(new_alpha , new_alpha_i)
                  }
                  for (s in 1:sample_size){
                      mu_sigma_beta_alpha_membership[[k,s,4]][t,]=new_alpha
                  }
           }
           ########################################################################################## beta update 
           for (t in 1:T){
               sum3=rep(0,6)
               sum4=rep(0,6)
               sum5=0
               new_beta=c()
               #######################################################################################
               for (s in 1:sample_size){
                    mu_1          =mu_sigma_beta_alpha_membership[[k-1,s,1]][t,1]
                    mu_2          =mu_sigma_beta_alpha_membership[[k-1,s,1]][t,2]
                    min_prime=mi+mu_2-mu_1
                    sum_membership_pro=rep(0,7)
                    ############################################################################
                    if (length(segment_reads_passing[[s,t]])>0){      
                          for (n in 1:length(segment_reads_passing[[s,t]])){
                              sample_insertion=segment_reads_passing[[s,t]][n]
                              p_z_jr=0
                              for (i in 1:6){            
                                    alpha_i=mu_sigma_beta_alpha_membership[[k-1,s,4]][t,i]
                                    beta_i =mu_sigma_beta_alpha_membership[[k-1,s,5]][t,i]
                                    row=round(((sample_insertion-mi)+1)/length_of_interval)
                                    d=dim(ProDensity)
                                    if ((row>=1)&&(row<=d[1])){
                                         p_1=ProDensity[row,3]
                                    }
                                    else{  
                                         p_1=0
                                    } 
                                    row=round(((sample_insertion-min_prime)+1)/length_of_interval) 
                                    if ((row>=1)&&(row<=d[1])){
                                         p_2=ProDensity[row,3]
                                    }
                                    else{  
                                         p_2=0
                                    }  
                                    if ((beta_i*p_1+(1-beta_i)*p_2)!=0){
                                         p_z_jr=p_z_jr + alpha_i*beta_i*p_1/(beta_i*p_1+(1-beta_i)*p_2) 
                                    }
                                    ##################################################################
                                    if ((beta_i*p_1+(1-beta_i)*p_2)!=0){
                                       sum_membership_pro[i] = sum_membership_pro[i] + beta_i*p_1/(beta_i*p_1+(1-beta_i)*p_2)
                                    }
                                    ##################################################################
                              }  
                              mu_sigma_beta_alpha_membership[[k,s,6]][t,n]=p_z_jr
                          }
                      }
                      #############################################################################
                      for (i in 1:6){
                              sum3[i] = sum3[i] + sum_membership_pro[i]*pro[[k,s]][t,i] 
                      }
                      for (i in 1:6){
                              sum4[i] = sum4[i] + length(segment_reads_passing[[s,t]])*pro[[k,s]][t,i] 
                      }
                      sum5 = sum5 + pro[[k,s]][t,7]
                          
               }
               ########################################################################################
               for (i in 1:6){
                   new_beta_i= (sum3[i] + (nu[i,1]-1)*sum5)/(sum4[i] + (nu[i,1]+nu[i,2]-2)*sum5)
                   new_beta=rbind(new_beta , min(1,new_beta_i))
               }
               for (s in 1:sample_size){
                      mu_sigma_beta_alpha_membership[[k,s,5]][t,]=new_beta
               } 
           }   
           ############################################################################################lambda update
           sum6=0
           sum7=0
           for (t in 1:T){
               for (s in 1:sample_size){
                   for (i in 1:6){
                          sum6 = sum6 + teta[i]*pro[[k,s]][t,i]     
                   }
                   sum7 = sum7 + length(segment_reads_passing[[s,t]])*pro[[k,s]][t,7]
               }
           }
           lambda_new = sum7/sum6
           lambda=rbind(lambda , lambda_new*rep(1,T))
           ########################################################################################### mu update
           for (t in 1:T){
               for (s in 1:sample_size){ 
                       observation=segment_reads_passing[[s,t]]
                       mu=Mu(Insertion_Size,k,s,t,observation, mu_sigma_beta_alpha_membership)$mu 
                       mu_sigma_beta_alpha_membership[[k,s,1]][t,]=mu
               }
           }
           ##########################################################################################
}
#####################################################################################################
iteration <- 2
CNVs <- matrix()
final_path=c()
for (s in 1:sample_size){
           path=c()
           for (t in 1:T){
                 I=which.max(pro[[iteration,s]][t,1:6])
                 C=max(pro[[iteration,s]][t,1:6])
                 path=cbind(path,I)
           }   
           final_path=rbind(final_path,path)
}
for (s in 1:sample_size){
        PREDICTIONS=final_path[s,]
        n=length(PREDICTIONS)
        i=1
        start_pos=1
        start_cnv=PREDICTIONS[i]
        CNV_num=0
        while (i<n){
              if    ((i<n)&&(start_cnv==PREDICTIONS[i+1])){
                     i=i+1
              }
              else{ 
                     end_pos=i*segment_size
                     if (start_cnv!=3){
                              CNV_num=CNV_num+1
                              if ((s==1)&&(CNV_num==1)){
                                  CNVs[1]=s
                                  CNVs[2]=start_pos
                                  CNVs[3]=end_pos
                                  CNVs[4]=start_cnv-1
                               }
                               else{
                                  CNVs=rbind(CNVs,c(s,start_pos,end_pos,start_cnv-1))
                               }
                     } 
                     i=i+1
                     start_pos=((i-1)*segment_size)+1
                     start_cnv=PREDICTIONS[i]
              } 
        } 
}
#####################################################################################################