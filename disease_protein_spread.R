#install.packages("network")
#library(network)
#help("network-package")

#===================================================================================
#	-NEIGH: is a network of neighbors
#	NET is the same saved as a network
#===================================================================================

L<-num_obs<-200  #number of cells (nodes)
NEIGH<-matrix(rbinom(L^2,1,.01),L)
for( i in 1:L)
  {
  NEIGH[i,]<-NEIGH[,i]
  }

diag(NEIGH)<-0

NET<-network(NEIGH, directed=FALSE)

#===================================================================================
#	Plotting functions
#	-COORDS lets it draw network once and saves the coordinates of each node
#	-COLS a vector to define colors in the plot 
#	-visualize: plots network, coloring the "high risk" nodes red
#===================================================================================

COORDS<-plot(NET)
COLS<-rep(1,L)

visualize_net<-function()
  {
  COLS[P>20]<-"red"
  COLS[P<=20]<-"green"
  COLS[A<1]<-"black"
  plot(NET,vertex.col=COLS, jitter=FALSE, coord=COORDS)
  }

visualize_nodes<-function()
	{
  COLS[P>20]<-"red"
  COLS[P<=20]<-"green"
  COLS[A<1]<-"black"
	points(COORDS,col=COLS,pch=19)
	}

#===================================================================================
#	-Set up Simulation 
#===================================================================================

W<<-sample(150:250,L,replace=TRUE)
Wi<<-W

setup<-function()
	{
  W <<- Wi
	P<<-rep(0,L)
  A<<-rep(1,L)
	}

#===================================================================================
#	Key functions of stochastic simulation
#===================================================================================

seed<-function(node)
  {
  P[node]<<-15
  }

death<-function(threshold)
  {
  for(n in 1:L)
    {
    if(P[n]>threshold)
      {
      A[n]<<-0
      P[n]<<-0
      W[n]<<-0
      }
    } 
  }

birthW<-function(b)
  {
  W<<-W+b
  }

deathW<-function(d)
  {
  W<<-W-d*W
  }

natagg<-function(gamma)
  {
  alive<<-which(A==1)
  for(z in alive)
  {
    P[z]<<-P[z]+gamma*W[z]
    W[z]<<-W[z]-gamma*W[z]
  }
  }

natdis<-function(alpha)
  {
  alive<<-which(A==1)
  for(m in alive)
    {
    if(P[m]-alpha[m]>=0)
      {
      P[m]<<-P[m]-alpha[m]
      W[m]<<-W[m]+alpha[m]
      }
    else
      {
      W[m]<<-P[m]+W[m]
      P[m]<<-0  
      }
    }
  }

agg<-function(beta)
  {
  alive<<-which(A==1)
  for (j in alive)
    {
    oldP<<-P[j]
    P[j]<<-P[j]+beta*P[j]*W[j]
    W[j]<<-W[j]-beta*oldP*W[j]
  }
  }

infect<-function(rho)
  {
  alive<<-which(A==1)
  for (j in alive)
    {
    index<<-c()
    for (k in alive)
      {
      if (NEIGH[k,j]==1)
        {
        index<<-c(index,k)
        }
      }
    P[j]<<-P[j]+rho*sum(P[index])
    P[index]<<-P[index]-rho*P[index]

    }
  }

first_death<-function(i)
{
  if (sum(A) > 199)
  {
    #print("No deaths")
  }
  if (sum(A) <= 199)
  {
    deaths <<- c(deaths,i)
    #cat("First death: ", min(deaths))
  }
  #return (min(deaths))
}

death60<-function(Death_percent)
{
  val<-Death_percent[reps]*0.6
  critical<-which(Death_percent>=val)
  return(min(critical))
}

#===================================================================================
#	Simulation Setup
#	-reps: define how many time steps
#	-setup()  set S to 1, I, E, R to 0
#	-PREV: an empty vector to keep track of prevalence at each timestep
#	-par sets up a plot window. par(ask=TRUE) requires a <Enter> between plots (slows sown our movie)
#
#	-b : probaility infection given exposure / time step (ie rate)
#
#	-Assign an index case
#	-infect it with infect()
#	-expose its neighbors with expose()
#===================================================================================
index_case<<-sample(L,1)
simulationsetup <- function(alpha,rho,threshold) #Wi*gamma,0.2,100
{
reps<<-150
setup()
PREV<<-rep(NA,reps)
par(ask = TRUE)		

# rate of change to P dependent on # existing P
beta<<-.0005 
# underlying rate, dependent on W
gamma<<-.0001
# just constant, not dependent on any population
alpha<<-alpha
# rate of migration, dependent on P of neighboring cells
rho<<-rho
# death rate of W, dependent on W
d<<-.0002
# birth rate of W, independent of W
b<<-Wi*d
#"death" threshold
threshold<<-threshold

deaths <<- c()

#Death Percent
Death_percent<<-rep(0,reps)

seed(index_case)
#visualize_net()
}

#===================================================================================
#	Simulation loop
#	-reps: define how many time steps
#	-outputs a movie
#	-after movie plots prevalence over time
#===================================================================================

main <- function(alpha,rho,threshold)
{
simulationsetup(alpha,rho,threshold)
for(i in 1:reps)
	{
    Sys.sleep(0.05)
    deathW(d)
    birthW(b)
    natagg(gamma)
    natdis(alpha)
    agg(beta)
    infect(rho)
    death(threshold)
    #visualize_net()
    first_death(i)
    Death_percent[i]<<-sum(A==0)/L
#     if(i%%10==0)
#       {
#       old.par <- par(mfrow=c(1, 2))
#       plot(Death_percent,typ ="o", col="red", ylim=c(0,1),ylab="Percent Dead",xlab="Time step")
#       visualize_net()
#       #par(old.par)
#       Sys.sleep(0)
#       }
  }
#   print(Death_percent)
  first_death_result <<- min(deaths)
  death60_result <<- death60(Death_percent)
  ret <<-c(mean(alpha),rho,threshold,first_death_result,death60_result,Death_percent)
#   print(ret)
  return(ret)
}

#RESULTS
plotdis <- function()
{
length <- reps+5
normal <<- main(Wi*gamma,0.2,100)
normalplot <<- c(normal[6:length])
normaldf <<- c(normal[1:5],normal[length])
stemcell <- main(Wi*gamma,0.2,175)
stemcellplot <- stemcell[6:length]
stemcelldf <- c(stemcell[1:5],stemcell[length])
antibody <- main(Wi*gamma,0.05,100)
antibodyplot <- antibody[6:length]
antibodydf <- c(antibody[1:5],antibody[length])
genetherapy <- main(rep(0.035,L),0.2,100)
genetherapyplot <- genetherapy[6:length]
genetherapydf <- c(genetherapy[1:5],genetherapy[length])
labelsdf <- c("alpha","rho","death threshold","first death","60% dead","final death %")
results <<- data.frame(labels=labelsdf,normal=normaldf,stemcell=stemcelldf,antibody=antibodydf,genetherapy=genetherapydf)
visualize_net()
plot(normalplot,typ ="o", col="red", ylim=c(0,1),ylab="Percent Dead",xlab="Time step")
lines(stemcellplot,typ ="o", col="blue")
lines(antibodyplot,typ ="o", col="green")
lines(genetherapyplot,typ ="o", col="black")
}
#visualize_net()
plotdis()
print(results)
