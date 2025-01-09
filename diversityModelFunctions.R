
###functions used in the modelling

#exclude zeros to find mean alpha diversity
alpha1per<-function(Q1)
{
    res<-matrix(NA,1,ncol(Q1))

    for(i in 1:ncol(Q1))
    {
        xy<-Q1[,i]
        res[,i]<-length(xy[xy>0.01])

    }

    return(res)
}
gamma1per<-function(Q1)
{
    res<-matrix(NA,nrow(Q1),1)

    for(i in 1:nrow(Q1))
    {
        xy<-Q1[i,]
        res[i,]<-length(xy[xy>0.01])

    }

    return(sum(res>0))
}


changing.values<-function(hval, wval)
{
	h<-hval
	maxtime <- 1000000

	#td <- params$td
	w_adjust <- wval
	params <- Parameters(aid)
	topology <- params$topology
	pinit <- 1 / n_sp

	#res<-matrix(NA,1,3)
	#colnames(res)<-c("alpha","beta","gamma")

		A <- ReadMat(filename, w_adjust)
		E_ls <- ReadEnvs(envFile)
		C <- ReproductionRate(n_sp, n_com, as.numeric(h), E_ls)
		Q <- Initial(n_sp, n_com, pinit)
		Preg <- rowSums(Q) / sum(Q)
		tolerance <- 0.0000001 * n_sp * n_com
		Qdiff <- numeric()
		biomass <- numeric()
		g <- 0

		while (g < maxtime)
		{
		  prevQ <- Q # update the prevQ using the previous Q
		  Q <- Update(Q, A, C, nu, Preg, n_com) # new Q
		  diff <- calcDiff(Q, prevQ)
		  Qdiff <- c(Qdiff, diff)
		  biomass <- c(biomass, sum(Q))
		  if (Qdiff[length(Qdiff)] <= tolerance) {
		    break
		  }
		  g <- g + 1
		}

		if (g == maxtime)
		{
		  print("runtime error")
		}
		alpha1<-mean(alpha1per(Q))
		gamma1<-gamma1per(Q)
		beta1<-gamma1/alpha1

	#res[1,]<-c(mean(alpha1),beta1,gamma1)

	return(c(alpha1,beta1,gamma1))
}


env.file.from.aid<-function(aid){
  if (aid == 5){
    envFile <- "N100Grid_Lowest"
  } else if (aid == 6){
    envFile <- "N100Grid_Mid"
  } else if (aid == 7){
    envFile <- "N100Grid_Highest"
  } else {
    stop("aid not recognised")
  }
  return(envFile)

}




changing.values.a<-function(hval, wval, aid, nu, maxtime = 100000)
{
  filename <- "Grid_undir"
  envFile <- env.file.from.aid(aid)

	h<-hval

	#td <- params$td
	w_adjust <- wval
	params <- Parameters(aid)
	topology <- params$topology
	pinit <- 1 / n_sp

	#res<-matrix(NA,1,3)
	#colnames(res)<-c("alpha","beta","gamma")

		A <- ReadMat(filename, w_adjust)

		E_ls <- ReadEnvs(envFile)
		C <- ReproductionRate(n_sp, n_com, as.numeric(h), E_ls)
		Q <- Initial(n_sp, n_com, pinit)
		Preg <- rowSums(Q) / sum(Q)
		tolerance <- 0.0000001 * n_sp * n_com
		Qdiff <- numeric()
		biomass <- numeric()
		g <- 0
#		browser()
		while (g < maxtime)
		{
		  prevQ <- Q # update the prevQ using the previous Q
		  Q <- Update(Q, A, C, nu, Preg, n_com) # new Q
		  diff <- calcDiff(Q, prevQ)
		  Qdiff <- c(Qdiff, diff)
		  biomass <- c(biomass, sum(Q))
		  if (is.na(Qdiff[length(Qdiff)] )){
		    g<-maxtime
		    break
		  }
		  if (Qdiff[length(Qdiff)] <= tolerance) {
		    break
		  }
		  g <- g + 1
		}


#		browser()
		alpha1<-(alpha1per(Q))
#		gamma1<-gamma1per(Q)
#		beta1<-gamma1/alpha1

		if (g == maxtime) # if the loop reaches the maxtime, then it will print an error message and return NA
		{
		  print("runtime error")
		  alpha1[!is.na(alpha1)] <- NA
		}

	#res[1,]<-c(mean(alpha1),beta1,gamma1)

	return(alpha1)
}



##from Susuki et al. 2021

#works for
# > nu
# [1] 0.5
# > n_com
# [1] 100
# > n_sp
#[1] 127

Parameters <- function(aid)
{

#w_adjust is a parameter to adjust the total dispersal rate by multiplying the matrix. ls_td[2]

	  ls_td <-
	    list(
	      c(0.00005, 0.0005),
	      c(0.0005, 0.005),
	      c(0.005, 0.05),
	      c(0.05, 0.5),
	      c(0.1, 1),
	      c(0.15, 1.5),
	      c(0.2, 2),
	      c(0.25, 2.5),
	      c(0.3, 3),
	      c(0.35, 3.5),
	      c(0.4, 4)
	    )

	# Grid seems the default optoin, then autocorrelation should be

	  tplenv_ls <-
	    list(
	      c("Complete", "Lowest"),
	      c("Linear", "Lowest"),
	      c("Linear", "Mid"),
	      c("Linear", "Highest"),
	      c("Grid", "Lowest"),
	      c("Grid", "Mid"),
	      c("Grid", "Highest"),
	      c("SmallWorld", "Lowest"),
	      c("SmallWorld", "Mid"),
	      c("SmallWorld", "Highest"),
	      c("Tree", "Lowest"),
	      c("Tree", "Mid"),
	      c("Tree", "Highest")
	    )

	  A <- expand.grid(tplenv_ls, ls_td)
	  tplenv <- unlist(A[aid, 1])
	  lstd <- unlist(A[aid, 2])
	  topology <- tplenv[1]
	  autocorr <- tplenv[2]
	  td <- lstd [1]
	  w_adjust <- lstd[2]
	  return(list(
	    topology = topology,
	    autocorr = autocorr,
	    td = td,
	    w_adjust = w_adjust
	  ))
}

Initial <- function(n_sp, n_com, pinit)
{
  # pinit: the initial Pik for all i and k. pinit*n_sp has to be <= 1.
  # The example they used in the paper was n_sp = 20, n_com = 20
  Q <- matrix(pinit, n_sp, n_com)
  return(Q)
}

ReadMat <- function(filename, w_adjust)
{
  # the imported matrix file is for total dispersal 0.1. w_adjust is a parameter to adjust the total dispersal rate by multiplying the matrix.
  mat <- as.matrix(read.table(filename, sep = ' ')) * w_adjust
  return(mat)
}

ReadEnvs <- function(filename)
{
  E <- t(read.table(filename))
  return(E)
}

ReproductionRate <- function(n_sp, n_com, h, E_ls)
{
  X <- matrix(0, n_sp, n_com)
  for (k in 1:n_com) {
    X[, k] <- seq(0, 1, length.out = n_sp)
  }
  E <- matrix(rep(E_ls, n_sp), n_sp, n_com,byrow=TRUE)
  D <- E - X
  C <- exp(-(D * D) / (2 * h)) / sqrt(2 * h * pi)
  return(C)
}

WrightFisher <- function(Q)
{
  Q <- matrix(0, nrow(Q), ncol(Q))
  return(Q)
}

Update <- function(Q, A, C, nu, Preg, n_com)
{
  # <Variables>
  #   Matrix Q (Qik): the fraction of sites occupied by species i in community k.
  # <Parameters>
  #   Matrix A (alk): dispersal rate from community l to k. diag(A) = 0.
  #   Matrix C (cik): per capita potential reproductive rate of species i in community k.

# browser()
  S <- colSums(A) # row sum of the matrix A, transposed into a row vector
  Z <-
    matrix(1, nrow(Q), ncol(Q)) - matrix(rep(S, each = n_sp), n_sp, n_com)
  Z[Z < 0] <- 0
  q <- (C * Q) %*% A + (Z * (C * Q))

  # Death Process
  Q <- WrightFisher(Q)

  # Stochastic sampling of new settlers
  #q <- q / matrix(colSums(q),nrow=1)
  q <- t(t(q)/colSums(q))
  # immigration process
  Vk <- 1
  if (nu > 0) {
    Q <- Q + (nu / n_sp) # Adding immigrants. sum_i(nu/n_sp) = nu. nu is fraction of immigrants per local community.
    Vk <- Vk - nu # fraction of empty sites per community
  }
  Q <- Q + Vk * q
  return(Q)
}

calcDiff <- function(Q, prevQ) {
  diff <- sum(abs(Q - prevQ))
  return(diff)
}



