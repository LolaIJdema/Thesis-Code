
clustMD_S=function (X, G, CnsIndx, OrdIndx, Nnorms, MaxIter, model, store.params = FALSE, 
          scale = FALSE, startCL = "hc_mclust", autoStop = FALSE, ma.band = 50, 
          stop.tol = NA) 
{
  Y <- as.matrix(X) #data
  N <- nrow(Y) #observations
  J <- ncol(Y) #variables
  if (scale) {
    if (CnsIndx > 0) #when there are cont. variables
      Y[, 1:CnsIndx] <- scale(Y[, 1:CnsIndx]) #scale them
  }
  K <- apply(Y, 2, max) #number of groups
  if (CnsIndx > 0) #except for continuous var
    K[1:CnsIndx] <- NA
  D <- J
  if (J > OrdIndx) 
    D <- OrdIndx + sum(K[(OrdIndx + 1):J] - 1)
  if (J > OrdIndx) {
    nom.ind.Z <- vector("list", J - OrdIndx)
    for (j in 1:(J - OrdIndx)) {
      if (j == 1) {
        start <- OrdIndx + 1
      }else {
        start <- OrdIndx + sum(K[(OrdIndx + 1):(OrdIndx + 
                                                  j - 1)] - 1) + 1
      }
      finish <- start + K[OrdIndx + j] - 2
      nom.ind.Z[[j]] <- c(start:finish)
    }
  }
  if ((model == "BD") & (OrdIndx > CnsIndx)) {
    if ((OrdIndx - CnsIndx) == 1) {
      patt.indx <- list()
      for (p in 1:max(Y[, OrdIndx])) {
        patt.indx[[p]] <- which(Y[, OrdIndx] == p)
      }
    }else {
      patt.tab <- data.frame(table(data.frame((Y[, (CnsIndx + 
                                                      1):OrdIndx]))))
      patt.tab <- patt.tab[patt.tab$Freq != 0, 1:(OrdIndx - 
                                                    CnsIndx)]
      patt.indx <- list()
      for (p in 1:nrow(patt.tab)) {
        patt.indx[[p]] <- which(apply(Y[, (CnsIndx + 
                                             1):OrdIndx], 1, clustMD::patt.equal, patt.tab[p, ]))
      }
    }
  }
  Ynames <- colnames(X)
  VarNames <- as.character(1:D)
  if (!is.null(Ynames)) {
    VarNames[1:OrdIndx] <- Ynames[1:OrdIndx]
    if (J > OrdIndx) {
      NomNames <- list()
      for (j in (OrdIndx + 1):J) {
        NomNames[[j - OrdIndx]] <- rep(NA, (K[j] - 1))
        for (k in 1:(K[j] - 1)) {
          NomNames[[j - OrdIndx]][k] <- paste(Ynames[j], 
                                              "_", k, sep = "")
        }
      }
      VarNames[(OrdIndx + 1):D] <- unlist(NomNames)
    }
  }else {
    for (j in 1:OrdIndx) {
      VarNames[j] <- paste("V", j, sep = "")
    }
    if (J > OrdIndx) {
      NomNames <- list()
      for (j in (OrdIndx + 1):J) {
        NomNames[[j - OrdIndx]] <- rep(paste("V", j, 
                                             sep = ""), (K[j] - 1))
        for (k in 1:(K[j] - 1)) {
          NomNames[[j - OrdIndx]][k] <- paste("V", j, 
                                              "_", k, sep = "")
        }
      }
      VarNames[(OrdIndx + 1):D] <- unlist(NomNames)
    }
  }
  VarNames_sht <- as.character(1:D)
  if (!is.null(Ynames)) {
    VarNames_sht[1:OrdIndx] <- substr(Ynames[1:OrdIndx], 
                                      1, 7)
    if (J > OrdIndx) {
      NomNames_sht <- list()
      for (j in (OrdIndx + 1):J) {
        NomNames_sht[[j - OrdIndx]] <- rep(NA, (K[j] - 
                                                  1))
        for (k in 1:(K[j] - 1)) {
          NomNames_sht[[j - OrdIndx]][k] <- paste(substr(Ynames[j], 
                                                         1, 7), "_", k, sep = "")
        }
      }
      VarNames_sht[(OrdIndx + 1):D] <- unlist(NomNames_sht)
    }
  }else {
    VarNames_sht <- VarNames
  }
  Ez <- array(NA, c(N, D, G))
  for (g in 1:G) Ez[, 1:J, g] <- Y
  if (OrdIndx > CnsIndx) {
    perc.cut <- perc.cutoffs(CnsIndx, OrdIndx, Y, N)
    zlimits <- array(NA, c(N, J, 2))
    zlimits[, 1:CnsIndx, 1] <- -Inf
    zlimits[, 1:CnsIndx, 2] <- Inf
    for (j in (CnsIndx + 1):OrdIndx) {
      for (k in 1:K[j]) {
        zlimits[Y[, j] == k, j, 1] <- perc.cut[[j]][k]
        zlimits[Y[, j] == k, j, 2] <- perc.cut[[j]][k + 
                                                      1]
      }
    }
  }else {
    perc.cut <- list()
    zlimits <- array(NA, c(N, J, 2))
  }
  Zstart <- function(Kj, y) {
    new.z <- rep(0, (Kj - 1))
    if (y == 1) {
      new.z <- msm::rtnorm((Kj - 1), mean = 0, sd = 1, 
                           upper = 0)
    }else {
      new.z[-(y - 1)] <- stats::rnorm((Kj - 2), mean = 0, 
                                      sd = 1)
      new.z[(y - 1)] <- msm::rtnorm(1, mean = 0, sd = 1, #truncated normal
                                    lower = max(new.z))
    }
    new.z
  }
  Zinit <- matrix(NA, N, D)
  Zinit[, 1:OrdIndx] <- Y[, 1:OrdIndx]
  if (J > OrdIndx) {
    for (j in (OrdIndx + 1):J) {
      for (i in 1:N) {
        Zinit[i, nom.ind.Z[[j - OrdIndx]]] <- Zstart(K[j], 
                                                     Y[i, j])
      }
    }
  }
  
  #####################################################################
  #Take sample of the data here
  #This sample replaces the data Y
  samp=sample(nrow(X), 0.2*nrow(X))
  X_samp=X[samp, ]
  Y_samp=Y[samp, ]
  
  Zinit_samp=Zinit[samp,]
  
  N_samp=as.integer(N*0.2)
  #####################################################################
  
  if (startCL == "kmeans") {
    if (CnsIndx > 0) {
      ind <- kmeans(Y_samp[, 1:CnsIndx], G)$cl
    }else {
      ind <- kmeans(Y_samp, G)$cl
    }
  }else if (startCL == "hclust") {
    temp <- hclust(dist(Y_samp))
    ind <- cutree(temp, G)
  }else if (startCL == "mclust") {
    if (CnsIndx > 0) {
      if (CnsIndx == 1) {
        ind <- mclust::Mclust(Y_samp[, 1:CnsIndx], G, "V")$cl
      }else {
        ind <- mclust::Mclust(Y_samp[, 1:CnsIndx], G, "VVV")$cl
      }
    }else {
      if (J == 1) {
        ind <- mclust::Mclust(Y_samp, G, "V")$cl
      }else {
        ind <- mclust::Mclust(Y_samp, G, "VVV")$cl
      }
    }
  }else if (startCL == "random") {
    ind <- sample(1:G, N_samp, replace = TRUE)
  }else if (startCL == "hc_mclust") {
    if (CnsIndx > 0) {
      ind <- mclust::hclass(mclust::hc(modelName = "VVV", 
                                       data = Y_samp[, 1:CnsIndx]), G)
    }else {
      stop("Cannot use hc_mclust since there are no continuous variables")
    }
  }else {
    stop("Unknown starting value algorithm chosen! \n \n        Choose from: kmeans, hclust, hc_mclust, mclust or random")
  }
  pi.vec <- table(ind)/N_samp
  mu <- matrix(NA, D, G)
  for (g in 1:G) mu[, g] <- colMeans(matrix(Zinit_samp[ind == g,], sum(ind == g), D))
  Sigma <- array(NA, c(D, D, G))
  for (g in 1:G) Sigma[, , g] <- diag(D)
  a <- matrix(1, G, D)
  likeStore <- rep(NA, MaxIter)
  if (store.params == TRUE) {
    ind.store <- matrix(NA, N_samp, MaxIter)
    Ez.store <- array(NA, c(N_samp, D, G, MaxIter))
    tau.store <- array(NA, c(N_samp, G, MaxIter))
    mu.store <- array(NA, c(D, G, MaxIter))
    lambda.store <- array(NA, c(G, D, MaxIter))
    a.store <- array(NA, c(G, D, MaxIter))
    Sigma.store <- array(NA, c(D, D, G, MaxIter))
    if (J > OrdIndx) 
      probs.nom.store <- array(NA, c(J - OrdIndx, max(K[(OrdIndx + 
                                                           1):J]), G, MaxIter))
  }
  #purely aesthetic
  pb <- txtProgressBar(style = 3)
  prog <- 1
  for (iter in 1:MaxIter) {
    if (iter > prog * MaxIter/10) 
      prog <- prog + 1
    setTxtProgressBar(pb, prog * 0.1)
    if (J > OrdIndx) 
      norms <- MASS::mvrnorm(Nnorms, mu = rep(0, max(
        K[(OrdIndx + 1):J]) - 1), 
        Sigma = diag(max(K[(OrdIndx + 1):J]) - 1))
    
    temp.E <- E.step(N_samp, G, D, CnsIndx, OrdIndx, zlimits, 
                     mu, Sigma, Y_samp, J, K, norms, nom.ind.Z, patt.indx, 
                     pi.vec, model, perc.cut)
    tau <- temp.E[[1]]
    sumTauEz <- temp.E[[2]]
    sumTauS <- temp.E[[3]]
    probs.nom <- temp.E[[4]]
    Ez <- temp.E[[5]]
    ind <- mclust::map(tau)
    temp.M <- M.step(tau, N_samp, sumTauEz, J, OrdIndx, D, G, 
                     Y_samp, CnsIndx, sumTauS, model, a, nom.ind.Z)
    pi.vec <- temp.M[[1]]
    mu <- temp.M[[2]]
    lambda <- temp.M[[3]]
    a <- temp.M[[4]]
    Sigma <- temp.M[[5]]
    likeStore[iter] <- ObsLogLikelihood(N_samp, CnsIndx, G, Y_samp, 
                                        mu, Sigma, pi.vec, patt.indx, zlimits, J, OrdIndx, 
                                        probs.nom, model, perc.cut, K)
    if (store.params == TRUE) {
      ind.store[, iter] <- ind
      tau.store[, , iter] <- tau
      mu.store[, , iter] <- mu
      lambda.store[, , iter] <- lambda
      a.store[, , iter] <- a
      Sigma.store[, , , iter] <- Sigma
    }
    if (autoStop & (iter > 5)) {
      autoStop.aitken <- (OrdIndx == J)
      autoStop.ma <- (OrdIndx < J) & (iter > (ma.band + 
                                                10)) & (iter%%10 == 0)
      if (autoStop.aitken) {
        if (is.na(stop.tol)) {
          stop("stop.tol not specified.")
        }
        check.diff <- (likeStore[iter] - likeStore[iter - 
                                                     1])/abs(1 + likeStore[iter])
        if (check.diff < stop.tol) {
          setTxtProgressBar(pb, 1)
          break
        }
      }
      if (autoStop.ma) {
        if (is.na(stop.tol)) {
          stop("stop.tol or ma.band not specified.")
        }
        MA.lag <- 10
        likeOld <- mean(likeStore[(iter - ma.band - MA.lag):(iter - 
                                                               MA.lag)])
        likeNew <- mean(likeStore[(iter - ma.band):iter])
        if (abs((likeNew - likeOld)/likeNew) < stop.tol) {
          setTxtProgressBar(pb, 1)
          break
        }
      }
      if (autoStop.ma) {
        if (is.na(stop.tol)) {
          stop("stop.tol or ma.band not specified.")
        }
        MA.lag <- 10
        likeOld <- mean(likeStore[(iter - ma.band - MA.lag):(iter - 
                                                               MA.lag)])
        likeNew <- mean(likeStore[(iter - ma.band):iter])
        if (abs((likeNew - likeOld)/likeNew) < stop.tol) {
          setTxtProgressBar(pb, 1)
          break
        }
      }
    }
  } #end for loop 
  
  ##################################################################
  #run one iteration over dataset here
  
  temp.E <- E.step(N, G, D, CnsIndx, OrdIndx, zlimits, 
                   mu, Sigma, Y, J, K, norms, nom.ind.Z, patt.indx, 
                   pi.vec, model, perc.cut)
  tau <- temp.E[[1]]
  sumTauEz <- temp.E[[2]]
  sumTauS <- temp.E[[3]]
  probs.nom <- temp.E[[4]]
  Ez <- temp.E[[5]]
  ind <- mclust::map(tau)
  temp.M <- M.step(tau, N, sumTauEz, J, OrdIndx, D, G, 
                   Y, CnsIndx, sumTauS, model, a, nom.ind.Z)
  pi.vec <- temp.M[[1]]
  mu <- temp.M[[2]]
  lambda <- temp.M[[3]]
  a <- temp.M[[4]]
  Sigma <- temp.M[[5]]
  like <- ObsLogLikelihood(N, CnsIndx, G, Y, 
                           mu, Sigma, pi.vec, patt.indx, zlimits, J, OrdIndx, 
                           probs.nom, model, perc.cut, K)
  if (model == "BD") {
    probs.nom <- z.moments(D, G, N, CnsIndx, OrdIndx, zlimits, 
                           mu, Sigma, Y, J, K, norms, nom.ind.Z, patt.indx)[[2]]
  }else {
    probs.nom <- z.moments_diag(D, G, N, CnsIndx, OrdIndx, 
                                zlimits, mu, Sigma, Y, J, K, norms, nom.ind.Z)[[2]]
  }
  
  obslike <- ObsLogLikelihood(N, CnsIndx, G, Y, mu, Sigma, 
                              pi.vec, patt.indx, zlimits, J, OrdIndx, probs.nom, model, 
                              perc.cut, K)
  BIChat <- 2 * obslike - npars_clustMD(model, D, G, J, CnsIndx, 
                                        OrdIndx, K) * log(N)
  
  close(pb)
  CompleteLike_i <- rep(NA, N)
  for (i in 1:N) {
    CompleteLike_i[i] <- log(pi.vec[ind[i]]) + mvtnorm::dmvnorm(Ez[i, 
                                                                   , ind[i]], mean = mu[, ind[i]], sigma = matrix(Sigma[, 
                                                                                                                        , ind[i]], nrow = D, ncol = D), log = TRUE)
  }
  ICLhat <- 2 * sum(CompleteLike_i) - npars_clustMD(model, 
                                                    D, G, J, CnsIndx, OrdIndx, K) * log(N)
  rownames(mu) <- VarNames
  rownames(Sigma) <- VarNames
  colnames(Sigma) <- VarNames
  if (store.params == TRUE) {
    params.store.list <- list(cl.store = ind.store, tau.store = tau.store, 
                              means.store = mu.store, A.store = a.store, lambda.store = lambda.store, 
                              Sigma.store = Sigma.store)
  }
  else {
    params.store.list <- NULL
  }
  
  out.clustMD_S <- list(model = model, G = G, Y = Y, cl = ind, 
                      tau = tau, means = mu, A = a, Lambda = lambda, Sigma = Sigma, 
                      BIChat = BIChat, ICLhat = ICLhat, paramlist = params.store.list, 
                      VarNames = VarNames, VarNames_sht = VarNames_sht, likelihood.store = like)
  class(out.clustMD_S) <- "clustMD_S"
  out.clustMD_S
}
