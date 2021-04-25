#' Best Model Discrimination (MD) Follow-Up Experiments
#'
#' Best follow-up experiments based on the MD criterion are suggested to
#' discriminate between competing models.
#'
#' @param X Matrix. Design matrix of the initial experiment.
#' @param y Vector. Response vector of the initial experiment.
#' @param Xcand Matrix. Candidate runs to be chosen for the follow-up design.
#' @param nMod Integer. Number of competing models.
#' @param p_mod Vector. Posterior probabilities of the competing models.
#' @param fac_mod Matrix. Active factors in the competing models.
#' @param nFDes Integer. Number of runs to consider in the follow-up experiment.
#' @param max_int Integer. Maximum order of interactions in the models.
#' @param g Numeric. Variance inflation factor for active effects.
#' @param Iter Integer. Maximum number of iterations for each search.
#' @param nStart Integer. Number of random starting designs.
#' @param top Integer. Highest MD follow-up designs recorded.
#'
#' @return A list with all the input and output parameters.
#' @return \item{X}{Matrix. The design matrix.}
#' @return \item{y}{Vector. The response vector.}
#' @return \item{Xcand}{Matrix. Candidate runs to be chosen for the follow-up design.}
#' @return \item{Runs}{Integer. Number of runs.}
#' @return \item{Fac}{Integer. Number of factors.}
#' @return \item{nMod}{Integer. Number of competing models.}
#' @return \item{p_mod}{Vector. Posterior probabilities of the competing models.}
#' @return \item{fac_mod}{Matrix. Active factors in the competing models.}
#' @return \item{nFDes}{Integer. Number of runs to consider in the follow-up experiment.}
#' @return \item{max_int}{Integer. Maximum order of the interactions in the models.}
#' @return \item{g}{Numeric. Variance inflation factor for active effects.}
#' @return \item{Iter}{Integer. Maximum number of iterations for each search.}
#' @return \item{nStart}{Integer. Number of random starting designs.}
#' @return \item{top}{Integer. Highest MD follow-up designs recorded.}
#' @return \item{MD}{Data frame. Designs points and MD for top designs.}
#' @return \item{MDtop}{Vector. MD for top designs.}
#' @return \item{DEStop}{Data frame. Top design points.}
#' @references {
#' Meyer, R. D., Steinberg, D. M. and Box, G. E. P. (1996). "Follow-Up Designs to Resolve
#' Confounding in Multifactor Experiments (with discussion)". Technometrics, Vol. 38,
#' No. 4, pp. 303{-}332.
#'
#' Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors in Fractionated
#' Screening Experiments". Journal of Quality Technology. Vol. 25. No. 2. pp. 94{–}105.
#'
#' }
#' @export
#'
#' @examples
#' #Example 1
#' library(BsMD2)
#' data(BM93e3)
#' X <- as.matrix(BM93e3[1:16,c(1,2,4,6,9)]) #matriz de diseño inicial
#' y <- as.vector(BM93e3[1:16,10]) #vector de respuesta
#' p_mod <- c(0.2356,0.2356,0.2356,0.2356,0.0566) #probabilidad posterior de los 5 modelos
#' fac_mod <- matrix(c(2,1,1,1,1,3,3,2,2,2,4,4,3,4,3,0,0,0,0,4),nrow=5,
#'                   dimnames=list(1:5,c("f1","f2","f3","f4")))
#' Xcand <- matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
#'                   -1,-1,-1,-1,1,1,1,1,-1,-1,-1,-1,1,1,1,1,
#'                   -1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,
#'                   -1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,
#'                   -1,1,1,-1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1),
#'                 nrow=16,dimnames=list(1:16,c("blk","f1","f2","f3","f4"))
#' )
#' injectionMolding <- MDopt(X = X, y = y, Xcand = Xcand, nMod = 5, p_mod = p_mod, fac_mod = fac_mod,
#' nStart = 25)
#'
#' #Example 2
#' data(M96e2,package="BsMD2")
#' X <- as.matrix(cbind(blk = rep(-1,8), M96e2[c(25,2,19,12,13,22,7,32), 1:5]))
#' y <- M96e2[c(25,2,19,12,13,22,7,32), 6]
#' pp <- BsProb1(X = X[,2:6], y = y, p = .25, gamma = .4, max_int = 3, max_fac = 5, top = 32)
#' p <- pp@p_mod
#' facs <- pp@fac_mod
#' Xcand <- as.matrix(cbind(blk = rep(+1,32), M96e2[,1:5]))
#' #e2 <- MDopt(X = X, y = y, Xcand = Xcand, nMod = 32, p_mod = p, fac_mod = facs, g = .4,
#' #Iter = 10, nStart = 25, top = 5)
MDopt <- function(X, y, Xcand, nMod, p_mod, fac_mod, nFDes = 4, max_int = 3, g = 2, Iter = 20, nStart = 10, top = 10){

  n <- length(y)
  fac <- ncol(X)-1

  Si <- list()
  Xi <- list()
  betai <- list()
  gammi <- list()
  efectos <- list()

  models <- matrix(0, nMod, fac)
  for (i in 1:nMod){
    models[i, fac_mod[i,][fac_mod[i,] != 0]] <- 1
  }

  Xfac <- X[ , -1]
  Xc <- Xcand[ , -1]

  if(max_int > 1){
    comb <- utils::combn(fac, 2)
    mat <- matrix(0, ncol = ncol(comb), nrow = nrow(models))
    for (j in 1:ncol(comb)) {
      fac1 <- comb[1,j]
      fac2 <- comb[2,j]
      mat[which(models[, fac1] + models[, fac2] == 2), j] <-  1
      Xfac <- cbind (Xfac, Xfac[, fac1]*Xfac[, fac2])
      Xc <- cbind (Xc, Xc[, fac1]*Xc[, fac2])
    }
    models <- cbind(models, mat)
  }

  if(max_int > 2) {
    comb <- utils::combn(fac, 3)
    mat <- matrix(0, ncol = ncol(comb), nrow = nrow(models))
    for (j in 1:ncol(comb)) {
      fac1 <- comb[1,j]
      fac2 <- comb[2,j]
      fac3 <- comb[3,j]
      mat[which(models[, fac1] + models[, fac2] + models[, fac3] == 3), j] <-  1
      Xfac <- cbind (Xfac, Xfac[ , fac1]*Xfac[ , fac2]*Xfac[ , fac3])
      Xc <- cbind (Xc, Xc[ , fac1]*Xc[ , fac2]*Xc[ , fac3])
    }
    models <- cbind(models, mat)
  }

  Xfac <- cbind(X[,1], Xfac)
  Xc <- cbind(Xcand[,1], Xc)
  models <- cbind(rep(1, nMod), models)

  for (i in 1:nMod){
    efectos[[i]] <- which(models[i,]==1)
    tam <- length(efectos[[i]])
    Xi[[i]] <-  cbind(rep(1, n), Xfac[ ,efectos[[i]]])

    mat <- matrix(0, tam+1, tam+1)
    if(nrow(mat)>1) mat[cbind(2:nrow(mat), 2:nrow(mat))] <- 1
    gammi[[i]] <- (1/g^2)*mat
    #gammi[1,1] <- 0.00001
    betai[[i]] <- solve(gammi[[i]]+t(Xi[[i]])%*%Xi[[i]])%*%t(Xi[[i]])%*%y
    Si[[i]] <- t(y-Xi[[i]]%*%betai[[i]])%*%(y-Xi[[i]]%*%betai[[i]]) + t(betai[[i]])%*%gammi[[i]]%*%betai[[i]]
  }

  MDr <- function(extra){
    nex <- length(extra)
    y_gorro_estrella <- list()
    V_estrella <- list()

    for (i in 1:nMod){
      Xiestrella <- cbind(rep(1, nex), Xc[extra, efectos[[i]]])
      y_gorro_estrella[[i]] <- Xiestrella%*%betai[[i]]
      V_estrella[[i]] <- diag(1, nex, nex) + Xiestrella%*%solve(gammi[[i]] + t(Xi[[i]])%*%Xi[[i]])%*%t(Xiestrella)
    }

    MD <- 0
    m <- 1:nMod

    for (i in m){
      for (j in m[m != i]){
        MD <- MD + p_mod[i]*p_mod[j]*(-nex + sum(diag(solve(as.matrix(V_estrella[[j]]))%*%as.matrix(V_estrella[[i]]))) + (n - 1)*((t(as.vector(y_gorro_estrella[[i]]-y_gorro_estrella[[j]]))%*%solve(as.matrix(V_estrella[[j]]))%*%(as.vector(y_gorro_estrella[[i]]-y_gorro_estrella[[j]])))/Si[[i]]))
      }
    }
    MD <- MD*.5
    return(MD)
  }

  df_MD <- data.frame(DesignPoints = 0, MD = 0)

  for (j in 1:nStart){
    extra <- sample(1:nrow(Xcand),nFDes, replace = TRUE)
    iter <- 1
    last_out <- 0
    last_in <- 1
    while(last_out != last_in && iter < Iter){

      dp <- paste(sort(extra), collapse = " ")
      if (length(which(df_MD[,1] == dp)) != 0) break

      df_MD <- rbind(df_MD, c(dp, round(MDr(extra), 2)))

      op <- numeric(0)
      for(i in 1:nFDes){
        op[i] <- MDr(extra[-i])
      }
      index <- which(op == max(op))[1]
      last_out <- extra[index]
      extra <- extra[-index]

      op <- numeric(0)
      for(i in 1:nrow(Xcand)){
        op[i] <- MDr(c(extra, i))
      }
      last_in <- which(op == max(op))[1]
      extra <- c(extra, last_in)

      iter <- iter+1
    }
  }

  df_MD <- df_MD[order(as.numeric(df_MD$MD), decreasing = TRUE),]
  df_MD <- df_MD[1:top,]
  rownames(df_MD) <- paste0("D",seq(top))

  colnames(X) <- c("blk", paste0("F", seq(fac)))
  rownames(X) <- paste0("r", seq(n))

  colnames(Xcand) <- c("blk", paste0("F", seq(fac)))
  rownames(Xcand) <- paste0("r", seq(nrow(Xcand)))

  names(p_mod) <- paste0("M", seq(nMod))

  colnames(fac_mod) <- NULL
  rownames(fac_mod) <- paste0("M", seq(nMod))

  MDtop <- as.numeric(df_MD$MD)

  DEStop <- as.data.frame(df_MD$DesignPoints)
  colnames(DEStop) <- "DesignPoints"
  rownames(DEStop) <- paste0("D", seq(top))

  MD <- methods::new("MDopt", X = X, y = y, Xcand = Xcand, Runs = n, Fac = fac, nMod = nMod, p_mod = p_mod, fac_mod = fac_mod, nFDes = nFDes, max_int = max_int, g = g, Iter = Iter, nStart = nStart, top = top, MD = df_MD, MDtop = MDtop, DEStop = DEStop)

  return(MD)
}

#' Print MDopt Object
#'
#' @param x MDopt object
#' @param ... additional parameters
#'
#' @return Prints MDopt object
#' @export
print.MDopt <- function(x, ...){
  cat("\nBase:\n")
  base <- c(x@Runs, x@Fac, x@max_int, x@g, x@nMod)
  names(base) <- c("Runs", "Fac", "max_int", "g", "nMod")
  print(base)
  cat("\nFollow up:\n")
  followup <- c(nrow(x@Xcand), x@nFDes, x@Iter, x@nStart)
  names(followup) <- c("nCand", "Runs", "Iter", "nStart")
  print(followup)
  cat("\nCompeting Models:\n")
  factors <- vector()
  for (i in 1:x@nMod) factors[i] <- paste(x@fac_mod[i,][x@fac_mod[i,]!=0], collapse = ", ")
  df_mod <- as.data.frame(cbind(x@p_mod, factors))
  colnames(df_mod) <- c("Prob", "Factors")
  rownames(df_mod) <- paste0("M", seq(x@nMod))
  print(df_mod)
  cat("\nCandidate runs:\n")
  print(x@Xcand)
  cat("\nTop runs:\n")
  print(x@MD)
}
#' Summary for MDopt Object
#'
#' @param object MDopt Object
#' @param ... additional parameters
#'
#' @return Summary for MDopt object.
#' @export
summary.MDopt <- function(object, ...){
  x <- object
  cat("\nBase:\n")
  base <- c(x@Runs, x@Fac, x@max_int, x@g, x@nMod)
  names(base) <- c("Runs", "Fac", "max_int", "g", "nMod")
  print(base)
  cat("\nFollow up:\n")
  followup <- c(nrow(x@Xcand), x@nFDes, x@Iter, x@nStart)
  names(followup) <- c("nCand", "Runs", "Iter", "nStart")
  print(followup)
  cat("\nTop runs:\n")
  print(x@MD)
}
