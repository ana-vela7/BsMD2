#' Posterior Probabilities from Bayesian Screening
#' Experiments
#'
#' Marginal factor posterior probabilities and model
#' posterior probabilities from designed screening
#' experiments are calculated according to Box and
#' Meyer's Bayesian procedure.
#'
#' @param X Matrix. The design matrix.
#' @param y Vector. The response vector.
#' @param p Numeric. Prior probability assigned to active factors.
#' @param gamma Numeric. Variance inflation associated to active factors.
#' @param max_int Integer <= 3. Maximum order of interactions considered in the models.
#' @param max_fac Integer. Maximum number of factors included in the models.
#' @param top Integer. Number of models to keep with the highest posterior probability.
#'
#' @return A list with all the input and output parameters.
#'
#' @return \item{X}{Matrix. The design matrix.}
#' @return \item{y}{Vector. The response vector.}
#' @return \item{n}{Integer. Number of runs.}
#' @return \item{col}{Integer. Number of columns in the design matrix.}
#' @return \item{max_int}{Integer <= 3. Maximum order of interactions considered in the models.}
#' @return \item{max_fac}{Integer. Maximum number of factors included in the models.}
#' @return \item{pi}{Numeric. Prior probability assigned to active factors.}
#' @return \item{gamma}{Numeric. Variance inflation associated with active factors.}
#' @return \item{top}{Integer. Number of models to keep with the highest posterior probability.}
#' @return \item{Prob_fac}{Data frame. Posterior probability for each factor.}
#' @return \item{Prob_mod}{Data frame. Posterior probability for each of the top models.}
#' @return \item{nfac_mod}{Vector. Number of active factors in each of the top models.}
#' @return \item{p_mod}{Vector. Posterior probability for each of the top models.}
#' @return \item{fac_mod}{Matrix. Active factors for each of the top models.}
#' @references {
#' Box, G. E. P and R. D. Meyer (1986). "An Analysis for Unreplicated Fractional Factorials".
#' Technometrics. Vol. 28. No. 1. pp. 11{–}18.
#'
#' Box, G. E. P and R. D. Meyer (1993). "Finding the Active Factors in Fractionated Screening Experiments".
#' Journal of Quality Technology. Vol. 25. No. 2. pp. 94{–}105.}
#' @export
#' @examples
#' #Example 1
#' library(BsMD2)
#' data("BM93e1")
#' X <- as.matrix(BM93e1[,2:6])
#' y <- BM93e1[,7]
#' drillAdvance.BsProb1 <- BsMD2::BsProb1(X, y, .25, 1.6, 3, 5)
#' plot(drillAdvance.BsProb1)
#' summary(drillAdvance.BsProb1)
#'
#' #Example 2
#' data("BM93e2")
#' X <- as.matrix(BM93e2[,1:7])
#' y <- BM93e2[,8]
#' pp <- BsMD2::BsProb1(X, y, .25, 1.5, 3, 7)
#' plot(pp)
#' summary(pp)
#'
#' #Example 3
#' data("BM93e3")
#' X16 <- as.matrix(BM93e3[1:16,2:9])
#' y16 <- BM93e3[1:16,10]
#' pp16 <- BsMD2::BsProb1(X16, y16, .25, 2, 3, 8)
#'
#' X <- as.matrix(BM93e3[,1:9])
#' y <- BM93e3[,10]
#' pp <- BsProb1(X, y, .25, 2, 3, 4)
#'
#'

BsProb1 <- function (X, y, p = 0.25, gamma = 2, max_int = 3, max_fac = ncol(X), top = 10) {

  fac <- ncol(X)

  models <- expand.grid(rep(list(0:1), fac))
  models <- models[rowSums(models) <= max_fac,]
  models_ini <- models

  if(max_int > 1){
    comb <- utils::combn(fac, 2)

    mat <- matrix(0, ncol = ncol(comb), nrow = nrow(models))
    for (j in 1:ncol(comb)) {
      fac1 <- comb[1,j]
      fac2 <- comb[2,j]

      mat[which(models[,fac1] + models[,fac2] == 2), j] <- 1
      X <- cbind (X, X[, fac1]*X[, fac2])
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
      mat[which(models[, fac1] + models[, fac2] + models[, fac3] == 3), j] <- 1
      X <- cbind (X, X[,fac1]*X[,fac2]*X[,fac3])
    }
    models <- cbind(models, mat)
  }

  fac_act <- rowSums(models_ini)
  prior <- (p / (1 - p))^fac_act

  X0 <- rep(1, length(y))
  Sbeta0 <- t(y-mean(y))%*%(y-mean(y))
  det0 <- (det(t(X0)%*%X0))^(.5)
  post_mod <- rep(0, nrow(models))

  for (i in 1:nrow(models)){
    efectos <- which(models[i,]==1)
    tam <- length(efectos)
    Xi <-  as.matrix(cbind(X0, X[ ,efectos]))
    mat <- matrix(0, tam+1, tam+1)
    if(nrow(mat)>1) mat[cbind(2:nrow(mat), 2:nrow(mat))] <- 1
    gammi <- (1/gamma^2)*mat
    betai <- solve(gammi+t(Xi)%*%Xi)%*%t(Xi)%*%y
    Sbetai <- t(y-Xi%*%betai)%*%(y-Xi%*%betai)
    p0 <- prior[i]*gamma^(-tam)
    p1 <- det0/det(gammi+t(Xi)%*%Xi)^(.5)
    p2 <- ((Sbetai + t(betai)%*%gammi%*%betai)/Sbeta0)^(-(length(y)-1)/2)
    post_mod[i] <- p0*p1*p2
  }

  post_fac <- t(as.matrix(post_mod))%*%as.matrix(models_ini)
  post_fac <- c(1, post_fac)

  post_fac <- post_fac/sum(post_mod)
  df_fac <- data.frame(Factor = c("NULL", paste0("F",seq(fac))), Prob = round(post_fac, 3))

  df_mod <- data.frame(indice = 1:nrow(models), Prob = round(post_mod/sum(post_mod),3))
  df_mod <- df_mod[order(-df_mod$Prob),][1:top, ]
  variables <- models_ini[df_mod$indice, ]
  lvar <- apply(variables, 1, function(x) which(x == 1))
  vars <- unlist(sapply(lvar, function(x) if(length(x)==0) "none"
                        else paste(x, collapse = ", ")))
  df_mod$Factors <- vars
  df_mod <- df_mod[,-1]
  rownames(df_mod) <- paste0("M", seq(top))

  fac_mod <- matrix(0, nrow = top, ncol = max_fac)
  for (i in 1:top){
    while (length(lvar[[i]]) < max_fac) lvar[[i]] <- c(lvar[[i]], 0)
    fac_mod[i, ] <- lvar[[i]]
  }
  colnames(fac_mod) <- NULL
  rownames(fac_mod) <- paste0("M", seq(top))

  nfac_mod <- rowSums(variables)
  names(nfac_mod) <- paste0("M", seq(top))

  rownames(X) <- paste0("r", seq(nrow(X)))
  colnames(X) <- paste0("F", seq(ncol(X)))

  p_mod <- df_mod$Prob
  names(p_mod) <- paste0("M", seq(top))

  #methods::setClass("BsProb1", methods::representation(X = "matrix", y = "vector", n = "numeric", col = "numeric", max_int = "numeric", max_fac = "numeric", pi = "numeric", gamma = "numeric", top = "numeric", Prob_fac = "data.frame", Prob_mod = "data.frame", nfac_mod = "vector", p_mod = "vector", fac_mod = "matrix"))

  pp <- methods::new("BsProb1", X = X[,1:fac], y = y, n = length(y), col = fac, max_int = max_int, max_fac = max_fac, pi = p, gamma = gamma, top = top, Prob_fac = df_fac, Prob_mod = df_mod, p_mod = p_mod, nfac_mod = nfac_mod, fac_mod = fac_mod)
  return(pp)
}

#' Print BsProb1 Object
#'
#' @param x BsProb1 object.
#' @param ... additional parameters
#'
#' @return Prints BsProb1 object
#' @export
print.BsProb1 <- function(x, ...){
  cat("Design Matrix:\n")
  print(x@X)
  cat("\nResponse vector:\n")
  print(x@y)
  cat("\nCalculations:\n")
  calc <- c(x@n, x@col, x@max_int, x@max_fac, x@pi, x@gamma)
  names(calc) <- c("Runs", "Factors", "max_int", "max_fac", "p", "g")
  print(calc)
  cat("\nFactor Probabilities:\n")
  print(x@Prob_fac)
  cat("\nModel Probabilities:\n")
  print(x@Prob_mod)
}

#' Summary for BsProb1 Object
#'
#' @param object BsProb1 Object
#' @param ... additional parameters
#'
#' @return Summary for BsProb1 object.
#' @export
summary.BsProb1 <- function(object, ...){
  x <- object
  cat("\nCalculations:\n")
  calc <- c(x@n, x@col, x@max_int, x@max_fac, x@pi, x@gamma)
  names(calc) <- c("Runs", "Factors", "max_int", "max_fac", "p", "g")
  print(calc)
  cat("\nFactor Probabilities:\n")
  print(x@Prob_fac)
  cat("\nModel Probabilities:\n")
  print(x@Prob_mod)
}

#' Plot BsProb1 object
#'
#' @param x Object of class BsProb1.
#' @param fac Logical. If true, the plot of each factor's posterior probability is displayed.
#' @param ... additional graphical parameters passed to plot.
#'
#' @return Plot of BsProb1 object
#' @export
plot.BsProb1 <- function(x, fac = TRUE, ...){
  if (fac) graphics::barplot(x@Prob_fac$Prob, names.arg = x@Prob_fac$Factor, main = "Factor Probabilities", xlab = "Factor", ylab = "Prob")
  else graphics::barplot(x@Prob_mod$Prob, names.arg = rownames(x@Prob_mod), main = "Model Probabilities", xlab = "Model", ylab = "Prob")
}
