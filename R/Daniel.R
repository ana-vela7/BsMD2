#' Normal Plot of Effects
#'
#' Normal plot of effects from a two level factorial experiment.
#'
#' @param fit Object of class [lm]. Fitted model from lm.
#' @param half Logical. If TRUE, half-normal plot of effects is displayed.
#'
#' @return The function returns the Normal Probability plot or Half-normal plot of the factor effects.
#' @export
#'
#' @examples
#' #Example 6.2 Douglas C. Montgomery
#' #A <- rep(c(-1,1),8)
#' #B <- rep(c(-1,-1,1,1),4)
#' #C <- rep(c(rep(-1,4),rep(1,4)),2)
#' #D <- c(rep(-1,8), rep(1,8))
#' #rate <- c(45,71,48,65,68,60,80,65,43,100,45,104,75,86,70,96)
#' #fit <- lm(rate~A*B*C*D)
#' #Daniel(fit)
#' #Daniel(fit, half = TRUE)
#'
#' #Example 8.6 Douglas C. Montgomery
#' #A <- rep(c(-1,1), 16)
#' #B <- rep(c(-1,-1,1,1), 8)
#' #C <- rep(c(rep(-1, 4), rep(1, 4)), 4)
#' #D <- rep(c(rep(-1,8),rep(1,8)),2)
#' #E <- c(rep(-1,16),rep(1,16))
#' #F <- A*B*C
#' #G <- A*B*D
#' #H <- B*C*D*E
#'
#' #data <- as.data.frame(cbind(A,B,C,D,E,F,G,H))
#' #Blocks <- c(3,2,4,1,1,4,2,3,1,4,2,3,3,2,4,1,2,3,1,4,4,1,3,2,4,1,3,2,2,3,1,4)
#' #s <- c(2.76,6.18,2.43,4.01,2.48,5.91,2.39,3.35,4.4,4.1,3.22,3.78,5.32,3.87,3.03,2.95,
#' #2.64,5.5,2.24,4.28,2.57,5.37,2.11,4.18,3.96,3.27,3.41,4.3,4.44,3.65,4.41,3.4)
#' #y <- log(s)
#'
#' #fit2 <- lm(y ~ Blocks + .^2, data = data)
#' #Daniel(fit2)
#' #Daniel(fit2, half = TRUE)

Daniel <- function(fit, half = FALSE){
  effects <- 2*stats::na.omit(stats::coef(fit)[-1])
  n <- length(effects)
  y <- 100*(1:n-0.5)/n
  if (half){
    x <- sort(abs(effects))
    main <- "Half-normal plot"
    }
  else{
    x <- sort(effects)
    main <- "Normal Probability plot"
  }
  graphics::plot(x, y, pch = 20, main = main, xlim = c(x[1], 2*x[n]-x[n-1]), xlab = "Effects", ylab = "Normal Score")
  graphics::text(x, y, names(x), cex = 0.6, pos = 4, col = "blue")
}
