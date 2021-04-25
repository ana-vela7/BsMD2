#' Lenth's Plot of Effects
#'
#' Plot of the factor effects with significance levels based on robust estimation of contrast standard errors.
#'
#' @param fit Object of class [lm]. Fitted model from lm.
#' @return The function returns Lenth's plot. Dashed lines correspond to the margin of error (ME) and dotted lines to the simultaneous margin of error (SME).
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
#' #Lenth(fit)

Lenth <- function(fit){
  effects <- 2*stats::na.omit(stats::coef(fit)[-1])
  n <- length(effects)
  med <- stats::median(abs(effects))
  s0 <- med*1.5
  PSE <- 1.5*stats::median(abs(effects[abs(effects)<2.5*s0]))
  m <- 2^4-1
  gamma <- 1-(1+.95^(1/m))/2
  d <- m/3
  ME <- stats::qt(.975,d)*PSE
  SME <- stats::qt(1-gamma,d)*PSE

  plot(-1, 0, xlim = c(0, n+1), ylim = c(-max(abs(effects)), max(abs(effects))), main = "Lenth's Plot", xlab = "", ylab = "Effects", xaxt = "none")
  graphics::axis(1, labels = names(effects), at = 1:n, las = 2)
  for (i in 1:n) graphics::lines(c(i, i), c(0, effects[i]), lwd = 2)
  graphics::abline(0, 0)
  graphics::abline(ME,0, lty = 2)
  graphics::abline(-ME,0, lty = 2)
  graphics::abline(SME,0, lty = 3)
  graphics::abline(-SME,0, lty = 3)
}
