rfe.plotcv <- function(fit, figure = "both", cex = 1)
{
        par(mar = c(5, 5, 5, 1))
        if(figure == "both") {
                par(mfrow = c(2, 1))
        }
        nc <- nrow(fit$individualErrors)
        if((figure == "both") | (figure == "total")) {
                plot(logb(fit$nFeatures, 2), fit$errorCV, ylim = c(-0.1, 0.8),
                        xlab = "log base 2 of Number of Genes", ylab =
                        "Overall Error", type = "n", yaxt = "n")
                lines(logb(fit$nFeatures, 2), fit$errorCV, col = 2)
                axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
                o <- fit$errorCV == min(fit$errorCV)
                points(logb(fit$nFeatures[o], 2), fit$errorCV[o], pch = "x")
                error.bars(logb(fit$nFeatures, 2), fit$errorCV - fit$errorSE,
                        fit$errorCV + fit$errorSE)
        }
        if((figure == "both") | (figure == "individual")) {
                plot(logb(fit$nFeatures, 2), fit$individualErrors[1,  ], ylim = c(
                        -0.1, 1.1), xlab = "log base 2 of Number of Genes",
                        ylab = "Class-Specific Errors", type = "n", yaxt = "n")
                axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
                for(i in 1:nrow(fit$individualErrors)) {
                        lines(logb(fit$nFeatures, 2), fit$individualErrors[i,  ],
                                lty = i)
                }
                legend(min(logb(fit$nFeatures, 2)) + ((max(logb(fit$nFeatures,
                        2)) - min(logb(fit$nFeatures, 2))) * 2)/3, 1, dimnames(
                        fit$individualErrors)[[1]], lty = (1:nc), xjust = 0.5, yjust
                         = 1, cex = cex, ncol = 1)
        }
        par(mfrow = c(1, 1))
}
