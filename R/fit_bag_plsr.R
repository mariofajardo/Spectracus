#' Fit bootstrap aggregated PLSR
#' 
#' @author Budiman Minasny
#' 
#' @param soilv soil variables of interest
#' @param spec matrix or dataframe containing spectra. Each row represents a spectrum, each column a wavelength.
#' @param nbag integer. Number of bootstrap parameters.
#' @param maxc integer. Maximum number of components for PLSR.
#' 
#' @return A list containing:
#' \item{model}{PLS models}
#' \item{oob_rmse}{mean out-of-bag RMSE}
#' \item{cal_rmse}{mean calibration RMSE}
#' 
#' @seealso \code{\link{predict_bag_plsr}}
#' 
#' @export



fit_bag_plsr<-function (soilv, spec, nbag, maxc) {
  # fit a bootstrap aggregated PLSR
  # soilv = soil variables of interest, spec= spectra, nbag, no. of bootstrap, maxc= max. no. of components used in PLSR
  # return variables: model = PLS models, oob_rmse = mean out-of-bag RMSE, cal_rmse = mean calibration RMSE
  # use pls package
  # Budiman March 2016
    nc <- maxc
    n <- length(soilv)
    v.pls <- vector(nbag, mode = "list")
    cal_rmse <- matrix(0, nrow = nbag, ncol = 1)
    oob_rmse <- matrix(0, nrow = nbag, ncol = 1)
    for (ibag in 1:nbag) {
      s <- sample.int(n, replace = TRUE)
      v.pls[[ibag]] <- plsr(soilv[s] ~ ., data= spec[s, ], maxc)
      pred.v <- predict(v.pls[[ibag]], ncomp = nc, newdata = spec[s, 
                                                                  ])
      err2 <- (soilv[s] - pred.v)^2
      cal_rmse[ibag] <- sqrt(mean(err2))
      pred.C <- predict(v.pls[[ibag]], ncomp = nc, newdata = spec[-s, 
                                                                  ])
      err2 <- (soilv[-s] - pred.C)^2
      oob_rmse[ibag] <- sqrt(mean(err2))
    }
    av_cal_rmse <- mean(cal_rmse)
    av_oob_rmse <- mean(oob_rmse)
    list(model.bpls = v.pls, oob_rmse = av_oob_rmse, cal_rmse = av_cal_rmse)
  }
  
