#' Converts an multiply imputed dataset obtained from DURR and IURR functions into a \code{mids} object
#'
#' This function converts imputed data obtained from \code{GMdurr} or \code{GMiurr} into an object of class \code{mids}. The original incomplete data set needs to be available so that we know where the missing data are. This function is useful to manipulate the imputed data sets with the help of functions from \code{mice} package.
#' @param data The original incomplete data set
#' @param imp The return value from \code{GMdurr} or \code{GMiurr} with the option of \code{included=FALSE}
#' @return An object of class \code{\link{mids}}
#' @export
#' @examples
#' data(GMdata)
#' GM.imp=GMdurr(data=GMdata)
#' imp=as.mice(data=GMdata,imp=GM.imp)


as.mice=function(data,imp){
  m=dim(imp$imp[[1]])[2]
  mice_temple=mice(data,m=m,maxit=0)
  mis_names=names(mice_temple$nmis)[mice_temple$nmis>0]
  gm_mice_temple=mice_temple
  for(i in 1:sum(mice_temple$nmis>0)){
    gm_mice_temple$imp[[mis_names[i]]]=imp$imp[[i]]
  }
  return(gm_mice_temple)
}



#' High-dimensional data set with general missing patterns.
#'
#' A data set containing 100 observations of 200 continuous variables (\code{z_1} to \code{z_200}). The first three variables have missing values with numbers of 28, 30 and 28, respectively.
#'
#' @format A data frame with 100 rows and 200 variables
#' @examples
#' data(GMdata)
#' @name GMdata
NULL