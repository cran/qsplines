#' @description Adapter to ensure minimal angles between two successive
#'   quaternions.
#' @noRd
.canonicalized <- function(quaternions){
  n_quaternions <- length(quaternions)
  out <- quaternion(length.out = n_quaternions)
  p <- H1
  for(i in seq_len(n_quaternions)){
    q <- quaternions[i]
    if(dotprod(p, q) < 0){
      q <- -q
    }
    out[i] <- q
    p <- q
  }
  out
}

.check_keyRotors <- function(keyRotors, closed){
  if(length(keyRotors) < 2L){
    stop("At least two keyRotors are required.")
  }
  if(closed){
    keyRotors <- c(keyRotors, keyRotors[1L])
  }
  .canonicalized(keyRotors)
}

.check_keyTimes <- function(keyTimes, n_quaternions){
  if(is.null(keyTimes)){
    return(seq_len(n_quaternions))
  }
  if(any(diff(keyTimes) <= 0)){
    stop("`keyTimes` must be an increasing vector of numbers.")
  }
  keyTimes
}

.check_time <- function(t, keyTimes, special = FALSE){
  n_keyTimes <- length(keyTimes)
  lastKeyTime <- keyTimes[n_keyTimes]
  if(t < keyTimes[1L] || t > lastKeyTime){
    stop("The interpolating times must be within the range of `keyTimes`.")
  }
  if(t < lastKeyTime){
    idx <- findInterval(t, keyTimes, left.open = FALSE, rightmost.closed = TRUE)
  }else{ # t = lastKeyTime
    if(special){
      idx <- n_keyTimes - 2L
    }else{
      idx <- n_keyTimes - 1L
    }
  }
  idx
}

.slerp <- function(q1, q2, t){
  (q2 * onion_inverse(q1))^t * q1
}

.isPositiveInteger <- function(x){
  is.numeric(x) && (length(x) == 1L) && (!is.na(x)) && (floor(x) == x)
}

#' @title Interpolate a vector of times
#' @description Linearly interpolate an increasing vector of times. This is
#'   useful to deal with the quaternions splines.
#'
#' @param times increasing vector of times
#' @param n integer, controls the number of interpolations: there will be
#'   \code{n-1} time values between two consecutive original times
#' @param last Boolean, whether to include or exclude the last element
#'
#' @return A vector, a refinement of the \code{times} vector.
#' @export
#'
#' @examples
#' library(qsplines)
#' interpolateTimes(1:4, n = 3)
#' interpolateTimes(c(1, 2, 4), n = 3)
interpolateTimes <- function(times, n, last = TRUE){
  stopifnot(.isPositiveInteger(n))
  n_times <- length(times)
  newtimes <- numeric(0L)
  for(i in seq_len(n_times-1L)){
    newtimes <- 
      c(newtimes, seq(times[i], times[i+1L], length.out = n + 1L)[-n-1L])
  }
  if(last){
    newtimes <- c(newtimes, times[n_times])
  }
  newtimes
}

.getQMatrix <- function(quaternions){
  stopifnot(is.quaternion(quaternions))
  as.matrix(quaternions)
}

.isNumericVector <- function(x){
  is.numeric(x) && !anyNA(x)
}
