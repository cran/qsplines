#' @title Barry-Goldman quaternions spline
#' @description Constructs a spline of unit quaternions by the Barry-Goldman
#'   method.
#'
#' @param keyRotors a vector of unit quaternions (rotors) to be interpolated;
#'   it is automatically appended with the first one to have a closed spline
#' @param keyTimes the times corresponding to the key rotors; must be an
#'   increasing vector of length \code{length(keyRotors)+1}; if \code{NULL},
#'   it is set to \code{c(1, 2, ..., length(keyRotors)+1)}
#' @param n_intertimes a positive integer used to linearly interpolate the 
#'   times given in \code{keyTimes} in order that there are
#'   \code{n_intertimes - 1} between two key times (so one gets the key
#'   times if \code{n_intertimes = 1}); if this argument is given, then 
#'   it has precedence over \code{times}
#' @param times the interpolating times, they must lie within the range of
#'   \code{keyTimes}; ignored if \code{n_intertimes} is given
#'
#' @return A vector of unit quaternions with the same length as \code{times}.
#'
#' @export
#'
#' @note The function does not check whether the quaternions given in
#'   \code{keyRotors} are unit quaternions.
#' 
#' @importFrom utils head
#' 
#' @examples 
#' library(qsplines)
#' # Using a Barry-Goldman quaternions spline to construct 
#' #   a spherical curve interpolating some key points on
#' #   the sphere of radius 5.
#'
#' # helper function: spherical to Cartesian coordinates
#' sph2cart <- function(rho, theta, phi){
#'   return(c(
#'     rho * cos(theta) * sin(phi),
#'     rho * sin(theta) * sin(phi),
#'     rho * cos(phi)
#'   ))
#' }
#'
#' # construction of the key points on the sphere
#' keyPoints <- matrix(nrow = 0L, ncol = 3L)
#' theta_ <- seq(0, 2*pi, length.out = 9L)[-1L]
#' phi <- 1
#' for(theta in theta_){
#'   keyPoints <- rbind(keyPoints, sph2cart(5, theta, phi))
#'   phi = pi - phi
#' }
#' n_keyPoints <- nrow(keyPoints)
#'
#' # construction of the key rotors; the first key rotor is the 
#' #   identity quaternion and rotor i sends the first key point 
#' #   to the key point i
#' keyRotors <- quaternion(length.out = n_keyPoints)
#' rotor <- keyRotors[1L] <- H1
#' for(i in seq_len(n_keyPoints - 1L)){
#'   keyRotors[i+1L] <- rotor <-
#'     quaternionFromTo(
#'       keyPoints[i, ]/5, keyPoints[i+1L, ]/5
#'     ) * rotor
#' }
#'
#' # Barry-Goldman quaternions spline
#' \donttest{rotors <- BarryGoldman(keyRotors, n_intertimes = 10L)
#'
#' # construction of the interpolating points on the sphere
#' points <- matrix(nrow = 0L, ncol = 3L)
#' keyPoint1 <- rbind(keyPoints[1L, ])
#' for(i in seq_along(rotors)){
#'   points <- rbind(points, rotate(keyPoint1, rotors[i]))
#' }
#'
#' # visualize the result with the 'rgl' package
#' library(rgl)
#' spheres3d(0, 0, 0, radius = 5, color = "lightgreen")
#' spheres3d(points, radius = 0.2, color = "midnightblue")
#' spheres3d(keyPoints, radius = 0.25, color = "red")}
BarryGoldman <- function(keyRotors, keyTimes = NULL, n_intertimes, times){
  stopifnot(is.quaternion(keyRotors))
  stopifnot(is.null(keyTimes) || .isNumericVector(keyTimes))
  stopifnot(missing(n_intertimes) || .isPositiveInteger(n_intertimes))
  stopifnot(missing(times) || .isNumericVector(times))
  keyRotors <- .check_keyRotors(keyRotors, closed = TRUE)
  n_keyRotors <- length(keyRotors)
  if(is.null(keyTimes) && missing(times)){
    keyTimes <- seq_len(n_keyRotors)
  }else if(length(keyTimes) != n_keyRotors){
    stop("Number of key times must be one more than number of key rotors.")
  }
  if(missing(times) && missing(n_intertimes)){
    stop(
      "You must supply either `n_intertimes` or `times`."
    )
  }
  if(!missing(n_intertimes)){
    times <- interpolateTimes(keyTimes, n_intertimes, FALSE)
  }
  keyRotors <- .getQMatrix(keyRotors)
  Q <- BarryGoldman_cpp(keyRotors, keyTimes, times)
  as.quaternion(Q)  
}
