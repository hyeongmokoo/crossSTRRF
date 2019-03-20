#' Log ratio of spatial densities
#' @references 
#' 
#' @author Hyeongmo Koo
#' @param x_y_t are x-y coordinates and time of the points in vector
#' @param spatial.bw is a spatial bandwith in a two column vector, first for x and second for y
#' @param temporal.bw is a temporal bandwidth in a scalar
#' @param t.size is a tempolar voxel size.
#' @param study.area is a boundary of a study area in the SpatialPolygonDataFrame class
#' @param nsim is the number of simulation. nsim is larger than 1, the significance test of the function is conducted.

cross.strrf <- function(x, y, t, spatial.bw, temporal.bw, t.size, study.area, nsim=0){
  require(spatstat) ##Required package
  require(smacpod)##Required package

  
  max.time <- floor(max(t)/t.size)
  
  sta.time <- 1
  t.grid <- max.time - sta.time + 1  
  
  res.strrf.lst <- vector("list", t.grid)
  res.strrf.sim.lst <- vector("list", t.grid)
  res.strrf.test.lst <- vector("list", t.grid)
  
  sample.dbf <- data.frame(x, y, t)
  sa.owin <- as.owin(study.area)
  
  for(i in sta.time:max.time){
    t.cri <- t.size * i
    sub.sample <- subset(sample.dbf, (t > t.cri - temporal.bw)&(t < t.cri +temporal.bw))
    
    time.vec <- sub.sample$t
    time.vec[which(sub.sample$t>=t.cri)]<-2
    time.vec[which(sub.sample$t<t.cri)]<-1
    m <- factor(time.vec, levels=1:2, labels = c("controls", "cases"))
    
    sample.ppp <- ppp(sub.sample$x, sub.sample$y, window = sa.owin, marks=m)
    
    k <- (sub.sample$t-(t.size*i))/temporal.bw
    wei.t <- (exp(-(k^2)/2))/(sqrt(2*pi))
    case <- subset(sample.ppp, marks=="cases")
    control <- subset(sample.ppp, marks=="controls")
    
    res.strrf <- logrrf.mod(sample.ppp, case=2, sigma=spatial.bw, edge=T, diggle=T, sigmacon=spatial.bw, weights = wei.t)
    
    if(nsim > 1){
      res.strrf.sim <- logrrf.mod(sample.ppp, case=2, level=0.95, sigma=spatial.bw, edge=T, diggle=T, sigmacon=spatial.bw, weights = wei.t, nsim=999)
      res.strrf.test <- logrr.test(res.strrf.sim)
      res.strrf.sim.lst[[i]] <- res.strrf.sim
      res.strrf.test.lst[[i]] <- res.strrf.test
      
    }
    res.strrf.lst[[i]] <- res.strrf
  }
  
  output <- list(res.strrf.lst, res.strrf.sim.lst, res.strrf.test.lst)
}



#' Log ratio of spatial densities
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.
#' 
#' @author Joshua French, the functions are adopted and modified from the "smacpod" package.

logrrf.mod <- function(x, sigma = NULL, sigmacon = NULL, case = 2, nsim = 0, level = 0.90, alternative = "two.sided", ..., bwargs = list(), weights=NULL, edge=TRUE, varcov=NULL, at="pixels", leaveoneout=TRUE, adjust=1, diggle=FALSE, nreport = 50)
{
  if(!is.element(alternative, c("two.sided", "greater", "lower"))) stop("alternative is not valid.")
  alpha = 1 - level

  # determine bandwidth if necessary
  if(is.function(sigma)) # use user-supplied function, if given
  {
    which_bwargs <- which(names(bwargs) %in% names(formals(sigma))[-1])
    if(length(which_bwargs) > 0 )
    {
      sigma = do.call(sigma, c(list(X = x, bwargs[which_bwargs])))
    }else
    {
      sigma = do.call(sigma, list(X = x))
    }
  }
  if(is.null(sigma)) # use spatstat::bw.relrisk if nothing given
  {
    which_bwargs <- names(bwargs) %in% names(formals(spatstat::bw.relrisk))[-1]
    if(length(which_bwargs) > 0 )
    {
      sigma = do.call(spatstat::bw.relrisk, c(list(X = x, bwargs[which_bwargs])))
    }else
    {
      sigma = do.call(spatstat::bw.relrisk, list(X = x))
    }
  }
  if(is.null(sigmacon)) sigmacon = sigma # determine sigmacon, if NULL
  
  cases = which(x$marks == levels(x$marks)[case])
  N1 = length(cases)
  r = spdensity(x = x[cases,], sigma = sigma, ..., weights = weights[cases],
                  edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                  adjust = adjust, diggle = diggle)
  g = spdensity(x = x[-cases,], sigma = sigmacon, ..., weights = weights[-cases],
                  edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                  adjust = adjust, diggle = diggle)
  r$v = log(r$v) - log(g$v)
  r$tolenv = NULL
  
  if(nsim > 0)
  {
    simr = array(0, dim = c(r$dim, nsim + 1))
    simr[,,1] = r$v
    if(nreport <= nsim) cat("Simulations completed: ")
    
    for(i in 1:nsim)
    {
      cases = sample(x$n, N1)
      fsim = spdensity(x = x[cases,], sigma = sigma, ..., weights = weights[cases],
                      edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                      adjust = adjust, diggle = diggle)
      gsim = spdensity(x = x[-cases,], sigma = sigmacon, ..., weights = weights[-cases],
                      edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                      adjust = adjust, diggle = diggle)

            simr[,,i + 1] = log(fsim$v) - log(gsim$v)
      if((i %% nreport) == 0){ cat(paste(i,"")); utils::flush.console() }
    }
    r$simr = simr
    r$tolenv = tolenv(r, level = level, alternative = alternative)
    class(r) = c("logrrenv", class(r))
  }
  
  r$window = x$window
 
  return(r)
}

tolenv = function(object, level = 0.90, alternative = "two.sided")
{
  alpha = 1 - level
  if(alternative == "two.sided")
  {
    tol = apply(object$simr, c(1, 2), quantile, 
                prob = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  }
  else if(alternative == "lower")
  {
    tol = apply(object$simr, c(1, 2), quantile, 
                prob = c(1 - level, 1), na.rm = TRUE)
  }
  else
  {
    tol = apply(object$simr, c(1, 2), quantile, 
                prob = c(0, level), na.rm = TRUE)
  }
  above = (object$simr[,,1] > tol[2,,]) + 0
  below = -1*(object$simr[,,1] < tol[1,,])
  both = above + below
  return(spatstat::im(mat = both, xcol = object$xcol, yrow = object$yrow))
}

