#' @title Initialise a finite element basis
#' 
#' @description This function initialises an object of class \code{FEBasis} which defines a set of 'radial `tent' basis functions over a pre-specified triangulation in 2-D
#'
#' @param p \code{n} \eqn{\times} 2 matrix of vertex locations.
#' @param t \code{m} \eqn{\times} 3 matrix of triangulations. Each row identifies which rows of \code{p} make up each triangle.
#' @param M \code{n} \eqn{\times} \code{n} mass matrix: \eqn{\langle \phi, \phi^T \rangle}.
#' @param K \code{n} \eqn{\times} \code{n} stiffness matrix: \eqn{\langle \nabla\phi, \nabla\phi^T \rangle}.
#' @return Object of class \code{FEBasis} 
#' @keywords finite elements, basis functions
#' @export
#' @examples
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
initFEbasis3d = function(p,t,M,K) {
  fn <- pars <- list()
  pars$p <- p
  pars$t <- t
  pars$M <- M
  pars$K <- K
  df <- data.frame(x = pars$p[,1],
                   y = pars$p[,2],
                   z = pars$p[,3]
                   n = 1:nrow(p))
  pars$vars <- df
  # Do tessellation
  Voronoi <- deldir(pars$p[,1],
                    pars$p[,2],
                    plotit='F',
                    sort=F,
                    rw=c(min(pars$p[,1])-0.00001,
                         max(pars$p[,1])+0.00001,
                         min(pars$p[,2])-0.00001,
                         max(pars$p[,2])+.00001))
  pars$pol <- PolygonfromVoronoi(Voronoi,pars$p)
  
  pars$vars$area_tess_km2 = rep(0,nrow(p))
  for (i in 1:nrow(p)) {
    pars$vars$area_tess_km2[i] <- area.poly(pars$pol[[i]])
  }
  this_basis <- new("FEBasis", pars=pars, n=nrow(p), fn=fn)
  return(this_basis)}




setMethod(".find_inc_matrix",signature(basis = "FEBasis"), function(basis,obs,mulfun = NULL, mask = NULL, n_grid=NULL, muldata = NULL, md5_wrapper=NULL) { 
  
  P <- Imat(nrow(obs))
  if (is(obs,"Obs"))
    if("P" %in% names(obs@args)) {
      P <- obs@args$P
    } 
  
  
  
  if (class(obs) == "data.frame") {
    C <- FindC(basis@pars$p,
               basis@pars$t,
               list(obs$x,obs$y),method="C")
  } else if(class(obs) ==  "Obs") {
    C <-  FindC(basis@pars$p,
                basis@pars$t,
                list(obs@df$x,obs@df$y),method="C")
    
  } else if(class(obs) ==  "Obs_poly") {
    if(is.null(n_grid)) {
      warning("n_grid not specified, defaulting to a grid of 400 points for integration over footprints")
      n_grid <- 400
    }
    
    
    if (is.null(mulfun))  mulfun <- 1
    
    if(is.null(md5_wrapper)) {
      C <- FindC_polyaverage(
        basis@pars$p,
        basis@pars$t,
        obs@pol,
        plotit=F,
        method="C",
        ds=n_grid,
        mulfun=mulfun,
        muldata=muldata)
    } else {
      stopifnot(class(md5_wrapper) == "function")
      C <- md5_wrapper(FindC_polyaverage,
                       basis@pars$p,
                       basis@pars$t,
                       obs@pol,
                       plotit=F,
                       method="C",
                       ds=n_grid,
                       mulfun=mulfun,
                       muldata=muldata)
    }
    if (!(is.null(mask))) {
      if(!(mask %in% names(getDf(basis)))) stop("Cannot find mask field in basis")
      C[,which(!basis[mask])] <- 0 
    }
    
    
  }
  
  C <- P %*% C
  
  return(C)
  
})


## Find C matrix when observations are isolated points
FindC <- function(p,tri,locs,method="R") {
  # p 
  if (length(locs[[1]]) > 0)  {
    t_num <- tsearch2(p[,1], p[,2], tri, locs[[1]], locs[[2]], bary = FALSE)
    z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
    b <- matrix(0,3,1)
    A <- matrix(0,3,3)
    A[,1] = 1
    
    
    if(method=="C") {
      nt <- length(t_num)
      X <- .C("element_interp",as.integer(as.integer(nt)),as.integer(t_num),as.integer(tri[,1]),
              as.integer(tri[,2]),as.integer(tri[,3]),as.double(p[,1]),as.double(p[,2]),
              as.double(locs[[1]]),as.double(locs[[2]]),i_ind1 = double(nt), i_ind2 = double(nt),
              i_ind3 = double(nt), j_ind1 = double(nt), j_ind2 = double(nt),
              j_ind3 = double(nt),z1 = double(nt), z2 = double(nt),z3 = double(nt))
      i_ind <- cbind(X$i_ind1,X$i_ind2,X$i_ind3)
      j_ind <- cbind(X$j_ind1,X$j_ind2,X$j_ind3)       
      z <- cbind(X$z1,X$z2,X$z3)       
      
    } else {
      
      for (i in 1:length(t_num)) {
        t_num_i <- t_num[i]
        this_tri <- tri[t_num_i,]
        this_p <- p[this_tri,]
        A[,2:3] <- this_p
        Ainv <- solve(A) 
        
        #        A <- matrix(c(1,1,1,p[tri[t_num_i,1],1],p[tri[t_num_i,2],1],
        #                 p[tri[t_num_i,3],1],p[tri[t_num_i,1],2],
        #                 p[tri[t_num_i,2],2],p[tri[t_num_i,3],2]),3,3)
        
        for (j in 1:3) {
          b[,] <- 0
          b[j] <- 1
          i_ind[i,j] <- i
          j_ind[i,j] <- this_tri[j]
          z[i,j] <- matrix(c(1,locs[[1]][i],locs[[2]][i]),1,3)%*%Ainv%*%b
        }
      }
    }
    
    C <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                      dims = c(length(locs[[1]]),dim(p)[1]))
    return(C)
  }  else {
    return(matrix(1,0,0))
  }
  
}