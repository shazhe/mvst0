#' @rdname pred_variance_large
#' @aliases pred_variance_large,list-Graph_2nodes-method
setMethod("pred_variance_large",signature(Results="list",G="Graph_2nodes"),function(Results,G) {
  
  Olist <- extractClass(G@v,"Obs")
  Glist <- extractClass(G@v,"process")
  
  C_full <- G@e[[1]]@Cmat
  Qobs <- getPrecision(Olist[[1]])
  
  x_mean <- Results$Post_GMRF@rep$x_mean
  Qtot <- Results$Post_GMRF@Q

  if (!("Partial_Cov" %in% names(Results))) {
      X <- cholPermute(Qtot)
      Partial_Cov <- Takahashi_Davis(Qtot,cholQp = X$Qpermchol,P = X$P)
  } else {
    Partial_Cov <- Results$Partial_Cov
  }
  Symbolic <- (t(C_full) %*% C_full)
  Symbolic@x = rep(1,length(Symbolic@x))
  Partial_Cov2 <- Partial_Cov* Symbolic # We only need these elements
  
  Olist[[1]]["pred_var"] <- rowSums((C_full %*% Partial_Cov2) * C_full)
  Olist[[1]]["pred_mean"] <- as.vector(C_full %*% x_mean)
  return(Olist[[1]])  
})

#' @rdname pred_variance_large
#' @aliases pred_variance_large,matrix-Graph_2nodes-method
setMethod("pred_variance_large",signature(Results="matrix",G="Graph_2nodes"),function(Results,G) {
  
  Olist <- extractClass(G@v,"Obs")
  Glist <- extractClass(G@v,"process")
  
  C_full <- G@e[[1]]@Cmat
  Qobs <- getPrecision(Olist[[1]])
  
  x_mean <- apply(x_samp,1,mean)
  Olist[[1]]["pred_mean"] <- as.vector(C_full %*% x_mean)
  Olist[[1]]["pred_var"] <- apply((C_full %*% Results),1,var)
  
  return(Olist[[1]])  
})


#' @rdname validate
#' @aliases validate,list-Graph_2nodes-method
setMethod("validate",signature(Results="list",G="Graph_2nodes"),function(Results,G,sim_obs=F) {
  
  Olist <- extractClass(G@v,"Obs")
  Glist <- extractClass(G@v,"process")
  
  C_full_test <- G@e[[1]]@Cmat
  Qobs_test <- getPrecision(Olist[[1]])
  
  x_mean <- Results$Post_GMRF@rep$x_mean
  Qtot <- Results$Post_GMRF@Q
  prec_delta <- last(Results$VB_results$prec_delta)
  eta <- last(Results$VB_results$eta)
  var_eta_post <- last(Results$VB_results$eta_var)
  b_test <- Olist[[1]]@df$fine_scale

  if("Qpermchol" %in% names(Results)) {
      X <- list(Qpermchol = Results$Qpermchol, P = Results$P)
  } else {
     cat("Doing Cholesky decomposition",sep="\n")
     X <- cholPermute(Qtot,matlab_server=matlab_server)
  }
  cat("Finding C %*% Sigma. This will take a while...",sep="\n")
  mm <- t(cholsolve(Qtot,t(C_full_test),perm=T,cholQp = X$Qpermchol, P = X$P))
  cat("Computing C %*% Sigma complete ...",sep="\n")
  var1 <- rowSums(mm*C_full_test)
  Sigma1 <- mm %*% t(C_full_test)   # Sigma1 is the variance due to the process
  
  # check
  #S <- chol2inv(chol(Qtot))
  #Sigma1 <- C_full_test %*% S %*% t(C_full_test)
  # CHECK PASSED
  
  if (!is.null(prec_delta)) {
    var2 <-  1/as.vector(prec_delta * exp(-b_test*eta + b_test^2*var_eta_post/2))
    Sigma2 <- diag(var2)  # Sigma2 is the variance due to the small-scale variation
  } else {
    var2 = 0
    Sigma2 <- Zeromat(nrow(Sigma1))
  }
  var3 <- (1/diag(Qobs_test))  # Sigma2 is the variance due to the observations
  Sigma3 <- sparsediag(var3)
  
  y_test <- getDf(Olist[[1]])
  if(sim_obs) {
   y_test$z <- as.numeric(C_full_test%*%x_mean + t(chol(forceSymmetric(Sigma1 + Sigma2 + Sigma3))) %*% rnorm(nrow(Qobs_test)) )
  }
  
  y_pred_mean <- C_full_test %*% matrix(x_mean)
  y_pred_var <- var1 + var2 + var3
  
  val_results <- y_test
  val_results$y_pred_mean <- as.numeric(y_pred_mean)
  val_results$var1 <- var1
  val_results$var2 <- var2
  val_results$var3 <- var3
  val_results$RMS <- as.vector((y_test$z - y_pred_mean)^2)*diag(Qobs_test)
  val_results$PMCC <- as.vector((y_test$z - y_pred_mean)^2 + (y_pred_var))
  val_results$CR1 <- as.vector((y_test$z - y_pred_mean)/sqrt(y_pred_var))   # CR1 is the standardised error
  val_results$CR2 <- sqrt(as.vector((y_test$z - y_pred_mean)^2/(y_pred_var)))
  
  
  # Compute Mahalnobis distance (See notes and Bastos and OHagan)
  val_results$Mahalanobis <- as.numeric(t(y_test$z - y_pred_mean) %*% chol2inv(chol(forceSymmetric(Sigma1 + Sigma2 + Sigma3))) %*% (y_test$z - y_pred_mean))
  
  # Eigen errors
  X <- eigen(Sigma1 + Sigma2 + Sigma3)
  G <- X$vectors %*% diag(sqrt(X$values))
  val_results$DE <- as.numeric(solve(G,(y_test$z - y_pred_mean)))
  
  # Pivoted Chol errors
#   X <- find_pivot_Bastos(Sigma1+Sigma2+Sigma3)
#   R <- X$R
#   pivot <- X$piv
#   P <- sparseMatrix(i=pivot,j=1:length(pivot),x=1)
  R <- chol(as(Sigma1+Sigma2+Sigma3,"matrix"),pivot=T)
  pivot <- attr(R,"pivot")
  P <- sparseMatrix(i=pivot,j=1:length(pivot),x=1)
  G <- P %*% t(R)
  val_results$DPC <- as.numeric(solve(G,(y_test$z - y_pred_mean)))
  return(val_results)
})

#' @rdname validate
#' @aliases validate,matrix-Graph_2nodes-method
setMethod("validate",signature(Results="matrix",G="Graph_2nodes"),function(Results,G,sim_obs=F) {
  
  Olist <- extractClass(G@v,"Obs")
  Glist <- extractClass(G@v,"process")
  
  C_full_test <- G@e[[1]]@Cmat
  Qobs_test <- getPrecision(Olist[[1]])
  
  x_mean <- apply(x_samp,1,mean)
  b_test <- Olist[[1]]@df$fine_scale
  
  cat("Finding C %*% Sigma. This will take a while...",sep="\n")
  Results_detrended <- Results - matrix(x_mean,nrow(Results),ncol(Results))
  Sigma1 <- (C_full_test %*% Results_detrended) %*% t(C_full_test %*% Results_detrended)/(ncol(Results)-1)
  var1 <-  diag(Sigma1)
  
  var2 = 0
  Sigma2 <- Zeromat(nrow(Sigma1))
  
  var3 <- (1/diag(Qobs_test))  # Sigma2 is the variance due to the observations
  Sigma3 <- sparsediag(var3)
  
  y_test <- getDf(Olist[[1]])
  if(sim_obs) {
    y_test$z <- as.numeric(C_full_test%*%x_mean + t(chol(forceSymmetric(Sigma1 + Sigma2 + Sigma3))) %*% rnorm(nrow(Qobs_test)) )
  }
  
  y_pred_mean <- C_full_test %*% matrix(x_mean)
  y_pred_var <- var1 + var2 + var3
  
  val_results <- y_test
  val_results$y_pred_mean <- as.numeric(y_pred_mean)
  val_results$var1 <- var1
  val_results$var2 <- var2
  val_results$var3 <- var3
  val_results$RMS <- as.vector((y_test$z - y_pred_mean)^2)*diag(Qobs_test)
  val_results$PMCC <- as.vector((y_test$z - y_pred_mean)^2 + (y_pred_var))
  val_results$CR1 <- as.vector((y_test$z - y_pred_mean)/sqrt(y_pred_var))   # CR1 is the standardised error
  val_results$CR2 <- sqrt(as.vector((y_test$z - y_pred_mean)^2/(y_pred_var)))
  
  
  # Compute Mahalnobis distance (See notes and Bastos and OHagan)
  val_results$Mahalanobis <- as.numeric(t(y_test$z - y_pred_mean) %*% chol2inv(chol(forceSymmetric(Sigma1 + Sigma2 + Sigma3))) %*% (y_test$z - y_pred_mean))
  
  # Eigen errors
  X <- eigen(Sigma1 + Sigma2 + Sigma3)
  G <- X$vectors %*% diag(sqrt(X$values))
  val_results$DE <- as.numeric(solve(G,(y_test$z - y_pred_mean)))
  
  # Pivoted Chol errors
  #   X <- find_pivot_Bastos(Sigma1+Sigma2+Sigma3)
  #   R <- X$R
  #   pivot <- X$piv
  #   P <- sparseMatrix(i=pivot,j=1:length(pivot),x=1)
  R <- chol(as(Sigma1+Sigma2+Sigma3,"matrix"),pivot=T)
  
  pivot <- attr(R,"pivot")
  P <- sparseMatrix(i=pivot,j=1:length(pivot),x=1)
  G <- P %*% t(R)
  val_results$DPC <- as.numeric(solve(G,(y_test$z - y_pred_mean)))
  return(val_results)
})



#' @rdname extractClass
#' @aliases extractClass,list-character-method
setMethod("extractClass",signature(L="list",Cl = "character"),function(L,Cl) {
  types <- lapply(L,function(df) {return(is(df,Cl))})
  Outputlist <-L[which(types==1)]
  return(Outputlist)
})
setMethod(".exist",signature(L="link_list",to="block",from="block"),function(L,to,from) {
  detect <- lapply(L,function(df) {
    if((df@from@uid == from@uid) & (df@to@uid == to@uid)) {
      return(1)
    } else {
      return(0)
    }
    
  })
  if(any(unlist(detect) == 1)) {
    return(which(detect==1))
  } else {
    return(0)
  }
})


#' @rdname compress
#' @aliases compress,Graph-method
setMethod("compress",signature(Graph="Graph"),function(Graph) {
  # Create ordering
  Olist <- extractClass(Graph@v,"Obs")
  Glist <- extractClass(Graph@v,"process")
  
  for(i in 1:length(Glist)) {
   if(class(Glist[[i]])=="GMRF_basis")
    Glist[[i]] <- Glist[[i]]@G 
  }
  
  # Produce big C Matrix
  Cmat <- Reduce("rBind",lapply(Olist,function(Out) { 
    Reduce("cBind",lapply(Glist,function(In) {
      if(l <- .exist(Graph@e,Out,In)) {
        return(Graph@e[[l]]@Cmat)
      } else {
        return(Zeromat(nrow(Out),nrow(In@Q))) 
      }
    }))
  }))
  
  O_reduced <- Reduce("concat",Olist)
  GMRF_reduced <- Reduce("concat",Glist)
  L <- link(GMRF_reduced, O_reduced, Cmat = Cmat)
  e_reduced <- new("link_list",list(L=L))
  v_reduced <- new("block_list",list(G=GMRF_reduced,O=O_reduced))
  Graph_reduced <- new("Graph_2nodes",e=e_reduced,v=v_reduced)
  return(Graph_reduced)
})


#' @rdname setGMRF
#' @aliases setGMRF,Graph_2nodes-method
setMethod("setGMRF",signature(Graph="Graph_2nodes", obj="GMRF"),
          function(Graph,obj) {
            for(i in 1:2) {
             if(is(Graph@v[[i]],"GMRF"))
               Graph@v[[i]] <- obj
            }
            return(Graph)
          })


#' @rdname Infer
#' @aliases Infer,Graph_2nodes-method
setMethod("Infer",signature(Graph="Graph_2nodes"),
          function(Graph,SW=0,Comb=NULL) {
            Olist <- extractClass(Graph@v,"Obs")
            Glist <- extractClass(Graph@v,"process")
            
            Obs <- Olist[[1]]
            Field <- Glist[[1]]
            t_axis <- Glist[[1]]@t_axis
            y_tot <- getDf(Obs)
            Qobs <- getPrecision(Obs)
            Q_full <- getPrecision(Field)
            x_prior <- getMean(Field)
            C_full <- Graph@e[[1]]@Cmat
            if(!SW) {
              ybar = t(C_full)%*%Qobs%*%y_tot$z + Q_full %*% x_prior
              Qtot <- t(C_full)%*%Qobs%*%C_full + Q_full 
              cat("Doing Cholesky and Takahashi",sep="\n")
              X <- cholPermute(Qtot,matlab_server=matlab_server)
              
              x_mean <- cholsolve(Qtot,ybar,perm=T,cholQp = X$Qpermchol, P = X$P)
              if(is.null(Comb)) {
                # Standard Gaussian update
                Partial_Cov <- Takahashi_Davis(Qtot,cholQp = X$Qpermchol,P = X$P)
                x_margvar <- diag(Partial_Cov)
              } else {
                Cov_lin <- cholsolveAQinvAT(Qtot,Comb,Lp = X$Qpermchol,P=X$P)
                x_mean_comb <- as.vector(Comb%*%x_mean)
                x_margvar_comb <- as.vector(diag(Cov_lin)) 
              }
            } else {
              if(is.null(Comb)) {
                stop("Will not attempt Sherman-Woodbury without linear combs.")
              } else {
                X <- cholPermute(Q_full,matlab_server=matlab_server)
                U1 <- cholsolve(Q_full,t(Comb),perm=T,cholQp = X$Qpermchol, P = X$P)
                Tmat <- chol2inv(chol(diag((1/(diag(Qobs)))) + cholsolveAQinvAT(Q_full,C_full,Lp=X$Qpermchol, P = X$P)))
                U2 <-   -cholsolve(Q_full,t(C_full) %*% (Tmat %*% (C_full %*% U1)),perm=T,cholQp = X$Qpermchol, P = X$P)
                Cov_lin <- Comb %*% (U1 + U2)
                mu_lin <- t(U1 + U2) %*% (t(C_full) %*% Qobs %*% y_tot$z + Q_full %*% x_prior)
                x_mean_comb <- as.vector(mu_lin)
                x_margvar_comb <- as.vector(diag(Cov_lin))
              }
            }
            
            if(is.null(Comb)) {
              rep <- cbind(Field@rep,data.frame(x_mean = as.numeric(x_mean), x_margvar=as.numeric(x_margvar)))
              Results_GMRF <- new("GMRF",mu=matrix(x_mean),Q = Qtot,rep=rep)
              for(i in 1:length(Graph@v)) {
                if(is(Graph@v[[i]],"Obs")) {
                  Graph@v[[i]]@df$residuals <- Graph@v[[i]]@df$z - as.numeric(Graph@e$L@Cmat %*% x_mean)
                }
              }
              Results_list <- list(Graph=Graph, Post_GMRF = Results_GMRF,Partial_Cov = Partial_Cov,
                                   Qpermchol = X$Qpermchol, P = X$P)
            } else {
              Comb_results <- list(mu=matrix(x_mean_comb),cov = Cov_lin)
              if(!SW) {
                Results_GMRF <- new("GMRF",mu=matrix(x_mean),Q = Qtot)
                Results_list <- list(Graph=Graph, Post_GMRF = Results_GMRF, Comb_results = Comb_results)
              } else {
                Results_list <- list(Graph=Graph, Comb_results = Comb_results)
              }
            }
            
            return(Results_list)
          })


#' @title link
#' @description This function takes am object of class \code{process} and an object of class \code{Obs} and creates a link between the two.
#' If the \code{process} is of class \code{GMRF_basis}, then an incidence matrix is constructed which maps the process to the observation. If
#' there is no basis set associated with the GMRF, then the incidence matrix needs to be specified manually.
#' @param Obj1 object of class \code{process}.
#' @param Obj2 object of class \code{Obs}.
#' @param Cmat the matrix mapping the process to the observation. Needs to be specified if no basis function set is assoicated with the process.
#' @param mul_factor a constant with which to scale the process, i.e. \eqn{y = mul_factor * C * x}
#' @param mulfun  as above, but a possibly spatially varying amplification of the process of the observation. This, for example, could be a spatially-varying
#' density function when converting height to mass. 
#' @param n_grid the number of grid points to use when integrating over the process with an observation of a large footprint. Defaults to a 400 x 400 grid.
#' @param mask the label of a process attribute. The observation is assumed to only be influenced by the process where this attribute is set to 1.
#' @return Object of class \code{linkGO}, a link between a GMRF and an observation.
#' @keywords link, incidence matrix
#' @export
#' @examples
#' \dontrun{
#' require(Matrix)
#' data(icesat)
#' data(surf_fe)
#'
#' ## First create observation object
#' icesat_obs <- Obs(df=icesat,
#'                  abs_lim = 5,
#'                  avr_method = "median",
#'                  box_size=100,
#'                  name="icesat")
#'
#' ## Now create GMRF defined over some FE basis
#' Mesh <- initFEbasis(p=surf_fe$p,
#'                     t=surf_fe$t,
#'                     M=surf_fe$M,
#'                     K=surf_fe$K)
#' 
#' mu <- matrix(0,nrow(Mesh),1)
#' Q <- sparseMatrix(i=1:nrow(surf_fe$p), j = 1:nrow(surf_fe$p), x = 1)
#'
#' my_GMRF <- GMRF(mu = mu, Q = Q,name="SURF",t_axis = 0:6)
#' SURF <-GMRF_basis(G = my_GMRF, Basis = Mesh)
#'
#' L1 <- link(SURF,icesat_obs)
#' }
link <- function(Obj1,Obj2,Cmat = NULL, mul_factor = NULL, mulfun = NULL, muldata=NULL,
                 n_grid = NULL, mask = NULL, md5_wrapper = NULL) {
  if (is(Obj1,"process") & is(Obj2,"Obs"))
  {
    .Object <- new("linkGO",from=Obj1,to=Obj2,Cmat = Cmat,
                   mul_factor = mul_factor,  mulfun = mulfun,muldata=muldata,
                   n_grid = n_grid, mask = mask,md5_wrapper=md5_wrapper)
  } else stop("Invalid object specification") 
  
  return(.Object)
}

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
                         method="R",
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
                     method="R",
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
  ## p -- the vertex locations
  ## tri -- triangulations, define the vertices of the triangle by identifying them as rows in p.
  ## locs -- list of xy coordinates of the observation points
  if (length(locs[[1]]) > 0)  {
    # Now for each observation point, find the triangle it falls in
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
## Like FindC_boxaverage2 butfor arbitrary polygons. Here ds is the number of points to use for integration
FindC_polyaverage  <- function(p,tri,polygons,plotit=F,mulfun = 0,muldata=NULL,method="R",ds=400)  {
  if (plotit == T) dev.new()
  n <- dim(p)[1]
  m <- length(polygons)
  C_j <- C_i <- C_z <- NULL
  
  ## transform p to be LongLat if it is on spere
  if(ncol(p) == 3){
   p0 <- p
   p <- do.call(cbind, Lxyz2ll(list(x = p[,1], y = p[,2], z = p[,3])))
   }
  
  
  # Find radius to consider
  max_tri_lengths <- apply(tri,1,function(x) max(rdist(p0[x,],p0[x,])))
  max_tri_id <- which(max_tri_lengths == max(max_tri_lengths))[1]
  max_tri_length <- max(rdist(p[tri[max_tri_id,],], p[tri[max_tri_id,],]))
  
  
  
  # For each box
  #cat("Constructing incidence matrix from supports ... ",sep="\n")
  #pb <- txtProgressBar(min = 0, max = m, style = 3)
#ttt1 <- proc.time()
  for (i in 1:m) {
    # Creates fine grid
    pv <- poly_xy(polygons[i])
    pnts_to_consider <- which(apply(rdist(p,pv),1,min) < max_tri_length)
    pnts <- p[pnts_to_consider,]
    tris_to_consider <- which(apply(tri,1,function(x) any(x %in% pnts_to_consider)))
    tris <- tri[tris_to_consider,]
    
    x_grid <- seq(min(pv[,1]),max(pv[,1]),length=ds)
    y_grid <- seq(min(pv[,2]),max(pv[,2]),length=ds)
    dx <- mean(diff(x_grid))
    dy <- mean(diff(y_grid))
    #GRID <- meshgrid(x_grid,y_grid)
    #x <- as.vector(GRID$x)
    #y <- as.vector(GRID$y)
    GRID <- expand.grid(y_grid,x_grid)
    x <- GRID[,2]
    y <- GRID[,1]

    xy_in_poly <- pnt.in.poly(cbind(x,y),pv)$pip
    x <- x[which(xy_in_poly == 1)]
    y <- y[which(xy_in_poly == 1)]
    
    # Now for each grid point find the triangle it falls in
    t_num <- tsearch2(p[,1], p[,2], tris, x, y, bary = FALSE)
    t_num <- t_num[!is.na(t_num)]   # Change this to cater for boundary effects later on!
    z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
    # Now see the magnitude for the three basis it touches
    if (class(mulfun) == 'function') {
      mul_vals <- cbind(mulfun(p[tris[t_num,1],]),
                        mulfun(p[tris[t_num,2],]),
                        mulfun(p[tris[t_num,3],]))
    } else {
      mul_vals = 1
    }
    
    
    if(method=="C") {
      nt <- length(t_num)
      X <- .C("element_interp",as.integer(as.integer(nt)),as.integer(t_num),as.integer(tris[,1]),
              as.integer(tris[,2]),as.integer(tris[,3]),as.double(p[,1]),as.double(p[,2]),
              as.double(x),as.double(y),i_ind1 = double(nt), i_ind2 = double(nt),
              i_ind3 = double(nt), j_ind1 = double(nt), j_ind2 = double(nt),
              j_ind3 = double(nt),z1 = double(nt), z2 = double(nt),z3 = double(nt))
      i_ind <- cbind(X$i_ind1,X$i_ind2,X$i_ind3)
      j_ind <- cbind(X$j_ind1,X$j_ind2,X$j_ind3)       
      z <- cbind(X$z1,X$z2,X$z3) 
    } else {
      for (j in 1:length(t_num))      {
        ## t_num[j] identify the triangle the j th grid point falls in
        ## A: first column = 1, 2-3 contains the x, y coords of the vertex of the triangles
        ## for 3d coords, A should be expend to be 3 by 4, with 2-4 cols for x,y,z coords
        A <- matrix(c(1,1,1,
                      p[tris[t_num[j],1],1],p[tris[t_num[j],2],1], p[tris[t_num[j],3],1],
                      p[tris[t_num[j],1],2],p[tris[t_num[j],2],2], p[tris[t_num[j],3],2]),
                    3,3)
        
        for (k in 1:3) {
          b <- matrix(0,3,1)
          b[k] <- 1
          ## def prev: z <- j_ind <- i_ind <- matrix(0,length(t_num),3)
          i_ind[j,k] <- j
          j_ind[j,k] <- tris[t_num[j],k]
          ## weight for interpolating the value at grid point from the triangle vertex
          z[j,k] <- matrix(c(1,x[j],y[j]),1,3)%*%solve(A)%*%b 
          if (class(mulfun) == 'function') {
            
          }
        }
      }   
    }
    ## simplify? sparse matrix
    z <- z*mul_vals              
    Mmat <- sparseMatrix(as.vector(i_ind),as.vector(j_ind),x=as.vector(z),
                         dims = c(length(x),n))
    Cm <-   colSums(Mmat)*dx*dy
    not_zero <- which(abs(Cm)>0)
    C_j <- c(C_j,not_zero)
    C_i <- c(C_i,rep(i,length(not_zero)))
    C_z <- c(C_z,Cm[not_zero])
    #setTxtProgressBar(pb, i)
    
  }
  C <- sparseMatrix(C_i,C_j,x = C_z,dims = c(m,n))
  #ttt2 <- proc.time()
  #ttt2 - ttt1
  return(C)
}

FindC_EOF <- function(X,Obs,roundobs=1) { # Requires x,y and n (n to maintin order after merge)
  Obs$x <- round(Obs$x/roundobs)*roundobs
  Obs$y <- round(Obs$y/roundobs)*roundobs
  C <- merge(Obs,X,all.x=T)
  C <- arrange(C,desc(-n))
  C <- C[,!(names(C) %in% c("x","y","n"))]
  return(C)
}

## Like FindC_boxaverage2 butfor arbitrary polygons
FindC_average_EOF  <- function(X,polygons,mulfun=1)  {
  ## x: coords and spatial process
  ## polygons: a list of polygons
  ## mulfun: weight used when taking the average
  n <- dim(X)[2] - 2
  m <- length(polygons)
  XY <- X[,1:2]
  S <- X[,-(1:2)]
  SS <- S * mulfun
  C <- matrix(0,m,n)
  
  # For each box
  for (i in 1:m) {
    # Find points which lie in box
    pv <- poly_xy(polygons[i])
    xy_in_poly <- pnt.in.poly(XY,pv)$pip
    C[i,] <- apply((SS[which(xy_in_poly==1),]),2,sum)
  }
  return(as(C,"dgCMatrix"))
}

