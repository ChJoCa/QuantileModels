#' MVMQ-CAViaR model estimation proposed by White et al. (2015).
#'
#' @param Y A matrix, xts or zoo object of multivariate series.
#' @param p Order of autoregressive quantile.
#' @param q Order of lag Y values.
#' @param tau A vector containing the quantiles of interest. Set default to 0.05.
#' @param band.hs Logical parameter passed to [quantreg::bandwidth.rq], if [TRUE] the Hall-Sheather bandwidth is computed, if [FALSE] it is computed the Bofinger bandwidth.
#' @param jac.method One of "Richardson" (default), "simple", or "complex". This determines the method to compute finite differences jacobian to as part of the standard errors calculations. See [numDeriv::jacobian] for more information.
#' @param sign.level The `alpha` parameter of [quantreg::bandwidth.rq] function to control the level of significance for intended confidence intervals.
#' @param quant.type One of the types available in [quantile]. This is used to initialize the quantile process.
#' @param optim.config A list containing the optimization-related parameters, see [nloptr.print.options()] and `details` for more information.
#' @param init.optim.config A list containing the univariate optimization-related parameters, see [nloptr.print.options()] for more information.
#' @param global.estim Should the optimization be performed in first place with global optimization with [GenSA::GenSA] in order to potentially achieve a better exploration of parameter space?
#' @param global.optim A list containing the global optimization-related parameters, see [GenSA::GenSA]  for more information.
#' @param jac.options A list passed to \code{method.args} in [numDeriv::jacobian].
#' @importFrom quantreg bandwidth.rq
#' @importFrom GenSA GenSA
#' @importFrom nloptr nloptr
#' @importFrom numDeriv jacobian
#' @importFrom zoo index
#' @importFrom xts as.xts
#' @returns A list containing different results from the estimation. Class `MVMQ_CAViaR`
#' @details
#' This implementation follows in essence the same optimization strategy as White et al. (2015), by first obtaining for each series their univariate CAViaR estimation with [CAViaR], and setting the rest of the parameters to 0. As mentioned above, staring from those univariate estimates, this function offers two different optimization strategies: i) feed this starting point to a local optimizer (by default `NLOPT_LN_SBPLX`) following the original work, in such case the user should set `global.estim`=FALSE (the default); ii) feed this starting point to [GenSA::GenSA] in order perform a global optimization with the intention to explore the parameter space, and the result from this global phase is then used as starting point for the same procedure as i), then the user should set `global.estim`=TRUE.
#'
#' Regarding the specification, at the moment, only the symmetric absolute value specification is available, having the following form:
#'  \deqn{\boldsymbol{f_t}(\boldsymbol{\theta}) = \boldsymbol{c} + \sum_{i=1}^p \boldsymbol{A_i} \boldsymbol{f_{t-i}(\boldsymbol{\theta})} + \sum_{j=1}^q \boldsymbol{B_j} |\boldsymbol{Y_{t-j}|}}
#' @references White, H., Kim, T. H., & Manganelli, S. (2015). VAR for VaR: Measuring tail dependence using multivariate regression quantiles. Journal of econometrics, 187(1), 169-188.
#' @export
#' @examples
#' \donttest{
#' Barclays <- MVMQ_CAViaR(MVMQ[,c(6,1)],tau =c(0.01,0.01),band.hs = TRUE)
#' summary(Barclays)
#' #or
#' Barclays
#' plot(Barclays,rows=2,columns=1)
#' }

MVMQ_CAViaR <- function(Y,p=1,q=1,tau=rep(0.05,ncol(Y)),band.hs=FALSE,jac.method="simple",jac.options=list(),sign.level=0.05,quant.type=7,
                   optim.config=list(),init.optim.config=list(),global.estim=FALSE,global.optim=list()){




  fechas <- as.Date(index(Y))
  nombres <- colnames(Y)
  serie <- unname(as.matrix(Y))
  T <- nrow(serie)
  N=ncol(serie)
  VaR=residuos=matrix(0.0,ncol = N,nrow = T)
  emp_quant <- diag(apply(serie[1:100,],MARGIN = 2,FUN = quantile,probs=tau,type=quant.type))

  dimensions=c(p,q)
  maximum=max(dimensions)

  if(jac.method=="Richardson"){
    options <- list(eps=1e-8, d=0.00001, zero.tol=sqrt(.Machine$double.eps/7e-4), r=6, v=4, show.details=FALSE)
    jacobian.method <- modifyList(options,jac.options)
  }else{
    options <- list(eps=1e-9)
    jacobian.method <- modifyList(options,jac.options)
  }
  optimization.config <- list(algorithm="NLOPT_LN_SBPLX",maxeval=50000,xtol_rel=1e-10,xtol_abs=0,ftol_rel=0)

  optimization.config <- modifyList(optimization.config,optim.config)


  global.optimization <- list(maxit=200,temperature=20000,visiting.param=6)
  global.optimization <- modifyList(global.optimization,global.optim)


  bandwithd <-bandwidth.rq(p=tau,n=T,hs=band.hs,alpha = sign.level)
  nsqr <- N*N
  #estimation
  lagged_part <- abs(embed(serie,dimension = 1+q)[,-(1:N)])
  for(j in 1:p){VaR[p,] <- emp_quant}

   init_constants <- vector(mode = "numeric",length = N)
  init_A_mats <- array(0.0,dim = c(N,N,p))
  init_B_mats <- array(0.0,dim = c(N,N,q))
  univariate_estimations <- matrix(0.0,ncol = N,nrow = 1+p+q)

  for (w in 1:N) {
    univariate_estimations[,w] <- CAViaR_init(serie[,w],p=p,q=q,tau = tau[w],quant.type=quant.type,optim.config=init.optim.config)
  }

  for (z in 1:p) {
    init_A_mats[,,z] <- diag(univariate_estimations[1+z,])
  }
  for (g in 1:q) {
    init_B_mats[,,g] <- diag(univariate_estimations[1+p+g,])
  }
  inital_A <- matrix(init_A_mats,nrow = N,ncol = N*p)
  inital_B <- matrix(init_B_mats,nrow = N,ncol = N*q)

  initial_point <- c(univariate_estimations[1,],as.vector(t(inital_A)),as.vector(inital_B))
  message("Begining optimization")
  if(global.estim==TRUE){
  estimacion <- GenSA(par = initial_point,fn=OBJ_F_MVMQ,lower = rep(-5,N+nsqr*p+nsqr*q),upper = rep(5,N+nsqr*p+nsqr*q),control = global.optimization,series=serie,tau=tau,lagged_abs_y=as.matrix(lagged_part),N=N,residuos=residuos,VaR_vector=VaR,p=p,q=q,maximum=maximum,time=T)
  estimacion <- list(solution=estimacion$par,objective=estimacion$value)

  estimacion <- nloptr(x0=estimacion$solution,eval_f =OBJ_F_MVMQ,lb = rep(-5,N+nsqr*p+nsqr*q) ,ub=rep(5,N+nsqr*p+nsqr*q),
                       opts = optimization.config,
                       series=serie,tau=tau,lagged_abs_y=as.matrix(lagged_part),
                       N=N,residuos=residuos,
                       VaR_vector=VaR,p=p,q=q,maximum=maximum,
                       time=T)
  } else
    estimacion <- nloptr(x0=initial_point,eval_f =OBJ_F_MVMQ,lb = rep(-5,N+nsqr*p+nsqr*q) ,ub=rep(5,N+nsqr*p+nsqr*q),
                         opts = optimization.config,
                         series=serie,tau=tau,lagged_abs_y=as.matrix(lagged_part),
                         N=N,residuos=residuos,
                         VaR_vector=VaR,p=p,q=q,maximum=maximum,
                         time=T)




  predict <- MVMQ_FILTER(estimacion$solution,lagged_abs_y=as.matrix(lagged_part),
                         N=N,
                         VaR_vector=VaR,p=p,q=q,maximum=maximum,
                         time=T)

  constants <- as.matrix(estimacion$solution[1:N])
  A_matrix <- matrix(estimacion$solution[(N+1):(N+nsqr*p)],ncol = N*p,nrow = N,byrow = TRUE)
  B_matrix <- matrix(estimacion$solution[((N+nsqr*p)+1):(N+nsqr*p+nsqr*q)],ncol = N*q,nrow = N,byrow =FALSE)


  residuoss <- serie-predict
  hits <- (residuoss<0)*1
  psi <- t(tau-t(hits))
  ka_te <- apply(abs(residuoss-quantile(residuoss,probs=0.5)),MARGIN=2,quantile,probs=0.5) # hacer quantile con APPLY como arriba
  ce_te <- as.vector(ka_te*(qnorm(tau+bandwithd)-qnorm(tau-bandwithd)))

  matrices_Q_and_V <- array(0.0,dim = c(N+nsqr*p+nsqr*q,N+nsqr*p+nsqr*q,2),dimnames =list(NULL,NULL,c("V","Q") ))
    message("Calculating Standard Errors")

    for (i in 1:N) {
      jaco_n <-  jacobian(SAV_FILTER_for_jaco,estimacion$solution,lagged_abs_y=as.matrix(lagged_part),
                          N=N,
                          VaR_vector=VaR,p=p,q=q,maximum=maximum,serie=(i-1),
                          time=T,method = jac.method,method.args = jacobian.method)
      V_jaco_n <- jaco_n*psi[,i]
      matrices_Q_and_V[,,1] <-matrices_Q_and_V[,,1]+crossprod(V_jaco_n)/T

      indicator <- (abs(residuoss[,i])<ce_te[i])*1
      Q_jaco_n  <- (jaco_n*indicator)
      matrices_Q_and_V[,,2] <- matrices_Q_and_V[,,2]+crossprod(Q_jaco_n)/(2*ce_te[i]*T)

    }


    Q_inv <- solve(matrices_Q_and_V[,,2])

    VcV <- (Q_inv%*%matrices_Q_and_V[,,1]%*%Q_inv)/T
    se <- sqrt(diag(VcV))


parameters=estimacion$solution
degrees=T-length(parameters)
p_values=2*(1-pt(abs(parameters/se),df=degrees))
upper_interval=parameters+qt(1-(sign.level/2),df=degrees)*se
lower_interval=parameters+qt(sign.level/2,df=degrees)*se

const_se <- as.matrix(se[1:N])
A_se <- matrix(se[(N+1):(N+nsqr*p)],ncol = N*p,nrow = N,byrow = TRUE)
B_se <- matrix(se[((N+nsqr*p)+1):(N+nsqr*p+nsqr*q)],ncol = N*q,nrow = N,byrow =FALSE)

const_pval <- as.matrix(p_values[1:N])
A_pval <- matrix(p_values[(N+1):(N+nsqr*p)],ncol = N*p,nrow = N,byrow = TRUE)
B_pval <- matrix(p_values[((N+nsqr*p)+1):(N+nsqr*p+nsqr*q)],ncol = N*q,nrow = N,byrow =FALSE)


const_li <- as.matrix(lower_interval[1:N])
A_li <- matrix(lower_interval[(N+1):(N+nsqr*p)],ncol = N*p,nrow = N,byrow = TRUE)
B_li <- matrix(lower_interval[((N+nsqr*p)+1):(N+nsqr*p+nsqr*q)],ncol = N*q,nrow = N,byrow =FALSE)

const_ui <- as.matrix(upper_interval[1:N])
A_ui <- matrix(upper_interval[(N+1):(N+nsqr*p)],ncol = N*p,nrow = N,byrow = TRUE)
B_ui <- matrix(upper_interval[((N+nsqr*p)+1):(N+nsqr*p+nsqr*q)],ncol = N*q,nrow = N,byrow =FALSE)

rownames(A_matrix)=rownames(B_matrix)=rownames(constants)=rownames(const_se)=rownames(A_se)=rownames(B_se)=rownames(const_pval)=rownames(A_pval)=rownames(B_pval)=rownames(const_li)=rownames(A_li)=rownames(B_li)=rownames(const_ui)=rownames(A_ui)=rownames(B_ui)=nombres
quantile_names <- paste0("q.",nombres,",",rep(1:p,each=N))
absolute_names <- paste0("|",nombres,"|",",",rep(1:q,each=N))
colnames(A_matrix)=quantile_names
colnames(B_matrix)=absolute_names
colnames(hits)=colnames(predict)=colnames(serie)=nombres
names(tau) <- nombres
output <- structure(list(estimacion,Cons=constants,A_mat=A_matrix,B_mat=B_matrix,
                         VaR=as.xts(predict,order.by=fechas),SE=list(constants=const_se,A_mat=A_se,B_mat=B_se),
                         CI_low=list(constants=const_li,A_mat=A_li,B_mat=B_li),
                         CI_up=list(constants=const_ui,A_mat=A_ui,B_mat=B_ui),
                         Pval=list(constants=const_pval,A_mat=A_pval,B_mat=B_pval),CI_level=sign.level,
                         tau=tau,data=as.xts(serie,order.by=fechas),loss=estimacion$objective,bandwithd=band.hs,hits_seq=hits),class="MVMQ_CAViaR")

return(output)

}
