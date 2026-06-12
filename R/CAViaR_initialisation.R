

#' Univariate caviar for MVMQ-CAViaR initialisation jiji
#' @noRd

CAViaR_init <- function(Y,p=1,q=1,tau=0.05,quant.type=7,optim.config=list()){

  serie <- unname(as.matrix(Y))
  T <- nrow(serie)
  VaR=residuos=matrix(0.0,ncol = 1,nrow = T)
  maximum=max(c(p,q))
  emp_quant <- quantile(serie[1:300,], probs=tau,type=quant.type)

  optimization.config <-  list(algorithm="NLOPT_GN_MLSL_LDS",
                               population=100,ranseed=0,maxeval=30000,xtol_rel=1e-8,
                               local_opts=list(algorithm="NLOPT_LN_NELDERMEAD",maxeval=1000,
                                               xtol_rel=1e-8,ftol_rel=0))

  optimization.config <- modifyList(optimization.config,optim.config)

  lagged_part <- abs(embed(Y,dimension = 1+q)[,-1])

  estimacion <- nloptr::nloptr(x0=rep(0.0,1+p+q),eval_f =OBJF_unc_caviar_SAV,lb = rep(-10,1+p+q) ,ub=rep(10,1+p+q),
                               opts = optimization.config,
                               yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                               emp_quantil=emp_quant,residuos=residuos,
                               VaR_vector=VaR,p=p,q=q,maximum=maximum,
                               time=T)
  return(estimacion$solution)

  }
