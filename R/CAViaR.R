
#' CAViaR model estimation
#'
#' This function allows the estimation for a general order the caviar model proposed by Engle & Manganelli (2004).
#'
#' @param Y A vector, matrix, zoo or xts object containing the univariate series.
#' @param p Order of autoregressive quantile.
#' @param q Order of lag Y values.
#' @param model.type The specification, one of the following:
#' - "SAV" (*Symmetric absolute value*)
#'
#' \deqn{f_t(\theta) = \beta_0 + \sum_{i=1}^p \beta_i f_{t-i}(\beta) + \sum_{j=1}^q \gamma_j |y_{t-j}|}
#'
#' - "AS" (*Asymmetric slope*)
#'
#' \deqn{f_t(\theta) = \beta_0 + \sum_{i=1}^p \beta_i f_{t-i}(\beta) + \sum_{j=1}^q \left( \gamma_{1,j} (y_{t-j})^+ + \gamma_{2,j} (y_{t-j})^- \right)}
#'
#'
#'   where \eqn{(y_{t-j})^+ = \max(y_{t-j}, 0)} and \eqn{(y_{t-j})^- = \min(y_{t-j}, 0)}
#'
#' - "INDGARCH" (*Indirect GARCH*)
#'
#' \deqn{f_t(\theta) = \left( \beta_0 + \sum_{i=1}^p \beta_i f_{t-i}^2(\beta) + \sum_{j=1}^q \gamma_j y_{t-j}^2 \right)^{1/2}}
#'
#' - "I-CAV" (*Improved CAViaR*)
#'
#'  \deqn{f_t(\theta) = \beta_0 + \beta_1 f_{t-1}(\theta) + (1 - \beta_1)\left(\frac{\nu}{1-\gamma_1} I(y_{t-1} > 0)+ \frac{\nu}{\gamma_1} I(y_{t-1} \leq 0)\right) |y_{t-1} - u|}
#'
#'
#'  where \eqn{\nu = \sqrt{\gamma_1^2 + (1 - \gamma_1)^2}}, \eqn{0<\gamma_1<1} and \eqn{u} is the sample mean.
#'
#' @param tau The quantile of interest. Set default to 0.05.
#' @param band.hs Logical parameter passed to [quantreg::bandwidth.rq], if [TRUE] the Hall-Sheather bandwidth is computed, if [FALSE] it is computed the Bofinger bandwidth.
#' @param jac.method One of "Richardson" (default), "simple", or "complex". This determines the method to compute finite differences jacobian to as part of the standard errors calculations. See [numDeriv::jacobian] for more information.
#' @param sign.level The `alpha` parameter of [quantreg::bandwidth.rq] function to control the level of significance for intended confidence intervals.
#' @param quant.type One of the types available in [quantile]. This is used to initialize the quantile process.
#' @param optim.config A list containing the optimization-related parameters, see [nloptr.print.options()] for more information. See `details` for more information.
#' @param refine.estim Should the result from the global optimization be refined with a gradient-based solver?
#' @param refinement.conf A list containing the refinement optimization-related parameters, see [nloptr.print.options()] for more information
#' @param jac.options A list passed to \code{method.args} in [numDeriv::jacobian].
#'
#' @importFrom stats qnorm quantile embed pt qt
#' @importFrom utils modifyList
#' @importFrom quantreg bandwidth.rq
#' @importFrom nloptr nloptr varmetric
#' @importFrom numDeriv jacobian
#' @importFrom xts as.xts
#' @importFrom zoo index
#' @references
#' Engle, R. F., & Manganelli, S. (2004). CAViaR: Conditional Autoregressive Value at Risk by Regression Quantiles. Journal of Business &amp; Economic Statistics, 22(4), 367–381.
#'
#' Huang, D., Yu, B., Fabozzi, F. J., & Fukushima, M. (2009). CAViaR-based forecast for oil price risk. Energy Economics, 31(4), 511-518.
#'
#' White, H., Kim, T. H., & Manganelli, S. (2015). VAR for VaR: Measuring tail dependence using multivariate regression quantiles. Journal of econometrics, 187(1), 169-188.
#' @details
#'The main difference in this implementation is the optimization procedure used to estimate the parameters. The original paper begins with a grid search and alternates between the Nelder-Mead and quasi-Newton optimization algorithms until convergence. In this function, however, a global optimization algorithm is used: the Multilevel Single Linkage ("NLOPT_GN_MLSL_LDS"). At each iteration, the sampled points are quasi-deterministic, using Sobol's low discrepancy sequences. Thus, different optimization runs are much less dependent on the choice of random seed.
#'
#'As local optimizer, it is set by default to the `NELDER-MEAD`. However, other oprion is `SUBPLEX` which uses the `NELDER-MEAD` is subspaces, and is claimed to be more efficient and robust than the later. One can change this as described in the [nloptr::nloptr] function.
#'
#'An alternative to MLSL is the Improved Stochastic Ranking Evolution Strategy ("NLOPT_GN_ISRES"). However, this algorithm is much more dependent on the choice of the random seed. Although it is theoretically the case that it should return the same estimates among different runs, in practice, the user should be careful when choosing this algorithm. Furthermore, the standard errors are computed using the numerical Jacobian from `numDeriv`, and the bandwidth is computed as described in White et al. (2015).
#' @returns A list containing different results from the estimation. Class `CAViaR_estim`.
#' @export
#' @examples
#' data=dataCAViaR
#' SAV <- CAViaR(Y=data$GM[1:2892],model.type = "SAV",
#' p=1,q=1,band.hs = TRUE,quant.type = 7,
#' tau=0.05,refine.estim = FALSE)
#'
#' summary(SAV)
#' #or
#' SAV
#'
#' plot(SAV)
#'
CAViaR <- function(Y,p=1,q=1,model.type="SAV",tau=0.05,band.hs=FALSE,jac.method="Richardson",jac.options=list(),sign.level=0.05,quant.type=7,
                   optim.config=list(),refine.estim=FALSE,refinement.conf=list()){


  if(jac.method=="Richardson"){
    options <- list(eps=1e-6, d=0.00001, zero.tol=sqrt(.Machine$double.eps/7e-5), r=4, v=2, show.details=FALSE)
    jacobian.method <- modifyList(options,jac.options)
  }else{
    options <- list(eps=1e-6)
    jacobian.method <- modifyList(options,jac.options)
    }

  fechas <- as.Date(index(Y))
  serie <- unname(as.matrix(Y))
  T <- nrow(serie)
  VaR=residuos=matrix(0.0,ncol = 1,nrow = T)
  dimensions=c(p,q)
  maximum=max(dimensions)
  emp_quant <- quantile(serie[1:300,], probs=tau,type=quant.type)
  bounds=rep(-5,1+p+q)

  bandwithd <- bandwidth.rq(p=tau,n=T,hs=band.hs,alpha = sign.level)

  optimization.config <- list(algorithm="NLOPT_GN_MLSL_LDS",
                              population=10,ranseed=0,maxeval=20000,xtol_rel=1e-8,
                              local_opts=list(algorithm="NLOPT_LN_NELDERMEAD",maxeval=800,
                                              xtol_rel=1e-8,ftol_rel=0))

  optimization.config <- modifyList(optimization.config,optim.config)

  refinement.config <- list(maxeval=100*maximum,xtol_rel=1e-13,ftol_rel=0)

  refinement.config <- modifyList(refinement.config,refinement.conf)


    if(model.type=="SAV"){

      param_names <- c("Beta 0",paste0("Beta ",1:p),paste0("Gamma ",1:q))
      lagged_part <- abs(embed(Y,dimension = 1+q)[,-1])

      message("Begining optimization")
      estimacion <- nloptr(x0=rep(0.0,1+p+q),eval_f =OBJF_unc_caviar_SAV,lb = rep(-10,1+p+q) ,ub=rep(10,1+p+q),
                                   opts = optimization.config,
                                   yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                   emp_quantil=emp_quant,residuos=residuos,
                                   VaR_vector=VaR,p=p,q=q,maximum=maximum,
                                   time=T)
         if(refine.estim==TRUE){
      estimacion <- varmetric(x0 =estimacion$solution,fn=OBJF_unc_caviar_SAV,lower =rep(-10,1+p+q)
                          ,upper=rep(10,1+p+q),
                          yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                          emp_quantil=emp_quant,residuos=residuos,
                          VaR_vector=VaR,p=p,q=q,maximum=maximum,
                          time=T,control=refinement.config)
      estimacion <- list(solution=estimacion$par,objective=estimacion$value)
      }
      predict <- unc_caviar_sav(estimacion$solution,lagged_y=as.matrix(lagged_part),
                                emp_quantil=emp_quant,
                                VaR_vector=VaR,p=p,q=q,maximum=maximum,
                                time=T)
      jaco <-  jacobian(unc_caviar_sav,estimacion$solution,lagged_y=as.matrix(lagged_part),
                        emp_quantil=emp_quant,
                        VaR_vector=VaR,p=p,q=q,maximum=maximum,
                        time=T,method = jac.method,method.args = jacobian.method)

      residuoss <- serie-predict
      ka_te <- quantile(abs(residuoss-quantile(residuoss,probs=0.5)),probs=0.5)
      ce_te <- ka_te*(qnorm(tau+bandwithd)-qnorm(tau-bandwithd))
      message("Calculating Standard Errors")
      A_mat=crossprod(jaco)/T
      D_mat=crossprod(jaco[abs(residuoss)<ce_te,])/(2*ce_te*T)
      D_inv=solve(D_mat)
      VCV=tau*(1-tau)*(D_inv%*%A_mat%*%D_inv)/T
      se=sqrt(diag(VCV))
      parameters=estimacion$solution
      degrees=T-length(parameters)
      p_values=2*(1-pt(abs(parameters/se),df=degrees))
      upper_interval=parameters+qt(1-(sign.level/2),df=degrees)*se
      lower_interval=parameters+qt(sign.level/2,df=degrees)*se
      }else

      if(model.type=="AS"){
        pos_neg <- paste0("Gamma",c("+","-"))
        gammas <- paste0(pos_neg,",",rep(1:q,each=2))
        param_names <- c("Beta 0",paste0("Beta ",1:p),gammas)
        lagged_ <- embed(Y,dimension = 1+q)[,-1]
        n <- nrow(lagged_)
        positive <- (lagged_>0)*lagged_
        negative <- (lagged_<=0)*lagged_
        todo <- cbind(positive,negative)
        sequencia <- as.vector(rbind(1:q,(q+1):(2*q)))
        lagged_part <- todo[,sequencia]

        message("Begining optimization")
        estimacion <- nloptr(x0=rep(0.0,1+p+2*q),eval_f =OBJF_unc_caviar_AS,lb = rep(-10,1+p+2*q) ,ub=rep(10,1+p+2*q),
                                     opts = optimization.config,
                                     yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                     emp_quantil=emp_quant,residuos=residuos,
                                     VaR_vector=VaR,p=p,q=2*q,maximum=maximum,
                                     time=T)
        if(refine.estim==TRUE){
          estimacion <- varmetric(x0 =estimacion$solution,fn=OBJF_unc_caviar_AS,lower =rep(-10,1+p+2*q)
                                      ,upper=rep(10,1+p+2*q),
                                      yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                      emp_quantil=emp_quant,residuos=residuos,
                                      VaR_vector=VaR,p=p,q=2*q,maximum=maximum,
                                      time=T,control=refinement.config)
          estimacion <- list(solution=estimacion$par,objective=estimacion$value)
        }
        predict <- unc_caviar_as(estimacion$solution,lagged_y=as.matrix(lagged_part),
                                       emp_quantil=emp_quant,
                                       VaR_vector=VaR,p=p,q=2*q,maximum=maximum,
                                       time=T)
        jaco <-  jacobian(unc_caviar_as,estimacion$solution,lagged_y=as.matrix(lagged_part),
                                    emp_quantil=emp_quant,
                                    VaR_vector=VaR,p=p,q=2*q,maximum=maximum,
                                    time=T,method = jac.method,method.args = jacobian.method)

        residuoss <- serie-predict
        ka_te <- quantile(abs(residuoss-quantile(residuoss,probs=0.5)),probs=0.5)
        ce_te <- ka_te*(qnorm(tau+bandwithd)-qnorm(tau-bandwithd))
        message("Calculating Standard Errors")
        A_mat=crossprod(jaco)/T
        D_mat=crossprod(jaco[abs(residuoss)<ce_te,])/(2*ce_te*T)
        D_inv=solve(D_mat)
        VCV=tau*(1-tau)*(D_inv%*%A_mat%*%D_inv)/T
        se=sqrt(diag(VCV))
        parameters=estimacion$solution
        degrees=T-length(parameters)
        p_values=2*(1-pt(abs(parameters/se),df=degrees))
        upper_interval=parameters+qt(1-(sign.level/2),df=degrees)*se
        lower_interval=parameters+qt(sign.level/2,df=degrees)*se

      }else

      if(model.type=="INDGARCH"){
        param_names <- c("Beta 0",paste0("Beta ",1:p),paste0("Gamma ",1:q))
        lagged_part <- (embed(Y,dimension = 1+q)[,-1])^2

        message("Begining optimization")
        estimacion <- nloptr(x0=rep(0.5,1+p+q),eval_f =OBJF_unc_caviar_INDGARCH,lb = rep(0,1+p+q) ,ub=rep(10,1+p+q),
                                     opts = optimization.config,
                                     yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                     emp_quantil=emp_quant,residuos=residuos,
                                     VaR_vector=VaR,p=p,q=q,maximum=maximum,
                                     time=T)
        if(refine.estim==TRUE){
          estimacion <- varmetric(x0 =estimacion$solution,fn=OBJF_unc_caviar_INDGARCH,lower =rep(0,1+p+q)
                                      ,upper =rep(10,1+p+q),
                                      yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                      emp_quantil=emp_quant,residuos=residuos,
                                      VaR_vector=VaR,p=p,q=q,maximum=maximum,
                                      time=T,control=refinement.config)
          estimacion <- list(solution=estimacion$par,objective=estimacion$value)
        }
        predict <- unc_caviar_indgarch(estimacion$solution,lagged_y=as.matrix(lagged_part),
                                       emp_quantil=emp_quant,
                                       VaR_vector=VaR,p=p,q=q,maximum=maximum,
                                       time=T)
        jaco <-  jacobian(unc_caviar_indgarch,estimacion$solution,lagged_y=as.matrix(lagged_part),
                                   emp_quantil=emp_quant,
                                   VaR_vector=VaR,p=p,q=q,maximum=maximum,
                                   time=T,method = jac.method,method.args = jacobian.method)

        residuoss <- serie-predict
        ka_te <- quantile(abs(residuoss-quantile(residuoss,probs=0.5)),probs=0.5)
        ce_te <- ka_te*(qnorm(tau+bandwithd)-qnorm(tau-bandwithd))
        message("Calculating Standard Errors")
        A_mat=crossprod(jaco)/T
        D_mat=crossprod(jaco[abs(residuoss)<ce_te,])/(2*ce_te*T)
        D_inv=solve(D_mat)
        VCV=tau*(1-tau)*(D_inv%*%A_mat%*%D_inv)/T
        se=sqrt(diag(VCV))
        parameters=estimacion$solution
        degrees=T-length(parameters)
        p_values=2*(1-pt(abs(parameters/se),df=degrees))
        upper_interval=parameters+qt(1-(sign.level/2),df=degrees)*se
        lower_interval=parameters+qt(sign.level/2,df=degrees)*se

      }else
        if(model.type=="I-CAV"){
          param_names <- c("Beta 0","Beta 1","Gamma 1")
          lagged=embed(serie,2)
          lagged_part=abs(lagged[,2]-mean(serie))
          positivo=(lagged[,2]>0)*1

          message("Begining optimization")
          estimacion <- nloptr(x0=c(0,0,0.1),eval_f =OBJF_cons_caviar_ICAV,lb = c(-10,-10,-0.9999999) ,ub=c(10,10,0.999999),
                                       opts = optimization.config,
                                       yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                       emp_quantil=emp_quant,residuos=residuos,
                                       VaR_vector=VaR,positive=positivo,
                                       time=T)

          if(refine.estim==TRUE){
            estimacion <- varmetric(x0 =estimacion$solution,fn=OBJF_cons_caviar_ICAV,lower =c(-10,-10,-0.9999999)
                                        ,upper =c(10,10,0.999999),
                                        yes=serie,tau=tau,lagged_y=as.matrix(lagged_part),
                                        emp_quantil=emp_quant,residuos=residuos,
                                        VaR_vector=VaR,positive=positivo,
                                        time=T,control=refinement.config)
            estimacion <- list(solution=estimacion$par,objective=estimacion$value)

          }
          predict <- unc_caviar_i_caviar(estimacion$solution,lagged_abs_y=as.matrix(lagged_part),
                                         emp_quantil=emp_quant,
                                         VaR_vector=VaR,positive=positivo,
                                         time=T)
          jaco <-  jacobian(unc_caviar_i_caviar,estimacion$solution,lagged_abs_y=as.matrix(lagged_part),
                                      emp_quantil=emp_quant,
                                      VaR_vector=VaR,positive=positivo,
                                      time=T,method = jac.method,method.args = jacobian.method)

          residuoss <- serie-predict
          ka_te <- quantile(abs(residuoss-quantile(residuoss,probs=0.5)),probs=0.5)
          ce_te <- ka_te*(qnorm(tau+bandwithd)-qnorm(tau-bandwithd))
          message("Calculating Standard Errors")
          A_mat=crossprod(jaco)/T
          D_mat=crossprod(jaco[abs(residuoss)<ce_te,])/(2*ce_te*T)
          D_inv=solve(D_mat)
          VCV=tau*(1-tau)*(D_inv%*%A_mat%*%D_inv)/T
          se=sqrt(diag(VCV))
          parameters=estimacion$solution
          degrees=T-length(parameters)
          p_values=2*(1-pt(abs(parameters/se),df=degrees))
          upper_interval=parameters+qt(1-(sign.level/2),df=degrees)*se
          lower_interval=parameters+qt(sign.level/2,df=degrees)*se
        }

        else{stop("Wrong model.type selected, please make sure to choose one supported model or review your syntax ;)")}

  result <- cbind(parameters,se,p_values,lower_interval,upper_interval)
  ci <- paste0(c(100*sign.level/2,100*(1-(sign.level/2))),"% CI")
  colnames(result) <- c("Coef.","S.E","P>|t|",ci)
  rownames(result) <- param_names
  resultados <- structure(list(Results=result,data=as.xts(serie,order.by=fechas),VaR=as.xts(predict,order.by=fechas),loss=estimacion$objective,set_up=list(tau=tau,spec=model.type,band=band.hs,p=p,q=q,optimization=optimization.config,refinement=refinement.config,refine=refine.estim,q.type=quant.type)),class="CAViaR_estim")
  return(resultados)
}
