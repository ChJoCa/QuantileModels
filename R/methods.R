#' Plot CAViaR
#'
#' @param x `CAViaR_estim` class object.
#' @param .by Frequency to display the dates in the plot, see [axis.Date].
#' @param .format Format of the displayed dates in the plot, see [axis.Date].
#' @param titl Optional title of the plot.
#' @param ... other arguments to plot
#'
#' @method plot CAViaR_estim
#' @importFrom graphics axis.Date lines par
#' @export
plot.CAViaR_estim <- function(x,.by="month",.format="%b-%Y",titl="VaR",...){
  fechas <- as.Date(zoo::index(x$data))
  ene <- length(fechas)
  limits <- ceiling(max(abs(cbind(x$VaR,x$data))))
  titulo <- paste0(titl)
  plot(fechas,as.vector(x$data),xaxt="n",xlab = "",ylab = "",ylim = c(-limits,limits),main = titulo,las=1,type = "l",cex.axis=0.8,...)
  lines(fechas,as.vector(x$VaR),col="red",type ="l")
  axis.Date(side=1,at=seq(fechas[1],fechas[ene],by=.by),format=.format,cex.axis=0.8,las=2)
}

#' Plot MVMQ_CAViaR
#'
#' @param x `MVMQ_CAViaR` class object.
#' @param rows Number of rows to display in the plot, passed to [par]
#' @param columns Number of columns to display in the plot, passed to [par]
#' @param .by Frequency to display the dates in the plot, see [axis.Date].
#' @param .format Format of the displayed dates in the plot, see [axis.Date].
#' @param titl Optional title of the plot.
#' @param ... other arguments to plot
#'
#' @method plot MVMQ_CAViaR
#' @importFrom graphics axis.Date lines par
#' @export
plot.MVMQ_CAViaR <- function(x,rows=2,columns=1,.by="month",.format="%b-%Y",titl="VaR at ",...){
  par(mfrow=c(rows,columns))
  fechas <- as.Date(zoo::index(x$data))
  ene <- length(fechas)
  equations_names <- colnames(x$VaR)
  for (j in equations_names) {
    limits <- ceiling(max(abs(cbind(x$VaR[,j],x$data[,j]))))
    titulo <- paste0(titl,x$tau[j]*100,"%",": ",j)
    plot(fechas,as.vector(x$data[,j]),xaxt="n",xlab = "",ylab = "",ylim = c(-limits,limits),main = titulo,las=1,type = "l",cex.axis=0.8,...)
    lines(fechas,as.vector(x$VaR[,j]),col="red",type ="l")
    axis.Date(side=1,at=seq(fechas[1],fechas[ene],by=.by),format=.format,cex.axis=0.8,las=2)

    }

par(mfrow=c(1,1))
  }

#' Summary CAViaR
#'
#' @param object `CAViaR_estim` class object.
#' @param conf.level Confidence level to the coverage test.
#' @param ... Other arguments passed to [print].
#' @param digits number of decimals to display, passed to [round].
#'
#' @method summary CAViaR_estim
#' @importFrom ufRisk covtest
#' @export
summary.CAViaR_estim <- function(object,digits=5,conf.level=0.95,...){
  cat("CAViaR estimation \n")
  cat("-------------------\n")
  cat("Model specification:",object$set_up$spec,"\n")
  cat("Quantile (tau):",object$set_up$tau,"\n")
  cat("Loss function value at estimates:",object$loss,"\n")
  cat("In-sample coverage:",round(sum(1*(object$data<object$VaR))/length(object$data),digits = digits),"\n")
  if(object$set_up$band){cat("Hall-Sheather bandwidth \n")}else{cat("Bofinger bandwidth \n")}
  cat("\n")
  cat("Estimation results: \n")
  cat(paste0(rep("=",(51+digits)),collapse = ""),"\n")
  rounded <- round(object$Results,digits =digits )
  print(rounded,...)
  cat(paste0(rep("=",(51+digits)),collapse = ""),"\n")
  cat("\n")
  cat("Coverage test \n")
  test <- covtest(obj = list(Loss=as.vector(object$data),VaR=as.vector(object$VaR),p=object$set_up$tau),conflvl = conf.level)
  cat("-------------------\n")
  cat("Kupiec conditional coverage test (LRcc), p-value:",round(test$p.cc,digits = digits),"\n")
  cat("Christoffersen independence test (LRind), p-value:",round(test$p.ind,digits = digits),"\n")
  cat("Christoffersen unconditional coverage test (LRuc), p-value:",round(test$p.uc,digits = digits),"\n")
  }


#' Summary MVMQ_CAViaR
#'
#' @param object `MVMQ_CAViaR` class object.
#' @param conf.level Confidence level to the coverage test.
#' @param ... Other arguments passed to [print].
#' @param digits number of decimals to display ,passed to [round].
#'
#' @method summary MVMQ_CAViaR
#' @importFrom ufRisk covtest
#' @export
summary.MVMQ_CAViaR <- function(object,digits=5,conf.level=0.95,...){
  cat("MVMQ CAViaR estimation \n")

  cat("Loss function at estimates:",round(object$loss,digits = digits),"\n")
  if(object$bandwithd){cat("Hall-Sheather bandwidth \n")}else{cat("Bofinger bandwidth \n")}
  equations_names <- rownames(object$A_mat)
  regressores_names <- c("Const.",colnames(object$A_mat),colnames(object$B_mat))
  ci <- paste0(c(100*object$CI_level/2,100*(1-(object$CI_level/2))),"% CI")
  results_names <- c("Coef.","S.E","P>|t|",ci)
  cat(paste0(rep("=",(51+digits)),collapse = ""),"\n")
  for (i in equations_names) {
    cat("Equation:", i,"\n")
    cat("Quantile (tau):",object$tau[i],"\n")
    cat("In sample coverage",round(sum(object$hits_seq[,i])/length(object$hits_seq[,i]),digits = digits),"\n")
    results <- cbind(c(object$Cons[i,],object$A_mat[i,],object$B_mat[i,]),
                     c(object$SE$constants[i,],object$SE$A_mat[i,],object$SE$B_mat[i,]),
                     c(object$Pval$constants[i,],object$Pval$A_mat[i,],object$Pval$B_mat[i,]),
                     c(object$CI_low$constants[i,],object$CI_low$A_mat[i,],object$CI_low$B_mat[i,]),
                     c(object$CI_up$constants[i,],object$CI_up$A_mat[i,],object$CI_up$B_mat[i,]))

    colnames(results) <- results_names
    rownames(results) <- regressores_names
    cat("Estimation results: \n")
    cat(paste0(rep("-",(51+digits)),collapse = ""),"\n")
    print(round(results,digits = digits),...)
    cat("\n")
    cat("Coverage test \n")
    test <- covtest(obj = list(Loss=as.vector(object$data[,i]),VaR=as.vector(object$VaR[,i]),p=object$tau[i]),conflvl = conf.level)
    cat("-------------------\n")
    cat("Kupiec conditional coverage test (LRcc), p-value:",round(test$p.cc,digits = digits),"\n")
    cat("Christoffersen independence test (LRind), p-value:",round(test$p.ind,digits = digits),"\n")
    cat("Christoffersen unconditional coverage test (LRuc), p-value:",round(test$p.uc,digits = digits),"\n")
    cat(paste0(rep("-",(51+digits)),collapse = ""),"\n")
    cat("\n")
  }

  }



#' Print MVMQ_CAViaR
#'
#' @param x `MVMQ_CAViaR` class object.
#' @param conf.level Confidence level to the coverage test.
#' @param ... Other arguments passed to [print].
#' @param digits number of decimals to display, passed to [round].
#'
#' @method print MVMQ_CAViaR
#' @importFrom ufRisk covtest
#' @export
print.MVMQ_CAViaR <- function(x,digits=5,conf.level=0.95,...){
  cat("MVMQ CAViaR estimation \n")

  cat("Loss function at estimates:",round(x$loss,digits = digits),"\n")
  if(x$bandwithd){cat("Hall-Sheather bandwidth \n")}else{cat("Bofinger bandwidth \n")}
  equations_names <- rownames(x$A_mat)
  regressores_names <- c("Const.",colnames(x$A_mat),colnames(x$B_mat))
  ci <- paste0(c(100*x$CI_level/2,100*(1-(x$CI_level/2))),"% CI")
  results_names <- c("Coef.","S.E","P>|t|",ci)
  cat(paste0(rep("=",(51+digits)),collapse = ""),"\n")
  for (i in equations_names) {
    cat("Equation:", i,"\n")
    cat("Quantile (tau):",x$tau[i],"\n")
    cat("In sample coverage",round(sum(x$hits_seq[,i])/length(x$hits_seq[,i]),digits = digits),"\n")
    results <- cbind(c(x$Cons[i,],x$A_mat[i,],x$B_mat[i,]),
                     c(x$SE$constants[i,],x$SE$A_mat[i,],x$SE$B_mat[i,]),
                     c(x$Pval$constants[i,],x$Pval$A_mat[i,],x$Pval$B_mat[i,]),
                     c(x$CI_low$constants[i,],x$CI_low$A_mat[i,],x$CI_low$B_mat[i,]),
                     c(x$CI_up$constants[i,],x$CI_up$A_mat[i,],x$CI_up$B_mat[i,]))

    colnames(results) <- results_names
    rownames(results) <- regressores_names
    cat("Estimation results: \n")
    cat(paste0(rep("-",(51+digits)),collapse = ""),"\n")
    print(round(results,digits = digits),...)
    cat("\n")
    cat("Coverage test \n")
    test <- covtest(obj = list(Loss=as.vector(x$data[,i]),VaR=as.vector(x$VaR[,i]),p=x$tau[i]),conflvl = conf.level)
    cat("-------------------\n")
    cat("Kupiec conditional coverage test (LRcc), p-value:",round(test$p.cc,digits = digits),"\n")
    cat("Christoffersen independence test (LRind), p-value:",round(test$p.ind,digits = digits),"\n")
    cat("Christoffersen unconditional coverage test (LRuc), p-value:",round(test$p.uc,digits = digits),"\n")
    cat(paste0(rep("-",(51+digits)),collapse = ""),"\n")
    cat("\n")
  }

}

#' Print CAViaR
#'
#' @param x `CAViaR_estim` class object.
#' @param conf.level Confidence level to the coverage test.
#' @param ... Other arguments passed to [print].
#' @param digits number of decimals to display, passed to [round].
#'
#' @method print CAViaR_estim
#' @importFrom ufRisk covtest
#' @export
print.CAViaR_estim <- function(x,digits=5,conf.level=0.95,...){
  cat("CAViaR estimation \n")
  cat("-------------------\n")
  cat("Model specification:",x$set_up$spec,"\n")
  cat("Quantile (tau):",x$set_up$tau,"\n")
  cat("Loss function value at estimates:",x$loss,"\n")
  cat("In-sample coverage:",round(sum(1*(x$data<x$VaR))/length(x$data),digits = digits),"\n")
  if(x$set_up$band){cat("Hall-Sheather bandwidth \n")}else{cat("Bofinger bandwidth \n")}
  cat("\n")
  cat("Estimation results: \n")
  cat(paste0(rep("=",(51+digits)),collapse = ""),"\n")
  rounded <- round(x$Results,digits =digits )
  print(rounded,...)
  cat(paste0(rep("=",(51+digits)),collapse = ""),"\n")
  cat("\n")
  cat("Coverage test \n")
  test <- covtest(obj = list(Loss=as.vector(x$data),VaR=as.vector(x$VaR),p=x$set_up$tau),conflvl = conf.level)
  cat("-------------------\n")
  cat("Kupiec conditional coverage test (LRcc), p-value:",round(test$p.cc,digits = digits),"\n")
  cat("Christoffersen independence test (LRind), p-value:",round(test$p.ind,digits = digits),"\n")
  cat("Christoffersen unconditional coverage test (LRuc), p-value:",round(test$p.uc,digits = digits),"\n")
}

