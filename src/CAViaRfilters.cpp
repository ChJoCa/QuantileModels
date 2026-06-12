//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

 //[[Rcpp::export]]
 arma::vec unc_caviar_indgarch(arma::rowvec coefficients,arma::uword p,arma::uword q, arma::uword maximum,
                          arma::mat lagged_y, double emp_quantil, arma::vec VaR_vector,arma::uword time){

   arma::rowvec ar_coef=coefficients.cols(1,p);

   arma::vec lagged_component=lagged_y*coefficients.cols(p+1,p+q).t();

   VaR_vector.rows(0,maximum-1)+=emp_quantil*emp_quantil;

   double ar_part=0.0;

   for(arma::uword t=maximum;t<time;t++){
      ar_part=0.0;
     for(arma::uword i=0;i<p;i++){
       ar_part+=ar_coef(i)*VaR_vector(t-1-i);
     }

     VaR_vector(t)= coefficients(0)+ar_part+lagged_component(t-q);

   }
   return arma::sign(emp_quantil)*arma::sqrt(VaR_vector);
 }


//[[Rcpp::export]]
double OBJF_unc_caviar_INDGARCH(arma::rowvec coefficients,const arma::vec& yes,double tau,arma::uword p,arma::uword q, arma::uword maximum,
                                arma::mat lagged_y, double emp_quantil,arma::vec residuos,arma::vec VaR_vector,arma::uword time){

  double loss=0.0;
  residuos=yes-unc_caviar_indgarch(coefficients,p,q,maximum,lagged_y,emp_quantil,VaR_vector,time);

  loss=arma::accu(tau*residuos.elem(arma::find(residuos>0.0)))+arma::accu((tau-1)*residuos.elem(arma::find(residuos<0.0)));

return loss;
}




//[[Rcpp::export]]
arma::vec unc_caviar_sav(arma::rowvec coefficients,arma::uword p,arma::uword q, arma::uword maximum,
                           arma::mat lagged_y, double emp_quantil, arma::vec VaR_vector,arma::uword time){

  arma::rowvec ar_coef=coefficients.cols(1,p);

  arma::vec lagged_component=lagged_y*coefficients.cols(p+1,p+q).t();

  VaR_vector.rows(0,maximum-1)+=emp_quantil;

  double ar_part=0.0;

  for(arma::uword t=maximum;t<time;t++){
    ar_part=0.0;
    for(arma::uword i=0;i<p;i++){
      ar_part+=ar_coef(i)*VaR_vector(t-1-i);
    }

    VaR_vector(t)= coefficients(0)+ar_part+lagged_component(t-q);

  }
  return VaR_vector;
}

//[[Rcpp::export]]
 double OBJF_unc_caviar_SAV(arma::rowvec coefficients,const arma::vec& yes,double tau,arma::uword p,arma::uword q, arma::uword maximum,
                                 arma::mat lagged_y, double emp_quantil,arma::vec residuos,arma::vec VaR_vector,arma::uword time){

   double loss=0.0;
   residuos=yes-unc_caviar_sav(coefficients,p,q,maximum,lagged_y,emp_quantil,VaR_vector,time);

   loss=arma::accu(tau*residuos.elem(arma::find(residuos>0.0)))+arma::accu((tau-1)*residuos.elem(arma::find(residuos<0.0)));

   return loss;
 }






 //[[Rcpp::export]]
 arma::vec unc_caviar_as(arma::rowvec coefficients,arma::uword p,arma::uword q, arma::uword maximum,
                           arma::mat lagged_y, double emp_quantil, arma::vec VaR_vector,arma::uword time){


   arma::rowvec ar_coef=coefficients.cols(1,p);
   arma::vec lagged_component=lagged_y*coefficients.cols(p+1,p+q).t();

   VaR_vector.rows(0,maximum-1)+=emp_quantil;
   double ar_part=0.0;
   arma::uword qu=q/2;
   for(arma::uword t=maximum;t<time;t++){
      ar_part=0.0;
     for(arma::uword i=0;i<p;i++){
       ar_part+=ar_coef(i)*VaR_vector(t-1-i);
     }

     VaR_vector(t)= coefficients(0)+ar_part+lagged_component(t-qu);

   }
   return VaR_vector;
 }

 //[[Rcpp::export]]
 double OBJF_unc_caviar_AS(arma::rowvec coefficients,const arma::vec& yes,double tau,arma::uword p,arma::uword q, arma::uword maximum,
                            arma::mat lagged_y, double emp_quantil,arma::vec residuos,arma::vec VaR_vector,arma::uword time){

   double loss=0.0;
   residuos=yes-unc_caviar_as(coefficients,p,q,maximum,lagged_y,emp_quantil,VaR_vector,time);

   loss=arma::accu(tau*residuos.elem(arma::find(residuos>0.0)))+arma::accu((tau-1)*residuos.elem(arma::find(residuos<0.0)));

   return loss;
 }



 //[[Rcpp::export]]
 arma::vec unc_caviar_i_caviar(arma::rowvec coefficients,
                               arma::vec lagged_abs_y, double emp_quantil,arma::vec positive, arma::vec VaR_vector,arma::uword time){

   double gamma=coefficients(2);

   double nu=arma::sign(emp_quantil)*std::sqrt((gamma*gamma)+((1.0-gamma)*(1.0-gamma)));

   VaR_vector(0)=emp_quantil;

   arma::vec determined_part=(nu/(1.0-gamma)*positive+(nu/gamma)*(-(positive-1.0)))%lagged_abs_y;

   for(arma::uword t=1;t<time;t++){

 VaR_vector(t)= coefficients(0)+coefficients(1)*VaR_vector(t-1)+(1.0-coefficients(1))*determined_part(t-1);

   }
   return VaR_vector;
 }
 //[[Rcpp::export]]
 double OBJF_cons_caviar_ICAV(arma::rowvec coefficients,const arma::vec& yes,double tau,arma::vec positive,
                            arma::mat lagged_y, double emp_quantil,arma::vec residuos,arma::vec VaR_vector,arma::uword time){

   double loss=0.0;

   residuos=yes-unc_caviar_i_caviar(coefficients,lagged_y,emp_quantil,positive,VaR_vector,time);

   loss=arma::accu(tau*residuos.elem(arma::find(residuos>0.0)))+arma::accu((tau-1)*residuos.elem(arma::find(residuos<0.0)));

   return loss;
 }




 // MVMQ-CAViaR

 //[[Rcpp::export]]
 arma::mat MVMQ_FILTER(arma::rowvec coefficients,
                               arma::mat lagged_abs_y,arma::uword N,arma::uword p,arma::uword q, arma::mat VaR_vector,arma::uword time,arma::uword maximum){

   arma::rowvec constants(coefficients.memptr(),N,false,true);
   arma::mat coef_quantiles(coefficients.memptr()+N,N*p,N,false,true);
   arma::mat coef_lagged(coefficients.memptr()+(N+N*N*p),N,N*q,false,true);

   arma::mat non_recursive_part=lagged_abs_y*coef_lagged.t();

   arma::mat cont_var;

   for(arma::uword t=maximum;t<time;t++){

      cont_var= arma::reverse(VaR_vector.rows(t-p,t-1),0).as_row();

     VaR_vector.row(t)= constants+(cont_var*coef_quantiles)+non_recursive_part.row(t-q);

   }
   return VaR_vector;
 }

 //[[Rcpp::export]]
 arma::vec SAV_FILTER_for_jaco(arma::rowvec coefficients, arma::uword serie,
                       arma::mat lagged_abs_y,arma::uword N,arma::uword p,arma::uword q, arma::mat VaR_vector,arma::uword time,arma::uword maximum){


   arma::rowvec constants(coefficients.memptr(),N,false,true);
   arma::mat coef_quantiles(coefficients.memptr()+N,N*p,N,false,true);
   arma::mat coef_lagged(coefficients.memptr()+(N+N*N*p),N,N*q,false,true);

   arma::mat non_recursive_part=lagged_abs_y*coef_lagged.t();

   arma::mat cont_var;

   for(arma::uword t=maximum;t<time;t++){

     cont_var= arma::reverse(VaR_vector.rows(t-p,t-1),0).as_row();

     VaR_vector.row(t)= constants+(cont_var*coef_quantiles)+non_recursive_part.row(t-q);

   }
   return VaR_vector.col(serie);
 }


 //[[Rcpp::export]]
 double OBJ_F_MVMQ(arma::rowvec coefficients,
                       arma::mat lagged_abs_y,arma::uword N,arma::uword p,arma::uword q, arma::mat VaR_vector,arma::uword time,arma::uword maximum,arma::vec tau,arma::mat residuos,arma::mat series){

   residuos=series-MVMQ_FILTER(coefficients,lagged_abs_y,N,p,q,VaR_vector,time,maximum);
   arma::mat hits=-arma::conv_to<arma::mat>::from(residuos<0.0);
   hits.each_row()+=tau.t();
   hits%=residuos;
   return arma::accu(hits);
 }
