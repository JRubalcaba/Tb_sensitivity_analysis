#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
  
// [[Rcpp::export]]
// HEAT TRANSFER

NumericVector heat_transfer(NumericVector Tb, double A, double C, double hc, double Te, double delta_T){
  NumericVector Tb_1(Tb.size());
  for(int i = 0; i < Tb.size(); i++){
    Tb_1[i] = Te + (Tb[i] - Te) * exp(-A/C * hc * delta_T);
  } 
  
  return Tb_1;
} 

// [[Rcpp::export]]
// PERFRORMANCE FUNCTION

NumericVector W_T(NumericVector Tb, double Topt) {
  NumericVector w(Tb.size());
  for(int i = 0; i < Tb.size(); i++){
    w[i] = abs(Tb[i] - Topt);
  }
  return w;
} 

// [[Rcpp::export]]
// BEHAVIORAL THERMOREGULATION

NumericVector pmultin(NumericVector W_i, NumericVector W_j, NumericVector q_j, double lambda){
  double Z = q_j[0] * exp(-lambda * W_j[0]);
  for(int i = 1; i < W_j.size(); i++){
    Z = Z + q_j[i] * exp(-lambda * W_j[i]);
  }
  NumericVector p(W_i.size());
  p[0] = q_j[0] * exp(-lambda * W_i[0]) / Z;
  for(int i = 1; i < W_i.size(); i++){
    p[i] = q_j[i] * exp(-lambda * W_j[i]) / Z;
  }
  return p;
} 
  
// [[Rcpp::export]]
// TRANSIENT-HEAT KERNEL

NumericMatrix kernel_q_i(NumericVector T_t, NumericVector T_t1, double Te, NumericVector mpar){
  double A = mpar["A"];
  double C = mpar["C"];
  double hc = mpar["hc"];
  double delta_T = mpar["delta_T"];
  double sigma = mpar["sigma"];

  NumericVector mu_i = heat_transfer(T_t, A, C, hc, Te, delta_T);

  NumericMatrix Q_i(T_t.size(), T_t1.size());
  for(int j=0; j < T_t.size(); j++){
    for(int k=0; k < T_t1.size(); k++){
      Q_i(j,k) = 1/(sigma*sqrt(2*3.1416))*exp(-(T_t1[k]-mu_i[j])*(T_t1[k]-mu_i[j])/(2*sigma*sigma));
    }
  }
  return Q_i;
}

// [[Rcpp::export]]
// THERMOREGULATION KERNEL

NumericMatrix kernel_p_i(NumericVector T_t, NumericVector Te, NumericVector q_j, NumericVector mpar){
  double A = mpar["A"];
  double C = mpar["C"];
  double hc = mpar["hc"];
  double delta_T = mpar["delta_T"];
  double Topt = mpar["Topt"];
  int n = Te.size();
  double lambda = mpar["lambda"];
  
  NumericMatrix w_i(T_t.size(), n);
  for(int i = 0; i < n; i++){
    double Te_i = Te[i];
    NumericVector T_t1 = heat_transfer(T_t, A, C, hc, Te_i, delta_T);
    w_i(_,i) = W_T(T_t1, Topt);
  }
  
  NumericMatrix P_i(T_t.size(), n);
  for(int i=0; i < T_t.size(); i++){
    P_i(i,_) = pmultin(w_i(i,_), w_i(i,_), q_j, lambda);
  }
  return P_i;
}

// [[Rcpp::export]]
// INTEGRATION

NumericMatrix matrix_operator(NumericMatrix P, NumericMatrix A, NumericVector v){
  NumericVector rows = A(_,1);
  NumericVector cols = A(1,_);
  
  NumericMatrix prod(rows.size(), cols.size());
  for(int i=0; i<rows.size(); i++){
    for(int j=0; j<cols.size(); j++){
      prod(i,j) = v[i] * A(i,j);
    }
  }
  NumericMatrix sum(rows.size(), cols.size());
  for(int i=0; i<rows.size(); i++){
    for(int j=0; j<cols.size(); j++){
      sum(i,j) = prod(i,j) + P(i,j);
    }
  }
  return sum;
}


// [[Rcpp::export]]

NumericMatrix kernel_P(NumericVector T_t, NumericVector T_t1, NumericVector Te, NumericVector q_j, NumericVector mpar){
  int n = Te.size();
  
  NumericMatrix P(T_t.size(),T_t1.size());
  for(int i=0; i<n; i++){
    double Te_i = Te[i];
    NumericMatrix Q_i = kernel_q_i(T_t, T_t1, Te_i, mpar);
    NumericMatrix P_i = kernel_p_i(T_t, Te, q_j, mpar);
    
    P = matrix_operator(P, transpose(Q_i), P_i(_,i));
  }
  return P;
}

// [[Rcpp::export]]
NumericVector GetEigenvector(NumericMatrix P){
  NumericVector cols = P(_,1);
  int n = cols.size();
  NumericVector v(n);
  NumericMatrix tempmat(n, n);
  for(int i=0; i<n; i++){
    v[i] = 1;
  }
  for(int k=0; k<90; k++){
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        tempmat(j,i) = v[i]*P(j,i);
      }
    }
    for(int i=0; i<n; i++){
      double sum = 0;
      for(int j=0; j<n; j++){
        sum = sum + tempmat(i,j);
      }
      v[i] = sum;
    }
  }
  double sum = 0;
  for(int i=0; i<n; i++){
    sum = sum + v[i];
  }
  for(int i=0; i<n; i++){
    v[i] = v[i]/sum;
  }
  return v;
}


// RUN METROPOLIS HASTINGS
////// [[Rcpp::export]]
//NumericMatrix Metrocpp(NumericVector T_t, NumericVector T_t1, NumericMatrix Te, NumericVector mpar, int iter, NumericVector data){
//  
//  NumericMatrix mcmc(iter, 3);
//  mcmc(0,0) = mpar["K"];
//  mcmc(0,1) = mpar["Topt"];
//  mcmc(0,2) = 0;
//  
//  NumericVector hour = Te(1,_);
//  NumericMatrix dailymat(T_t.size(), hour.size());
//  
//  for(int i=0; i<hour.size(); i++){
//    NumericVector Te_t = Te(_,i);
//    NumericMatrix P(T_t.size(), T_t1.size());
//    P = kernel_P(T_t, T_t1, Te_t, mpar);
//    NumericVector v = GetEigenvector(P);
//    
//    dailymat(_,i) = v;
//  }
//  
//  NumericVector null_dist = dailymat(_,0);
//  for(int i=1; i<hour.size(); i++){
//    for(int j=0; j<T_t.size(); j++){
//      null_dist[j] = null_dist[j]+dailymat(j,i);
//    }
//  }
//  for(int i=0; i<T_t.size(); i++){
//    null_dist[i] = null_dist[i]/14;
//  }
//  
//  for(int i=0; i<data.size(); i++){
//    int x_i = data[i];
//    mcmc(0,2) = mcmc(0,2) + log(null_dist[x_i-1]);
//  }
//  //mcmc(0,2) = mcmc(0,2) + log(R::rnorm(priorTopt, 4));

  /// RUN CHAIN
  //  for(int i=1; i<iter; i++){
  //    mcmc(i,0) = mcmc(i-1,0) + R::rnorm(0, 0.05);
  //  mcmc(i,1) = mcmc(i-1,1) + R::rnorm(0, 0.05);
  //  mcmc(i,2) = 0;
  //  
  //  mpar["K"] =  mcmc(i,0);
  //  mpar["Topt"] =  mcmc(i,1);
  //  
  //  for(int j=0; j<hour.size(); j++){
  //    NumericVector Te_t = Te(_,j);
  //    NumericMatrix P(T_t.size(), T_t1.size());
  //    P = kernel_P(T_t, T_t1, Te_t, mpar);
  //    NumericVector v = GetEigenvector(P);
  //    
  //    dailymat(_,j) = v;
  //  }
  //  
  //  NumericVector null_dist = dailymat(_,0);
  //  for(int j=1; j<hour.size(); j++){
  //    for(int k=0; k<T_t.size(); k++){
  ////      null_dist[k] = null_dist[k]+dailymat(k,j);
  //  }
  //  }
  ////  for(int j=0; j<T_t.size(); j++){
  //  null_dist[j] = null_dist[j]/14;
  //  }
  //  
  //  for(int j=0; j<data.size(); j++){
  //    int x_i = data[j];
  //    mcmc(i,2) = mcmc(i,2) + log(null_dist[x_i-1]);
  //  }
    //mcmc(0,2) = mcmc(0,2) + log(R::rnorm(priorTopt, 2));
    
    //  double prob = exp(mcmc(i-1,2) - mcmc(i,2));
    
    //  if(R::runif(0,1) < prob){
    //  mcmc(i,0) = mcmc(i-1,0);
    //  mcmc(i,1) = mcmc(i-1,1);
    //}else{
    //  mcmc(i,0) = mcmc(i,0);
    //  mcmc(i,1) = mcmc(i,1);
    //}
    //Rcout << "Iteration: " << i << " of " << iter << std::endl;
    //  }
  
  //return mcmc;
  //}