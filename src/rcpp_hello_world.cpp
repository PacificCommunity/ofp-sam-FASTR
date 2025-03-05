
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}

// [[Rcpp::export]]
NumericMatrix YPR(NumericVector age, NumericVector sel, NumericVector m, NumericVector mat, NumericVector wgt, double alpha, double beta, double sigma, double inc) {

  int n = age.size();
  NumericVector N(n);
  NumericVector EqS(400);
  NumericVector EqY(400);
  NumericMatrix res(400,4);
  
  double F;
  double Z;
  double SPR;
  double YPR;
  double Yield;
  double rec;
  N[0]= 1;
  F   = 0;
  
  for(int i = 0; i < 400; i++){
    SPR = 1 * mat[0]*wgt[0];
    YPR = (F*sel[0])/(F*sel[0]+m[0])*N[0]*(1-exp(-F*sel[0]-m[0]));
    
    for(int j=1; j<n; j++){
      N[j] = N[j-1]*exp(-(F*sel[j-1]+m[j-1]));
      if(j==n-1){
        Z = F*sel[j] + m[j];
        N[j] = N[j]/(1.0-exp(-Z));    // plus group
      }
      SPR += (N[j]*mat[j]*wgt[j]);
      YPR += (F*sel[j])/(F*sel[j]+m[j])*N[j]*(1-exp(-F*sel[j]-m[j]));
    }
//    Rcout << "SPR " << SPR <<std::endl;
    
    EqS[i] = alpha*(exp(0.5*sigma))*SPR-beta;

    if(EqS[i] <=0)
      EqS[i] = 0;

    rec = (alpha*EqS[i])/(beta+EqS[i])*exp(0.5*sigma);

    Yield = 0;
    for(int j=0; j<=n; j++) {
	Yield += (F*sel[j]/(F*sel[j]+m[j])*N[j]*rec*(1-exp(-F*sel[j]-m[j]))) * wgt[j];
    }

    F = F+inc;

    EqY[i]  =Yield;
    res(i,0)=EqS[i];
    res(i,1)=EqY[i];
    res(i,2)=SPR;
    res(i,3)=YPR;
  }
//  return EqY;
  return res;
}



