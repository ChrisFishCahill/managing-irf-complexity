#include <TMB.hpp>
/*
 * omniscient manager for spasmodic recruitment
 *                cahill oct 2022
 */
template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace Eigen;
  
  // Data
  DATA_INTEGER(n_year); 
  DATA_INTEGER(n_age);
  DATA_SCALAR(vbk);
  DATA_SCALAR(s);
  DATA_SCALAR(cr);
  DATA_VECTOR(ages);
  DATA_VECTOR(recmult);
  DATA_INTEGER(obj_ctl);
  
  // transformed data -- set up leading vectors
  vector<Type> ninit(n_age);
  vector<Type> w(n_age);
  vector<Type> vul(n_age);
  vector<Type> mwt(n_age);

  Type ro = 1;
  ninit(0) = 2;
  w(0) = pow((1 - exp(-vbk)),3);
  vul(0) = 1/(1 + exp(-0.5*(1-5)));
  mwt(0) = w(0) / (1 + exp(-0.5*(1 - 7)));

  for(int a = 1; a < n_age; a ++){
    ninit(a) = ninit(a-1)*s;
    w(a) = pow((1 - exp(-vbk*(a+1))), 3);
    vul(a) = 1 / (1 + exp(-0.5*(ages(a) - 5)));
    mwt(a) = w(a) / (1 + exp(-0.5*(ages(a) - 7)));
  }
  ninit(n_age-1) = ninit(n_age - 1) / ( 1 -s ); // plus group
  Type spro = 0;
  for(int a = 0; a < n_age; a++){spro += ninit(a)*mwt(a);}

  Type reca = cr / spro;
  Type recb = (cr - 1) / (ro * spro);
  
  // parameters to solve 
  PARAMETER_VECTOR(Ut);
  
  vector<Type> abar(n_year);
  vector<Type> yield(n_year);
  vector<Type> utility(n_year);
  vector<Type> ssb(n_year);
  vector<Type> n(n_age);
  abar.setZero();
  yield.setZero();
  utility.setZero();
  ssb.setZero();
  n.setZero();

  n = ninit;
  for(int t = 0; t < n_year; t++){
    vector<Type> sumterms(3);
    sumterms.setZero();
    for(int a = 0; a < n_age; a++){
      sumterms(0) += vul(a)*n(a)*w(a);                    // sumproduct(vul*n*w) across a
      sumterms(1) += mwt(a)*n(a);                         // sumproduct(mwt * n)
      sumterms(2) += ages(a)*n(a);
    }
    yield(t) = Ut(t)*sumterms(0);
    utility(t) = pow(yield(t),0.6);
    ssb(t) = sumterms(1);
    abar(t) = sumterms(2) / sum(n);                       //sumproduct(ages*n) / sum(n)
    for(int a = 0; a < n_age; a++){
      n(a) = s*n(a)*(1-vul(a)*Ut(t));
    }
    n(n_age - 1) = n(n_age - 1) + n(n_age - 2);           // plus group
    for(int a = (n_age-2); a --> 1; ){ n(a) = n(a - 1); } // advance fish one a
    n(0) = reca*ssb(t) / (1 + recb*ssb(t))*recmult(t);    // recruits
  }
  
  REPORT(ssb);
  REPORT(yield);
  REPORT(abar);
  REPORT(utility);
  REPORT(ro); 
  REPORT(spro); 
  REPORT(reca); 
  REPORT(recb); 
  REPORT(mwt); 
  REPORT(w); 
  REPORT(vul); 

  // objective function
  Type obj = 0;
  if(obj_ctl == 0){// yield objective
    obj -= sum(yield);
  }
  if(obj_ctl == 1){// hara utility objective
    obj -= sum(utility);
  }
  return obj; 
}
