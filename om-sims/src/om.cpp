#include <TMB.hpp>
/*
 * omniscient manager for spasmodic recruitment
 *          cahill & walters oct 2022
 */
template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace Eigen;
  
  DATA_INTEGER(n_year); 
  DATA_INTEGER(n_age);
  DATA_SCALAR(vbk);      // von bertalanffy k 
  DATA_SCALAR(s);        // avg survival
  DATA_SCALAR(cr);       // compensation ratio
  DATA_SCALAR(rinit);    // initial number age 1 recruits
  DATA_SCALAR(uo);       // average historical exploitation rate
  DATA_SCALAR(asl);      // vul parameter 1
  DATA_SCALAR(ahv);      // vul parameter 2
  DATA_SCALAR(ahm);      // age half mature 
  DATA_SCALAR(upow);     // utility power
  DATA_VECTOR(ages); 
  DATA_VECTOR(recmult);  // recruitment sequence
  DATA_INTEGER(obj_ctl); // 0 = MAY, 1 = HARA utility
  
  // transformed data -- set up leading vectors
  vector<Type> n(n_age);
  vector<Type> w(n_age);
  vector<Type> vul(n_age);
  vector<Type> mwt(n_age);
  n.setZero(); w.setZero(); vul.setZero(); mwt.setZero(); 
  
  Type ro = rinit;  
  Type sso = 1; // avg survivorship to a
  n(0) = rinit;
  w(0) = pow((1 - exp(-vbk)),3);
  vul(0) = 1/(1 + exp(-asl*(ages(0) - ahv)));
  mwt(0) = w(0) / (1 + exp(-asl*(ages(0) - ahm)));

  for(int a = 1; a < n_age; a ++){
    sso *= s*(1-vul(a-1)*uo); 
    n(a) = rinit*sso;
    w(a) = pow((1 - exp(-vbk*(ages(a)))), 3);
    vul(a) = 1 / (1 + exp(-asl*(ages(a) - ahv)));
    mwt(a) = w(a) / (1 + exp(-0.5*(ages(a) - ahm)));
  }
  n(n_age-1) = n(n_age - 1) / ( 1 - s * vul(n_age - 1) * uo); // plus group
  Type spro = 0;
  for(int a = 0; a < n_age; a++){spro += n(a)*mwt(a);}
  Type reca = cr / spro;
  Type recb = (cr - 1) / (ro * spro);
  
  // parameters to solve 
  PARAMETER_VECTOR(Ut);
  
  // vectors for results
  vector<Type> abar(n_year);
  vector<Type> yield(n_year);
  vector<Type> utility(n_year);
  vector<Type> ssb(n_year);
  vector<Type> vulb(n_year); 
  abar.setZero(); yield.setZero(); utility.setZero(); ssb.setZero(); vulb.setZero(); 
  for(int t = 0; t < n_year; t++){
    vector<Type> sumterms(1);
    sumterms.setZero();
    for(int a = 0; a < n_age; a++){
      vulb(t) += vul(a)*n(a)*w(a);                          // sumproduct(vul*n*w) across a
      ssb(t) += mwt(a)*n(a);                                // sumproduct(mwt * n)
      sumterms(0) += ages(a)*n(a);                          // sumproduct(n,a)
    }
    yield(t) = Ut(t)*vulb(t);
    utility(t) = pow(yield(t), upow);
    abar(t) = sumterms(0) / sum(n);                         // sumproduct(ages*n) / sum(n)
    for(int a = 0; a < n_age; a++){
      n(a) = s*n(a)*(1-vul(a)*Ut(t));
    }
    n(n_age - 1) = n(n_age - 1) + n(n_age - 2);             // plus group
    for(int a = (n_age - 2); a > 0; a--){n(a) = n(a - 1);}  // advance fish one a
    n(0) = reca*ssb(t) / (1 + recb*ssb(t))*recmult(t);      // recruits
  }
  
  REPORT(ssb);
  REPORT(yield);
  REPORT(vulb); 
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
  if(obj_ctl == 0){  // yield objective
    obj -= sum(yield);
  }
  if(obj_ctl == 1){  // hara utility objective
    obj -= sum(utility);
  }
  return obj; 
}
