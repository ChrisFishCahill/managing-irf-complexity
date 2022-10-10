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
  
  // set up leading vectors
  vector<Type> n(n_age);
  vector<Type> vul(n_age);
  vector<Type> wt(n_age);
  vector<Type> mat(n_age);
  vector<Type> Lo(n_age);   
  vector<Type> mwt(n_age);
  vector<Type> Lf(n_age);
  n.setZero(); vul.setZero(); wt.setZero(); mat.setZero();
  Lo.setZero(); mwt.setZero(); Lf.setZero(); 
  Type sbro = 0;  
  Type ro = rinit;  
  wt(0) = pow((1 - exp(-vbk)),3);
  vul(0) = 1/(1 + exp(-asl*(ages(0) - ahv)));
  mwt(0) = wt(0) / (1 + exp(-asl*(ages(0) - ahm)));
  
  for(int a = 0; a < n_age; a ++){
    vul(a) =1 /( 1 + exp(-asl*(ages(a) - ahv))); 
    wt(a) = pow((1 - exp(-vbk*(ages(a)))), 3); 
    mat(a) = 1/(1 + exp(-asl*(ages(a) - ahm))); 
    if(a == 0){ 
      Lo(a) = 1;
      Lf(a) = 1; 
      n(a) = rinit*Lf(a);
    } 
    if(a > 0){
      Lo(a) = Lo(a-1)*s;
      Lf(a) = Lf(a-1)*s*(1 - vul(a-1)*uo);
      if(a == (n_age - 1)){ // plus group 
        Lo(a) = Lo(a) / (1 - s); 
        Lf(a) = Lf(a - 1) * s * (1 - vul(a-1)*uo) / (s * (1 - vul(a-1)*uo)); 
      }
      n(a) = n(a-1)*Lf(a);
    }
    mwt(a) = mat(a)*wt(a); 
    sbro += Lo(a)*mwt(a); 
  } 
  Type reca = cr/sbro; 
  Type recb = (cr - 1) / (ro*sbro); 
  
  // parameters to solve 
  PARAMETER_VECTOR(Ut);
  
  // vectors for results
  vector<Type> abar(n_year);
  vector<Type> yield(n_year);
  vector<Type> utility(n_year);
  vector<Type> ssb(n_year);
  vector<Type> vulb(n_year); 
  abar.setZero(); yield.setZero(); utility.setZero(); 
  ssb.setZero(); vulb.setZero(); 

  for(int t = 0; t < n_year; t++){
    vector<Type> sumterms(1);
    sumterms.setZero();
    for(int a = 0; a < n_age; a++){
      vulb(t) += vul(a)*n(a)*wt(a);                         // sumproduct(vul*n*w) across a
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
  REPORT(sbro); 
  REPORT(reca); 
  REPORT(recb); 
  REPORT(mwt); 
  REPORT(wt); 
  REPORT(vul); 
  REPORT(n);
  REPORT(Lf); 
  REPORT(Lo); 

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