data {
  int<lower = 0> n_year; 
  real<lower = 0> vbk; 
  real<lower = 0, upper = 1> surv; 
  real<lower = 0> cost; 
  vector<lower = 0, upper = 1>[n_year] R;
  int obj_ctl; // 0 = yield, 1 = utility
}
parameters {
  vector<lower = 0, upper = 1>[n_year] Ut; 
}
transformed parameters {
  vector<lower=0>[n_year] Adom;   // age of Dominant cohort
  vector<lower=0>[n_year] Wdom;   // weight of Dominant cohort
  vector<lower=0>[n_year] Ndom;   // relative numbers Dom cohort
  vector<lower=0>[n_year] Nresid; // number from Dom cohort left over
  vector<lower=0>[n_year] B;      // biommass 
  vector<lower=0>[n_year] yield;       
  vector<lower=0>[n_year] utility;      

  real Ut_sum = sum(Ut); 
  
  // initialize
  Adom[1] = 1; 
  Ndom[1] = 1; 
  Wdom[1] = (1-exp(-vbk * Adom[1]))^3; 
  Nresid[1] = 0; 
  B[1] = Ndom[1] * Wdom[1] + Nresid[1];
  yield[1] = Ut[1]*B[1]; 
  utility[1] = yield[1]^0.6;
  
  for(t in 2:n_year){
    if(R[t-1] == 0){
      Adom[t] = Adom[t - 1] + 1;                          // increment age by one year
      Ndom[t] = Ndom[t - 1] * surv * (1 - Ut[t - 1]);
      Nresid[t] = Nresid[t - 1] * surv * (1 - Ut[t - 1]); // carry over residual 
    }
    else {
      Adom[t] = 1; 
      Ndom[t] = R[t - 1]; 
      Nresid[t] = Ndom[t - 1] * surv * (1 - Ut[t - 1]); 
    }
    Wdom[t] = (1 - exp(-vbk * Adom[t]))^3; 
    B[t] = Ndom[t] * Wdom[t] + Nresid[t];
    yield[t] = Ut[t] * B[t]; 
    utility[t] = yield[t]^0.6; 
  }
}
model {
 // implicit Ut ~ uniform(0,1); 
 // loop through years and add up yield or utility
 for(t in 1:n_year){
   if(obj_ctl == 0){
     target += sum(yield) - cost * Ut_sum; 
   }
   if(obj_ctl == 1){
     target += sum(utility) - cost * Ut_sum; 
   }
   // print("target = ", target()); 
 }
}
