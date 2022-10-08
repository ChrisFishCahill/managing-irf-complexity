/*
* Omniscient manager for spasmodic recruitment 
         Cahill and Walters Oct 2022
*/
data {
  int<lower = 0> n_year; 
  int<lower = 0> n_age; 
  real<lower = 0> vbk;                    // von Bertalanffy k 
  real<lower = 0, upper = 1> s;           // survival 
  real<lower = 0> cr;                     // goodyear compensation ratio
  vector[n_year] recmult;                 // recruitment multiplier
  int obj_ctl;                            // 0 = yield, 1 = utility
}
transformed data {
  vector<lower=0>[n_age] ninit;           // initial numbers at age
  vector<lower=0>[n_age] w;               // weight at age
  vector<lower=0>[n_age] vul; 
  vector<lower=0>[n_age] mwt;             // maturity x weight for SSB calc 
  vector<lower=0>[n_age] ages;               
  real spro = 0; 
  real ro = 1; 
  
  // initialize leading vectors
  ninit[1] = 2;                              
  w[1] = (1 - exp(-vbk))^3; 
  vul[1] = 1/(1 + exp(-0.5*(1 - 5))); 
  mwt[1] = w[1] /(1 + exp(-0.5*(1 - 7))); 
  ages[1] = 1; 
  for(a in 2:n_age){
    ninit[a] = ninit[a - 1] * s; 
    w[a] = (1 - exp(-vbk*a))^3; 
    vul[a] = 1 / (1 + exp(-0.5*(a - 5))); 
    mwt[a] = w[a] / (1 + exp(-0.5*(a - 7)));  
    ages[a] = a;
  }
  
  ninit[n_age] = ninit[n_age] / (1-s); 
  
  for(a in 1:n_age){spro += ninit[a]*mwt[a];} 
}
parameters {
  vector<lower = 0, upper = 1>[n_year] Ut; 
}
transformed parameters {
  real reca = cr / spro; 
  real recb = (cr - 1) / (ro*spro); 
  vector<lower=0>[n_year] abar;      // average age of population
  vector<lower=0>[n_year] yield;       
  vector<lower=0>[n_year] ssb; 
  vector<lower=0>[n_age] n = ninit; 
  
  for(t in 1:n_year){
    yield[t] = Ut[t]*sum(vul.*n.*w); 
    ssb[t] = sum(mwt.*n); 
    abar[t] = sum(ages.*n) / sum(n); 
    n = s*n.*(1-vul*Ut[t]);
    n[20] += n[19];                                  // plus group
    for(a in (n_age-1):2){n[a] = n[a - 1];}          // advance fish one age
    n[1] = reca*ssb[t] / (1+recb*ssb[t])*recmult[t]; // recruits
  }
}
model {
 // loop through years and add up yield or utility
 for(t in 1:n_year){
   if(obj_ctl == 0){ // yield
     target += sum(yield);
   }
   if(obj_ctl == 1){ // utility
     target += sum(yield^0.6);
   }
 }
}
