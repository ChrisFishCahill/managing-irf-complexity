/*
* BERTA_single_lake
* Bayesian Estimation of Recruitment Trends in Alberta
* Cahill 30 December 2021
*/
data {
  int <lower=0> n_surveys;             // number of surveys, i.e., rows of caa matrix
  int <lower=0> n_ages;                // number of ages 
  int <lower=0> n_obs;                 // n_ages*n_surveys 
  int <lower=0> n_years;               // number of years 
  int <lower=0> caa[n_surveys, n_ages];// catch at age matrix 
  vector[n_surveys] prop_aged;         // proportion critters aged in survey k,t
  vector[n_surveys] effort;            // survey effort in k,t
  int year[n_surveys];                 // year indicator 
  int <lower=0> ages[n_ages];          // helper vector
  vector[2] survey_yrs;                // years for which to calculate rbar
  int<lower=1> which_year;             // deliniate early vs. late
  real<lower=0> v_prior_early;         // prior, early F
  real<lower=0> v_prior_late;          // prior, late F
  vector<lower=0>[2] prior_sigma_v;    // prior variance for F's (early and late)
  real<lower=0> Ro_mean;               // prior mean, Ro
  real<lower=0> Ro_sd;                 // prior sd, Ro
  real<lower=0> ln_ar_sd;              // prior sd, ln_ar
  real<lower=0> prior_mean_w;          // prior sd
  real<lower=0> prior_sigma_w;         // prior sd
  real<lower=0> vbk;                   // lake specific vbk  
  real<lower=0> linf;                  // lake specific linf 
  real<lower=0> a50;                   // lake specific female a50
  real<lower=0> wl_beta;               // lake specific wl beta
  real<lower=0> lbar;                  // global average linf all lakes
  real<lower=0> M;                     // Instantaneous natural mortality
  real<lower=0> theta;                 // Lorenzen exponent
  real<lower=0> phi;                   // vulnerability parameter (nets)
  real<lower=0> psi;                   // vulnerability parameter (angling)
  vector[2] G_bound;                   // bounds for G
  real<lower=0> get_SSB_obs;           // logical get SSB observed
  real<lower=0> prior_sigma_G;         // sd of G
  int Rinit_ctl;                       // G*Ro (0) or a more complicated form (1)
  int<lower=0> length_Fseq;            // length of Fseq for generated quantities calcs
  vector[length_Fseq] Fseq;            // sequence of values to iterate across for Fmsy
  int<lower=0> rec_ctl;                // Ricker (0), BH (1)
  real<lower=0> cr_prior;              // compensation ratio prior 
  int<lower=0> vul_ctl;              // length2 (0) or logistic (1) fishing vul
}
transformed data {
  int<lower=0> caa_obs[n_obs];         // rowsums for survey k,t
  vector<lower=0>[n_ages] v_a;         // net vulnerability age a
  vector<lower=0>[n_ages] v_f_a;       // angling vulnerability age a
  vector<lower=0>[n_ages] l_a;         // length at age a 
  vector<lower=0>[n_ages] M_a;         // M at age a 
  vector<lower=0>[n_ages] f_a;         // fec at age a
  vector<lower=0>[n_ages] w_a;         // Weight at age a
  real<lower=0> sbro;                  // spawning biomass per recruit unfished 
  vector[n_surveys] SSB_Cn;            // SSB obs numerator
  vector[n_surveys] SSB_Cd;            // SSB obs denominator
  vector[n_surveys] SSB_C;             // SSB obs
  vector<lower=0>[n_ages] Lo;          // survivorship unfished (F=0) 
  real ln_ar_mean;                     // ln_ar mean
  int counter = 0;                     // counter for caa_obs
  real ah_vul = 3.5;                   // logistic selectivity parameter
  real sd_vul = 0.5;                   // logistic selectivity parameter
  // calculate vul, length-age, M-age, fec-age 
  sbro = 0; // initialize
  for(a in 1:n_ages){ 
    v_a[a] = ((linf/lbar)*(1 - exp(-vbk * ages[a])))^phi; 
    if(vul_ctl == 0){
      //v_f_a[a] = ((linf/lbar)*(1 - exp(-vbk * ages[a])))^psi; // used in paper - length^2
      v_f_a[a] = (1 - exp(-vbk * ages[a]))^psi;  // length^2 selectivity 
    }
    if(vul_ctl == 1){
      v_f_a[a] = 1 / (1 + exp(-(ages[a] - ah_vul) / sd_vul)); // logistic selectivity 
    
    }
    l_a[a] = (linf/lbar)*(1 - exp(-vbk * ages[a])); 
    M_a[a] = M/l_a[a]^theta;
    w_a[a] = 0.00001*(linf*(1 - exp(-vbk * ages[a])))^wl_beta;
    if(a < a50){
        f_a[a] = 0; 
      } else { 
        // relative weight at age assumed to follow vb
        f_a[a] = fmax(0, (l_a[a]^wl_beta)); 
      }
      if(a == 1){ 
        Lo[a] = 1;
      } else{
        Lo[a] = Lo[a-1]*exp(-M_a[a-1]); 
      }
      sbro += f_a[a]*Lo[a];
  }
  ln_ar_mean = log(cr_prior/sbro); 
  
  // calculate the rowsums for each survey, observed SSB
  // SSB_C(t)=sum over a of fec(a)*C(a,t)/[vage(a)Nnet(t)Paged(t)]
  for(i in 1:n_surveys){
    SSB_Cn[i] = 0; 
    SSB_Cd[i] = 0; 
    SSB_C[i] = 0; 
    for(a in 1:n_ages){
      counter = counter + 1; 
      SSB_Cn[i] += f_a[a]*caa[i, a]*(1/v_a[a]); 
      caa_obs[counter] = caa[i,a]; 
    }
    SSB_Cd[i] = prop_aged[i]*effort[i]; 
    SSB_C[i] = SSB_Cn[i] / SSB_Cd[i];
  }
}
parameters {
  vector<lower=0>[2] v;                              // early period F v[1] or late F v[2]
  real<lower=0> Ro;                                  // average unfished recruitment 
  real ln_ar;                                        // ln stock-recruit a
  vector[n_years-2] w;                               // recruitment anomalies--first 2 yrs for initiazation
  real<lower=G_bound[1], upper=G_bound[2]> G;        // is population at equilibrium (1) or declining (<1)
  //real<lower=0> phi;                               // negative binomial parameter
  //vector<lower=0, upper=1>[n_lakes] su_stock;      // estimate stocked fish survival
}
transformed parameters {
  vector<lower=0>[n_years] F_vec;                    // Insantaneous fishing mortality
  real Nat_array[n_ages, n_years];                   // Numbers at age  
  vector [n_years] SSB;                              // Spawning stock biomass
  vector<lower=0>[n_obs] caa_pred;                   // predicted catch at age
  real<lower=0> br;                                  // derived stock-recruit b from sbro, Ro
  vector<lower=0>[n_years] R2;                       // kick out recruits rather than N(a,t)
  vector[n_ages] Su_Fearly;                          // survivorship including fishing--early
  vector[n_ages] Su_Flate;                           // survivorship including fishing--late    
  real sbrf_early;                                   // spawning biomass per recruit fished
  real sbrf_late;                                    // spawning biomass per recruit fished
  real<lower=0> pinit;                               // how much depletion from Ro
  real<lower=0> Rinit;                               // initial recruitment
  vector[n_surveys] SSB_obs;                         // SSB observed
  vector[n_years] pred_N_catch;                      // predicted catch N/ha, *not vulN*
  vector[n_years] pred_B_catch;                      // predicted catch biomass/ha
  real<lower=0> SPR;                                 // spawning potential ratio
  real<lower=0> SSB_bar;                             // average ssb survey years
  real<lower=0> SBR;                                 // spawning biomass ratio
  real<lower=0> counter_SSB;                         // hack for SBR calcs
  real<lower=0> cr;                                  // compensation ratio kick out to check math
  real<lower=0> sbro_report;                         // report sbro  
  real ln_ar_mean_report;                            // report the ar mean 
  vector<lower=0>[n_ages] l_a_report;                // report length at age a 
  vector<lower=0>[n_ages] Lo_report;                 // report survivorship unfished (F=0) 
  vector<lower=0>[n_ages] v_a_report;                // report net vulnerability age a
  vector<lower=0>[n_ages] v_f_a_report;              // report angling vulnerability age a
  vector<lower=0>[n_ages] f_a_report;                // report fec at age a
  vector<lower=0>[n_ages] w_a_report;                // report weight at age a
  vector<lower=0>[n_ages] M_a_report;                // report M at age a
  
  // calculate sbrf
  sbro_report = sbro; 
  ln_ar_mean_report = ln_ar_mean; 
  sbrf_early = 0; 
  sbrf_late = 0; 
  for(a in 1:n_ages){
    if(a == 1){
      Su_Fearly[a] = 1; 
      Su_Flate[a] = 1; 
      } else {
        Su_Fearly[a] = Su_Fearly[a-1]*exp(-M_a[a-1] - v_f_a[a-1]*v[1]); 
        Su_Flate[a] = Su_Flate[a-1]*exp(-M_a[a-1] - v_f_a[a-1]*v[2]); 
      }
    sbrf_early += f_a[a]*Su_Fearly[a]; 
    sbrf_late += f_a[a]*Su_Flate[a]; 
      
    // fill out report vectors
    l_a_report[a] = l_a[a]; 
    Lo_report[a] = Lo[a];
    v_a_report[a] = v_a[a]; 
    v_f_a_report[a] = v_f_a[a]; 
    f_a_report[a] = f_a[a]; 
    w_a_report[a] = w_a[a]; 
    M_a_report[a] = M_a[a]; 
  }
  SPR = sbrf_late / sbro; 

  // Calculate recruitment b's, Rinit's
  if(rec_ctl==0){ // ricker b
    br = (ln_ar + log(sbro)) / (Ro*sbro); 
  }
  if(rec_ctl==1){ // bev-holt b
    br = (exp(ln_ar)*sbro - 1)/ (Ro*sbro); 
  }
  if(Rinit_ctl == 0){
    Rinit = G*Ro; 
  }
  if(Rinit_ctl == 1){
    if(rec_ctl==0){ // ricker Req
      Rinit = (ln_ar + log(sbrf_early)-log(G)) / (br*sbrf_early);
    }
    if(rec_ctl==1){ // bev-holt Req
      Rinit = (exp(ln_ar)*(sbrf_early-1)) / (br*sbrf_early); 
    }
  }
  pinit = Rinit / Ro; 

  // Initialize F(t) fishing rate vector from start year to 2018
  for(t in 1:n_years){
    if(t < which_year){
        F_vec[t] = v[1];
      }
      if(t >= which_year){
        F_vec[t] = v[2];
      }
      SSB[t] = 0;  
      pred_N_catch[t] = 0; 
      pred_B_catch[t] = 0;
  }
  
  // Initialize the N(at) array age structure for t = 1,2
  cr = exp(ln_ar)*sbro; 
  for(t in 1:2){
    Nat_array[1, t] =  Rinit;
    if(t == 1){
      for(a in 2:n_ages){
        Nat_array[a, t] = Nat_array[a-1, t]*exp(-M_a[a-1] - v_f_a[a-1]*F_vec[1]); 
        SSB[t] += Nat_array[a ,t]*f_a[a];
        pred_N_catch[t] += Nat_array[a ,t]*v_f_a[a];  
        pred_B_catch[t] += Nat_array[a ,t]*v_f_a[a]*w_a[a];  
        }
        pred_N_catch[t] = pred_N_catch[t]*(1-exp(-F_vec[1]));
        pred_B_catch[t] = pred_B_catch[t]*(1-exp(-F_vec[1]));
      }
      if(t == 2){
        for(a in 2:n_ages){
          Nat_array[a, t] = Nat_array[a, 1]; 
        }
      }
    SSB[t] = SSB[1]; 
    pred_N_catch[t] = pred_N_catch[1]; 
    pred_B_catch[t] = pred_B_catch[1]; 
  }
  
  // Calculate the N(at) array and derived outputs
  R2[1] = Nat_array[1, 1];
  R2[2] = Nat_array[1, 2];
    
  SSB_bar = 0;
  counter_SSB = 0; 

  for(t in 3:n_years){
    if(rec_ctl==0){// ricker
      Nat_array[1, t] = SSB[t-2]*exp(ln_ar - br*SSB[t-2] + w[t-2]);
    }
    if(rec_ctl==1){// beverton-holt
      Nat_array[1, t] = SSB[t-2]*exp(ln_ar + w[t-2]) / (1 + br*SSB[t-2]);
    }
    for(a in 2:n_ages){
      Nat_array[a, t] = Nat_array[a-1, t-1]*exp(-M_a[a-1] - v_f_a[a-1]*F_vec[t-1]);  
      SSB[t] += Nat_array[a ,t]*f_a[a];
      pred_N_catch[t] += Nat_array[a ,t]*v_f_a[a]; 
      pred_B_catch[t] += Nat_array[a ,t]*v_f_a[a]*w_a[a]; 
    }
    pred_N_catch[t] = pred_N_catch[t]*(1-exp(-F_vec[t]));
    pred_B_catch[t] = pred_B_catch[t]*(1-exp(-F_vec[t]));
    R2[t] = Nat_array[1, t];
    // calculate mean SSB across survey years
    if(t >= survey_yrs[1] && t <= survey_yrs[2]){
      counter_SSB += 1; 
      SSB_bar += SSB[t];
     }
  }
  SSB_bar = SSB_bar / counter_SSB;
  SBR = SSB_bar / (Ro*sbro);     
  
  // Calculate the preds vector
  // C(k,a,t)=N(a,t)*Nnet(t)Paged(t)*v_a(a) 
  for(hack in 1:1){// hack to initialize j = 0
    int j = 0; 
    for(i in 1:n_surveys){
      for(a in 1:n_ages){
        j += 1;
        caa_pred[j] = 0;                        // initialize
        caa_pred[j] = Nat_array[a, year[i]]*    // numbers(a,t)
        prop_aged[i]*                           // prop_aged
        effort[i]*                              // survey effort
        v_a[a];                                 // vulnerability to survey gear
      }
      if(get_SSB_obs==1){
        SSB_obs[i] = SSB_C[i];  
      }
    }
  }
}
model {
  // priors:
  v[1] ~ normal(v_prior_early, prior_sigma_v[1]); 
  v[2] ~ normal(v_prior_late, prior_sigma_v[2]); 
  Ro ~ lognormal(Ro_mean, Ro_sd); 
  ln_ar ~ normal(ln_ar_mean, ln_ar_sd);
  G ~ normal(0,prior_sigma_G); 
  to_vector(w) ~ normal(prior_mean_w, prior_sigma_w); 
  
  // likelihood
  caa_obs ~ poisson(caa_pred); 
  
  // NOTE: Stan is creating log(posterior) = log(likelikhood) + log(priors)
  
  // extra stuff for negative binomial model and stocking survival 
  // phi ~ cauchy(0,3);
  // su_stock ~ beta(2,2); 
  // caa_obs ~ neg_binomial_2(caa_pred, phi); 
}
generated quantities{
  real Fmsy;           // instantaneous fishing mortality @ MSY                
  real MSY;            // weight yield kg/ha
  real F_ratio;        // Flate / F_msy  
  real F_early_ratio;  // Fearly / F_msy                        
  real b_ratio;        // average ssb survey years / pristine ssb
  
  // Fmsy, MSY subroutine
  Fmsy = 0; 
  MSY = 0; 
  for(i in 1:length_Fseq){
    real sbrf = 0; 
    real ypr = 0; 
    real su = 1; 
    real Req = 0; 
    real Yeq = 0; 
    for(a in 1:n_ages){
      sbrf += su*f_a[a]; // accumulate spawning biomass per recruit
      ypr += su*(1-exp(-Fseq[i]*v_f_a[a]))*w_a[a]; 
      su = su*exp(-M_a[a] - Fseq[i]*v_f_a[a]); 
      }
      if(rec_ctl == 0){ // ricker
        Req = log( exp(ln_ar)*(sbrf) ) / (br*sbrf); // Botsford predction of Req for F[i]
      }
      if(rec_ctl == 1){ // bh
        Req = ( exp(ln_ar)*sbrf-1.0 ) / (br*sbrf); 
      }
      Yeq = Req*ypr; // predicted equilibrium yield
      if(Yeq > MSY){
        MSY = Yeq; 
        // jitter values to make a purdy histogram:
        Fmsy = Fseq[i] + 0.01*(uniform_rng(0,1)-0.5); 
      } else {
        // Yeq for this F is lower than highest value already found,
        // so can exit the subroutine
        continue; 
      }
   }
   // Kobe plot calculations 
   F_early_ratio = v[1] / Fmsy;
   F_ratio = v[2] / Fmsy;
   b_ratio = SSB_bar / (Ro*sbro);
}
