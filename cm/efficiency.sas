%macro coh_a_cc_b;
%* mle for alpha with case-control beta and with true beta;
%* Inputs: cohort betas_cl;
%* Outputs: cohsamp cohbeta0

* check efficiency;
* bhat used for zbeta is from the case-control sample;

  data coh1;
  merge cohort betas_cl;
  by trial;
  zbeta = z*bhat;
  zbeta0 = z*log(&or_z);
  run;

title2 'Cohort logodds with sampled beta hat';
 proc printto log=nul; run;
proc logist data=coh1 descending outest=ests1 covout noprint;
  by trial;
  model d= / offset=zbeta;
  run;
proc printto log=log; run;

data estparms;
  set ests1; where _TYPE_ eq 'PARMS';
  run;
data estparms;
  merge estparms betas_cl (keep=trial bhat vbhat);
  by trial;
  run;

proc summary data=estparms n mean var min max;
  where vbhat lt 10;
  output out=cohsamp mean(intercept)=logodds var(intercept)=vlogodds
  mean(bhat)=bhat var(bhat)=vbhat;
  run;
proc datasets; delete ests1 estparms; quit;
title2;

title2 'Cohort logodds with beta_0';
 proc printto log=nul; run;
proc logist data=coh1 descending outest=ests1 covout noprint;
  by trial;
  model d= / offset=zbeta0;
  run;
proc printto log=log; run;

data estparms;
  set ests1; where _TYPE_ eq 'PARMS';
  run;
data estparms;
  set estparms;
  bhat = log(&or_z);
  vbhat = 0;
  by trial;
  run;

proc summary data=estparms n mean var min max;
  where vbhat lt 10;
  output out=cohbeta0 mean(intercept)=logodds var(intercept)=vlogodds
  mean(bhat)=bhat var(bhat)=vbhat;
  run;
proc datasets; delete ests1 estparms; quit;
%mend;

%macro co_cohbeta;
%*co using cohort beta - to gage alpha variability due to case-control beta;
%*Inputs: cohort betas_cohort;
%*Outputs: co_cohbeta co_beta0

* Use cohort beta in case-only estimator;
data co;
  merge cohort betas_cohort ;
  by trial;
  invphi = d*exp(-z*bhat);
 run;

proc summary data=co; 
  by trial;
  output out=cocluster sum(invphi)= suminvphi sum(d)=totcases;
  run;
data logodds;
  merge cocluster betas_cohort;
  logodds = log(suminvphi / (&n-totcases));
  run;

title2 'Case-only estimator with cohort beta hat';
 proc summary data=logodds n mean var min max;
 where vbhat lt 10;
  var logodds bhat;
  output out=co_cohbeta mean(logodds)=logodds var(logodds)=vlogodds
        mean(bhat)=bhat var(bhat)=vbhat;
  run;
  proc datasets; delete co cocluster logodds; quit; 
 title2;

* Use beta_0 in case-only estimator;
data co;
  merge cohort betas_cohort ;
  by trial;
  invphi = d/(&or_z**z);
 run;

proc summary data=co; 
  by trial;
  output out=cocluster sum(invphi)= suminvphi sum(d)=totcases;
  run;
data logodds;
  merge cocluster;
  logodds = log(suminvphi / (&n-totcases));
  bhat = log(&or_z);
  vbhat = 0;
  run;

title2 'Case-only estimator with beta_0';
 proc summary data=logodds n mean var min max;
 where vbhat lt 10;
  var logodds bhat;
  output out=co_beta0 mean(logodds)=logodds var(logodds)=vlogodds
        mean(bhat)=bhat var(bhat)=vbhat;
  run;
  proc datasets; delete co cocluster logodds; quit; 
 title2;

%mend;

