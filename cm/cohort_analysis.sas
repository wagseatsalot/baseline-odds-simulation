%macro cohort_analysis;
%* betas_cohort - cohort mles for each trial;
%* cohortests - simulation results;
title2 'Cohort mle';
 proc printto log=nul; run;
proc logist data=cohort descending outest=ests1 covout noprint;
  by trial;
  model d= z;
  run;
proc printto log=log; run;

data betas_cohort (keep=trial bhat vbhat logodds cov_ab_mle);
  set ests1;
  by trial; 
  retain bhat logodds cov_ab_mle;
  if _TYPE_ eq 'PARMS' then do; bhat = z; logodds= intercept; delete; end;
  if _TYPE_ eq 'COV' and _NAME_ eq 'Intercept' then do; cov_ab_mle = z; output ; end;
  if _type_ eq 'COV' and _name_ eq 'z' then vbhat = z;
  run;
proc means data=betas_cohort n mean var min max;
 where vbhat lt 10;
 var logodds bhat cov_ab_mle;
  output out=cohortests mean(logodds)=logodds var(logodds)=vlogodds
    mean(bhat)=bhat var(bhat)=vbhat mean(cov_ab_mle)=cov_ab 
   var(cov_ab_mle)=vcov_ab;
  run;

proc datasets; delete betas ests1; quit;
title2;
%mend;
