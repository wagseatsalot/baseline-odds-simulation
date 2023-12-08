%macro output_results;
%* Inputs: simulation results from each method;
%* Outputs: reslog.case_only - accumulates case only simulation info;
%*          reslog.results - accumulates info from all methods;
%*          sampmethod comparisons.lst - txt output to current directory;

data allests;
  set case_only (keep=logodds vlogodds bhat vbhat in=in_co)
	D_rb  (keep=logodds vlogodds bhat vbhat in=in_D_rb)
        uncond(keep=logodds vlogodds bhat vbhat in=in_uncond)
        survey(keep=logodds vlogodds bhat vbhat in=in_survey)
        cohortests(keep=logodds vlogodds bhat vbhat in=in_cohort)
	co_cohbeta(keep=logodds vlogodds bhat vbhat in=in_co_cohbeta)
	cohsamp(keep=logodds vlogodds bhat vbhat in=in_cohsamp)
;
  retain lambda0 &lambda0 or &or_z n &n clustsize &clustsize m &m 
          meancases &meancases maxsetsize &maxsetsize ntrials &ntrials;
  length method $ 30 simparms $ 50;


  simparms = trim(left("&simparms"));
/*
  lambda0 = input("&lambda0",6.2);
  or      = input("&or_z",5.1); 
  n       = input("&n",7.);
  clustsize= input("&clustsize",7.); 
  m=        input("&m",4.);
  meancases = input("&meancases",7.1);
  maxsetsize=input("&maxsetsize",7.1);
  ntrials=  input("&ntrials",7.1);
*/
  if in_co then method = 'Case only';
  if in_D_rb then method = 'Rao-Blackwell R-D conditional';
  if in_uncond then method='Uncond logistic with weights';
  if in_survey then method='Survey logistic';
  if in_cohort then method='Full cohort';
  if in_co_cohbeta then method='CO with cohort beta';
  if in_cohsamp then method='Coh with samp beta';
label logodds='Mean log baseline odds';
label vlogodds='Empirical variance of log odds' ;
label bhat='Mean beta hat';
label vbhat='Empirical variance of beta hat';
label method='Method';
run;

proc append force data=allests base=reslog.results; run;

data case_only;
  set case_only;
  retain lambda0 &lambda0 or &or_z n &n clustsize &clustsize m &m 
          meancases &meancases maxsetsize &maxsetsize ntrials &ntrials;
  length method $ 12 simparms $ 50;
  method = 'Case only';
  simparms = trim(left("&simparms"));
/*  lambda0 = input("&lambda0",6.2);
  or      = input("&or_z",5.1); 
  n       = input("&n",7.);
  clustsize= input("&clustsize",7.); 
  m=        input("&m",4.);
  meancases = input("&meancases",7.1);
  maxsetsize=input("&maxsetsize",7.1);
  ntrials=  input("&ntrials",7.1);
*/
run;

proc append force data=case_only base=reslog.case_only; run;

proc printto; run;
proc printto file="&sampmethod comparison.lst";
title2 "lambda0=&lambda0 or=&or_z n=&n clustsize=&clustsize m=&m meancases=&meancases 
    maxsetsize=&maxsetsize ntrials=&ntrials";
proc print data=allests label;
  var method logodds vlogodds bhat vbhat;
  format logodds vlogodds bhat vbhat 7.3;
run;
proc printto; run;
title;

%mend;