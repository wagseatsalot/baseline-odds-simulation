%macro shift_survey_cm;
%* Inputs: cm, case_cont;
%* Outputs: uncond, survey;

data cm2 (keep=trial d offset z surveyw);
  merge cm (in=in_cm) case_cont;
  by trial cluster;
  if not in_cm then delete;
  retain w;
  w0 = .00001;
  w1 = .00001;
  if (_msamp0-_ncases0 gt 0) then w0 = (_ntotss0-_ncases0)/(_msamp0-_ncases0);
  if (_msamp1-_ncases1 ne 0) then w1 = (_ntotss1-_ncases1)/(_msamp1-_ncases1);
  if c eq 0 then w = w0; else w = w1; 
    offset = log(w);
	surveyw = d + (1-d)*w;  * for use with survey estimator;
  run;
proc printto log=nul; run;
proc logist data=cm2 descending noprint outest=ests covout noprint;
  by trial;
  model d= z / offset=offset;
  run;
proc printto log=log; run;
title2 'Shifted baseline odds parameter';
data estparms;
  set ests; where _TYPE_ eq 'PARMS';
  bhat = z;
  run;

proc means data=estparms n mean var min max;
  var intercept bhat;
  output out=uncond mean(intercept)=logodds var(intercept)=vlogodds
      mean(bhat)=bhat var(bhat)=vbhat;
  run;
proc datasets; delete estparms ests; quit;

proc printto log=nul;
proc surveylogistic data=cm2;
  by trial;
  model d (descending) = z;
  weight surveyw;
  ods output parameterestimates=ests;
  run;
proc printto; run;
proc printto log=log; run;
data interc (keep=trial intercept);
  set ests; where variable eq 'Intercept';
  intercept = estimate;
  run;
data betas (keep=trial bhat);
  set ests; where variable eq 'z';
  bhat = estimate;
  run;
data estparms;
   merge interc betas;
   by trial;
   run;

title2 'Survey estimator';
proc means data=estparms n mean var min max;
  var intercept bhat;
  output out=survey mean(intercept)=logodds var(intercept)=vlogodds
              mean(bhat)=bhat var(bhat)=vbhat;
  run;
proc datasets; delete ests intercept bhat; quit;
%mend;