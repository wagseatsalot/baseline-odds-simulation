%macro make_cm;
%* 1:1 counter-match sample (within cluster strata);
%* Input: cohort;
%* Output: cm, betas_cl;

proc sort data=cohort; by trial cluster descending d rand; run;
data cm (keep=trial cluster triclust d c z totalss);
  merge cohort clustinfo;
  by trial cluster;
  array _ntotss[0:1] _ntotss0 _ntotss1;
  array _msamp[0:1] _msamp0 _msamp1;
  array _ii[0:1] _temporary_;
  retain _msamp0 _msamp1 _ntotss0 _ntotss1;

  if first.trial or first.cluster then do;
* goal is to get 2*_totcases in the sample, even if there are fewer
* than _totcases in the stratum;
     dif0 =  _ntotss0 - _totcases; dif1 = _ntotss1 - _totcases; 
    _msamp0 = min(_ntotss0,_totcases-min(0,dif1)); 
    _msamp1 = min(_ntotss1,_totcases-min(0,dif0));
    _ii[0] = 0; _ii[1] = 0; 
    end;
   _ii[c] = _ii[c] + 1; 
   totalss = _ntotss[c];
    if _ii[c] le _msamp[c] then output; else; delete;
 run;

%make_counter_matched(cm,triclust,d,c,0,1,totalss,z,setvars=trial cluster,outdsn=case_cont)

ods select parameterestimates;
ods output parameterestimates=temp;
%cm_analysis ;
ods output;
ods select all;

data betas_cl; 
  set temp;
  by trial;
  retain bhat vbhat;
  if parameter eq 'beta' then bhat=estimate;
  if parameter eq 'beta' then vbhat=standarderror**2;
  if LAST.TRIAL then output; else delete;
  keep trial bhat vbhat;
run;
 proc means data=betas_cl n mean var min max; 
   var bhat vbhat;
   run;
 proc datasets; delete temp; quit; 
%mend;