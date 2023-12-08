%macro cm_normal_probs;
options noxwait;


***************************************************
* Read in parameters for a given run
* from work.parms
***************************************************
 ;
proc sql noprint;
  select count(*)
    into :nruns
    from parms;
quit;

%put There will be &nruns runs;

%do run_nm = 1 %to &nruns;
data _null_;
  
  run_nm = &run_nm;
  set parms point=run_nm ;
* lambda0 - baseline odds;
  call symput('lambda0',trim(left(put(lambda0,12.3))));

* Correlation between the c and z;
  call symput('rho2',trim(left(put(rho2,12.3))));

* Z Rate ratio;
  call symput('or_z',trim(left(put(or_z,12.3))));
  stop;

run;
  
* Create simulation parameter variable;
%let simparms = &rho2 &lambda0 &or_z;

options noxwait;

* counter-matching variable c indicates if z >=0;
data cohort (keep=trial cluster c z d rand);
  retain lambda0 &lambda0 or_z &or_z;
* or_z is the or over upper 5th percentile to lower 5th percentile;
  or_z_sd =  exp(log(or_z)/4);

  do trial = 1 to &ntrials;
  do cluster = 1 to &n/&clustsize;
     clustmean = rand('normal')*sqrt(&rho2);
   do i = 1 to &clustsize;
     z = clustmean + rand('normal')*sqrt(1-&rho2);
	 c = z ge 0;
     lambda = lambda0*or_z_sd**z;
     p = lambda/(1+lambda);
     d = ranuni(0) le p;
     rand = ranuni(0);
    output;
    end;
    end;
    end;
  length d 3;
run;

proc summary data=cohort;
  where c eq 0;
  var d cluster;
  output out=cohcstr0 (drop=_freq_ _type_) sum(d)=_ncases0 n(cluster)=_ntotss0;
  by trial cluster;
  run;
proc summary data=cohort;
  where c eq 1;
  var d cluster;
  output out=cohcstr1 (drop=_freq_ _type_) sum(d)=_ncases1 n(cluster)=_ntotss1;
  by trial cluster;
  run;
data clustinfo (keep=trial cluster triclust _ncases0 _ncases1 _ntotss0 _ntotss1  _totcases _ntot);
  merge cohcstr0 (in=in0) cohcstr1 (in=in1);
  by trial cluster;
  triclust = trial + cluster/1000;
  if not in0 then do; _ncases0 = 0; _ntotss0 = 0; end;
  if not in1 then do; _ncases1 = 0; _ntotss1 = 0; end;
  _ntot = _ntotss0+_ntotss1;
  _totcases = _ncases0+ _ncases1;
  run;

proc summary data=cohort; 
  by trial;
  output out=cases sum(d)=totcases;
  run;

* put mean number of cases in a macro variable;
proc sql noprint;
  select mean(totcases)
    into :meancases
    from cases;
quit;
proc datasets; delete cases; quit;


/*
%cohort_analysis;
* Input: cohort;
* Output: betas_cohort, cohortests;
*/
%make_cm;
* Input: cohort;
* Output: cm, betas_cl;

%co_est(cm);
* Input: cm, betas_cl;
* Output: cotrial, case_only;

%cm_co_var;
* Inputs: cotrial, case_cont case-only;
* Output: cm_co_var case-only (adds var_a and cov_ab stats);

%est_probs;

* clean up a few things;
/* proc datasets; delete cm_co_var cmco logodds; quit; */

* add simulation info and output probs results;
data probs;
  set probs;
  retain lambda0 &lambda0 or &or_z n &n clustsize &clustsize m &m 
          meancases &meancases maxsetsize &maxsetsize ntrials &ntrials;
  length method $ 12 simparms $ 50;
  method = 'Case only';
  simparms = trim(left("&simparms"));
run;

proc append force data=probs base=reslog.probs; run;

 
%end;
%mend;
