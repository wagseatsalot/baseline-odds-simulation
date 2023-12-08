%macro cm_normal;
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

%cohort_analysis;
* Input: cohort;
* Output: betas_cohort, cohortests;

%make_cm;
* Input: cohort;
* Output: cm, betas_cl;

%co_est(cm);
* Input: cm, betas_cl;
* Output: cotrial, case_only;

%cm_co_var;
* Inputs: cotrial, case_cont case-only;
* Output: case-only (adds var_a and cov_ab stats);

* clean up a few things;
proc datasets; delete cmco logodds; quit;

%cm_rb; * rao-blackwellized alpha estimate;
* Inputs: betas_cl case_cont
* Output: D_rb (simulation summary of rb estimator);


%coh_a_cc_b; * use mle for alpha with case-control beta;
* Inputs: cohort betas_cl;
* Outputs: cohsamp;

%co_cohbeta; * use mle beta hat in case-only alpha;
* Inputs: cohort betas_cohort;
* Outputs: co_cohbeta;

%shift_survey_cm;
* Inputs: cm, case_cont;
* Outputs: uncond, survey;

%output_results;
* Inputs: simulation results from each method;
* Outputs: reslog.case_only - accumulates case only simulation info;
*          reslog.results - accumulates info from all methods;
*          sampmethod comparisons.lst - txt output to current directory.;
 
%end;
%mend;
