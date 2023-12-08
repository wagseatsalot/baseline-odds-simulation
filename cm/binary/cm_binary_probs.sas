%macro cm_binary_probs;
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
* lam	bda0 - baseline odds;
  call symput('lambda0',trim(left(put(lambda0,12.3))));

* Prob of C=1;
  call symput('pz',trim(left(put(pz,12.3))));

* Prob of Z=1;
  call symput('pc',trim(left(put(pc,12.3))));

* Correlation between the c and z;
  call symput('theta',trim(left(put(theta,12.3))));

* Z Rate ratio;
  call symput('or_z',trim(left(put(or_z,12.3))));

* C Rate ratio ;
  call symput('or_c',trim(left(put(or_c,12.3))));


* calculate pi parameterization and put into macro variables;

* pi11 is the solution to a quadratic equation;

   a = theta -1;
   b = -(1 + (pz + pc)*(theta-1));
   c = theta*pz*pc;

   if a = 0 then pi11 = pz*pc;
   else do;
* solve quadratic;

     disc = b*b - 4*a*c;
     if disc lt 0 then do; put 'This combo of parameters doesn''t work, discriminant <0 ' 
      theta= pz= pc= disc=;  abort; end;

     plusroot = (-b + sqrt(disc))/(2*a);
     minusroot = (-b - sqrt(disc))/(2*a);
     if 0 le plusroot le 1 then pi11 = plusroot;
     if 0 le minusroot le 1 then pi11 = minusroot;
     if not ((0 le plusroot le 1) or (0 le minusroot le 1)) then do;
     put 'Neither solution is between zero and 1 ' theta= pz= pc= plusroot= minusroot=;
	 abort;
	 end;  
   end; 
* valid pi11 so compute others;
   pi10 = pz - pi11;
   pi01 = pc - pi11;
   pi00 = theta*pi01*pi10/pi11;

* output into macro variables;
  call symput('pi11',trim(left(put(pi11,12.3))));
  call symput('pi10',trim(left(put(pi10,12.3))));
  call symput('pi01',trim(left(put(pi01,12.3))));
  call symput('pi00',trim(left(put(pi00,12.3))));

stop;
run;
  %put pi00=&pi00  pi10=&pi10  pi01=&pi01  pi11=&pi11;
  
* Create simulation parameter variable;
%let simparms =&lambda0 &pz &pc &theta &or_z &or_c;


* dist of z in cluster depends on value of c;
data cohort (keep=trial cluster clustval c z d rand);
  retain lambda0 &lambda0 or_z &or_z pi00 &pi00  
      pi10 &pi10  pi01 &pi01  pi11 &pi11;

 * counter-matching variable set to have c-z correlation;

* set sensitivity and specificity equal to .9;
if _n_ eq 1 then do;
  pc1 = .9;  * if z=1, then p(c=1) = pc1;
  pc0 = .1; * if z=0 , then p(c=1) = pc0;
  end;

  do trial = 1 to &ntrials;
  do cluster = 1 to &n/&clustsize;
     clustval = ranuni(0) le &pc;
   do i = 1 to &clustsize;
     if clustval eq 0 then z = ranuni(0) le pi10/(pi10 + pi00); 
     else z = ranuni(0) le pi11/(pi11 + pi01);
* counter-matching variable;
     if z eq 1 then c = ranuni(0) le pc1;
     else c = ranuni(0) le pc0;
     lambda = lambda0*or_z**z;
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
