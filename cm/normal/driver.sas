*****************************
*****************************
* driver program.sas
******************************;
* computes asymptotic quantities for variance calculations
* two-stage: cluster samples with replacement then 1:1 cm sample
* of controls from clusters;
******************************;
ods html close;
ods listing;
proc datasets kill; quit;
%include "../../deleteallmacrovars.sas";
%deleteAllMacroVars();
%include "cm_normal.sas";
%include "../../cohort_analysis.sas";
%include "../make_cm.sas";
%include "../cm_analysis.sas";
%include "../make_counter_matched.sas";
%include "../../co_est.sas";
%include "../cm_co_var.sas";
%include "../cm_rb.sas";
%include "../../efficiency.sas";
%include "../cm_shift_survey.sas";
%include "../../output_results.sas";


* other parameters;
%let n=1000;
%let clustsize=5;
%let m=2;
%let ntrials = 1000;
%let sampmethod= CM;
libname reslog '.';

* parameters input;
%include "../../parms_normal.sas";

options nomprint;
%cm_normal;

/*
%let maxsetsize=20;
%let maxcaseset=20;
%let n=500;
*/
/*
ods html close;
ods listing file='case_only.lst';
proc sort data=reslog.case_only out=co; by or  descending clustsize descending meancases; run;
proc print data=co;
  where meancases lt 40 and clustsize ne 1000;
  var meancases or clustsize logodds vlogodds var_a_R cov_ab cov_ab_R;
  format logodds bhat vbhat evbhat 7.2 vlogodds var_a_R cov_ab cov_ab_R 7.3;
run;
ods listing;
*/


