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
%include "cm_binary.sas";
%include "cm_binary_probs.sas";
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
%include "../../est_probs.sas";


* other parameters;
%let n=1000;
%let clustsize=5;
%let m=2;
%let ntrials = 1000;
%let sampmethod= CM;
libname reslog '.';

* parameters input;
%include "parms_binary.sas";

%let z0=0;
%let z1=1;

options nomprint;
%cm_binary_probs;

/*
%let maxsetsize=20;
%let maxcaseset=20;
%let n=500;
*/
/*

ods listing file='probabilities.lst';
proc sort data=reslog.probs out=probs; by or descending clustsize; run;
proc print data=probs;
  var meancases clustsize or phat0 inphat0CI inCI0 phat1 inphat1CI inCI1;
  format phat0 inphat0CI inCI0 phat1 inphat1CI inCI1 percent7.1 ;
run;
ods listing;

*/
