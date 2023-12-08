************************************************************
* make_counter_matched.sas - macro to arrange data
*   with one line per counter-matched case-control set 
*   for CLR analysis using NLMIXED
*   lines of data have cases first, with samstr indicator
************************************************************
* MACRO INPUT VARIABLES
* Required:
*  indsn_rs - input risk set or case-control data set file
*  setno - set number variable
*  samstr - sampling stratum (consecutive integers)
*  low - lowest samstr value
*  high - high samstr value 
*  totalss - number in cohort in sampling stratum
*  cc - case-control variable cc=1 case, cc=0 control
*  vlist - variables to output, separated by spaces 
* Optional:
*  setvars - set level variables to be included in output data
*  multcases - yes = multiple cases per set. Numbers o
*   cases and total in set are included in the output data
*      (default=no)
*  outdsn - name of output data set (default=case_cont)
************************************************************
* OUTPUT DATA SET
* setno - set number
* setvars - set level variables
* _ncases - total number of cases (multicase = yes)
* _ntot - total number in set (multicase = yes)
* _ntotsslow-_ntotsshigh - cohort number in each sampling strata
* _msamplow-_msamphigh - number sampled from each sampling strata
* _ncaseslow-_ncaseshigh - number of cases in each sampling strata
* _z- covariates starting with var1 over all subjects, then var2 over all subjects, etc
* _cc - case-control status indicators
* _samstr - sampling stratum for each subject
************************************************************
* GLOBAL MACRO VARIABLES 
* declared global so it can be used in array sizing in analysis
* maxsetsize - size of largest case-control set
* maxcaseset - size of largest case set
************************************************************ 
* Information needed about the data, input data set name,
* variables, and the maximum set size are included in
* the data set label.
************************************************************
;
%macro make_counter_matched(indsn_rs,setno,cc,samstr,low,high,totalss,vlist,setvars=,multcases=yes,outdsn=case_cont);

%global maxsetsize; %global maxcaseset;

%* get the number of variables in vlist;

%let nvars = 0;
%do %until (%scan(&vlist,%eval(&nvars+1)) lt 0);
  %let nvars = %eval(&nvars+1);
%end; 

* get maximum set size and put it in a macro variable;
proc freq data=&indsn_rs;
  table &setno /noprint out=_ntot (rename=(count=_ntot) drop=percent);
  run;
proc sql noprint;
  select max(_ntot)
    into :maxsetsize
    from _ntot;
quit;

* Get maximum number of cases and put into a macro variable;
proc freq data=&indsn_rs;
  table &setno /noprint out=_totcases (rename=(count=_totcases) drop=percent);
  where &cc;
  run;
proc sql noprint;
  select max(_totcases)
    into :maxcaseset
    from _totcases;
quit;

* get rid of leading blanks;
%let maxsetsize= &maxsetsize;
%let maxcaseset= &maxcaseset;

* Get number of cases in each sampling stratum;
proc freq data=&indsn_rs;
  table &setno * &samstr /noprint out=_ncases (rename=(count=_cases) drop=percent) sparse;
  where &cc;
  run;
* Get number sampled from each sampling stratum;
proc freq data=&indsn_rs;
  table &setno * &samstr /noprint out=_msamp (rename=(count=_samp) drop=percent) sparse;
  run;

%do isam = &low %to &high;
data _ncases&isam;
  merge _ncases _msamp;
  by &setno &samstr;
  _ncases&isam = _cases;
  _msamp&isam = _samp;
  where &samstr eq &isam;
  drop _cases _samp;
  run;
%end;

* get rid of leading blank;
%let maxsetsize= &maxsetsize;

%* size of array to hold covariate values for case-control set;
%let arraysize=%eval(&nvars*&maxsetsize);

%* sort to get cases first;
proc sort data=&indsn_rs out=_rsdata; by &setno &samstr descending &cc; run;

%* for multiple case case-control sets, add in number in set and
%* number of cases;

%let ncases=;
%let ntot=;
%let nsamst=%eval(&high-&low+1);
%if %upcase(multcases) ne NO %then %do;
%let ncases=_ncases;
%let ntot=_ntot;
%end;


data &outdsn (keep= &setno &setvars _totcases &ntot _ncases&low-_ncases&high _msamp&low-_msamp&high
_ntotss&low-_ntotss&high _z1-_z&arraysize _cc1-_cc&maxsetsize
label="input data set: &indsn_rs, variable list: &vlist, 
maximum set size: &maxsetsize, maximum case set size: &maxcaseset");
merge %do isam=&low %to &high; &ncases&isam %end; _totcases &ntot _rsdata;
by &setno;
array vars{&nvars} &vlist;
retain nsamst &nsamst _ntotss&low-_ntotss&high;
array _ntotss{&low:&high} _ntotss&low-_ntotss&high;
array _z{&nvars,&maxsetsize} ;
array _samstr{&maxsetsize} _samstr1-_samstr&maxsetsize;
array _cc{&maxsetsize} _cc1-_cc&maxsetsize;
retain _i _z1-_z&arraysize  _ntotss&low-_ntotss&high _cc1-_cc&maxsetsize 
  _samstr1-_samstr&maxsetsize ;

* initialize covariate arrays to zero;
 if first.&setno then do;
   _i = 0;
   do _ii=1 to &maxsetsize;
    do _ivar = 1 to &nvars;
    _z(_ivar,_ii)=.; 
    end;
    _cc{_ii} = .;
    _samstr{_ii} = .;
   end;
* initialize sampling strata;
   do _ii = &low to &high;
     _ntotss{_ii} = 0;
     end;
  end;
  
* get total number in the sampling strata;
  _ntotss{&samstr} = &totalss;
 _i = _i + 1;
  
* sampling stratum for each subject;
  _samstr{_i} = &samstr;
  _cc{_i} = &cc;
* put covariates for each subject into z array;
* loop through all variables ;
  do _ivar = 1 to &nvars;
   _z{_ivar,_i} = vars{_ivar}; 
   end;
  
 label _ntot='_ntot: number of case-control subjects';
  label _totcases='_totcases: total number of cases';
%do istr = &low %to &high;
  label _ncases&istr = "_ncases&istr: cases in &samstr stratum &istr";
  label _ntotss&istr = "_ntotss&istr: total in &samstr stratum &istr";
  label _msamp&istr = "_msamp&istr: sampled in &samstr stratum &istr";
  %end;

* when through the data for this cm set, output;
 if last.&setno then output;
 run;

%put input data set: &indsn_rs;
%put variable list: &vlist; 
%put maximum set size: &maxsetsize;
%put maximum case set size: &maxcaseset;
%put output data set: &outdsn;
%put Above information is in the output data set label;
* clean up;

proc datasets;
  delete %do isam=&low %to &high; &ncases&isam %end; _totcases _ntot _rsdata;
  quit;

%mend;
