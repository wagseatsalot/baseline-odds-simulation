%macro cm_analysis;

* For 1:1 cm data;

%*let maxsetsize=16;
%*let maxcaseset=8;
proc nlmixed data=case_cont; 
  by trial;
  parms beta=0;
*  profile beta /alpha=.05;
  array z{&maxsetsize} _z1-_z&maxsetsize; 
  array _cc{&maxsetsize} _cc1-_cc&maxsetsize;
  array phi[&maxsetsize] phi1-phi&maxsetsize;
  array _ntotss{2} _ntotss0 _ntotss1;
  array _ncases{2} _ncases0 _ncases1;
  array _msamp{2} _msamp0 _msamp1;

* initialize subject-or;
 ntot=_ntot; casesetor = 1; 
 do i = 1 to _ntot;
    phi[i] = exp(z[i]*beta);
    casesetor = casesetor * phi[i]**_cc[i];
    end;
*    put casesetor=;


*******************************************************
* RECURSION ARRAYS
* recursion arrays, set size to at least the max number of cases;
* a suffix is previous recursion level, b is current level;
* b0 prefix - sum of products
* NLMIXED does not allow for negative array indices so
*  start arrays at 1 and shift all array references by one
*******************************************************;

  array b0a[%eval(&maxcaseset+1)] _temporary_;
  array b0b[2,%eval(&maxcaseset+1)] _temporary_;

  one = 1;
  
* calculate sums of products using the recursive formula;



* initialize recursion array;
 do isam = 1 to 2;
 b0b[isam,0+one] = 1;
  do j = 1 to &maxcaseset; 
    b0b[isam,j+one] = 0; 
  end; 
  end;

**********************************************************
* calculate sums for each sampling stratum using recursion
**********************************************************;

 low = 0; high = 0;
 do isam = 1 to 2; * loop over sampling strata;

* get indices for members of the sampling stratum;
 if _msamp[isam] ne 0 then do; 
  low = high + 1;
  high = high + _msamp[isam];

* initialize a arrays (last level in the recursion);
  b0a[0+one] = 1; 
  do l = 1 to _totcases;
    b0a[l+one] = 0; 
    end; 

* do the recursion;
  do k = low to high;
    do l = 1 to min(k,_totcases);
      b0b[isam,l+one] = b0a[l+one] + b0a[l-1+one]*phi[k];
      end;

* re-initialize the m-1 step array;
   do l = 0 to min(k,_totcases);
     b0a[l+one]=b0b[isam,l+one]; 
     end; /* reinitialization loop */
    end;  /* recursion loop */
   end; /* _msamp ne 0 */
 end; /* sampling strata loop */


* Multiply contributions by the appropriate weights;
   do isam = 1 to 2;
    weight = 1;
    b0b[isam,0+one] =b0b[isam,0+one]*weight;
    do l = 1 to min(_msamp[isam],_totcases);
     weight =  weight*(_ntotss[isam]-l+1)/(_msamp[isam]-l+1);
     b0b[isam,l+one] = b0b[isam,l+one]*weight;
 *    put cluster= isam= l= weight= b0b[isam,l+1]=;
     end;
 end; /* sampling strata loop */

* sum of sets quantaties;
  
  den = 0;
  do k = max(0,_totcases-_msamp1) to min(_totcases,_msamp0);
   den = den + b0b[1,k+one]*b0b[2,(_totcases-k)+one];
*   put k= b0b[1,k+one]= b0b[2,(_totcases-k)+one]= den=;
   end;

    
*  put casesetor= den=;
*  CLL contribution from the counter-matched set;   
  dum = 1; 
  L = log(casesetor) - log(den);
  
  model dum ~ general(L);
  run;
 %mend;
