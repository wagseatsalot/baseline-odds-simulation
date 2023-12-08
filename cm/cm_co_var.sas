%macro cm_co_var;
%* calculate the variance of case-only alpha hat from cluster cm data;
%* Inputs
%*  case_control - one line per cluster cm data;
%*  cotrial - case-only alphas with clr beta hats and vbeta hats;


data cm_co_var (keep=trial logodds bhat vbhat sigma2_D sigma2_R var_a_R var_b cov_ab_R);
  merge cotrial case_cont (in=in_case_cont);
  by trial;
  if not in_case_cont then delete;
******************************************
* input arrays:
    z - covariate values (will expand to multiple covariates later);
  array zt[1,&maxsetsize] _z1-_z%eval(1*&maxsetsize);
  array cct[&maxsetsize] _cc1-_cc%eval(1*&maxsetsize); * case-control status;


******************************************
* INPUT values:

* lambda - odds of disease with hat alpha, hat beta;
  array lambda[&maxsetsize] _temporary_;

* p_est - probability of disease with hat alpha, hat beta;
  array p_est[&maxsetsize] _temporary_;

* Summary values for the set;
  array _ntotss[2] _ntotss0 _ntotss1;
  array _ncases[2] _ncases0 _ncases1;
  array _msamp[2] _msamp0 _msamp1;

*******************************************************
* RECURSION ARRAYS
* recursion arrays, set size to at least the max number of cases;
* a suffix is previous recursion level, b is current level;
* b0 prefix - sum of products
* b1 prefix - sum of sum of single var times products
* b2 prefix - sum of product of two sums times products
*******************************************************;

  array b0a[0:&maxsetsize] _temporary_;
  array b0b[2,0:&maxsetsize] _temporary_;

* b1z - b1 with the covariate z as the summed value;
  array b1za[1,0:&maxsetsize] _temporary_;
  array b1zb[1,2,0:&maxsetsize] _temporary_;

* b1p - b1 with the 1/p as the summed value;
  array b1pa[0:&maxsetsize] _temporary_;
  array b1pb[2,0:&maxsetsize] _temporary_;

* b1qp2 - b1 with the q/p^2 as the summed value;
  array b1qp2a[0:&maxsetsize] _temporary_;
  array b1qp2b[2,0:&maxsetsize] _temporary_;

* b2zp - b2 with Z and 1/p as the summed values;
  array b2Zpa[0:&maxsetsize] _temporary_;
  array b2Zpb[2,0:&maxsetsize] _temporary_;

  retain sigma2_D Z_lambdainv_D totcases U_cov sigma2_R;
* initialize sums;
  if first.trial then do;
    sigma2_D=0; Z_lambdainv_D=0; totcases=0; U_cov=0; sigma2_R=0; 
	end;

* estimated quantities plugging in hat alpha and hat beta;
  if in_case_cont then do;
  do i = 1 to _ntot;
    lambda[i] = exp(logodds+zt[1,i]*bhat);
    p_est[i] = lambda[i] / (1+lambda[i]);
    end;
    end;

* initialize recursion array;
 do isam = 1 to 2;
  do j = 0 to &maxcaseset; 
    b0b[isam,j] = 0; b1zb[1,isam,j] = 0; b1pb[isam,j]=0; b1qp2b[isam,j]=0; b2zpb[isam,j]=0;
  end; 
  b0b[isam,0] = 1;
  end;

* calculate weights and variance pieces using the recursive formula;
/*  do i = 1 to 2; do j = 1 to &maxsetsize; put zt[i,j]= ; end; end; */

* calculate sums for each sampling stratum;
 low = 0; high = 0;
 do isam = 1 to 2;
 if _msamp[isam] ne 0 then do; 
  low = high + 1;
  high = high + _msamp[isam];

* initialize a arrays;
do l = 0 to _totcases;
   b0a[l] = 0; b1za[1,l]= 0; b1pa[l]=0; b1qp2a[l] = 0; b2zpa[l]= 0;
   end; 
   b0a[0] = 1;

* do the recursion;
  do k = low to high;
    do l = 1 to min(k,_totcases);
      b0b[isam,l] = b0a[l] + b0a[l-1]*lambda[k];
      b1zb[1,isam,l] = b1za[1,l] + b1za[1,l-1]*lambda[k] + b0a[l-1]*zt[1,k]*lambda[k];
      b1pb[isam,l] = b1pa[l] + b1pa[l-1]*lambda[k] + b0a[l-1]*(1/p_est[k])*lambda[k];
      b1qp2b[isam,l] = b1qp2a[l] + b1qp2a[l-1]*lambda[k] + b0a[l-1]*((1-p_est[k])/p_est[k]**2)*lambda[k];
      b2zpb[isam,l] = b2zpa[l] + b2zpa[l-1]*lambda[k] + lambda[k]*zt[1,k]*b1pa[l-1]
        + (lambda[k]/p_est[k])*b1za[1,l-1] + (lambda[k]*zt[1,k]/p_est[k])*b0a[l-1];
      end;

* re-initialize the m-1 step array;
   do l = 1 to min(k,_totcases);
     b0a[l]=b0b[isam,l]; b1za[1,l]=b1zb[1,isam,l]; b1pa[l]=b1pb[isam,l]; 
     b1qp2a[l]=b1qp2b[isam,l]; b2zpa[l]=b2zpb[isam,l];
     end; /* reinitialization loop */
    end;  /* recursion loop */
   end; /* _msamp ne 0 */
 end; /* sampling strata loop */


* Multiply contributions by the appropriate weights;
*  put cluster= _msamp[1]= _ntotss[1]= _msamp[2]= _ntotss[2]= _totcases= _z1= _z2=;
   do isam = 1 to 2;
    weight = 1;
*    put cluster= isam= "l=0 " weight= b0b[isam,0]= b1zb[1,isam,0]= b1pb[isam,0]= b1qp2b[isam,0]= b2zpb[isam,0]=;
    do l = 1 to min(_msamp[isam],_totcases);
     weight =  weight*(_ntotss[isam]-l+1)/(_msamp[isam]-l+1);

 *    put isam= l= weight= b0b[isam,l]= b1zb[1,isam,l]= b1pb[isam,l]= b1qp2b[isam,l]= b2zpb[isam,l]=;

     b0b[isam,l] = b0b[isam,l]*weight;
     b1zb[1,isam,l] = b1zb[1,isam,l]*weight;
     b1pb[isam,l] = b1pb[isam,l]*weight;
     b1qp2b[isam,l] = b1qp2b[isam,l]*weight;
     b2zpb[isam,l] = b2zpb[isam,l]*weight;
     end;
 end; /* sampling strata loop */


* cases only quantities;

  do i = 1 to _ntot;
	Z_lambdainv_D = Z_lambdainv_D + cct[i]*zt[1,i]/lambda[i];
	sigma2_D = sigma2_D + cct[i]*(1-p_est[i])/p_est[i]**2;
    end;

* sum of sets quantaties;
  
  b0sum = 0;
  do k = max(0,_totcases-_msamp1) to min(_totcases,_msamp0);
   b0sum = b0sum + b0b[1,k]*b0b[2,_totcases-k];
   end;
*   put cluster= b0sum=;

  b1zsum = 0; b1psum = 0; b1qp2sum = 0;
  do k = max(1,_totcases-_msamp1) to min(_totcases,_msamp0);
   b1zsum = b1zsum + b1zb[1,1,k]*b0b[2,_totcases-k]; 
   b1psum = b1psum + b1pb[1,k]*b0b[2,_totcases-k]; 
   b1qp2sum = b1qp2sum + b1qp2b[1,k]*b0b[2,_totcases-k]; 
   end;
  do k = max(1,_totcases-_msamp0) to min(_totcases,_msamp1);
   b1zsum = b1zsum + b1zb[1,2,k]*b0b[1,_totcases-k]; 
   b1psum = b1psum + b1pb[2,k]*b0b[1,_totcases-k]; 
   b1qp2sum = b1qp2sum + b1qp2b[2,k]*b0b[1,_totcases-k]; 
   end;

  b2zpsum = 0;
  do k = max(1,_totcases-_msamp1) to min(_totcases,_msamp0);
   b2zpsum = b2zpsum + b2zpb[1,k]*b0b[2,_totcases-k]
    + b1zb[1,1,k]*b1pb[2,_totcases-k]; 
*   put k= _totcases= b2zpb[1,k]= b0b[2,_totcases-k]= b1zb[1,1,k]= b1pb[2,_totcases-k]= b2zpsum= ;
   end;
  do k = max(1,_totcases-_msamp0) to min(_totcases,_msamp1);
   b2zpsum = b2zpsum + b2zpb[2,k]*b0b[1,_totcases-k]
    + b1zb[1,2,k]*b1pb[1,_totcases-k]; 
*   put k= _totcases= b2zpb[1,k]= b0b[2,_totcases-k]= b1zb[1,1,k]= b1pb[2,_totcases-k]= b2zpsum= ;
   end;

  U_cov = U_cov + b2zpsum/b0sum - b1zsum*b1psum/b0sum**2;
  sigma2_R = sigma2_R + b1qp2sum/b0sum;
*  put cluster= b0sum= b2zpsum= b1zsum= b1psum= U_cov= sigma2_R=;

  totcases = totcases + _totcases;

 * put trial= cluster= U_cov= sigma2_R= sigma2_D= ;
 
  if last.trial then do;
    var_a_R = (Z_lambdainv_D**2*vbhat - 2*Z_lambdainv_D*vbhat*U_cov + sigma2_D)/(&n-totcases)**2;
	var_b = vbhat;
	cov_ab_R = (U_cov - Z_lambdainv_D)* vbhat /(&n-totcases);
	output;
	end;
run;

proc corr data=cotrial cov outp=covoutput;
  where vbhat lt 10;
  var logodds bhat;
  run;
data cov_ab (keep=cov_ab);
  set covoutput;
  cov_ab = bhat;
  if _n_ eq 1 then output; else delete;
  label cov_ab="cov_ab=empirical alpha,beta covariance";
  run;

proc means data=cm_co_var mean var stderr min max;
  where vbhat lt 10;
  var sigma2_D sigma2_R var_a_R cov_ab_R;
  output out=CO_var_est mean= var(var_a_R)=vvar_a_R var(cov_ab_R)=vcov_ab_R;
  run;
data case_only;
  merge case_only cov_ab CO_var_est;
  run;

proc datasets; delete cov_ab covoutput; quit; 

%mend;
