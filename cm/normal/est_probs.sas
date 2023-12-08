%macro est_probs;
* estimated probabilities at macro z0 and z1;

data trialprobs; 

  set &sampmethod._co_var;

  z0= &z0; z1=&z1;
  or_z_sd =  exp(log(&or_z)/4);


  phat0_true = &lambda0*exp(log(&or_z)/4)**z0 /(1+&lambda0*exp(log(&or_z)/4)**z0);
  phat0 = exp(logodds+z0*bhat)/(1+exp(logodds+z0*bhat));
  vphat0 = (var_a_R + 2*z0*cov_ab_R + (z0**2)*vbhat)* ( phat0*(1-phat0) )**2;
  CIlow0 = phat0 - 1.96*sqrt(vphat0);
  CIhigh0 = phat0 + 1.96*sqrt(vphat0);
  inphat0CI = CIlow0 le phat0_true le CIhigh0;

  phat1_true = &lambda0*exp(log(&or_z)/4)**z1 /(1+&lambda0*exp(log(&or_z)/4)**z1);
  phat1 = exp(logodds+&z1*bhat)/(1+exp(logodds+&z1*bhat));
  vphat1 = (var_a_R + 2*z1*cov_ab_R + (z1**2)*vbhat)* ( phat1*(1-phat1) )**2;
  CIlow1 = phat1 - 1.96*sqrt(vphat1);
  CIhigh1 = phat1 + 1.96*sqrt(vphat1);
  inphat1CI = CIlow1 le phat1_true le CIhigh1;

  loglogphat0 = log(-log(phat0));
  varloglog0 = (var_a_R + 2*z0*cov_ab_R + (z0**2)*vbhat)* ( phat0*(1-phat0)/(-phat0*log(phat0)) )**2;
  if varloglog0 le 0 then put z0= var_a_R= cov_ab_R= vbhat= varloglog0=;
  CIlow0= exp(-exp( log(-log(phat0)) + 1.96*sqrt(varloglog0) ));
  CIhigh0 = exp(-exp( log(-log(phat0)) - 1.96*sqrt(varloglog0) ));
  inCI0 = CIlow0 le phat0_true le CIhigh0;

  loglogphat1 = log(-log(phat1));
  varloglog1 = (var_a_R + 2*z1*cov_ab_R + (z1**2)*vbhat)* ( phat1*(1-phat1)/(-phat1*log(phat1)) )**2;
  if varloglog1 le 0 then put z1= var_a_R= cov_ab_R= vbhat= varloglog1=;
  CIlow1= exp(-exp( log(-log(phat1)) + 1.96*sqrt(varloglog1) ));
  CIhigh1 = exp(-exp( log(-log(phat1)) - 1.96*sqrt(varloglog1) ));
  inCI1 = CIlow1 le phat1_true le CIhigh1;

  drop CIlow0 CIhigh0 CIlow1 CIhigh1;

  output;
  run;

 proc summary data=trialprobs;
   where vbhat lt 10;
   var bhat vbhat phat0 phat1 vphat0 vphat1 loglogphat0 loglogphat1 varloglog0 varloglog1;

   output out=probs mean(bhat)=bhat var(bhat)=vbhat mean(vbhat)=evbhat
   mean(phat0)=phat0  mean(phat1)=phat1 var(phat0)=vphat0  var(phat1)=vphat1 
   mean(vphat0)=evphat0 mean(vphat1)=evphat1
   mean(loglogphat0)=loglogphat0 mean(loglogphat1)=loglogphat1 var(loglogphat0)=vloglogphat0 var(loglogphat1)=vloglogphat1
   mean(varloglog0)=evarloglog0 mean(varloglog1)=evarloglog1
   mean(inphat0CI)=inphat0CI mean(inphat1CI)=inphat1CI mean(inCI0)=inCI0 mean(inCI1)=inCI1; 
   run;
%mend;
