data parms;
  input pz pc theta lambda0 or_z or_c;
/*  pz = 1-pz; pc = 1-pc; theta = 1/theta; lambda_0 = lambda_0*or_z; or_z = 1/or_z; */
  cards;
 .5 .3 2 0.4925 1 1
 .5 .3 2 0.3411 2 1
 .5 .3 2 0.2265 4 1
run;
/*
  .3 .3 2 .1110 1 1
  .3 .3 2 .0865 2 1
  .3 .3 2 .0617 4 1
  .3 .3 2 .0309 1 1
  .3 .3 2 .0239 2 1
  .3 .3 2 .0165 4 1

*/
