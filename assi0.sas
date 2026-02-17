/*Assignment 0*/
proc iml;

mu = 10;
sigma = 3;
n=25;
n_sim = 1000;


*Simulation;

X = randfun(n_sim || n ,"Normal",mu,sigma);

/*efficiently generates a matrix or vector of random numbers from a specified statistical distribution*/

/*Syntax: x = RANDFUN(N, "Distribution" <,param1> <,param2> <,param3>) */

/*compute sample means*/

Xbar = X[,:]; *row means;

*print(X);
*print(Xbar);

/*comparison*/

emp_mean= mean(Xbar);
emp_std = std(Xbar);

print  emp_mean mu ;
print  emp_std (sigma/sqrt(n));

z=1.96;

/* Compute CI lower and upper bounds for each sample */
Lower = Xbar - z * sigma / sqrt(n);
Upper = Xbar + z * sigma / sqrt(n);

/* Check if CI contains the true mean */
ContainsMu = (Lower <= mu) & (Upper >= mu);

print ContainsMu;

/* Compute proportion of CIs that contain mu */
coverage = mean(ContainsMu);

print coverage;

/* Create SAS dataset */
create simdata from Xbar[colname={"Xbar"}];
append from Xbar;
close simdata;



quit;

/* Now plot histogram */
proc sgplot data=simdata;
   histogram Xbar;
   density Xbar / type=normal;
run;


proc sort data=simdata; by Xbar; run;

data ecdf;
   set simdata;
   n + 1;
   Fn = n / 1000;
run;


data cdfplot;
   set ecdf;
   theor_cdf = probnorm((Xbar - 10)/0.6); /* Z = (x - mu)/sigma */
run;



proc sgplot data=cdfplot;
   series x=Xbar y=Fn        / lineattrs=(color=blue thickness=2) legendlabel="Empirical CDF";
   series x=Xbar y=theor_cdf / lineattrs=(color=red  thickness=2 pattern=shortdash) legendlabel="Theoretical CDF";
   xaxis label="Sample Mean (\bar{X})";
   yaxis label="CDF";
   keylegend / location=inside position=topright across=1;
run;





run;
