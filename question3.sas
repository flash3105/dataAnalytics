proc iml;


pi = {0.3333333, 0.6666667};
mu = {10, 30};
sigma = {2, 8};


start normalPDF(x, mu, sigma);
   z = 1/(sqrt(2*constant("pi"))*sigma) * exp(-(x-mu)**2 / (2*sigma**2));
   return(z);
finish;

/* Grid */
x = do(-10, 60, 0.2);
f = j(ncol(x),1,0);

do i = 1 to ncol(x);
   f[i] = pi[1]*normalPDF(x[i], mu[1], sigma[1]) +
          pi[2]*normalPDF(x[i], mu[2], sigma[2]);
end;


print mu;
/* Quick check
print f[1:10];     should see small numbers increasing near first peak
print f[150:160];  near second peak (~30) */



/* Create dataset for plotting */
create mixPDF var {"x" "f"};
do i = 1 to ncol(x);
   xVal = x[i];
   fVal = f[i];
   append;
end;
close mixPDF;

minF = min(f);
print minF;   *check non-negativity;



dx = 0.2;  /* step size */
totalArea = sum(f) * dx; /* simple Riemann sum approximation */
print totalArea;

/* Indices where x <= 40 */
idx = loc(x <= 40);

/* Approximate CDF at 40 */
F40 = sum(f[idx]) * dx;
print F40;



/* Compute cumulative sum for the full CDF */
CDF = cusum(f) * dx;   /* cusum = cumulative sum in IML */

/* Find index where CDF first exceeds 0.2 */
idx20 = loc(CDF >= 0.2)[1];  /* take first match */
x20 = x[idx20];              /* 20th percentile */
print x20;


quit;

proc sgplot data=mixPDF;
   series x=x y=f / lineattrs=(thickness=2 color=blue);
   xaxis label="x";
   yaxis label="f(x)";
run;
