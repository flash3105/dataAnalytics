

/* ============================================================
   Question 3: Kolmogorov-Smirnov Statistic — Bootstrap
   ============================================================ */

proc iml;

/* ============================================================
   STEP 1: Generate n=200 from Exponential(mean=100)
   Using Probability Integral Transformation (PIT):
       U ~ Uniform(0,1)
       X = -100 * ln(1 - U)  [inverse CDF of Exponential]
   ============================================================ */

call randseed(12345);          /* set seed for reproducibility  */
n      = 200;
lambda = 1/100;                /* rate parameter = 1/mean       */

U = j(n, 1, 0);
call randgen(U, "Uniform");    /* generate U ~ Uniform(0,1)     */

/* Inverse CDF of Exponential: F^{-1}(u) = -1/lambda * ln(1-u) */
x = (-1/lambda) * log(1 - U);

print "First 5 generated values:";
print (x[1:5]);

/* ============================================================
   STEP 2: KS Statistic assuming F* = Normal
   KS statistic: T = max|S(x) - F*(x)|
   where S(x) = empirical CDF, F*(x) = N(mu, sigma^2)
   ============================================================ */

/* Estimate Normal parameters from data (MLE) */
mu_hat    = mean(x);
sigma_hat = std(x);
print mu_hat sigma_hat;

/* Sort data for empirical CDF */


/* Compute KS statistic */
start ks_normal(x, mu, sigma);
    call sort(x, {1});
    n   = nrow(x);
    /* Empirical CDF values at each sorted point */
    S   = (1:n)` / n;             /* S(x_i) = i/n                  */
   * S_minus = (0:n-1)` / n;       /* S(x_i^-) = (i-1)/n            */

    /* Theoretical Normal CDF at each sorted point */
    F_star = cdf("Normal", x, mu, sigma);
    /* KS = max of |S(x_i) - F*(x_i)| and |S(x_i^-) - F*(x_i)| */
    D = max( abs(S - F_star) );
    return(D);
finish;

T_obs_normal = ks_normal(x, mu_hat, sigma_hat);
print "Observed KS Statistic (Normal):";
print T_obs_normal;

/* ============================================================
   STEP 3: Bootstrap estimate of empirical distribution of T
   1000 bootstrap resamples
   ============================================================ */
B = 1000;
T_boot    = j(B, 1, 0);          /* store bootstrap KS statistics */

do i = 1 to B;
    /* Resample n obs WITH replacement from x */
    idx =T(sample(1:n, n, "replace"));

    x_boot = x[idx];

    call sort(x_boot, {1});

    /* Re-estimate parameters from bootstrap sample */
    mu_b    = mean(x_boot);
    sigma_b = std(x_boot);

    /* Compute KS statistic on bootstrap sample */
    T_boot[i] = ks_normal(x_boot, mu_b, sigma_b);
end;

/* Sort bootstrap distribution for CI */
call sort(T_boot, {1});

/* ============================================================
   STEP 4: 95% Confidence Interval for KS statistic
   Using bootstrap percentiles: [2.5%, 97.5%]
   ============================================================ */
lower = T_boot[ ceil(0.025 * B) ];
upper = T_boot[ ceil(0.975 * B) ];

print "95% Bootstrap CI for KS Statistic (Normal):";
print lower upper;
*print T_boot;



/* ============================================================
   STEP 5: Test H0: data is from a Normal distribution
   Reject H0 if T_obs > upper bound of CI
   ============================================================ */
print "Observed KS Statistic vs 95% CI upper bound:";
print T_obs_normal upper;

if T_obs_normal > upper then
    print "REJECT H0: Data is NOT from a Normal distribution";
else
    print "FAIL TO REJECT H0: Data is consistent with Normal distribution";

/* ============================================================
   STEP 6: Repeat with F* = Exponential
   ============================================================ */

/* KS statistic for Exponential */
start ks_exponential(x, lambda);
    n      = nrow(x);
    S      = (1:n)` / n;
    S_minus= (0:n-1)` / n;
    /* Theoretical Exponential CDF: F(x) = 1 - exp(-lambda*x) */
    F_star = 1 - exp(-lambda # x);
    D      = max( abs(S - F_star) || abs(S_minus - F_star) );
    return(D);
finish;

/* MLE for exponential: lambda = 1/mean */
lambda_hat     = 1 / mean(x);
T_obs_exp      = ks_exponential(x, lambda_hat);

print "Observed KS Statistic (Exponential):";
print T_obs_exp;

/* Bootstrap for Exponential */
T_boot_exp = j(B, 1, 0);

do i = 1 to B;
    idx        = ceil( n * randfun(n, "Uniform") );
    x_boot     = x[idx];
    call sort(x_boot, {1});
    lambda_b        = 1 / mean(x_boot);
    T_boot_exp[i]   = ks_exponential(x_boot, lambda_b);
end;

call sort(T_boot_exp, {1});

lower_exp = T_boot_exp[ ceil(0.025 * B) ];
upper_exp = T_boot_exp[ ceil(0.975 * B) ];

print "95% Bootstrap CI for KS Statistic (Exponential):";
print lower_exp upper_exp;
print T_boot_exp;

print "Observed KS Statistic vs 95% CI upper bound (Exponential):";
print T_obs_exp upper_exp;

if T_obs_exp > upper_exp then
    print "REJECT H0: Data is NOT from an Exponential distribution";
else
    print "FAIL TO REJECT H0: Data is consistent with Exponential distribution";

/* ============================================================
   Compare both models
   ============================================================ */
print "--- Summary ---";
print "KS Normal:      " T_obs_normal;
print "KS Exponential: " T_obs_exp;

if T_obs_normal < T_obs_exp then
    print "Normal is a BETTER fit (smaller KS statistic)";
else
    print "Exponential is a BETTER fit (smaller KS statistic)";

/* ============================================================
   Output bootstrap distributions for plotting
   ============================================================ */
create boot_results var {"T_boot" "T_boot_exp"};
append from T_boot T_boot_exp;
close boot_results;

quit;

/* ============================================================
   Plot bootstrap distributions of T
   ============================================================ */
proc sgplot data=boot_results;
    density T_boot /
            lineattrs=(color=red thickness=2)
            legendlabel="Bootstrap KS — Normal";
    density T_boot_exp /
            lineattrs=(color=blue thickness=2 pattern=dash)
            legendlabel="Bootstrap KS — Exponential";
    xaxis label="KS Statistic T";
    yaxis label="Density";
    title "Bootstrap Distribution of KS Statistic — Normal vs Exponential";
    keylegend / location=inside position=topright;
run;
