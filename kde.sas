dm 'odsresults;clear';
title;
title1;
title2;
footnote;

goptions reset=global gunit=pct border
          ftext=swissb htitle=6 htext=3;

symbol1 interpol=stepj width=4
        color=blue
        value=dot
        height=.3;

PROC IMPORT OUT= WORK.kde_data
            DATAFILE="C:\Users\flash\Documents\2026\Masters\MVA\assi4\q2data.xls"
            DBMS=XLS REPLACE;
     SHEET="Sheet1";
     GETNAMES=NO;
RUN;


proc sgplot data=kde_data;
        histogram A;
run;

proc kde data=kde_data;
    univar A /  ngrid=200 unistats out=kde_sj;
run;


proc iml;

use kde_data;
read all var {"A"} into x1;
close kde_data;

ns = 300;
x  = T( sample(x1, ns, "NoRep") );   /* FIX: transpose to column vector */
n  = nrow(x);
call sort(x, {1});

/* Now mu, sigma, lambda will all be scalars */
mu     = mean(x);
sigma  = std(x);
lambda = 1 / mu;
print mu sigma lambda;

start kde(x, h, ngrid, mu, sigma, lambda);
    n     = nrow(x);
    min_x = x[1] - 3*h;
    max_x = x[n] + 3*h;
    grid  = T( do(min_x, max_x, (max_x-min_x)/ngrid) );

    /* Normal density — now works because mu and sigma are scalars */
    norm_dens = exp(-0.5 * ((grid - mu)/sigma)##2)
                / (sigma * sqrt(2*constant('PI')));

    /* Exponential density */
    exp_dens = lambda # exp(-lambda # grid);

    /* Safe loc() check */
    neg_idx = loc(grid < 0);
    if ncol(neg_idx) > 0 then exp_dens[neg_idx] = 0;

    start gaussian_kernel(u);
        return( exp(-0.5 * u##2) / sqrt(2*constant('PI')) );
    finish;

    density = j(nrow(grid), 1, 0);
    do i = 1 to nrow(grid);
        u          = (grid[i] - x) / h;
        density[i] = sum( gaussian_kernel(u) ) / (n*h);
    end;

    return( grid || density || norm_dens || exp_dens );
finish;

*h  = 3;
*h = 1.06 * std(x) * nrow(x)**(-1/5);  *silverman's rule of thumb;
*h=std(x);
*h=1;
h=0.5;
print h;
kde_result = kde(x, h, 200, mu, sigma, lambda);

create kde_plot var {"grid" "density" "normal_dens" "exp_dens"};
append from kde_result;
close kde_plot;

quit;

proc sgplot data=kde_plot;
    series x=grid y=density /
           lineattrs=(color=blue thickness=2 pattern=solid)
           legendlabel="KDE Estimate";
    series x=grid y=normal_dens /
           lineattrs=(color=red thickness=2 pattern=dash)
           legendlabel="Normal Density";
    series x=grid y=exp_dens /
           lineattrs=(color=green thickness=2 pattern=dot)
           legendlabel="Exponential Density";
    xaxis label="X values";
    yaxis label="Density";
    title "KDE vs Normal vs Exponential Density";
    keylegend / location=inside position=topright;
run;
