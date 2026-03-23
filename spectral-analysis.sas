proc import out = sasuser.wolfer datafile ="C:\Users\flash\Documents\2026\Masters\TRA880\Unit2\wolfer.csv" dbms = csv replace;
        getnames = yes;
        datarow = 2;
run;

goptions reset = all i = join;
axis1 label = (angle = 90 'Wölfer sunspot numbers');
axis2 label = ('Year') minor = (number = 9) order = 1750 to 1900 by 10;
symbol1 color=black value = dot;
title1 'Time plot of the annual Wölfer sunspot numbers';
proc gplot data = sasuser.wolfer;
        plot sunspots * year / vaxis = axis1 haxis = axis2;
run;

goptions reset = all;
title1 'Annual Wölfer sunspot numbers';
proc arima data = sasuser.wolfer plots(only) = series(acf pacf);
        identify var = sunspots nlag = 25;
run;

goptions reset = all;
title1 'Spectral analysis';
data sunspot;
set sasuser.wolfer;
run;

proc spectra data = sunspot out = spectral p s adjmean whitetest;
        var sunspots;
    weights 1 2 3 4 3 2 1;
run;

proc print data = spectral;
run;

goptions reset=all i=join;
axis1 label = ('Frequency') minor = (number = 1) order = 0 to 3.1416 by 0.3927;
symbol1 color = black value = dot;
title1 'Spectral analysis';
proc gplot data = spectral;
        plot (p_01 s_01) * freq / haxis = axis1;
run;

goptions reset = all i = join;
axis1 label = ('Period') minor = (number = 4) order = 0 to 80 by 5;
symbol1 color = black value = dot;
title1 'Spectral analysis';
proc gplot data = spectral;
        plot (p_01 s_01) * period / haxis = axis1 href = 11;
run;


PROC ARIMA data = sasuser.wolfer;
    IDENTIFY var = sunspots nlag = 25;          /* Step 1: ID */
    ESTIMATE p = 2 method = CLS;                /* Step 2: Fit AR(2) CLS */
    ESTIMATE p = 1 method = CLS;                /* Compare AR(1) */
    ESTIMATE p = 2 q = 1 method = CLS;          /* Compare ARMA(2,1) */
RUN;


goptions reset = all;

title1 'AR(2) Forecast of Wolfer Sunspot Numbers (1901-1925)';


proc arima data = sasuser.wolfer;

    identify var = sunspots nlag = 25;
    estimate p = 2 method = CLS;

    forecast lead = 25

    alpha = 0.05

    out   = sasuser.forecasts

   id    =year;


run;

quit;

title1 'Forecast Values:1901-1925';
proc print data = sasuser.forecasts noobs;
        where year>=1901;
        var year forecast std l95 u95;
        format forecast std l95 u95 8.2;
run;


data observed;
        set sasuser.wolfer;
        type='Observed';
        value=sunspots;
        lower=.;
        upper=.;
        keep year value type lower upper;
run;


data fitted_forecast;
        set sasuser.forecasts;
        if year<=1900 then do;
           type = 'Fitted';
           value= forecast;
           lower=.;
           upper=.;
        end;

        if year >=1901 then do;
           type = 'Forecast';
           value = forecast;
           lower =l95;
           upper =u95;
        end;
        keep year value type lower upper;
run;

data combined;
        set observed fitted_forecast;
run;

proc sort data = combined;
        by year type;
run;

proc sgplot data=combined noautolegend;

   band x=year
        lower=lower
        upper=upper /
        transparency=0.5
        fillattrs=(color=lightcoral)
        legendlabel='95% confidence interval';

   series x=year y=value /
        lineattrs=(color=black thickness=1.5)
        markers
        markerattrs=(symbol=circlefilled size=4 color=black)
        legendlabel='Observed (1750-1900)'
        ;
   where type='Observed';

   series x=year y=value /
        lineattrs=(color=steelblue thickness=1.5 pattern=shortdash)
        legendlabel='Fitted values'
        ;
   where type='Fitted';

   series x=year y=value /
        lineattrs=(color=tomato thickness=2.5)
        markers
        markerattrs=(symbol=circlefilled size=5 color=tomato)
        legendlabel='Forecast (1901-1925)'
        ;
   where type='Foreca';

   refline 1900 /
        axis=x
        lineattrs=(color=gray thickness=1.5 pattern=dash)
        label='End of observed data (1900)'
        labelloc=inside;

   xaxis label='Year'
         values=(1750 to 1930 by 10)
         grid;

   yaxis label='Wolfer Sunspot Number'
         min=0
         grid;

   keylegend / location=inside
               position=topleft
               across=1;

run;




proc sgplot data=combined noautolegend;

   band x=year
        lower=lower
        upper=upper /
        transparency=0.5
        fillattrs=(color=lightcoral)
        name='ci'
        legendlabel='95% Confidence Interval';

   series x=year y=value /
        group=type
        lineattrs=(thickness=1.5)
        name='series';

   refline 1900 /
        axis=x
        lineattrs=(color=gray thickness=1.5 pattern=dash)
        label='End of observed data (1900)'
        labelloc=inside;

   xaxis label='Year' values=(1750 to 1930 by 10) grid;
   yaxis label='Wolfer Sunspot Number' min=0 grid;

   keylegend 'series' 'ci' /
        location=inside
        position=topleft
        across=1;

run;


 /* Create trigonometric variables */
data sasuser.wolfer_trig;
    set sasuser.wolfer;
    pi    = constant('pi');
    omega = 2 * pi / 11;
    cos_t = cos(omega * year);
    sin_t = sin(omega * year);
run;

/* Fit harmonic regression by CLS */
proc reg data=sasuser.wolfer_trig;
    model sunspots = cos_t sin_t;
    output out=sasuser.trig_resid
           predicted=yhat
           residual=resid;
run;
quit;


ods graphics on;
proc arima data=sasuser.trig_resid;
    identify var=resid(0);
run;
quit;
ods graphics off;



proc univariate data=sasuser.trig_resid normal;
    var resid;
    histogram resid / normal kernel;       /* separate options with spaces, not parentheses */
    qqplot resid / normal(mu=est sigma=est);
run;
quit;


/* Create forecast period */
data future;
    do year = 1901 to 1925;
        pi    = constant('pi');
        omega = 2 * pi / 11;
        cos_t = cos(omega * year);
        sin_t = sin(omega * year);
        output;
    end;
run;

data combined_all;
    set sasuser.wolfer_trig future;
run;


proc reg data=sasuser.wolfer_trig outest=estimates noprint;
    model sunspots = cos_t sin_t;
run;

/* Score (predict) for all years */
proc score data=combined_all score=estimates out=forecasted type=parms;
    var cos_t sin_t;
run;


data final_plot;
    set forecasted;

    if year <= 1900 then do;
        type='Observed';
        observed = sunspots;
        fitted   = model1;
    end;
    else do;
        type='Forecast';
        forecast = model1;
    end;
run;


proc sgplot data=final_plot;

   /* Observed */
   series x=year y=observed /
        markers
        lineattrs=(color=black)
        legendlabel='Observed';

   /* Fitted */
   series x=year y=fitted /
        lineattrs=(color=steelblue pattern=shortdash)
        legendlabel='Fitted';

   /* Forecast */
   series x=year y=forecast /
        lineattrs=(color=tomato thickness=2)
        legendlabel='Forecast';

   refline 1900 / axis=x
        lineattrs=(pattern=dash);

   xaxis label='Year';
   yaxis label='Sunspot Number';

run;
