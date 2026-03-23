title1 'Question 3 - First Observations';

data sasuser.question3;
    infile "C:\Users\flash\Documents\2026\Masters\TRA880\Unit2\kilowatt.txt"
           dlm=' '
           missover;

    input kwh;

    t = _n_;
    year = 2011 + int((t-1)/12);
    month = mod(t-1,12) + 1;
    date = mdy(month,1,year);

    format date monyy7.;
run;

proc print data=sasuser.question3(obs=10);
run;

/* --- Verify import --- */
title1 'Monthly Electricity Consumption - First and Last 5 obs';

proc print data=sasuser.question3 (obs=5) noobs;
    var t year month date kwh;
run;

proc print data=sasuser.question3 noobs;
    where t >= 92;
    var t year month date kwh;
run;

/* ============================================================
   STEP 2: Summary Statistics
   ============================================================ */
title1 'Descriptive Statistics - Monthly Electricity Consumption (kWh)';

proc means data=sasuser.question3
           n mean std min p25 median p75 max skewness kurtosis;
    var kwh;
run;

/* ============================================================
   STEP 3: Time Plot
   ============================================================ */
goptions reset=all;

title1 'Monthly Household Electricity Consumption - Virginia (2011-2018)';

axis1 label=(angle=90 'Electricity Consumption (kWh)');
axis2 label=('Date');

symbol1 color=steelblue value=dot i=join;

proc gplot data=sasuser.question3;
    plot kwh*date / vaxis=axis1 haxis=axis2;
run;
quit;

/* ============================================================
   STEP 4: Distribution - Histogram
   ============================================================ */
title1 'Distribution of Monthly Electricity Consumption';

proc univariate data=sasuser.question3 normal;
    var kwh;
    histogram kwh / normal kernel;
    qqplot kwh / normal(mu=est sigma=est);
run;

/* ============================================================
   STEP 5: ACF and PACF - up to lag 36 (3 years)
   ============================================================ */
goptions reset=all;

title1 'ACF and PACF - Monthly Electricity Consumption';

proc arima data=sasuser.question3
           plots(only)=series(acf pacf);
    identify var=kwh nlag=36;
run;
quit;

/* ============================================================
   STEP 6: Spectral Analysis
   ============================================================ */
goptions reset=all;

title1 'Spectral Analysis - Monthly Electricity Consumption';

proc spectra data=sasuser.question3
             out=sasuser.spectral
             p s
             adjmean
             whitetest;
    var kwh;
    weights 1 2 3 4 3 2 1;
run;

/* --- Print spectral output --- */
title1 'Spectral Output Table';

proc print data=sasuser.spectral noobs;
    var freq period p_01 s_01;
    format freq 8.4 period 8.2 p_01 s_01 12.2;
run;

/* ============================================================
   STEP 7: Plot spectrum vs FREQUENCY
   ============================================================ */
goptions reset=all;

title1 'Spectrum vs Frequency - Monthly Electricity Consumption';

axis1 label=('Frequency (radians)') minor=(number=1)
      order=0 to 3.1416 by 0.5236;

symbol1 color=steelblue value=dot i=join;
symbol2 color=tomato value=dot i=join;

proc gplot data=sasuser.spectral;
    plot (p_01 s_01)*freq / haxis=axis1 overlay legend;
run;
quit;

/* ============================================================
   STEP 8: Plot spectrum vs PERIOD
   ============================================================ */
goptions reset=all;

title1 'Spectrum vs Period - Monthly Electricity Consumption';

axis1 label=('Period (months)') minor=(number=1)
      order=0 to 50 by 5;

symbol1 color=steelblue value=dot i=join;
symbol2 color=tomato value=dot i=join;

proc gplot data=sasuser.spectral;
    plot (p_01 s_01)*period / haxis=axis1
                                 href=6 12
                                 overlay;
run;
quit;

/* ============================================================
   STEP 9: Identify the dominant frequency
   ============================================================ */
title1 'Dominant Frequency - Peak of Smoothed Spectrum';

proc sql;
    select freq,
           period,
           s_01 as smoothed_spectrum
    from sasuser.spectral
    having s_01 = max(s_01);
quit;




/* ============================================================
   Question 3(c): Fit Harmonic Regression Model using CLS (OLS)

   Model:
   X_t = µ
         + A1 cos(2pt/6) + B1 sin(2pt/6)
         + A2 cos(2pt/12) + B2 sin(2pt/12)
         + e_t

   where t = 1, 2, ..., 96
   (Jan 2011 = t=1, ..., Dec 2018 = t=96)
   ============================================================ */

/* ------------------------------------------------------------
   STEP 1: Construct trigonometric regressors
   ------------------------------------------------------------ */
data sasuser.kilowatt_trig;
    set sasuser.question3;

    pi = constant('pi');

    /* Semi-annual cycle (period = 6 months) */
    cos_6  = cos(2*pi*t/6);
    sin_6  = sin(2*pi*t/6);

    /* Annual cycle (period = 12 months) */
    cos_12 = cos(2*pi*t/12);
    sin_12 = sin(2*pi*t/12);
run;

/* --- Inspect first 13 observations --- */
title1 'Trigonometric Regressors (First 13 Observations)';

proc print data=sasuser.kilowatt_trig(obs=13) noobs;
    var t kwh cos_6 sin_6 cos_12 sin_12;
    format cos_6 sin_6 cos_12 sin_12 8.4;
run;

/* ------------------------------------------------------------
   STEP 2: Fit harmonic regression model (CLS = OLS)
   ------------------------------------------------------------ */
goptions reset=all;

title1 'Harmonic Regression Fit - Monthly Electricity Consumption';

proc reg data=sasuser.kilowatt_trig
         plots(only)=(diagnostics residuals fitplot);

    model kwh = cos_6 sin_6 cos_12 sin_12 / clb;

    output out=sasuser.trig_resid3
           predicted=yhat
           residual=resid;
run;
quit;

/* ------------------------------------------------------------
   STEP 3: Parameter estimates
   ------------------------------------------------------------ */
title1 'Estimated Harmonic Coefficients';

proc reg data=sasuser.kilowatt_trig;
    model kwh = cos_6 sin_6 cos_12 sin_12;
run;
quit;

/* ------------------------------------------------------------
   STEP 4: Observed vs Fitted values
   ------------------------------------------------------------ */
goptions reset=all;

title1 'Observed vs Fitted Values - Harmonic Model';

axis1 label=(angle=90 'Electricity Consumption (kWh)');
axis2 label=('Time Index (t)') order=1 to 96 by 12;

symbol1 color=black  value=dot  i=join;
symbol2 color=red    value=none i=join line=2;

proc gplot data=sasuser.trig_resid3;
    plot kwh*t=1
         yhat*t=2
         / overlay vaxis=axis1 haxis=axis2 legend;
run;
quit;

/* ------------------------------------------------------------
   STEP 5: Residual time plot
   ------------------------------------------------------------ */
goptions reset=all;

title1 'Residuals over Time - Harmonic Model';

axis1 label=(angle=90 'Residuals (kWh)');
axis2 label=('Time Index (t)') order=1 to 96 by 12;

symbol1 color=steelblue value=dot i=join;

proc gplot data=sasuser.trig_resid3;
    plot resid*t / vaxis=axis1
                   haxis=axis2
                   vref=0;
run;
quit;

/* ------------------------------------------------------------
   STEP 6: Residual autocorrelation diagnostics
   ------------------------------------------------------------ */
goptions reset=all;

title1 'Residual ACF and PACF - Harmonic Model';

ods graphics on;

proc arima data=sasuser.trig_resid3;
    identify var=resid nlag=36;
run;
quit;

ods graphics off;

/* ------------------------------------------------------------
   STEP 7: Residual normality diagnostics
   ------------------------------------------------------------ */
title1 'Residual Normality Diagnostics - Harmonic Model';

proc univariate data=sasuser.trig_resid3 normal;
    var resid;
    histogram resid / normal kernel;
    qqplot resid / normal(mu=est sigma=est);
run;




/* ------------------------------------------------------------
   STEP 1: Create future dataset (2019)
   ------------------------------------------------------------ */
data sasuser.future2019;
    do t = 97 to 108;

        year  = 2011 + int((t-1)/12);
        month = mod(t-1,12) + 1;
        date  = mdy(month,1,year);
        format date monyy7.;

        pi = constant('pi');

        /* Harmonic terms */
        cos_6  = cos(2*pi*t/6);
        sin_6  = sin(2*pi*t/6);
        cos_12 = cos(2*pi*t/12);
        sin_12 = sin(2*pi*t/12);

        output;
    end;
run;



data sasuser.kilowatt_all;
    set sasuser.kilowatt_trig sasuser.future2019;
run;

/* ------------------------------------------------------------
   STEP 3: Fit model and forecast
   ------------------------------------------------------------ */
proc reg data=sasuser.kilowatt_all noprint;
    model kwh = cos_6 sin_6 cos_12 sin_12;
    output out=sasuser.forecast_all predicted=yhat;
run;
quit;


/* ------------------------------------------------------------
   STEP 4: Time plot (Observed + Fitted + Forecast)
   ------------------------------------------------------------ */
goptions reset=all;

title1 'Observed, Fitted and Forecasted Electricity Consumption (2011–2019)';

axis1 label=(angle=90 'Electricity Consumption (kWh)');
axis2 label=('Date');

/* --- Define symbols --- */
symbol1 color=black value=dot  i=join;        /* Observed */
symbol2 color=blue  value=none i=join;        /* Fitted */
symbol3 color=red   value=none i=join line=2; /* Forecast */

/* --- Define legend --- */
legend1 label=none
        value=('Observed' 'Fitted / Forecast')
        position=(top right inside)
        mode=share;

/* --- Plot --- */
proc gplot data=sasuser.forecast_all;
    plot kwh*date=1
         yhat*date=2
         / overlay
           vaxis=axis1
           haxis=axis2
           legend=legend1;
run;
quit;
