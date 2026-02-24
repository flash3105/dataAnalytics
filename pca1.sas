*PCA;
*Tshemollo Rapolai;


proc import datafile="C:\Users\flash\Documents\2026\Masters\MVA\Assi1\pcomp.xlsx"
    out=mva
    dbms=xlsx
    replace;
    /* change if your sheet has a different name */
    getnames=yes;   /* use first row as variable names */
run;



proc iml;

use Mva;
read all  into X  ;


*print X;

cov1 = cov(X); *covariance matrix;

print cov1;

*standardise , by getting a correlation matrix;

cor = corr(X);

print cor;

*eigenvalues ;

call eigen(v,u,cor);

print v u;     *v contains the eigenvalues and u contains eigenvectors;


*v : Variance explained by that principal component;


pc1 = v[1]/7;
pc2 = v[2]/7;     *proportion of the explained variation;


print pc1 pc2;


/* Compute principal component scores */
Z = X * u;

/* Keep first two PCs */
PC12 = Z[,1:2];

create PCData from PC12[colname={"PC1" "PC2"}];
append from PC12;
close PCData;


quit;

run;


proc sgplot data=PCData;
   scatter x=PC1 y=PC2;
   xaxis label="Principal Component 1";
   yaxis label="Principal Component 2";
   title "Scatter Plot of First Two Principal Components";
run;




libname mylib "C:\Users\flash\Documents\2026\Masters\MVA\Assi1";  /* folder where your SAS dataset is */

proc iml;

*use Mva;
*use mylib.dat2;
use mylib.dat3;


read all  into X  ;


*print X;

cov1 = cov(X); *covariance matrix;

print cov1;

*standardise , by getting a correlation matrix;

cor = corr(X);

print cor;

*eigenvalues ;

call eigen(v,u,cor);

print v u;     *v contains the eigenvalues and u contains eigenvectors;


*v : Variance explained by that principal component;

total = trace(cor);
print total;

pc1 = v[1]/total;
pc2 = v[2]/total;
pc3 = v[3]/total;
pc4 =v[4]/total;
     *proportion of the explained variation;


print pc1 pc2 pc3 pc4;


/* Compute principal component scores */
Z = X * u;

/* Keep first two PCs */
PC12 = Z[,1:2];

create PCData from PC12[colname={"PC1" "PC2"}];
append from PC12;
close PCData;


quit;

run;


proc sgplot data=PCData;
   scatter x=PC1 y=PC2;
   xaxis label="Principal Component 1";
   yaxis label="Principal Component 2";
   title "Scatter Plot of First Two Principal Components";
run;
