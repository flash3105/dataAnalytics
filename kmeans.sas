
quit ;
dm 'odsresults;clear';
title ;

libname q1 "C:\Users\flash\Documents\2026\Masters\MVA\assi2" ;
proc gplot data=q1.q1 ;
plot x1*x2 ;
run ;


proc iml ;


use q1.q1 ;
 read all into x ;

/*x=x[1:20,] ;*/

print "first 20 obs in x" (x[1:20,]) ;


n=nrow(x) ;

s3={1,2,3};
cen = x[s3,] ;

cen={10 30,
     20 40,
       40 40} ;


*cen={20 30,
     20 40,
       40 40} ;



icen=cen ;

print icen ;* x ;

con=1010101 ;

do tt = 1 to 20 while (con>0);

print "Iteration:" tt ;

print  (cen[1,])  x ;

v1=cen[1,] // x ;
print v1 ;

d1=(distance(v1,"L2"))[,1] ;
print d1 ;

print v1 d1 ;



v2=cen[2,] // x ;
d2=(distance(v2,"L2"))[,1] ;



v3=cen[3,] // x ;
d3=(distance(v3,"L2"))[,1] ;

d= (d1 || d2 || d3)[2:nrow(v1),] ;

/*print d ;*/


md = d[,><];
mdi= d[,>:<];

print v1 d1 v2 d2 v3 d3 , d md mdi;

gr=    (d=(md)) ;
print gr ;
print gr  x;


ncen1=(gr[,1]#x)[+,]  /     (gr[,1])[+,] ;
/*print ncen1 ;*/

ncen2=(gr[,2]#x)[+,]/(gr[,2])[+,] ;

ncen3=(gr[,3]#x)[+,]/(gr[,3])[+,] ;


print ncen1 ncen2 ncen3 ;

cenres = cenres // (i || ncen1 || ncen2 || ncen3 );




ocen=cen ;
cen = ncen1 // ncen2 // ncen3 ;
print cen ;


con = (abs(ocen-cen))[<>] ;
print con ;

*%%print cen ;

*%%print "end it" tt ;


con1= con1 // ( tt  || con  );
end ;



print con1 ;

print cenres ;


*cluster evaluation measures;

pdat = x || mdi ;
call sort(pdat,{3});
print pdat ;
x=pdat ;
cnt=J(1,nrow(gr),1)*gr ;
print cnt ;
cnt=cusum(cnt) ;
print cnt ;


tt=(x[,1:2] - repeat(J(1,n,1/n)*x[,1:2],n))##2 ;


* just to check something ;
/*
what1=x[,1:2] ;
what2=repeat(J(1,n,1/n)*x[,1:2],n) ;
print what1 what2 ;
*/

t = sum((x[,1:2] - repeat(J(1,n,1/n)*x[,1:2],n))##2)  ;
print tt t ;

wk = J(n,5,.) ;

cdat = pdat[,1:2] ;
mdat=J(n,2,.);
mdat = cen[pdat[,3],] ;


www=(cdat-mdat)##2 ;
w = sum((cdat-mdat)##2) ;
print cdat mdat www w ;

r2=1-w/t ;
psF = ((t-w) / (nrow(cen)-1)  )/ (w / (n-nrow(cen)) ) ;
print t w r2 psf;


/* Save cluster assignments to a dataset for plotting */
create work.pdat from pdat[colname={"x1" "x2" "group"}];
append from pdat;



*plotting data;



/* Number of iterations */
niter = nrow(cenres);

/* Create iteration column */
iter = T(1:niter);

/* Combine */
cenres2 = iter || cenres;

/* Add initial centroids */
init = {10 30,
        20 40,
        40 40};

init_row = init[1,] || init[2,] || init[3,];
init_full = 0 || init_row;

cenres2 = init_full // cenres2;

/* Save to dataset */
cnames = {"iter"
          "c1x" "c1y"
          "c2x" "c2y"
          "c3x" "c3y"};



create cent_path from cenres2[colname=cnames];
append from cenres2;
close cent_path;
quit;


data cent_long;
    set cent_path;

    centroid=1; x=c1x; y=c1y; output;
    centroid=2; x=c2x; y=c2y; output;
    centroid=3; x=c3x; y=c3y; output;

    keep iter centroid x y;
run;






data clustered;
    merge q1.q1 work.pdat;
run;

proc sgplot data=cent_long;
    styleattrs datacontrastcolors=(red green blue);
    series x=x y=y / group=centroid markers;
    xaxis label="X1 coordinate";
    yaxis label="X2 coordinate";
    title "Centroid Convergence Path";
run;







symbol1 color=red value=dot height=1;
symbol2 color=blue value=circle height=1;
symbol3 color=green value=triangle height=1;

proc gplot data=clustered;
    plot x1*x2=group;
run;
quit;
