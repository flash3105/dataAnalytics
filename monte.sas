
quit ;
dm 'odsresults;clear';

options ls=85 nodate pageno=1 nosource ;
title "Monte Carlo Integration" ;

%macro mcint(a,b,np) ;

proc iml ;

start eval(x) ;

y=5160*(1- 0.8*exp(-0.1*x));

return(y) ;
finish eval ;

x=((1:2000)/100)`;

y=eval(x) ; *evaluate the function in the points x ;
            *calculate the function value only - no integration;


xy = x || y ;
nm = {"x" "y"} ;

create a from xy[colname=nm] ;
append from xy ;


maxy = round((xy[,2])[<>]+1) ;
print maxy ;


*Monte Carlo integration section ;

a=&a ;
b=&b ;

npoint=&np ;
ds = J(npoint,1,0);

xc = ranuni(ds)*(b-a)+a ;
yc = ranuni(ds)*maxy ;

yy = eval(xc) ;

cprop = (yc<yy)[+]/npoint ;
area = cprop*((b-a)*maxy) ;

*print  xc yc yy cprop area;

datall = (xy || J(nrow(x),1,5) )  //
          (xc || yc || (yc<yy) ) ;

print area ;

nm1 = {"x" "y" "gr"} ;
create dall from datall[colname=nm1] ;
append from datall ;

create area from area[colname="area"] ;
append from area;

quit ;
%mend ;


%mcint(2,20,100000) ;
* lower boundary , upper boundary , number of points to generate ;


symbol1 interpol=none width=1
                color=blue
        value=dot
        height=.3;

symbol2 interpol=none width=1
                color=red
        value=dot
        height=.3;



proc gplot data=dall ;
plot y*x / href=(3 9) chref=purple vref=20 cvref=green;
where gr=5 ;

run ;
quit ;

proc gplot data=dall ;
plot y*x=gr ;
run ;
title ;
quit ;




proc import datafile="C:\Users\flash\Documents\2026\Masters\MVA\Assignments\q2.xlsx"
    out=ekt.q2
    dbms=xlsx
    replace;
    getnames=yes;   /* use first row as variable names */
run;
quit ;
dm 'odsresults;clear';
title ;

symbol1 value=dot
        height=2 i=none;
symbol2 value=dot
        height=2 i=none;
symbol3 value=dot
        height=2 i=none;

proc gplot data=ekt.q2 ;
plot x2*x1=y ;
title "k-nn example" ;

run ;

proc iml ;
 use ekt.q2 ;
 read all  into xy ;
/* xy=xy[1:25,];*/
 n=nrow(xy) ;
 print n xy ;

 y=xy[,3] ;
 x=xy[,1:2] ;
 *x=(x-x[:,])/std(x) ;
 n=nrow(xy) ;
 k=25 ;

 print y x  n;

 dist= distance(x) ;
 print dist ;


 do i = 1 to n ;
  *print "Observation:" i ;
  check = dist[,i] || y ;
  *print check ;
  call sort(check,{1}) ;
  *print check ;
  check = check[1:k,] ;
  *print k check ;

  cnt = sum((check[,2]=1)) || (check[,2]=2)[+]
                           || (check[,2]=3)[+]   ;
  *print cnt ;

  yh = cnt[<:>] ;
/*  print yh ; */

  yh_comb = yh_comb // yh ;
 end ;

  cor = (y=yh_comb)[:] ;

  print y yh_comb  cor;


  res = x || y || yh_comb ;
  nm={"x1" "x2" "y" "yh"} ;

  create result from res[colname=nm] ;
  append from res ;


quit ;

title ;

*new;





proc gplot data=result ;
plot x2*x1=yh ;
plot x2*x1=y ;
run ;
