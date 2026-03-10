  /*PROC IMPORT DATAFILE="/home/u64159743/sasuser.v94/q1.xlsx"
            OUT=a
            DBMS=XLSX
            REPLACE;
            GETNAMES=YES;
RUN;*/

%web_drop_table(a);


FILENAME REFFILE "C:\Users\flash\Documents\2026\Masters\MVA\assi3\cluster.XLS"  ;

PROC IMPORT DATAFILE=REFFILE
      DBMS=XLS
      OUT=a;
      GETNAMES=YES;
RUN;

PROC CONTENTS DATA=a; RUN;


libname q1 "C:\Users\flash\Documents\2026\Masters\MVA\assi3" ;

proc gplot data=q1.X ;
plot  x2*x1 ;
run ;

%web_open_table(a);

/*data a ;
set a  (firstobs=9 obs=15) ;
run ;*/

proc print data=a ;
run ;

proc gplot data=q1.X ;
plot  x2*x1 ;
run ;

proc iml ;

********************* discuss at the end ;
/*This is a func that calculates the row and column of the clusters*/
start rowcol(pos,r,c);
if (pos/c - round(pos/c)) = 0 then do ; /*checking if it is a whole number */
 *print (pos/c) (round(pos/c)) ;
 row = round(pos/c) ;
 col = c ;
end ;
else do ; /*if not a whole number*/
  row = floor(pos/c) +1;
  col = pos-(row-1)*c ;
  *print row col ;
end ;
********************;
*print row col ;
rc = row || col ;
return rc;
finish rowcol ; /*close the function*/


use a  ;
 read all into x ;

 p=ncol(x) ;
 n=nrow(x);
 print x n ;
 nm =  (char(1):char(n))` ;/*store indices*/
 print nm ;
 nt = n ;

 dist=distance(x) + round(I(n)*(max(distance(x))+100)) ; /*adding 289 to the main diagonal*/
 print dist ;

 nm1=nm[1:nrow(x)] ;
 print dist[colname=nm1 rowname=nm1] ;

do j =  1 to  (nt-1) ;/*iteration to nrow -1*/

  print "********* START ***************";
  print  j ;

 r=nrow(dist) ;/*get the remaining rows in each iteration*/
 c=ncol(dist);/*get the remaining rows in each iteration*/
 n=r ;                  /*set n to the number of rows remaining*/
 print dist ;

 pos= (loc(dist=min(dist)))[1] ; /*look for the position of the minimum value per row and bring back the pos in the entire matrix*/

 print pos ;

 rc=(rowcol(pos,r,c))`;
 rct = rct // (j||n||rc`) ; /*continuosly add to the table as you iterate , add iteration number , rows and indices of the smallest*/
 call sort(rc,{1}) ; /*sort ascending?*/
 print rc ;

 dd1=dist[,rc[1]] ; /* rc has row and column of the smallest , since distance matrix is a square ... 3-1 and 1-3 will have the same value*/
 dd2=dist[,rc[2]] ;
 print dd1 dd2 ;

 ndist = (dist[,rc[1]] || dist[,rc[2]])[,><] ;
 *ndist = (dist[,rc[1]] || dist[,rc[2]])[,<>] ; /*get the biggest value per row between two column , returns a vector , 7x1...complete linkage*/

 ndist[rc[1]]=9090909 ; *replace the distance between a data point and its cluster ...it should be zero;
 print ndist ;

 dist[,rc[1]]=ndist ; *replace the smallest col with new cluster distance;
 dist[rc[1],]=ndist` ;*replace also the row with new distances to the cluster , have to use the transpose,new rep;

 print dist ;

 rkeep = loc((1:n)^=rc[2]) ;*columns indices where it does equal to the other col , because now it is represented by 1;
 print rkeep ;

 dist =dist[rkeep,rkeep] ;
 print dist ;

 print  "************* End ********************";
 end ;

 print rct ;

 nm1 = nm  // "*********************************************************8**********************************************" ;
 nmjj = nm1 ;
 print nm1 ;
 do jj = 1 to  nrow(rct) ;
   nmjj[rct[jj,3],]=rowcatc(nmjj[rct[jj,3],]  ||"+"||
                            nmjj[rct[jj,4],]);
   rkpos =loc((1:nt-jj+1)^=rct[jj,4]) ;
   nmjj=nmjj[rkpos,] ;
   nm1 = nm1 || (nmjj // J(jj+1,1,"-o-")  )  ;
   print jj nmjj nm1 ;
 end ;

 fclus = nm1[1:nrow(nm1)-1,] ;
 print fclus ;


 nclus=3 ;

 print nt ;
 fclus_s = fclus[1:nclus,nt-nclus+1] ;

 print fclus_s ;

 do i = 1 to nrow(fclus_s) ;
    print "Extract" i ;
    csel = fclus_s[i,] ;
    print  i csel ;
    delims = ' +';
    n = countw(csel, delims);
    words = scan(csel, 1:n, delims);  /* pass parameter vector: create vector of words */
    print words;
    words_n  = num(words)` ;
    print n words_n ;
    clusdat = clusdat // (J(nrow(words_n) ,1,i) || words_n) ;
end ;
print clusdat ;

  xdat1 = x ;
  mean = J(1,nrow(xdat1),1/nrow(xdat1))   *    xdat1 ;
  print mean ;

  xc = xdat1-mean ;
  T = trace(xc`*xc) ;
  print xdat1 xc T;

  do i = 1 to nclus ;
      cind = cind || (clusdat[,1]=i) ;
      xdatc =   (clusdat[,1]=i)#xdat1 ;
      nc = sum(clusdat[,1]=i) ;
      meanc = meanc // xdatc[+,]/nc ;
      print xdatc nc meanc ;
  end ;

  clmean = meanc[clusdat[,1],];
  print  cind meanc clmean ;

  w = trace((xdat1-clmean)`*(xdat1-clmean)) ;
  print nclus t w ;

    mm=repeat(mean,nclus) ; *print mm ;
    nk = (repeat((cind[+,]),p))` ;
    b = sum(nk#(meanc-mm)##2) ;
    bbbbb = (nk#(meanc-mm)##2) ;
    bk = (nk#(meanc-mm)##2)[,+] ;
    bs = (nk#(meanc-mm)##2)[+,] ;
    print meanc mm b bk bs bbbbb nk;

    t1= w+b ;
    print t t1 ;

    r2 = 1 - w/t ;
    r2k = bk/t ;
    r2s = bs/t ;
    print r2 r2k r2s;

    fc = ((t-w)/(nclus-1)) / (w/(nt-nclus)) ;
    print t w nclus nc, fc ;
  quit ;
quit ;
dm 'odsresults;clear';
title ;
  /*PROC IMPORT DATAFILE="/home/u64159743/sasuser.v94/q1.xlsx"
            OUT=a
            DBMS=XLSX
            REPLACE;
            GETNAMES=YES;
RUN;*/


%web_drop_table(a);

/*
FILENAME REFFILE "C:\Users\flash\Documents\2026\Masters\MVA\assi3\cluster.XLS"  ;

PROC IMPORT DATAFILE=REFFILE
      DBMS=XLS
      OUT=a;
      GETNAMES=YES;
RUN; */



libname a "C:\Users\flash\Documents\2026\Masters\MVA\assi2" ;

PROC CONTENTS DATA=a; RUN;


%web_open_table(a);

/*data a ;
set a  (firstobs=9 obs=15) ;
run ;*/

proc print data=a ;
run ;

proc gplot data=a ;
plot  x2*x1 ;
run ;

proc iml ;

********************* discuss at the end ;
/*This is a func that calculates the row and column of the clusters*/
start rowcol(pos,r,c);
if (pos/c - round(pos/c)) = 0 then do ; /*checking if it is a whole number */
 *print (pos/c) (round(pos/c)) ;
 row = round(pos/c) ;
 col = c ;
end ;
else do ; /*if not a whole number*/
  row = floor(pos/c) +1;
  col = pos-(row-1)*c ;
  *print row col ;
end ;
********************;
*print row col ;
rc = row || col ;
return rc;
finish rowcol ; /*close the function*/


use a  ;
 read all into x ;

 p=ncol(x) ;
 n=nrow(x);
 print x n ;
 nm =  (char(1):char(n))` ;/*store indices*/
 print nm ;
 nt = n ;

 dist=distance(x) + round(I(n)*(max(distance(x))+100)) ; /*adding 289 to the main diagonal*/
 print dist ;

 nm1=nm[1:nrow(x)] ;
 print dist[colname=nm1 rowname=nm1] ;

do j =  1 to  (nt-1) ;/*iteration to nrow -1*/

  print "********* START ***************";
  print  j ;

 r=nrow(dist) ;/*get the remaining rows in each iteration*/
 c=ncol(dist);/*get the remaining rows in each iteration*/
 n=r ;                  /*set n to the number of rows remaining*/
 print dist ;

 pos= (loc(dist=min(dist)))[1] ; /*look for the position of the minimum value per row and bring back the pos in the entire matrix*/

 print pos ;

 rc=(rowcol(pos,r,c))`;
 rct = rct // (j||n||rc`) ; /*continuosly add to the table as you iterate , add iteration number , rows and indices of the smallest*/
 call sort(rc,{1}) ; /*sort ascending?*/
 print rc ;

 dd1=dist[,rc[1]] ; /* rc has row and column of the smallest , since distance matrix is a square ... 3-1 and 1-3 will have the same value*/
 dd2=dist[,rc[2]] ;
 print dd1 dd2 ;

 ndist = (dist[,rc[1]] || dist[,rc[2]])[,><] ;
 *ndist = (dist[,rc[1]] || dist[,rc[2]])[,<>] ; /*get the biggest value per row between two column , returns a vector , 7x1...complete linkage*/

 ndist[rc[1]]=9090909 ; *replace the distance between a data point and its cluster ...it should be zero;
 print ndist ;

 dist[,rc[1]]=ndist ; *replace the smallest col with new cluster distance;
 dist[rc[1],]=ndist` ;*replace also the row with new distances to the cluster , have to use the transpose,new rep;

 print dist ;

 rkeep = loc((1:n)^=rc[2]) ;*columns indices where it does equal to the other col , because now it is represented by 1;
 print rkeep ;

 dist =dist[rkeep,rkeep] ;
 print dist ;

 print  "************* End ********************";
 end ;

 print rct ;

 nm1 = nm  // "*********************************************************8**********************************************" ;
 nmjj = nm1 ;
 print nm1 ;
 do jj = 1 to  nrow(rct) ;
   nmjj[rct[jj,3],]=rowcatc(nmjj[rct[jj,3],]  ||"+"||
                            nmjj[rct[jj,4],]);
   rkpos =loc((1:nt-jj+1)^=rct[jj,4]) ;
   nmjj=nmjj[rkpos,] ;
   nm1 = nm1 || (nmjj // J(jj+1,1,"-o-")  )  ;
   print jj nmjj nm1 ;
 end ;

 fclus = nm1[1:nrow(nm1)-1,] ;
 print fclus ;


 nclus=3 ;

 print nt ;
 fclus_s = fclus[1:nclus,nt-nclus+1] ;

 print fclus_s ;

 do i = 1 to nrow(fclus_s) ;
    print "Extract" i ;
    csel = fclus_s[i,] ;
    print  i csel ;
    delims = ' +';
    n = countw(csel, delims);
    words = scan(csel, 1:n, delims);  /* pass parameter vector: create vector of words */
    print words;
    words_n  = num(words)` ;
    print n words_n ;
    clusdat = clusdat // (J(nrow(words_n) ,1,i) || words_n) ;
end ;
print clusdat ;

  xdat1 = x ;
  mean = J(1,nrow(xdat1),1/nrow(xdat1))   *    xdat1 ;
  print mean ;

  xc = xdat1-mean ;
  T = trace(xc`*xc) ;
  print xdat1 xc T;

  do i = 1 to nclus ;
      cind = cind || (clusdat[,1]=i) ;
      xdatc =   (clusdat[,1]=i)#xdat1 ;
      nc = sum(clusdat[,1]=i) ;
      meanc = meanc // xdatc[+,]/nc ;
      print xdatc nc meanc ;
  end ;

  clmean = meanc[clusdat[,1],];
  print  cind meanc clmean ;

  w = trace((xdat1-clmean)`*(xdat1-clmean)) ;
  print nclus t w ;

    mm=repeat(mean,nclus) ; *print mm ;
    nk = (repeat((cind[+,]),p))` ;
    b = sum(nk#(meanc-mm)##2) ;
    bbbbb = (nk#(meanc-mm)##2) ;
    bk = (nk#(meanc-mm)##2)[,+] ;
    bs = (nk#(meanc-mm)##2)[+,] ;
    print meanc mm b bk bs bbbbb nk;

    t1= w+b ;
    print t t1 ;

    r2 = 1 - w/t ;
    r2k = bk/t ;
    r2s = bs/t ;
    print r2 r2k r2s;

    fc = ((t-w)/(nclus-1)) / (w/(nt-nclus)) ;
    print t w nclus nc, fc ;
  quit ;
quit ;
dm 'odsresults;clear';
title ;
  /*PROC IMPORT DATAFILE="/home/u64159743/sasuser.v94/q1.xlsx"
            OUT=a
            DBMS=XLSX
            REPLACE;
            GETNAMES=YES;
RUN;*/


%web_drop_table(a);

/*
FILENAME REFFILE "C:\Users\flash\Documents\2026\Masters\MVA\assi3\cluster.XLS"  ;

PROC IMPORT DATAFILE=REFFILE
      DBMS=XLS
      OUT=a;
      GETNAMES=YES;
RUN; */



libname a "C:\Users\flash\Documents\2026\Masters\MVA\assi3" ;

PROC CONTENTS DATA=a.X; RUN;


*%web_open_table(a);


/*data a ;
set a  (firstobs=9 obs=15) ;
run ;*/

*proc print data=a ;

*run ;

proc gplot data=a.X ;
plot  x2*x1 ;
run ;

proc iml ;

********************* discuss at the end ;
/*This is a func that calculates the row and column of the clusters*/
start rowcol(pos,r,c);
if (pos/c - round(pos/c)) = 0 then do ; /*checking if it is a whole number */
 *print (pos/c) (round(pos/c)) ;
 row = round(pos/c) ;
 col = c ;
end ;
else do ; /*if not a whole number*/
  row = floor(pos/c) +1;
  col = pos-(row-1)*c ;
  *print row col ;
end ;
********************;
*print row col ;
rc = row || col ;
return rc;
finish rowcol ; /*close the function*/


use a.X  ;
 read all into x ;

 p=ncol(x) ;
 n=nrow(x);
 print x n ;
 nm =  (char(1):char(n))` ;/*store indices*/
 print nm ;
 nt = n ;

 dist=distance(x) + round(I(n)*(max(distance(x))+100)) ; /*adding 289 to the main diagonal*/
 print dist ;

 nm1=nm[1:nrow(x)] ;
 print dist[colname=nm1 rowname=nm1] ;

do j =  1 to  (nt-1) ;/*iteration to nrow -1*/

  print "********* START ***************";
  print  j ;

 r=nrow(dist) ;/*get the remaining rows in each iteration*/
 c=ncol(dist);/*get the remaining rows in each iteration*/
 n=r ;                  /*set n to the number of rows remaining*/
 print dist ;

 pos= (loc(dist=min(dist)))[1] ; /*look for the position of the minimum value per row and bring back the pos in the entire matrix*/

 print pos ;

 rc=(rowcol(pos,r,c))`;
 rct = rct // (j||n||rc`) ; /*continuosly add to the table as you iterate , add iteration number , rows and indices of the smallest*/
 call sort(rc,{1}) ; /*sort ascending?*/
 print rc ;

 dd1=dist[,rc[1]] ; /* rc has row and column of the smallest , since distance matrix is a square ... 3-1 and 1-3 will have the same value*/
 dd2=dist[,rc[2]] ;
 print dd1 dd2 ;

 ndist = (dist[,rc[1]] || dist[,rc[2]])[,><] ;
 *ndist = (dist[,rc[1]] || dist[,rc[2]])[,<>] ; /*get the biggest value per row between two column , returns a vector , 7x1...complete linkage*/

 ndist[rc[1]]=9090909 ; *replace the distance between a data point and its cluster ...it should be zero;
 print ndist ;

 dist[,rc[1]]=ndist ; *replace the smallest col with new cluster distance;
 dist[rc[1],]=ndist` ;*replace also the row with new distances to the cluster , have to use the transpose,new rep;

 print dist ;

 rkeep = loc((1:n)^=rc[2]) ;*columns indices where it does equal to the other col , because now it is represented by 1;
 print rkeep ;

 dist =dist[rkeep,rkeep] ;
 print dist ;

 print  "************* End ********************";
 end ;

 print rct ;

 nm1 = nm  // "*********************************************************8**********************************************" ;
 nmjj = nm1 ;
 print nm1 ;
 do jj = 1 to  nrow(rct) ;
   nmjj[rct[jj,3],]=rowcatc(nmjj[rct[jj,3],]  ||"+"||
                            nmjj[rct[jj,4],]);
   rkpos =loc((1:nt-jj+1)^=rct[jj,4]) ;
   nmjj=nmjj[rkpos,] ;
   nm1 = nm1 || (nmjj // J(jj+1,1,"-o-")  )  ;
   print jj nmjj nm1 ;
 end ;

 fclus = nm1[1:nrow(nm1)-1,] ;
 print fclus ;


 nclus=3 ;

 print nt ;
 fclus_s = fclus[1:nclus,nt-nclus+1] ;

 print fclus_s ;

 do i = 1 to nrow(fclus_s) ;
    print "Extract" i ;
    csel = fclus_s[i,] ;
    print  i csel ;
    delims = ' +';
    n = countw(csel, delims);
    words = scan(csel, 1:n, delims);  /* pass parameter vector: create vector of words */
    print words;
    words_n  = num(words)` ;
    print n words_n ;
    clusdat = clusdat // (J(nrow(words_n) ,1,i) || words_n) ;
end ;
print clusdat ;

  xdat1 = x ;
  mean = J(1,nrow(xdat1),1/nrow(xdat1))   *    xdat1 ;
  print mean ;

  xc = xdat1-mean ;
  T = trace(xc`*xc) ;
  print xdat1 xc T;

  do i = 1 to nclus ;
      cind = cind || (clusdat[,1]=i) ;
      xdatc =   (clusdat[,1]=i)#xdat1 ;
      nc = sum(clusdat[,1]=i) ;
      meanc = meanc // xdatc[+,]/nc ;
      print xdatc nc meanc ;
  end ;

  clmean = meanc[clusdat[,1],];
  print  cind meanc clmean ;

  w = trace((xdat1-clmean)`*(xdat1-clmean)) ;
  print nclus t w ;

    mm=repeat(mean,nclus) ; *print mm ;
    nk = (repeat((cind[+,]),p))` ;
    b = sum(nk#(meanc-mm)##2) ;
    bbbbb = (nk#(meanc-mm)##2) ;
    bk = (nk#(meanc-mm)##2)[,+] ;
    bs = (nk#(meanc-mm)##2)[+,] ;
    print meanc mm b bk bs bbbbb nk;

    t1= w+b ;
    print t t1 ;

    r2 = 1 - w/t ;
    r2k = bk/t ;
    r2s = bs/t ;
    print r2 r2k r2s;

    fc = ((t-w)/(nclus-1)) / (w/(nt-nclus)) ;
    print t w nclus nc, fc ;
  quit ;
quit ;
dm 'odsresults;clear';
title ;
