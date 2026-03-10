proc iml;

  /* Define the data matrix: 15 obs x 6 variables */
  x = {1 0 1 0 1 0,
       1 0 1 0 1 1,
       1 0 1 0 0 0,
       0 0 1 0 1 0,
       1 0 1 1 1 1,
       0 1 0 1 0 1,
       0 1 0 1 0 0,
       0 1 0 1 1 0,
       0 1 0 0 0 1,
       1 1 0 1 0 1,
       0 0 0 0 0 0,
       1 0 0 0 0 0,
       0 0 0 0 0 1,
       0 0 0 1 0 0,
       0 0 1 0 0 0};

  n = nrow(x);   /* 15 observations */
  p = ncol(x);   /* 6 variables     */

  /* Initialise all three similarity matrices */
  RR = J(n, n, 0);   /* Russell and Rao   */
  JC = J(n, n, 0);   /* Jaccard           */
  SM = J(n, n, 0);   /* Sokal-Michener    */

  do i = 1 to n;
    do j = 1 to n;

      cp = sum(x[i,]  # x[j,]);              /* concordant presences  : both = 1 */
      ca = sum((1-x[i,]) # (1-x[j,]));       /* concordant absences   : both = 0 */
      pa = sum(x[i,]  # (1-x[j,]));          /* presence-absence mismatch        */
      ap = sum((1-x[i,]) # x[j,]);           /* absence-presence mismatch        */

      /* Russell and Rao: cp / p */
      RR[i,j] = cp / p;

      /* Jaccard: cp / (cp + pa + ap) */
      denom_JC = cp + pa + ap;
      if denom_JC = 0 then JC[i,j] = 1;     /* both all-zero vectors -> identical */
      else JC[i,j] = cp / denom_JC;

      /* Sokal-Michener: (cp + ca) / p */
      SM[i,j] = (cp + ca) / p;

    end;
  end;

  /* Row/column labels */
  labels = char(1:n, 2);

  print RR[rowname=labels colname=labels
           label="Russell and Rao Similarity Matrix"];

  print JC[rowname=labels colname=labels
           label="Jaccard Similarity Matrix"];

  print SM[rowname=labels colname=labels
           label="Sokal-Michener Similarity Matrix"];

quit;
