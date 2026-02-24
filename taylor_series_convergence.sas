/* ============================================================
   TRUNCATED TAYLOR SERIES - Derivative Approximation
   Watching convergence term by term for f(x) = sin(x) at x=5
   True value: f'(5) = cos(5) = 0.28366...
   ============================================================ */

proc iml;

   /* --- Target function and true derivative --- */
   x0     = 5;
   h      = 0.5;          /* step size (large so we can see convergence) */
   true_d = cos(x0);      /* exact first derivative for comparison       */

   /* ============================================================
      PART 1: First Derivative via Taylor Series
      
      We expand f(x0+h) in a Taylor series:
        f(x0+h) = f(x0) + h*f'   + (h^2/2!)*f''  + (h^3/3!)*f'''  + ...
        f(x0-h) = f(x0) - h*f'   + (h^2/2!)*f''  - (h^3/3!)*f'''  + ...
      
      Subtracting:
        f(x0+h) - f(x0-h) = 2h*f' + 2*(h^3/3!)*f''' + ...
      
      Rearranging for f':
        f'  = [f(x0+h) - f(x0-h)] / 2h
            - (h^2/6)*f'''
            - (h^4/120)*f'''''  ...  (truncation error terms)
      
      We simulate "adding terms back in" iteratively below.
   ============================================================ */

   print "=== PART 1: First Derivative of sin(x) at x=5 ===";
   print "True value: cos(5) = " true_d;
   print " ";

   /* Precompute needed function values */
   f0   = sin(x0);
   fp1  = sin(x0 + h);
   fm1  = sin(x0 - h);
   fp2  = sin(x0 + 2*h);
   fm2  = sin(x0 - 2*h);
   fp3  = sin(x0 + 3*h);
   fm3  = sin(x0 - 3*h);

   /* --- Stage 1: Basic central difference (2-point, O(h^2)) --- */
   approx1 = (fp1 - fm1) / (2*h);
   err1    = abs(approx1 - true_d);
   term    = 1;
   print "Stage" term "  [2-point central diff, O(h^2)]" approx1 "  Error=" err1;

   /* --- Stage 2: 4-point formula (adds h^4 correction, O(h^4)) --- 
      Uses:  (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / 12h            */
   approx2 = (-fp2 + 8*fp1 - 8*fm1 + fm2) / (12*h);
   err2    = abs(approx2 - true_d);
   term    = 2;
   print "Stage" term "  [4-point formula,    O(h^4)]" approx2 "  Error=" err2;

   /* --- Stage 3: 6-point formula (adds h^6 correction, O(h^6)) ---
      Uses: (f(x+3h) - 9f(x+2h) + 45f(x+h) - 45f(x-h) + 9f(x-2h) - f(x-3h)) / 60h */
   approx3 = (fp3 - 9*fp2 + 45*fp1 - 45*fm1 + 9*fm2 - fm3) / (60*h);
   err3    = abs(approx3 - true_d);
   term    = 3;
   print "Stage" term "  [6-point formula,    O(h^6)]" approx3 "  Error=" err3;

   print " ";
   print "--- Convergence summary: error shrinks as we include more Taylor terms ---";
   results = approx1 // approx2 // approx3;
   errors  = err1    // err2    // err3;
   stages  = {1, 2, 3};
   print stages results errors;


   /* ============================================================
      PART 2: Richardson Extrapolation (iterative refinement)
      
      This is the most elegant iterative version of Taylor series.
      We compute the central difference with step h, h/2, h/4 ...
      and combine them to cancel higher order error terms one by one.
      
      If D(h) = f' + c2*h^2 + c4*h^4 + ...
      Then: D_better = [4*D(h/2) - D(h)] / 3   (cancels h^2 term)
      Then: D_best   = [16*D_better(h/2) - D_better(h)] / 15  (cancels h^4)
      
      This is Richardson's table, and it converges fast!
   ============================================================ */

   print " ";
   print "=== PART 2: Richardson Extrapolation (iterative Taylor refinement) ===";
   print "Starting step size h=0.5, halving each iteration";
   print " ";

   h_cur  = 0.5;
   nsteps = 6;

   /* Build first column of Richardson table (raw central differences) */
   D      = j(nsteps, nsteps, 0);   /* Richardson table */
   h_vec  = j(nsteps, 1, 0);

   do i = 1 to nsteps;
      h_vec[i] = h_cur;
      D[i,1]   = (sin(x0 + h_cur) - sin(x0 - h_cur)) / (2*h_cur);
      h_cur    = h_cur / 2;
   end;

   /* Fill Richardson table column by column */
   do j = 2 to nsteps;
      do i = j to nsteps;
         power    = 4**(j-1);
         D[i,j]   = (power * D[i,j-1] - D[i-1,j-1]) / (power - 1);
      end;
   end;

   /* Print diagonal (best estimate at each stage) */
   print "Iter   StepSize        Approx              Error";
   do i = 1 to nsteps;
      best_est = D[i,i];
      err_i    = abs(best_est - true_d);
      print i   h_vec[i]   best_est   err_i;
   end;

   print " ";
   print "Notice how the error drops dramatically each iteration!";
   print "This is the power of iterative Taylor series refinement.";


   /* ============================================================
      PART 3: Same idea for f(x,y) = 3xy + x^2 + y^2
      Partial derivative d/dx at (5,3) using Richardson
      True value: df/dx = 3y + 2x = 3(3) + 2(5) = 19
   ============================================================ */

   print " ";
   print "=== PART 3: Richardson on partial derivative of z=3xy+x^2+y^2 at (5,3) ===";
   print "True dz/dx = 3y + 2x = 19";
   print " ";

   x0    = 5;
   y0    = 3;
   true_p = 19;

   start zFun(x,y);
      return(3*x*y + x##2 + y##2);
   finish;

   h_cur  = 0.5;
   Dp     = j(nsteps, nsteps, 0);
   h_vec2 = j(nsteps, 1, 0);

   do i = 1 to nsteps;
      h_vec2[i] = h_cur;
      Dp[i,1]   = (zFun(x0+h_cur, y0) - zFun(x0-h_cur, y0)) / (2*h_cur);
      h_cur     = h_cur / 2;
   end;

   do j = 2 to nsteps;
      do i = j to nsteps;
         power    = 4**(j-1);
         Dp[i,j]  = (power * Dp[i,j-1] - Dp[i-1,j-1]) / (power - 1);
      end;
   end;

   print "Iter   StepSize        Approx              Error";
   do i = 1 to nsteps;
      best_est = Dp[i,i];
      err_i    = abs(best_est - true_p);
      print i   h_vec2[i]   best_est   err_i;
   end;

   print " ";
   print "For a polynomial, Richardson converges to machine precision very quickly";
   print "because the Taylor series terminates (finite number of non-zero terms).";

quit;
