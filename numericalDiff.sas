proc iml;

/* Step size */
h = 0.0001;

/* Function f(x) = sin(x) */
start f(x);
   return(sin(x));
finish;

/* Point to evaluate */
x0 = 5;

/* First derivative (central difference) */
f_prime = (f(x0 + h) - f(x0 - h)) / (2*h);

/* Second derivative (central difference) */
f_doubleprime = (f(x0 + h) - 2*f(x0) + f(x0 - h)) / (h##2);

/* Print results */
print "f'(5)" f_prime;
print "f''(5)" f_doubleprime;

quit;


proc iml;

/* Define function */
start zFun(x,y);
   return(3*x*y + x##2 + y##2);
finish;

/* Evaluation point */
x0 = 5;
y0 = 3;

/* Small step size */
h = 0.0001;

/* First partial derivatives */
dz_dx = (zFun(x0+h, y0) - zFun(x0-h, y0)) / (2*h);
dz_dy = (zFun(x0, y0+h) - zFun(x0, y0-h)) / (2*h);

/* Second partial derivatives */
d2z_dx2 = (zFun(x0+h, y0) - 2*zFun(x0,y0) + zFun(x0-h, y0)) / (h##2);

/* Mixed partial derivative */
d2z_dxdy = ( zFun(x0+h,y0+h)
           - zFun(x0+h,y0-h)
           - zFun(x0-h,y0+h)
           + zFun(x0-h,y0-h) ) / (4*h##2);

/* Print results */
print dz_dx dz_dy d2z_dx2 d2z_dxdy;

quit;
