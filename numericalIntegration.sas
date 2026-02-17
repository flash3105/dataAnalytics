

proc iml;

mu = {0 0};                 /* row vector */
sigma = I(2);

invSigma = inv(sigma);
detSigma = det(sigma);

start bivnorm(x) global(mu, invSigma, detSigma);
   k = ncol(mu);
   diff = x - mu;
   exponent = -0.5 * diff * invSigma * diff`;
   const = 1 / ((2*constant("pi"))##(k/2) * sqrt(detSigma));
   return( const * exp(exponent) );
finish;

/* Create grid */
xGrid = do(-4, 4, 0.05);
yGrid = do(-4, 4, 0.05);
zMat  = j(ncol(xGrid), ncol(yGrid), .);

do i = 1 to ncol(xGrid);
   do j = 1 to ncol(yGrid);
      point = xGrid[i] || yGrid[j];
      zMat[i,j] = bivnorm(point);
   end;
end;

/*numerical integration*/

h = 0.05;                       /* grid spacing */
totalProb = sum(zMat) * h##2;  /* double integral approximation */

print totalProb;



posSum = 0;

do i = 1 to ncol(xGrid);
   do j = 1 to ncol(yGrid);
      if xGrid[i] > 0 & yGrid[j] > 0 then
         posSum = posSum + zMat[i,j];
   end;
end;

volumePositive = posSum * h##2;


print volumePositive;

/* Convert to dataset */
create density var {"x" "y" "z"};

do i = 1 to ncol(xGrid);
   do j = 1 to ncol(yGrid);
      x = xGrid[i];
      y = yGrid[j];
      z = zMat[i,j];
      append;
   end;
end;

close density;



/* Monte Carlo Integration */

call randseed(123);

a = -4;
b = 4;
areaXY = (b - a)##2;       /* 8^2 = 64 */

maxz = max(zMat);          /* maximum density height */

npoint = 100000;

start cal ;
/* Generate random (x,y) in square */
xRand = j(npoint,1);
yRand = j(npoint,1);
zRand = j(npoint,1);

call randgen(xRand, "Uniform", a, b);
call randgen(yRand, "Uniform", a, b);
call randgen(zRand, "Uniform", 0, maxz);

/* Evaluate density at random (x,y) */
count = 0;

do i = 1 to npoint;
   pt = xRand[i] || yRand[i];
   if zRand[i] < bivnorm(pt) then
      count = count + 1;
end;

prop = count / npoint;

/* Volume of 3D box */
boxVolume = areaXY * maxz;

/* Monte Carlo estimate */
mcVolume = prop * boxVolume;

finish;

call cal;

print mcVolume;




quit;
run;





proc g3d data=density;
   plot x*y=z;
run;
quit;



/* Direct Monte Carlo for positive quadrant */
proc iml;

call randseed(123);
npoint = 100000;

/* Generate X and Y from standard normal */
xRand = j(npoint,1);
yRand = j(npoint,1);

call randgen(xRand, "Normal", 0, 1);
call randgen(yRand, "Normal", 0, 1);

/* Count points where x>0 and y>0 */
countPos = sum( (xRand > 0) & (yRand > 0) );

/* Probability estimate for positive quadrant */
posProb = countPos / npoint;

print posProb;

quit;
