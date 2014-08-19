#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gaussianCurveFit.h"

/* Function which fits a Gaussian curve to the histogrammed data */
fitData(double y_dat[], int len, double p[NUM_PARAMS], int print, int plot)
{
  int i,j;
  int iteration = 0;
  final = 0;
  lambda = lambda_0;
  /* Procedure for each iteration */
  while (iteration < MAX_ITERATIONS && final == 0)
    {
      for (i=0;i<len;i++)
	{
	  func[i] = p[0] * expf(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
	  r[i] = y_dat[i] - func[i];
	}
      //printf("\n\nIteration %d\n", iteration+1);
      updateJac(len, p);
      solveDelta(y_dat, len);
      checkX2(y_dat, len, p, iteration);
      // printf("Delta: %g, Amplitude: %g\n", delta[0], p[0]);
      // printf("Delta: %g, Mean: %g\n", delta[1], p[1]);
      //printf("Delta: %g, Sigma: %g\n", delta[2], p[2]);
      iteration++;
    }
  if (print)
    {
      printf("\nAfter %d iterations the fit converged.\n\n", iteration);
      computeErrors(y_dat, len);
      printf("Degrees of Freedom: %d\n", len-NUM_PARAMS-1);
      printf("Final sum of squares of residuals: %g\n", ssr);
      printf("Reduced Chi-square: %g\n", ssr/ (double) (len-NUM_PARAMS-1));
      
      printf("Amplitude: %g\n", p[0]);
      printf("Mean: %g\n", p[1]);
      printf("Sigma: %g\n", p[2]);
    }
  
  if (plot)
    {
      FILE *pipe = popen("gnuplot -persist", "w");
      fprintf(pipe, "reset\n");
      //fprintf(pipe, "set terminal png font '/usr/share/fonts/webcore/arial.ttf' 11\n");
      // fprintf(pipe, "set output '/home/acheeseman/Desktop/co60_spectra/spectrum_new_%d.png'\n", z++);
      fprintf(pipe, "set style fill solid 0.5\n");
      fprintf(pipe, "set style line 2 lt 1 lw 2 pt 3 lc rgb 'blue'\n");
      //fprintf(pipe, "set xrange [1800:2200]\n");
      fprintf(pipe, "set xtics out nomirror\n");
      fprintf(pipe, "set mxtics 10\n");
      fprintf(pipe, "plot '-' w boxes notitle, '-' w l ls 2 t 'Gaussian Fit'\n");
      for (i=0;i<len;i++)
	{
	  fprintf(pipe, "%d\t%g\n", i, y_dat[i]);
	}
       fprintf(pipe, "e\n");
       for (i=0;i<len;i++)
	{
	  fprintf(pipe, "%d\t%g\n", i, func[i]);
	}
	fprintf(pipe, "e\n");
      //fprintf(pipe, "pause mouse\n");
      fclose(pipe);
    }
  return 0;
}

/* Updates the Jacobian, its transpose, their product, and its
   diagonal for a Gaussian function */
updateJac(int len, double p[NUM_PARAMS])
{
  int i,j,k;
  for (i=0;i<len;i++)
    {
      jac[i][0] = exp(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
      jac[i][1] = p[0]*((i-p[1])/(p[2]*p[2]))*exp(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
      jac[i][2] = p[0]*(((i-p[1])*(i-p[1]))/(p[2]*p[2]*p[2]))*exp(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
      
      jacT[0][i] = exp(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
      jacT[1][i] = p[0]*((i-p[1])/(p[2]*p[2]))*exp(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
      jacT[2][i] = p[0]*(((i-p[1])*(i-p[1]))/(p[2]*p[2]*p[2]))*exp(-((i-p[1])*(i-p[1]))/(2.*p[2]*p[2]));
    }

  // compute Hessian
  memset(H, 0, sizeof(H));
  for (i=0;i<NUM_PARAMS;i++)
    {
      for (j=0;j<NUM_PARAMS;j++)
	{
	  for (k=0;k<len;k++)
	    {
	      H[i][j] += jacT[i][k]*jac[k][j];
	    }
	}
    }
  // diag(H)
  for (i=0;i<NUM_PARAMS;i++)
    {
      for (j=0;j<NUM_PARAMS;j++)
	{
	  if ( i==j )
	    {
	      d[i][j] = H[i][j];
	    }
	  else
	    {
	      d[i][j] = 0.0;
	    }
	}
    }

  return 0;
}

/* Solves for delta (amount to increment the parameter by) from:
   [jTj + lambda*d]delta = jacT[y_dat - f(p)]                   */
solveDelta(double y_dat[], int len)
{
  int i,j,k;
  double A[NUM_PARAMS][NUM_PARAMS];
  double b[NUM_PARAMS];
  double C[NUM_PARAMS][NUM_PARAMS+1];

  memset(A, 0, sizeof(A));
  memset(b, 0, sizeof(b));
  // Calculate A & b
   for (i=0;i<NUM_PARAMS;i++)
    {
      for (j=0;j<NUM_PARAMS;j++)
	{
	  A[i][j] = H[i][j] + (lambda*d[i][j]);
	}
    }
  for (i=0;i<NUM_PARAMS;i++)
    {
      for (j=0;j<len;j++)
	{
	  b[i] += jacT[i][j]*r[j];
	}
    }

  // Gaussian Elimination
  int max;
  double t;
  for (i=0;i<NUM_PARAMS;i++)
    {
      C[i][NUM_PARAMS] = b[i];
      for (j=0;j<NUM_PARAMS;j++)
	{
	  C[i][j] = A[i][j];
	}
    }
  for (i=0;i<NUM_PARAMS;i++) {
    max = i;
    for (j=i+1;j<NUM_PARAMS;j++)
      if (C[j][i]>C[max][i])
	max = j;
    for (j=0;j<NUM_PARAMS+1;j++) {
      t = C[max][j];
      C[max][j] = C[i][j];
      C[i][j] = t;
    }
    for (j=NUM_PARAMS;j>=i;j--)
      for (k=i+1;k<NUM_PARAMS;k++)
	C[k][j] -= C[k][i]/C[i][i] * C[i][j];
  }
  for (i=NUM_PARAMS-1;i>=0;--i)
    {
      C[i][NUM_PARAMS] = C[i][NUM_PARAMS]/C[i][i];
      C[i][i] = 1;
      for (j=i-1;j>=0;j--)
	{
	  C[j][NUM_PARAMS] -= C[j][i]*C[i][NUM_PARAMS];
	  C[j][i] = 0;
	}
    }
  // set delta
  for (i=0;i<NUM_PARAMS;i++)
    {
      delta[i] = C[i][NUM_PARAMS];
     }
  return 0;
}

/*Check for convergence, then check if X2(p)-X2(p+delta) > epsilon_3*deltaT(lambda*delta + jT(y-f(p)),
and update lambda accordingly                                                                         */
checkX2(double y_dat[], int len, double p[NUM_PARAMS], int iteration)
{
  const double epsilon_1 = pow(10,-7);	// convergence tolerance for gradient
  const double epsilon_2 = pow(10,-7);	// convergence tolerance for parameters
  const double epsilon_3 = pow(10,-6);	// convergence tolerance for Chi-square

  int i,j;
  double p_new[NUM_PARAMS];
  double func_new[len];
  double x2_old, x2_new;
  for (i=0;i<NUM_PARAMS;i++)
    {
      p_new[i] = p[i] + delta[i];
    }
   for (i=0;i<len;i++)
    {
      func_new[i] = p_new[0] * expf(-((i-p_new[1])*(i-p_new[1]))/(2.*p_new[2]*p_new[2]));
    }

   double temp2[NUM_PARAMS] = {0};
   for (i=0;i<NUM_PARAMS;i++)
     {
       for (j=0;j<len;j++)
	 {
	   temp2[i] += jacT[i][j]*r[j];
	 }
     }
   
   // calculate Chi-squared values
   x2_old = (1./2)*dotProd(y_dat,y_dat,len) - dotProd(y_dat,func,len) + (1./2)*dotProd(func,func,len);
   x2_new = (1./2)*dotProd(y_dat,y_dat,len) - dotProd(y_dat,func_new,len) + (1./2)*dotProd(func_new,func_new,len);

   // test for convergence
   double max1 = 0;
   double max2 = 0;

   for (i=0;i<NUM_PARAMS;i++)
     {
       if (fabs(temp2[i]) > max1)
	 {
	   max1 = fabs(temp2[i]);   
	 }
     }
   for (i=0;i<NUM_PARAMS;i++)
     {
        if ((fabs(delta[i]/p[i])) > max2)
	 {
	   max2 = fabs((delta[i]/p[i]));
     	 }
     }

   if (max1 < epsilon_1 | max2 < epsilon_2 | x2_old/len < epsilon_3)
     {
       final = 1;
     }
   else
     {
       // test if new parameters are sufficiently better
       double t;
       t = epsilon_3 * dotProd(delta,temp2,NUM_PARAMS);
      
       
       if ((x2_old-x2_new) > t)
     {
       for (i=0;i<NUM_PARAMS;i++)
	 {
	   p[i] = p_new[i];
	 }
       lambda /= 10;
     }
       else
	 {
	   lambda *= 10;
	 }
     }
 return 0;
}

/* Computes the errors associated with the final fit to the data */
computeErrors(double y_dat[], int len)
{
  int i;
  ssr = 0;
  for (i=0;i<len;i++)
    {
      ssr += y_dat[i] - func[i];
    }
  return 0;
}

/* Computes the dot product of two vectors of length 'len' */
double dotProd(double a[], double b[], int len)
{
  double p = 0;
  if (sizeof(a) == sizeof(b))
    {
      int i;
      for (i=0;i<len;i++)
	{
	  p+= (a[i]*b[i]);
	}
      return p;
    }
  else
    {
      return 0;
    }
}
