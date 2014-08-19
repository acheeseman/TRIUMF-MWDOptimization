/* Functions to process waveforms input from a txt file and evaluate energy via MWD method.
Waveforms contained in a .mid file can be written to a txt file using Tig-Replay code.
Loops through values of the K and M parameters and calculates the energy resolution (FWHM) each iteration.

To compile: gcc -lm -g -o km_optimize_14 main.c gnuplot.c gaussianCurveFit.c
To run: ./km_optimize_14 ../inputfile.txt ../outputfile.txt */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gnuplot.h"
#include "gaussianCurveFit.h"

// waveform properties here - change these based on the input file
// values are in units of 10 ns samples
#define NUM_SAMPLES 1000	      // number of samples in each waveform
#define MAX_WAVEFORMS 150000          // maximum number of waveforms to process
#define TIME_THRESHOLD 500            // threshold crossing time of pulse
#define PEAK_ENERGY 975.65            // energy of peak to calibrate spectrum and calculate energy resolution
#define PEAK_NUM 3                    // which peak in the spectrum we want to fit

// signal processing constants
#define WINDOW_WIDTH 500              // window width parameter (L) for MWD 
#define INTEGRATION_DELAY 100         // time after threshold crossing to start integration
#define K_MIN 50
#define K_MAX WINDOW_WIDTH - INTEGRATION_DELAY
#define K_STEP 50
#define M_MIN 1000
#define M_MAX 10000
#define M_STEP 500

// curve fitting
#define GAUSSIAN_PEAK_WIDTH 25         // half the number of samples to include in the peak fit   

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define length(x) sizeof(x)/sizeof(x[0])

static float v[NUM_SAMPLES] = {0};
static float y[NUM_SAMPLES] = {0};
static int pulseHeights[MAX_WAVEFORMS] = {0};

/* Function Declarations */
int deconvolutionFilterFloat(int L, int M, float input[NUM_SAMPLES], float output[NUM_SAMPLES]);
float pulseHeightFloat(float input[NUM_SAMPLES], int K);

int main(int argc, char *argv[])
{
  int i,j,n;
  FILE *inputfile,*outputfile;
  
  // mwd parameters
  int k;                   // integration time
  int l = WINDOW_WIDTH;    // l is fixed
  int m;                   // decay constant

  int min;
  int max;
  int chan;

  double fwhm;                 // FWHM = 2.35 * sigma
  double fwhm_min = 100000;    // to keep track of the minimum FWHM
  int k_min;                   // to keep track of the k param for which fwhm is minimized
  int m_min;                   // to keep track of the m param for which fwhm is minimized

  int num_waveforms;
  int good_waveforms = 0;

  double p[NUM_PARAMS];
  double *histData = NULL;
  double peakData[2*GAUSSIAN_PEAK_WIDTH+1];

  int peak_centre, peak_start;
  double peak_amp;

  /* Read in waveform data from input file*/
  if ((inputfile = fopen(argv[1],"r")) == NULL)
    {
      fprintf(stdout, "Input file does not exist.\n");
    }
  else
    {
      // check if output file already exists
      if ((outputfile = fopen(argv[2],"wx")) == NULL)
	{
	  fprintf(stdout, "Could not open output file or file already exists.\n");
	}
      else
	{
	  fclose(outputfile);
	  // count the total number of waveforms in file
	  char c;
	  long total_waveforms;
	  long lines;
	  
	  while ((c=getc(inputfile)) != EOF)
	    {
	      if (c=='\n')
		{
		  lines++;
		}
	    }
	  total_waveforms = (long) (lines/NUM_SAMPLES*1.0 + 0.5);
	  num_waveforms = MIN(MAX_WAVEFORMS, total_waveforms);
	  for (k=K_MIN;k<=K_MAX;k+=K_STEP)
	    {
	      outputfile = fopen(argv[2],"a");
	      fprintf(outputfile, "\n");
	      for (m=M_MIN;m<=M_MAX;m+=M_STEP)
		{
		  printf("k = %d, m = %d\n", k, m);
		  good_waveforms = 0;
		  fseek(inputfile, 0, SEEK_SET);
		  for (i=0;i<num_waveforms;i++)
		    {
		      for (j=0;j<NUM_SAMPLES;j++)
			{
			  fscanf(inputfile, "%f\n", &v[j]);
			}
		      // filter waveform with MWD and evaluate pulse height
		      differenceFilterFloat(l, v, y);
		      deconvolutionFilterFloat(l, m, y, y);
		      pulseHeights[i] = (int) (pulseHeightFloat(y, k) + 0.5);
		      //plotLine(v, NUM_SAMPLES, 0);
		      // plotLine(y, NUM_SAMPLES, 0);
		    }
		  
		  /* Plot Energy Spectrum */
		  // Find minimum and maximum pulse heights   
		  max = 0;
		  min = 100000;
		 
		  for (i=0;i<num_waveforms;i++)
		    {
		      if (pulseHeights[i] > max)
			{
			  max = pulseHeights[i];
			}
		    }
		  for (i=0;i<num_waveforms;i++)
		    {
		      if (pulseHeights[i] > 0 && pulseHeights[i] < min)
			{
			  min = pulseHeights[i];
			}
		    }
		  histData = (double*)calloc(max+1, sizeof(double));

		  for (i=0;i<num_waveforms;i++)
		    {
		      if (pulseHeights[i] > 0)
		  	{
		  	  good_waveforms++;
		  	  chan = pulseHeights[i];
		  	  histData[chan] += 1.0;
		  	}
		    }
		  plotHisto(histData, max, 0, 1500);
	  
		  /* Fit peak to find FWHM */
		  // we only want to fit the peak, not the whole spectrum 	  
		  // find the peak we want to fit
		  peak_start = peakFind(histData, max, PEAK_NUM, 300, 50);
		  peak_amp = 0;
		  for (i=peak_start;i<peak_start+GAUSSIAN_PEAK_WIDTH+1;i++)
		    {
		      if (histData[i] > peak_amp)
		  	{
		  	  peak_amp = histData[i];
		  	  peak_centre = i;
		  	}
		    }
		  free(histData);
		  for (i=0;i<GAUSSIAN_PEAK_WIDTH*2+1;i++)
		    {
		      peakData[i] = histData[(peak_centre-GAUSSIAN_PEAK_WIDTH)+i];
		  }
		  p[0] = peak_amp*5.5;
		  p[1] = GAUSSIAN_PEAK_WIDTH;
		  p[2] = 0.5;
		  fitData(peakData, GAUSSIAN_PEAK_WIDTH*2+1, p, 0, 0);
		  // calculate the fwhm from fit parameters and keep track of the minimum value
		  if (2.35*p[2]<fwhm_min)
		    {
		      fwhm_min = 2.35*p[2];
		      k_min = k;
		      m_min = m;
		    }
		  fprintf(outputfile, "%d\t%d\t%g\n", k, m, 2.35*p[2]*PEAK_ENERGY/(p[1]-GAUSSIAN_PEAK_WIDTH+peak_centre));
		  printf("Processed %d waveforms, %d counts in spectrum\n", num_waveforms, good_waveforms);
		}
	      fclose(outputfile);
	    }
	  fclose(inputfile);
	  //print summary to screen
	  printf("Finished processing file %s\n", argv[1]);
	  printf("Mminimum energy resolution: %g channels and %g keV at k = %d and m = %d\n", fwhm_min, fwhm_min*PEAK_ENERGY/(p[1]-GAUSSIAN_PEAK_WIDTH+peak_centre), k_min, m_min);
	}
    }
  return 0;
}

/*** SIGNAL PROCESSING FUNCTIONS ***/
int differenceFilterFloat(int L, float input[NUM_SAMPLES], float output[NUM_SAMPLES])
{
  int i;
  for (i=L;i<NUM_SAMPLES;i++)
    {
      output[i] = input[i] - input[i-L];
    }
  return 0;
}
/* Difference and Decay Correction - floating point */
deconvolutionFilterFloat(int L, int M, float input[NUM_SAMPLES], float output[NUM_SAMPLES])
{ 
  int i,j;
  float sum;
  float waveform_temp[NUM_SAMPLES];
  for (i=L;i<NUM_SAMPLES;i++)
    {
      sum = 0.0;
      for (j=i-L;j<=i-1;j++)
	{
	  sum += input[j];
	}
      waveform_temp[i] = input[i] + (1/ ((float) M) ) * sum ;
    }
   for (i=L;i<NUM_SAMPLES;i++)
     {
       output[i] = waveform_temp[i];
     }
  return 0;
}
/* Returns pulse height of the filtered waveform */
float pulseHeightFloat(float input[NUM_SAMPLES], int K)
{
    //evaluate the energy for one waveform
  float pulseHeight = 0;
  int j;
  for (j=0;j<K;j++)
    {
      pulseHeight += (1/ ((float) K) ) * (input[TIME_THRESHOLD+INTEGRATION_DELAY+j]);
     }
  return pulseHeight;
}

/* Finds the nth  peak of an input array */
int peakFind(double input[], int len, int peakNum, int thresh, int width)
{
  int peakList[100] = {0};
  int peakWidthHalf = (int) MAX(1,floor(width*0.5));
  int max;
  int minIndx, maxIndx, peakIndx;
  int i, j, n;
  n = 0;
  
  for (i=0;i<len;i++)
    {
      if ( (int) (input[i]+0.5) > thresh)
	{
	  minIndx = i - peakWidthHalf;
	  maxIndx = i + peakWidthHalf;
	  if ( minIndx < 0 ) { minIndx = 0; }
	  if ( maxIndx > len ) { maxIndx = len; }
	  max = 0;
	  for (j=minIndx;j<maxIndx;j++)
	    {
	      if ( input[j] > max ) { 
		max = input[j];
		peakIndx = j; }
	    }
	  if ( peakIndx == i ) { peakList[n++] = i; }
	}
    }
      if (n < peakNum)
	{
	  printf("Could not find peak number %d\n", peakNum);
	  return -1;
	}
  return peakList[peakNum-1];
}
