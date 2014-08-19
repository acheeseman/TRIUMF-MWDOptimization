TRIUMF-MWDOptimization
======================

#####Author: Alison Cheeseman, 2014

C code which processes waveforms contained in a .txt file and evaluates their energy using the Moving Window Deconvolution
algorithm. The algorithm consists of three filters, a differentiation, a decay correction, and an integration. These filters are characterized by three parameters which influence the energy response of the system. The parameters are referred to by:

L: The differentiation width of the first filter.

M: The decay constant of the preamplifier signal - used in the decay correction

K: The integration time, used in the third filter (K must be less than L + the risetime of the signal)

This code loops through values of the K and M parameters for a fixed value of L, producing a spectrum and calculating the
energy resolution (FWHM) of a prespecified peak by fitting a Gaussian function to the peak using the Levengberg-Marquardt
algorithm.

The code can be modified to run on different sets of waveforms and to fit peaks of different energies by modifying constants
defined in main.c

To compile: gcc -lm -g -o km_optimize_14 main.c gnuplot.c gaussianCurveFit.c

To run: ./km_optimize_14 .../inputfile.txt .../output file.txt
