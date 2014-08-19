#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gnuplot.h"

/* Plots the input signal as a line plot in gnuplot */
int plotLine(float input[], int len, int start_sample)
{
  int i;
  FILE *pipe = popen("gnuplot -persist", "w");
  //fprintf(pipe, "set terminal png\n");
  //fprintf(pipe, "set output '/home/acheeseman/high_rate/decay_correction/waveform_%d.png'\n", i);
  fprintf(pipe, "set style line 1 lw 1 lc rgb 'black'\n");
  fprintf(pipe, "set xlabel 'Time (10 ns samples)'\n");
  fprintf(pipe, "set xrange [%d:]\n", start_sample);
  fprintf(pipe, "plot '-' w l ls 1 notitle\n");
  for(i=0;i<len;i++)
    {
      fprintf(pipe, "%.2f\n", input[i]);
    }
  fprintf(pipe, "e\n");
  fprintf(pipe, "pause mouse\n");
  fclose(pipe);
  return 0;
}

/* Plots two input signals on the same axes in gnuplot */
int plotLine2(float input1[], float input2[], int len1, int len2, int start_sample)
{
  int i;
  FILE *pipe = popen("gnuplot -persist", "w");
  //fprintf(pipe, "set terminal png\n");
  //fprintf(pipe, "set output '/home/acheeseman/high_rate/decay_correction/waveform_%d.png'\n", i);
  fprintf(pipe, "set style line 1 lw 1 lc rgb 'black'\n");
  fprintf(pipe, "set style line 2 lw 1 lc rgb 'red'\n");
  fprintf(pipe, "set xlabel 'Time (10 ns samples)'\n");
  fprintf(pipe, "set xrange [%d:]\n", start_sample);
  fprintf(pipe, "plot '-' w l ls 1 notitle, '-' w l ls 2 notitle\n");
  for(i=0;i<len1;i++)
    {
      fprintf(pipe, "%.2f\n", input1[i]);
    }
  fprintf(pipe, "e\n");
  for(i=0;i<len2;i++)
    {
      fprintf(pipe, "%.2f\n", input2[i]);
    }
  fprintf(pipe, "e\n");
  fprintf(pipe, "pause mouse\n");
  fclose(pipe);
  return 0;
}

/* Plots the input data as a box plot */
int plotHisto(double input[], int maximum, int x_start, int x_end)
{
  int i;
  FILE *pipe = popen("gnuplot -persist", "w");
  fprintf(pipe, "reset\n");
  //fprintf(pipe, "set terminal png size 1280,480\n");
  //fprintf(pipe, "set output '/home/acheeseman/high_rate_test_data/spectrum-%d.png'\n", k);
  fprintf(pipe, "set style fill solid 0.5\n");
  fprintf(pipe, "set yrange [0:]\n");
  fprintf(pipe, "set xrange [%d:%d]\n", x_start, x_end);
  fprintf(pipe, "set xtics nomirror out\n");
  fprintf(pipe, "set mxtics\n");
  fprintf(pipe, "set xlabel 'Pulse Height'\n");
  fprintf(pipe, "set ylabel 'Counts'\n");
  fprintf(pipe, "plot '-' w boxes notitle\n");
  for (i=0;i<=maximum;i++)
    {
      fprintf(pipe, "%d\t%g\n", i, input[i]);
    }
  fprintf(pipe, "e\n"); 
  fclose(pipe);
  return 0;
}
