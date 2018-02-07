//
//  draw_pop_mass_mean.h
//  
//
//  Created by Nichole Carlson on 2/6/18.
//

#ifndef draw_pop_mass_mean_h
#define draw_pop_mass_mean_h

#include jt_deconvolution_main.h

int rmvnorm(double *, double **, int, double *, unsigned long *, int);
//need a mvnormcdf calculation
double rnorm(double, double, unsigned long *);
double kiss(unsigned long *);
int cholesky_decomp(double **, int);
double **cholesky_invert(int, double **); */
#endif /* draw_pop_mass_mean_h */
