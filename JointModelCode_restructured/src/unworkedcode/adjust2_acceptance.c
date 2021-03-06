/*******************************************************************/
/*************************adjust2_acceptance.c *********************/
/*******************************************************************/

/**fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"



/*********************************************************************/
/*START OF adjust2_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust2_acceptance: this adjusts the proposal variance-covriance matrix for
					 baseline and halflife based on the inputted acceptance rate
					 and proposal matrix. If the acceptance rate is too high or
					 too low, proposal variance will be adjusted.
					 ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
					 acceptance counter divided by the attempt counter
					 double **X; the current proposal variance-covariance matrix
					 double corr; the correlation between the two proposal variances
					 RETURNS: None; update to the proposal variance is made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new diagonal elements of proposal variance-covariance matrix based on inputs

 SUBROUTINES USED
 None
**************************************************************************/

void adjust2_acceptance(double x, double **X, double corr)
{
	double y;

	y = 1. + 1000.*(x - .25)*(x - .25)*(x - .25);
	if (y < .90) {
		y = .90;
		X[0][0] *= y;
		X[1][1] *= y;
		X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
	}
	if (y > 1.1) {
		y = 1.1;
		X[0][0] *= y;
		X[1][1] *= y;
		X[0][1] = X[1][0] = corr * sqrt(X[0][0] * X[1][1]);
	}
}

