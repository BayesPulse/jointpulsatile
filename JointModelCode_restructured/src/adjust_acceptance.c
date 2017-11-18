/*******************************************************************/
/*************************adjust_acceptance.c *****************************/
/*******************************************************************/

/***FIX HERE***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"



/*********************************************************************/
/*START OF adjust_acceptance SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*adjust_acceptance: this adjusts the proposal variances based on the inputted
					acceptance rate and proposal variance. If the acceptance rate is
					too high or too low, proposal variance will be adjusted.
					ARGUMENTS: double x; the inputted acceptance rate; usually inputted as the
					acceptance counter divided by the attempt counter
					double *X; the current proposal variance
					RETURNS: None; update to the proposal variance is made internally
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 y: new proposal variance based on inputs

 SUBROUTINES USED
 None
 **************************************************************************/

void adjust_acceptance(double x, double *X)
{
	double y;

	y = 1. + 1000.*(x - .35)*(x - .35)*(x - .35);
	if (y < .9)
		y = .9;
	if (y > 1.1)
		y = 1.1;

	*X *= y;
}

