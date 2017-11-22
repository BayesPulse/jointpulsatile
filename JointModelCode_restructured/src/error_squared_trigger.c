/*******************************************************************/
/*************************error_squared_trigger.c *****************************/
/*******************************************************************/

/**fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/*********************************************************************/
/*START OF error_squared SUBROUTINE*/
/*********************************************************************/
/*********************************************************************/
/*error_squared: This subroutine calculates the sum of squared error. It uses
				the list of pulses and the common parameters to evaluate the
				likelihood and calculate the sum of the squared error between
				observed concentration and expected under the current parameters
				ARGUMENTS: double **ts; this is the matrix of observed data (a column of
				times and S columns of log(concentration)
				Subject_type *sublist; this is the current list of subjects that exist
				Common_parms *parms; the current values of the common parameters
				int N; the number of observations in **ts;
				RETURNS: double ssq; this value is needed for one of the parameters
				of the posterior distribution of model error variance
*********************************************************************/
/*********************************************************************/
/*VARIABLE DEFINITIONS
 i,j: generic counters
 *mean: mean concentration under current parameters and pulses
 ssq: sum of squared differences between observed log(concentration) and expected
 values
 *subject: current list of subjects and their qualities

 SUBROUTINES USED
 mean_concentration: found in birthdeath.c; sums across mean_contrib vectors for
 each pulse; returns this vector
 **************************************************************************/

double error_squared_l(double **ts, Subject_type *sublist, Common_parms *parms, int N) {
	int i, j;
	double *mean, ssq;
	Subject_type *subject;
	double *mean_concentration(Node_type *, Common_parms *, int, Node_type *, double);

	subject = sublist->succ;
	ssq = 0;
	for (j = 0; j < parms->numsub; j++){

		mean = mean_concentration(subject->driver, parms, N, subject->driver, subject->basehalf_l[0]);
		for (i = 0; i < N; i++)
			ssq += (ts[i][j + 1] - mean[i])*(ts[i][j + 1] - mean[i]);

		subject = subject->succ;
		free(mean);
	}

	return ssq;
}
