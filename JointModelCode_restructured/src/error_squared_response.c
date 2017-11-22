/*******************************************************************/
/*************************error_squared_response.c *****************************/
/*******************************************************************/

/**fix here**/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


double error_squared_response(double **ts, Subject_type *sublist, Common_parms *parms, int N) {
	int i, j;
	double *mean, ssq;
	Subject_type *subject;
	double *mean_concentration(Node_type *, Common_parms *, int, Node_type *, double);

	subject = sublist->succ;
	ssq = 0;
	for (j = 0; j < parms->numsub; j++){

		mean = mean_concentration(subject->response, parms, N, subject->response, subject->basehalf_f[0]);
		for (i = 0; i < N; i++)
			ssq += (ts[i][j + 1] - mean[i])*(ts[i][j + 1] - mean[i]);

		subject = subject->succ;
		free(mean);
	}

	return ssq;
}
