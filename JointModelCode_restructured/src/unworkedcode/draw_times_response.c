/*******************************************************************/
/*************************draw_times_response.c *****************************/
/*******************************************************************/

/***fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"



void draw_times_response(Subject_type *sublist, Common_parms *parms, Common_parms *parms_cox, double **ts,
	int N, unsigned long *seed, double v) {
	int i, ninter1,flag,k;
	Subject_type *subject;
	Node_type *subnode, *dnode,*node;
	double alpha, plikelihood, ptime, like_ratio, prior_ratio, tmp, current_time, *tmp_mean_contrib, condintnew, condintold, prop_ratio;
	double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);
	double phi(double, double, double);
	double skernel_1(Node_type *, Node_type *, Common_parms *);

	/*Allocate Memory*/
	tmp_mean_contrib = (double *)calloc(N, sizeof(double));

	/*Go to start of subject list*/
	subject = sublist->succ;

	/*Index used for likelihood2*/
	parms->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){
		k = 0;
		/*Go to first pulse for this subject*/
		subnode = subject->response->succ;

		/*Loop thru pulses for this subject*/
		while (subnode != NULL) {

			/*Increase denominator of acceptance rate for time*/
			ntime_f++;

			/*Find current likelihood for this subject*/
			current_like = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

			/*Compute proposal time*/
			ninter1 = 0;
			flag = 1;
			/*simulate a pulse position uniformly */
			/*if the new pulse location is too close to the old one, try it again*/
			/*like an hard core process*/
			while (flag == 1)
			{

				ptime = rnorm(subnode->time, v, seed);
				node = subject->response->succ;
				ninter1 = 0;
				while (node != NULL)
				{

					if (fabs(node->time - ptime) <= 0) {
						ninter1++;
					}
					node = node->succ;
				}
				if (ninter1 == 0) {
					flag = 0;


				}

			}

			/*Only proceed if our proposed time is reasonable */
			if (ptime <= fitend && ptime > fitstart) {

				prior_ratio = 0;

				condintnew = 0;
				condintold = 0;

				dnode = subject->driver->succ;
				while (dnode != NULL) {
					condintnew += exp(parms_cox->lrho) / sqrt(2 * 3.14159*parms_cox->nu) * exp(-1 / (2 * parms_cox->nu) * (ptime - dnode->time) * (ptime - dnode->time));
					condintold += exp(parms_cox->lrho) / sqrt(2 * 3.14159*parms_cox->nu) * exp(-1 / (2 * parms_cox->nu) * (subnode->time - dnode->time) * (subnode->time - dnode->time));
					dnode = dnode->succ;
				}

				if (condintnew == 0)
				{
					condintnew = 1e-300;
				}
				if (condintold == 0)
				{
					condintold = 1e-300;
				}

				/*Combine it all for the prior ratio*/
				condintnew = log(condintnew);
				condintold = log(condintold);

				prior_ratio = condintnew - condintold;
				prop_ratio = 0;
				/*prop_ratio = log(1 - phi(0, ptime, sqrt(v))) - log(1 - phi(0, subnode->time, sqrt(v)));*/

				/*Save current time*/
				current_time = subnode->time;

				/*Set pulse's time to the proposal value*/
				subnode->time = ptime;

				/*Save mean_contrib of this pulse*/
				for (i = 0; i < N; i++)
					tmp_mean_contrib[i] = subnode->mean_contrib[i];

				/*Recalculate mean_contrib of this pulse assuming proposed value*/
				mean_contribution(subnode, ts, parms, N, subject->basehalf_f[1]);

				/*Calculate the likelihood for this subject under the proposed value*/
				plikelihood = likelihood2(subject->response, ts, parms, N, subject->response, subject->basehalf_f[0]);

				/*Calculate the likelihood ratio*/
				like_ratio = plikelihood - current_like;

				/*Calculate log rho; set alpha equal to min(0,log rho)*/
				alpha = (0 < (tmp = (prior_ratio + like_ratio))) ? 0 : tmp;

				/*If log U < log rho, we accept proposed value. Increase acceptance
				count by one*/
				if (log(kiss(seed)) < alpha) {
					atime_f++;
					subnode->lambda = skernel_1(subnode, subject->driver, parms_cox);

				}

				/*Otherwise, we reject, and we have to revert back to current values*/
				else {

					/*Set pulse's time back to current time*/
					subnode->time = current_time;

					/*Set mean_contrib of this pulse back to current value*/
					for (i = 0; i < N; i++)
						subnode->mean_contrib[i] = tmp_mean_contrib[i];
				}
			} /*End of if time is feasible statement*/

			/*Advance one pulse*/
			subnode = subnode->succ;
			k++;
		}/*End of loop through pulses*/

		/*Go to next subject*/
		subject = subject->succ;
		parms->subindex++;

	}/*End of loop through subjects*/

	/*Free Memory*/
	free(tmp_mean_contrib);

}
