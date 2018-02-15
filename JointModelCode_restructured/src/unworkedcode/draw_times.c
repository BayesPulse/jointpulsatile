/*******************************************************************/
/*************************draw_times_trigger.c *********************/
/*******************************************************************/

/**fix here***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


void draw_times_trigger(Subject_type *sublist, Common_parms *parms_cox, double **ts,
	int N, unsigned long *seed, double v) {
	int i, j, k, l, ninterdiff1, ninterdiff2, ninterdiff3, num_parents, npulse;
	Subject_type *subject;
	Node_type *subnode, *node, *node2, *fnode;
	double alpha, plikelihood, new_loc, like_ratio, prior_ratio, tmp, current_time, *tmp_mean_contrib, lratio, lprop_ratio, tempdenomsum;
	double current_like, time_diff1, time_diff1_new, time_diff2, time_diff2_new, lkparent, kparent, lknew, knew;
	double rnorm(double, double, unsigned long *);
	double kiss(unsigned long *);
	void mean_contribution(Node_type *, double **, Common_parms *, int, double);
	double likelihood2(Node_type *, double **, Common_parms *, int, Node_type *, double);

	/*Allocate Memory*/

	tmp_mean_contrib = (double *)calloc(N, sizeof(double));

	/*Go to start of subject list*/
	subject = sublist->succ;
	npulse = subject->numnode_l;

	/*Index used for likelihood2*/
	parms_cox->subindex = 0;

	/*Loop thru subjects*/
	while (subject != NULL){
		k = 0;
		/*Go to first pulse for this subject*/
		subnode = subject->driver->succ;

		/*Loop thru pulses for this subject*/
		while (subnode != NULL) {

			/*Increase denominator of acceptance rate for time*/
			ntime_l++;

			/*Find current likelihood for this subject*/
			current_like = likelihood2(subject->driver, ts, parms_cox, N, subject->driver, subject->basehalf_l[0]);

			/*Compute proposal time*/
			new_loc = rnorm(subnode->time, sqrt(v), seed);

			/*Only proceed if our proposed time is reasonable */
			if (new_loc <= fitend && new_loc > fitstart) {

				ninterdiff1 = 0;
				ninterdiff2 = 0;
				ninterdiff3 = 0;
				node = subject->driver->succ;

				while (node != NULL && ninterdiff1 == 0){
					/*        printf("node->time %lf new_loc %lf\n",node->time,new_loc);*/
					if (subnode->time != node->time && fabs(node->time - new_loc) <= parms_cox->Rinter1) {
						ninterdiff1++;
					}
					else {
						if (subnode->time != node->time && fabs(node->time - new_loc) <= parms_cox->Rinter2) {
							ninterdiff2++;
						}
						if (subnode->time != node->time && fabs(subnode->time - node->time) <= parms_cox->Rinter2) {
							ninterdiff3++;
						}
					}
					node = node->succ;
				}
				/* printf("new_loc %lf node2->time %lf\n",new_loc,node2->time);
				 printf("ninterdiff1 %d\n",ninterdiff1);
				 printf("ninter2 %d\n",ninterdiff2);*/
				if (ninterdiff1 == 0) {
					/*conditional intensity for strauss/hardcore priors_cox*/
					lprop_ratio = (double)(ninterdiff2 - ninterdiff3) * log(parms_cox->gamma);

					/*component 1 of the conditional intensity for the parent given the data*/
					//lprop_ratio += exp(parms_cox->lrho) * (phi(fitend, subnode->time, sqrt(parms_cox->nu)) - phi(fitstart, subnode->time, sqrt(parms_cox->nu)) - phi(fitend, new_loc, sqrt(parms_cox->nu)) + phi(fitstart, new_loc, sqrt(parms_cox->nu)));

					/*ratio of proposal densities*/
					/*  lprop_ratio += log(phi(fitend, new_loc, pparents_sd) - phi(fitstart, new_loc, pparents_sd)) - log(phi(fitend, node2->time, pparents_sd) - phi(fitstart, node2->time, pparents_sd));*/

					/*product in the posterior cond. intensity */
					lratio = 0;
					k = 0;
					fnode = subject->response->succ;
					while (fnode != NULL){
						lkparent = -1 / (2 * parms_cox->nu) * ((fnode->time - subnode->time) * (fnode->time - subnode->time));
						kparent = exp(parms_cox->lrho) / sqrt(2 * 3.14159 * parms_cox->nu) * exp(lkparent);
						lknew = -1 / (2 * parms_cox->nu) * ((fnode->time - new_loc) * (fnode->time - new_loc));
						knew = exp(parms_cox->lrho) / sqrt(2 * 3.14159 * parms_cox->nu) * exp(lknew);
						tempdenomsum = fnode->lambda - kparent + knew;
						if (tempdenomsum == 0)
							tempdenomsum = 1e-300;
						lratio += log(tempdenomsum) - log(fnode->lambda);
						fnode = fnode->succ;
						k++;

					}            /*Save current time*/
					lprop_ratio += lratio;

					current_time = subnode->time;

					/*Set pulse's time to the proposal value*/
					subnode->time = new_loc;

					/*Save mean_contrib of this pulse*/
					for (i = 0; i < N; i++)
						tmp_mean_contrib[i] = subnode->mean_contrib[i];

					/*Recalculate mean_contrib of this pulse assuming proposed value*/
					mean_contribution(subnode, ts, parms_cox, N, subject->basehalf_l[1]);

					/*Calculate the likelihood for this subject under the proposed value*/
				plikelihood = likelihood2(subject->driver, ts, parms_cox, N, subject->driver, subject->basehalf_l[0]);
				

					like_ratio = plikelihood - current_like;

					/*Calculate log rho; set alpha equal to min(0,log rho)*/
					alpha = (0 < (tmp = (like_ratio + lprop_ratio))) ? 0 : tmp;

					/*If log U < log rho, we accept proposed value. Increase acceptance
					count by one*/
					if (log(kiss(seed)) < alpha) {
						atime_l++;
						fnode = subject->response->succ;
						while (fnode != NULL){

							fnode->lambda = skernel_1(fnode, subject->driver, parms_cox);
							fnode = fnode->succ;
						}
					}

					/*Otherwise, we reject, and we have to revert back to current values*/
					else {

						/*Set pulse's time back to current time*/
						subnode->time = current_time;

						/*Set mean_contrib of this pulse back to current value*/
						for (i = 0; i < N; i++)
							subnode->mean_contrib[i] = tmp_mean_contrib[i];
					}
				}
			}		/*End of if time is feasible statement*/

			/*Advance one pulse*/
			subnode = subnode->succ;
			k++;
		}/*End of loop through pulses*/

		/*Go to next subject*/
		subject = subject->succ;
		parms_cox->subindex++;

	}/*End of loop through subjects*/

	/*Free Memory*/
	free(tmp_mean_contrib);

}
