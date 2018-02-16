/*******************************************************************/
/*************************draw_fe_prior_a_var.c *****************************/
/*******************************************************************/

#include "draw_pop_mass_mean_SD.h"

void draw_pop_mass_mean_SD(Subject_type *sublist, Priors *priors, Common_parms *parms_f, unsigned long *seed, Hyper_priors *hyper) {
	/*declare varibles */
	int i, j, nsubj;
	double diff[2], **sigma_a_w, **sigma_a_w_inv;
	Subject_type *subject;
	/*declare functions */
	

	/* allocate memory*/
	sigma_a_w_inv = (double **)calloc(2, sizeof(double *));
	sigma_a_w = (double **)calloc(2, sizeof(double *));

	for (j = 0; j < 2; j++) {
		sigma_a_w_inv[j] = (double *)calloc(2, sizeof(double));
		sigma_a_w[j] = (double *)calloc(2, sizeof(double));
	}

	/* create the input matrix for the wishart */
	subject = sublist->succ;
	nsubj = parms_f->numsub;
	/*  priors->fe_mean_l[0] = 20.0;
	  priors->fe_mean_f[0] = 25.0;
	  subject->theta_l[0] = 21.0;
	  subject->theta_f[1] = 25.184;*/

	while (subject != NULL) {
		diff[0] = subject->theta_l[0] - priors->fe_mean_l[0];
		diff[1] = subject->theta_f[0] - priors->fe_mean_f[0];
		for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			sigma_a_w_inv[i][j] += diff[i] * diff[j];
		subject = subject->succ;

	}

	for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
		sigma_a_w_inv[i][j] += hyper->sig_a_inv[i][j];

	if (!cholesky_decomp(sigma_a_w_inv, 2)) {
		printf("S_inv matrix for full cond of sigma_a_inv is not PSD \n");
		exit(0);
	}

	sigma_a_w = cholesky_invert(2, sigma_a_w_inv);
	if (!cholesky_decomp(sigma_a_w, 2)) {
		printf("S matrix for full cond of sigma_a_inv is not PSK \n");
		exit(0);
	}
	/* acutually should not be decomposed*/
	rwishart(priors->fe_precision, sigma_a_w, 2, 4 + nsubj, seed,1);
	
	for (i = 0; i < 2; i++) {
		free(sigma_a_w[i]);
		free(sigma_a_w_inv[i]);
	}
	free(sigma_a_w);
	free(sigma_a_w_inv);
}

