/*******************************************************************/
/*************************write_mcmc_output.c *****************************/
/*******************************************************************/

/***FIX HERE***/
#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


/**NEEDS A FUNCTION WRAPPER**/
/**CODE WILL NOT RUN AT THIS POINT**/
		/***************************************************************/
	
		/* Every 50th iteration, print estimates to files*/
		if (!(i%NN)) {
			subject = sublist->succ;
			while (subject != NULL){
				num_node = 0;
				new_node = subject->driver->succ;
				while (new_node != NULL){
					num_node++;
					new_node = new_node->succ;
				}
				/*Print subject specific parameters to c1sk.out files*/

				fprintf(subject->csub_l, "%d %d %lf %lf %lf %lf %lf %lf \n", i, num_node, subject->theta_l[0], subject->theta_l[1], subject->basehalf_l[0], subject->basehalf_l[1], parms_l->rho, parms_l->nu);
				fflush(subject->csub_l);
				new_node = subject->driver->succ;
				num_node2 = 0;
				/*Print pulse specific parameters to p1sk.out files*/
				while (new_node != NULL) {
					fprintf(subject->psub_l, "%d %d %d %lf %lf %lf %lf %lf\n", i, num_node, num_node2, new_node->theta[0], new_node->theta[1], new_node->time, new_node->eta[0], new_node->eta[1]);
					num_node2++;
					new_node = new_node->succ;
				}

				fflush(subject->psub_l);
				num_node = 0;
				new_node = subject->response->succ;
				while (new_node != NULL){
					num_node++;
					new_node = new_node->succ;
				}
				/*Print subject specific parameters to c1sk.out files*/
				fprintf(subject->csub_f, "%d %d %lf %lf %lf %lf\n", i, num_node, subject->theta_f[0], subject->theta_f[1], subject->basehalf_f[0], subject->basehalf_f[1]);
				fflush(subject->csub_f);
				new_node = subject->response->succ;
				num_node2 = 0;
				/*Print pulse specific parameters to p1sk.out files*/
				while (new_node != NULL) {
					fprintf(subject->psub_f, "%d %d %d %lf %lf %lf %lf %lf\n", i, num_node, num_node2, new_node->theta[0], new_node->theta[1], new_node->time, new_node->eta[0], new_node->eta[1]);
					num_node2++;
					new_node = new_node->succ;
				}

				fflush(subject->psub_f);

				
				subject = subject->succ;


			}
			for (k = 0; k < 2; k++)
			for (j = 0; j < 2; j++)
				sigma_p[k][j] = priors->fe_precision[k][j];

			if (!cholesky_decomp(sigma_p, 2)){
				printf("not PSD matrix A\n");
				exit(0);
			}
			priors->re_var = cholesky_invert(2, sigma_p);
			/*Print common parameters to c1.out file*/
			fprintf(common_l, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf %lf %lf   %lf %lf\n", i, parms_l->nu,parms_l->rho,priors->fe_mean_l[0], priors->fe_mean_l[1], priors->re_var[0][0], priors->re_var[0][1], priors->fe_precision_wl, parms_l->re_precision[0], parms_l->re_precision[1],priors->meanbh_l[0], priors->meanbh_l[1], priors->varbh_l[0], priors->varbh_l[1], parms_l->sigma);
			fflush(common_l);

			fprintf(common_f, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, priors->fe_mean_f[0], priors->fe_mean_f[1], priors->re_var[1][1], priors->re_var[0][1],priors->fe_precision_wl, parms_f->re_precision[0], parms_f->re_precision[1], priors->meanbh_f[0], priors->meanbh_f[1], priors->varbh_f[0], priors->varbh_f[1], parms_f->sigma);
			fflush(common_f);
			fflush(stdout);
		}

		/*Print acceptance rates*/
		if (!(i%NNN)) {
			printf("\n\n");
			printf("iter = %d \n", i);
			printf("lh: mu_b %.2lf mu_h %.4lf mu_a %.2lf mu_w %.2lf  v %.4le\n",
				priors->meanbh_l[0], priors->meanbh_l[1], priors->fe_mean_l[0], priors->fe_mean_l[1], parms_l->sigma);
			printf("fsh: mu_b %.2lf mu_h %.4lf mu_a %.2lf mu_w %.2lf  v %.4le\n",
				priors->meanbh_f[0], priors->meanbh_f[1], priors->fe_mean_f[0], priors->fe_mean_f[1], parms_f->sigma);
			subject = sublist->succ;
			printf("lh: pct s_ma = %.2lf  pct A_ki = %.2lf , eta_m=%.2lf, eta_w=%.2lf\n",
				(double)afepmv_l / (double)nfepmv_l, (double)arem_l / (double)nrem_l, (double)ae_ml / (double)ne_ml, (double)ae_wl / (double)ne_wl);
			printf("pct s_mw = %.2lf pct s_w = %.2lf pct s2p_ki = %.2lf \n",
				(double)afepwv_l / (double)nfepwv_l, (double)afewv_l / (double)nfewv_l, (double)arew_l / (double)nrew_l);
			printf("pct s_b = %.2lf pct s_h = %.2lf pct B-HL = %.2lf pct time = %.2lf \n",



				(double)afebv_l / (double)nfebv_l, (double)afehv_l / (double)nfehv_l, (double)adelta_l / (double)ndelta_l, (double)atime_l / (double)ntime_l);
			printf("fsh: pct s_ma = %.2lf  pct A_ki = %.2lf , eta_m=%.2lf, eta_w=%.2lf\n",
				(double)afepmv_f / (double)nfepmv_f, (double)arem_f / (double)nrem_f, (double)ae_mf / (double)ne_mf, (double)ae_wf / (double)ne_wf);
			printf("pct s_mw = %.2lf pct s_w = %.2lf pct s2p_ki = %.2lf \n",
				(double)afepwv_f / (double)nfepwv_f, (double)afewv_f / (double)nfewv_f, (double)arew_f / (double)nrew_f);
			printf("pct s_b = %.2lf pct s_h = %.2lf pct B-HL = %.2lf pct time = %.2lf \n \n \n ",
				(double)afebv_f / (double)nfebv_f, (double)afehv_f / (double)nfehv_f, (double)adelta_f / (double)ndelta_f, (double)atime_f / (double)ntime_f);
			
			printf("%d    \n", sublist->succ->numnode_l);
			printf("%d    \n", sublist->succ->numnode_f);

			tdriver = sublist->succ->driver->succ;
			while (tdriver != NULL){
			printf("%lf  %lf  %lf   \n", tdriver->time,tdriver->theta[0],tdriver->theta[1]);
			tdriver = tdriver->succ;
		}

	/*	tresponse = tsublist->succ->response->succ;
			while (tresponse != NULL){
				printf("%lf  %lf  %lf   \n", tresponse->time, tresponse->theta[0], tresponse->theta[1]);		
				tresponse = tresponse->succ;
			}*/


			fflush(stdout);
		}

