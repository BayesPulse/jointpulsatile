/*******************************************************************/
/*************************gamma_code.c *****************************/
/*******************************************************************/

#include "jt_deconvolution_main.h"
#include "jt_birthdeath.h"
#include "jt_mcmc.h"


// the gamma distribution PDF
double gamm(double x)
{
	double ret = (1.000000000190015 +
		76.18009172947146 / (x + 1) +
		-86.50532032941677 / (x + 2) +
		24.01409824083091 / (x + 3) +
		-1.231739572450155 / (x + 4) +
		1.208650973866179e-3 / (x + 5) +
		-5.395239384953e-6 / (x + 6));

	return ret * sqrt(2 * M_PI) / x * pow(x + 5.5, x + .5) * exp(-x - 5.5);
}




// the gamma distribution PDF
double gammaPdf(double x, double a, double b)
{
	if (x <= 0 || a <= 0 || b <= 0)
		return 0.0;
	else
		return exp(-x * b) * pow(x, a - 1.0) * pow(b, a) / gamm(a);
}
