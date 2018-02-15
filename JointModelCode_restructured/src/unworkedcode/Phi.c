
double Phi(double y)
{
	/* Returns standard normal distribution function evaluated at y */
	/* Has absolute error of order 10^(-7) */
	/* Uses approximation to erfc from Numerical Recipes in C */
	double t, z, ans, x, inv_sqrt2 = 0.7071067811865475;
	x = y*inv_sqrt2;
	z = fabs(x);
	t = 1.0 / (1.0 + .5*z);
	ans = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
		t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
		t*(-0.82215223 + t*0.17087277)))))))));
	if (!(x < 0))
		return 1.0 - 0.5*ans;
	else
		return 0.5*ans;
}

double phi(double y, double mu, double s)
{
    /* Returns normal distribution function evaluated at y */
    /* Returns standard normal distribution function evaluated at y */
    /* Has absolute error of order 10^(-7) */
    /* Uses approximation to erfc from Numerical Recipes in C */
    double t, z, ans, x;
    
    x = (y - mu) / (1.4142135623730951*s);
    z = fabs(x);
    t = 1.0 / (1.0 + .5*z);
    ans = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(0.37409196 + t*(0.09678418 +
                                                                       t*(-0.18628806 + t*(0.27886807 + t*(-1.13520398 + t*(1.48851587 +
                                                                                                                            t*(-0.82215223 + t*0.17087277)))))))));
    if (x >= 0){
        return 1.0 - 0.5*ans;
    }
    else{
        return 0.5*ans;
    }
}
