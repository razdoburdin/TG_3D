#include <stdio.h>
#include <math.h>

double rho_bg(double z, double n)
{
	return pow(0.5/(n+1)*(1-z*z),n);
}

double a2_bg(double z, double n)
{
	return 0.5/n*(1-z*z);
}

double dp_0__dz__rho_0_bg(double z)
{
	return -z;
}

double drho_0__dz__rho_0_bg(double z,double n)
{
	return -2*n*z/(1-z*z);
}

void get_bg_type(char** s)
{
	sprintf((*s),"polytropic");
}
