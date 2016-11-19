#include <stdio.h>
#include <math.h>

double rho_bg(double z, double n)
{
	return exp(-0.5*z*z);
}

double a2_bg(double z, double n)
{
	return 1;
}

double dp_0__dz__rho_0_bg(double z)
{
	return -z;
}

double drho_0__dz__rho_0_bg(double z,double n)
{
	return -z;
}

void get_bg_type(char** s)
{
	sprintf((*s),"isothermal");
}
