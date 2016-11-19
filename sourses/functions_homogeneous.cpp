#include <stdio.h>
#include <math.h>

double rho_bg(double z, double n)
{
	return pow(0.5/(n+1),n);
}

double a2_bg(double z, double n)
{
	return 0.5/n;
}

double dp_0__dz__rho_0_bg(double z)
{
	return 0;
}

double drho_0__dz__rho_0_bg(double z,double n)
{
	return 0;
}

void get_bg_type(char** s)
{
	sprintf((*s),"homogeneous");
}
