#include "functions.h"
int perturbation::boundary(parameters data, background bg)
{
	double x=-10;
	if (bg.rho_0[1]*Re_w[1]<pow(10,x) && bg.rho_0[data.Z-2]*Re_w[data.Z-2]<pow(10,x))
	{
		Re_w[0]=0;
		Re_w[data.Z-1]=0;
		return 0;
	}
	else
	{
		fprintf(stderr,"\nPerturbations come to the boundary! Lz=%.1lf\n",data.Lz);
		return -1;

	}
}

void background::get_boundary_type(char** s)
{
	sprintf((*s),"infinite");
}
