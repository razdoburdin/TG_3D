#include "functions.h"

int perturbation::boundary(parameters data, background bg)
{
	Re_w[0]=Re_w[1];
	Re_w[data.Z-1]=Re_w[data.Z-2];
	return 0;
}

void background::get_boundary_type(char** s)
{
	sprintf((*s),"second");
}
