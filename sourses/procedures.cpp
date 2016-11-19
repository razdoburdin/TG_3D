#include <stdio.h>
#include "classes.h"

void G_kx(parameters data, background bg, double kx_min, double kx_max,double d_kx)
{
	data.kx=kx_min;
	optimal* singular_vectors;
	singular_vectors=(optimal*) malloc (1*sizeof(optimal));
	
	while (data.kx<=kx_max)
	{
		fprintf(stderr,"kx=%lf\n",data.kx);
		singular_vectors[0]=optimal(data,bg,singular_vectors,0);
		singular_vectors[0].write_G(data,bg,1);
		data.kx+=d_kx;
	}
	free(singular_vectors);
	singular_vectors=NULL;
}
