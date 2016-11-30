//Copyright 2016 Dmitry N. Razdoburdin.

/*This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
TG_3D is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with TG_3D; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA*/

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
