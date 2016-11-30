//Copyright 2016 Dmitry N. Razdoburdin.

/*This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
TG_3D is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with TG_3D; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA*/
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
