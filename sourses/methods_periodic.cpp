//Copyright 2016 Dmitry N. Razdoburdin.

/*This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
TG_3D is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
License for more details. You should have received a copy of the GNU General Public License along with TG_3D; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA*/
#include <stdio.h>

#include "functions.h"

int perturbation::boundary(parameters data, background bg)
{
	Re_w[0]=(Re_w[1]+Re_w[data.Z-2])*0.5;
	Re_w[data.Z-1]=(Re_w[1]+Re_w[data.Z-2])*0.5;
	return 0;
}

void background::get_boundary_type(char** s)
{
	sprintf((*s),"pereodic");
}
