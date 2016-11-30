//Copyright 2016 Dmitry N. Razdoburdin.
//
//This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
//
//TG_3D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with TG_3D. If not, see <http://www.gnu.org/licenses/>.

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
