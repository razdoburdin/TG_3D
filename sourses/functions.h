//Copyright 2016 Dmitry N. Razdoburdin.

//This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
//TG_3D is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with TG_3D; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "classes.h"
#include <stdio.h>
#include <math.h>

double coordinate(double j, double dz, double Lz);

double rho_bg(double z, double n);

double a2_bg(double z, double n);

double dp_0__dz__rho_0_bg(double z);

double drho_0__dz__rho_0_bg(double z,double n);

void chek_result_floder();

void timeformat(struct timeval time_start,struct timeval time_end,char **out);

void progress(int persent);

void get_bg_type(char** s);
