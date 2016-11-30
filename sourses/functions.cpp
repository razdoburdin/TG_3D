//Copyright 2016 Dmitry N. Razdoburdin.

//This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
//TG_3D is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with TG_3D; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <sys/stat.h>
#include <dirent.h>
#include <sys/time.h>

#include "classes.h"

double coordinate(double j, double dz, double Lz)//calculation of coordinate 
{
	return -Lz+j*dz;
}

void chek_result_floder()
{
	DIR* result;
	char* name;
	name=new char[7];

	snprintf(name,7,"result");
	if ((result=opendir(name))==NULL)
	{
		mkdir(name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	else
	{
		closedir(result);
	}
	delete [] name;
}

void timeformat(struct timeval time_start,struct timeval time_end,char **out)
{
	short houres,minutes;
	double seconds;
	double time_seconds;

	time_seconds=time_end.tv_sec-time_start.tv_sec+0.000001*(time_end.tv_usec-time_start.tv_usec);
	minutes=(int) (time_seconds/60);
	houres=minutes/60;
	minutes=minutes-60*houres;
	seconds=time_seconds-3600*houres-60*minutes;

	short size=13;
	(*out)=new char[size];
	short i=0;
	if (houres<10)
	{
		i+=snprintf((*out)+i,size-i,"0");
	}
	i+=snprintf((*out)+i,size-i,"%d:",houres);
	if (minutes<10)
	{
		i+=snprintf((*out)+i,size-i,"0");
	}
	i+=snprintf((*out)+i,size-i,"%d:",minutes);
	if (seconds<10)
	{
		i+=snprintf((*out)+i,size-i,"0");
	}
	i+=snprintf((*out)+i,size-i,"%.2lf",seconds);
}

void progress(int persent)
{
	fprintf(stderr,"\b\b\b\b");
	if (persent<10)
	{
		fprintf(stderr,"  %d",persent);
		fprintf(stderr,"%%");
	}
	if ((persent>=10) && (persent<100))
	{
		fprintf(stderr," %d",persent);
		fprintf(stderr,"%%");
	}
	if (persent==100)
	{
		fprintf(stderr,"%d",persent);
		fprintf(stderr,"%%");
	}
}
