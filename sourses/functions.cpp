#include <math.h>
#include <stdio.h>
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
