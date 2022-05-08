/* File: crelib.cpp
 * Author: CRE
 * Last Edited: Fri Oct 28 16:18:03 2016
 */

#include "crelib.h"
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>

void cre::die(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vfprintf(stderr,Format,argptr);
	fprintf(stderr, "\n");
	exit(1);
}

void cre::warn(const char * Format, ...)
{
	va_list argptr;
	va_start(argptr, Format);
	vfprintf(stderr,Format,argptr);
	fprintf(stderr, "\n");
}

unsigned long long LastTime=0;

void cre::updateTime(const char * What, const char * OtherMessage, FILE* LogFile)
{
	unsigned long long ThisTime;
	if (LastTime==0)
	{
		LastTime=getTimeInMuse();
		fprintf(LogFile, "Performed time initial.\n");
	}
	ThisTime=getTimeInMuse();
	fprintf(LogFile,"%s cost %lld ms.",What, (ThisTime-LastTime)/1000);
	LastTime=ThisTime;
	if (OtherMessage!=NULL) fprintf(LogFile," %s",OtherMessage);
	fprintf(LogFile,"\n");
}

unsigned long long cre::getTimeInMuse()
{
	timeval VTime;
	gettimeofday(&VTime, NULL);
	return VTime.tv_sec*1000000ull+VTime.tv_usec;
}
