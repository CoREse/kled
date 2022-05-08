#pragma once
#ifndef CRELIB_H
#define CRELIB_H
#include <stdio.h>
#include <stdlib.h>
#include "defines.h"
namespace cre
{
	void die(const char * Format=NULL, ...);
	void updateTime(const char * What, const char * OtherMessage=NULL, FILE * LogFile=stderr);
	unsigned long long getTimeInMuse();//get time in microseconds. Muse no minasan, daisuki!^.^

	void warn(const char * Format, ...);
}
#endif
