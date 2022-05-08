#pragma once
#ifndef CRELIB_DEFINE_H
#define CRELIB_DEFINE_H

#ifndef myalloc
#define myalloc(n,T) (T*)malloc(n*sizeof(T))
#endif

#ifndef STRING_HELPER
#define STRING_HELPER(X) #X
#endif

#ifndef STRING
#define STRING(X) STRING_HELPER(X)
#endif

#define MIN(a,b) (a>b?b:a)
#define MAX(a,b) (a>b?a:b)

namespace cre
{
	typedef unsigned int uint;
	typedef unsigned long long uint64;
	typedef long long int64;
	typedef unsigned char uchar;
}
#endif
