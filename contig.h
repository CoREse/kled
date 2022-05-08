#ifndef KLED_CONTIG_H
#define KLED_CONTIG_H

#include <string>

struct Contig
{
	const unsigned int Size;
	std::string Name;
	Contig(const char *Name,unsigned int Size);
};

#endif