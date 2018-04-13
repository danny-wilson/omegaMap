/************************************/
/*	main.cpp 17th November 2005		*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*	www.danielwilson.me.uk			*/
/************************************/

#include "omegaMap.h"
#include <random.h>

int main(const int argc, const char* argv[]) {
	if(argc<2) error("SYNTAX: ini-filename [options...]");
	omegaMap oM;
	Random ran;
	oM.go(argc-1,argv+1,ran,argv[1]);
	return 0;
}

