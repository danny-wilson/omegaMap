/********************************************/
/*	myerror.h 17th November 2005			*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_ERROR_H
#define _MYUTILS_ERROR_H

#include <stdio.h>
#include <stdlib.h>

namespace myutils
{
	inline void error(const char* error_text)
	{
		printf("ERROR: ");
		printf("%s\n", error_text);
		exit(13);
	}

	inline void warning(const char* warning_text)
	{
		printf("WARNING: ");
		printf("%s\n", warning_text);
		return;
	}

};

#endif