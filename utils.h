/********************************************/
/*	utils.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _MYUTILS_UTILS_H_
#define _MYUTILS_UTILS_H_

namespace myutils {
template<typename T>
void SWAP(T &a, T &b) {
	T c = a;
	a = b;
	b = c;
}

template<typename T>
T MIN(T a, T b) {
	return (a<b) ? a : b;
}

template<typename T>
T MAX(T a, T b) {
	return (a<b) ? b : a;
}
};

#endif // _MYUTILS_UTILS_H_