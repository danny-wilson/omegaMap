/************************************/
/*	decode.h 23rd February 2005		*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*	www.danielwilson.me.uk			*/
/************************************/

#ifndef _OMEGAMAP_ANALYSE_H_
#define _OMEGAMAP_ANALYSE_H_

#include <omegaMap.h>

class omegaMapAnalyse : public omegaMapBase {
public:
	bool analyseInitialized;

	omegaMapAnalyse() {
		analyseInitialized = false;
	}

	omegaMapAnalyse& go(const char* datafile, const char* inifile);
	omegaMapAnalyse& open(const char* datafile);	
	omegaMapAnalyse& initialize(const char* inifile);
	omegaMapAnalyse& initialize(const char* outfile_in, Vector<int> &fields_in, const int thinning_in);
	omegaMapAnalyse& recycle();
	omegaMapAnalyse& extract();
//	omegaMapAnalyse& extract(void (*f)(void*),void *farg);
	omegaMapAnalyse& iterate(const int iter_in);

	/* moves */

	omegaMapAnalyse& change_oBlock();
	omegaMapAnalyse& extend_oBlock();
	omegaMapAnalyse& split_oBlock();
	omegaMapAnalyse& merge_oBlock();
	omegaMapAnalyse& change_rBlock();
	omegaMapAnalyse& extend_rBlock();
	omegaMapAnalyse& split_rBlock();
	omegaMapAnalyse& merge_rBlock();
	omegaMapAnalyse& change_mu();
	omegaMapAnalyse& change_kappa();
	omegaMapAnalyse& change_indelLambda();
};

#endif // _OMEGAMAP_ANALYSE_H_