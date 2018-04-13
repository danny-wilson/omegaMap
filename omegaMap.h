/************************************/
/*	omegaMap.h 23rd February 2005	*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*  www.danielwilson.me.uk			*/
/************************************/

#ifndef _OMEGA_MAP_H_
#define _OMEGA_MAP_H_

#include <myutils.h>
#include <lotri_matrix.h>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <limits>

using namespace myutils;
using namespace std;

class oMatrix {
/* observation matrix associated with an omega block */
public:
	LowerTriangularMatrix< Vector<double> > *gamma;	/* probability of observing pair (i,j) for k haplotypes */
	LowerTriangularMatrix< Vector<double> > *gamma2;/* probability of observing pair (i,j) for k haplotypes under indel model */
	Matrix<double> R;								/* eigenvectors for the mutation rate matrix */
	Vector<double> lambda;							/* eigenvalues for the mutation rate matrix */
	double mu, kappa, omega;						/* parameters of the mutation rate matrix */
};

class oBlock {
/* omega block: defines the continguous sites that share the same omega */
public:
	int start, end;
	oBlock *_5prime, *_3prime;
	oMatrix *oMat;
};

class rBlock {
/* recombination block: defines the contiguous sites that share the same rho */
public:
	int start, end;
	rBlock *_5prime, *_3prime;
	double rho;
};

/* oEvent and oEventStart are derived from base classes called oEventBinary and oEventStartBinary
   so that istream& operator>>() and ostream& operator<<() functions can be overloaded separately for
   each type. When it comes to reading/writing, the class can then be typecast when the <</>> operator
   is called in order to specify which version of the overloaded function to call					*/
class oEventBinary {
/* records accepted and rejected moves in the MCMC */
public:
	enum types {DEFAULT,START,CHANGE_OBLOCK,EXTEND_OBLOCK,SPLIT_OBLOCK,MERGE_OBLOCK,
		CHANGE_MU,CHANGE_KAPPA,CHANGE_RBLOCK,EXTEND_RBLOCK,SPLIT_RBLOCK,MERGE_RBLOCK,CHANGE_INDEL,CHANGE_ORDER} type;
	bool accepted;
	double param[2];
	double likelihood;
	double alpha;
};
class oEvent : public oEventBinary {};

class oEventStartBinary : public oEvent {
/* records the initial state of the MCMC */
public:
	int nblocks;
	int nrblocks;
	double mu,kappa,indelLambda;
	Vector<int> oBlockStart;
	Vector<int> rBlockStart;
	Vector<double> omega;
	Vector<double> rho;
};
class oEventStart : public oEventStartBinary {};

class omegaMapBase {
/* a skeleton version of class omegaMap for reconstructing the chain from a datafile */
/*public:
	string omegaVersion;
	int n;
	int L;
	int niter;
	int seed;
	Vector<double> pi;
	int norders;
	double muPrior,kappaPrior,omegaPrior,rhoPrior;
	double oBlockPrior;
	double rBlockPrior;
	double oSplitMergeProbability;
	double rSplitMergeProbability;
	oEventStart eventStart;
	Vector<oEvent> event;
public:
	omegaMapData& inputFromText(const char* filename);*/

public:	/* Member variables */

	/* Settings */
	string omegaVersion;			// current version of omegaMap
	bool initialized;				// indicates whether the class has been initialized
	bool coutput;					// print counter to screen or not
	string hapfile;					// filename for FASTA data
	bool efficientLikelihood;		// turn on efficient likelihood calculation
	int stickyDebugEval;			// site at which inefficient likelihood is evaluated
//	double muPrior;					// rate parameter for the exponential prior on mu
//	double kappaPrior;				//  "       "      "   "       "         "  on kappa
//	double omegaPrior;				//  "       "      "   "       "         "  on omega
//	double rhoPrior;				//  "       "      "   "       "         "  on rho
//	double indelLambdaPrior;		//  "       "      "   "       "         "  on indelLambda
	double oBlockPrior;				// probability parameter for the prior on omega block length
	double oBlockPriorRatio;
	double rBlockPrior;				// probability parameter for the prior on recombination block length
	double rBlockPriorRatio;
	bool writeBinary;
	Vector<double> prMOVE;			// MCMC move proposal probabilities. Indexed by:
									//		2:	CHANGE_OBLOCK
									//		3:	EXTEND_OBLOCK
									//		4:	SPLITMERGE_OBLOCK
									//		5:	-----------------
									//		6:	CHANGE_MU
									//		7:	CHANGE_KAPPA
									//		8:	CHANGE_RBLOCK
									//		9:	EXTEND_RBLOCK
									//		10:	SPLITMERGE_RBLOCK
									//		11:	-----------------
									//		12:	CHANGE_INDEL
									//		13:	CHANGE_ORDER
	string omegaModel;				// constant, variable or independent
	string rhoModel;				// constant, variable or independent

	enum Priors {
		none,fixed,uniform,improper_uniform,exponential,gamma,exponential_ratio,inverse,improper_inverse,log_normal
	} muPrior,kaPrior,omPrior,rhPrior,inPrior;
	double muPriorParam[2],kaPriorParam[2],omPriorParam[2],rhPriorParam[2],inPriorParam[2];

	/* Parameters */
	int n;							// number of haplotypes
	int L;							// length (in codons) of haplotypes
	int norders;					// number of orderings for PAC likelihood
	Vector<double> pi;				// equilibrium frequencies of the codons
	int niter;						// total number of steps in the MCMC

	/* Storage */
	Random *ran;					// pointer to random number generator
	map<char,int> baseToInt;		// converts TUCAG- to 112345
	map<int,char> intToBase;		// converts 012345 to NTCAG-

	DNA dna;						// an internal copy of the haplotypes
	Vector< Vector<int> > codon;	// an internal conversion of haplotypes from nucleotides into codons
	Matrix< Vector<int>* > H;		// use this to access haplotypes *(H[ordering][hap])

	Vector<double> sqrtPi;			// square root of the codon equilibrium frequencies
	Matrix<double> pamlWork;		// workspace for paml's diagonalization routine
	Vector<oBlock> _block;			// each position has a potential _block
	Vector<oBlock*> block;			// for each position points to the active _block
	Vector<oMatrix*> oMatTemp;		// storage space for rate matrices for proposed moves
	Vector<rBlock> _rblock;			// each inter-position has a potential _rblock
	Vector<rBlock*> rblock;			// for each inter-position points to the active _rblock

	int nblocks;					// current number of omega blocks
	int nrblocks;					// current number of recombination blocks
	int alphaMargin;				// stores the 5'-most position at which alpha is not up to date
	int betaMargin;					// stores the 3'-most position at which  beta is not up to date
	double ****alpha;				// storage for the forward algorithm
	double ****beta;				// storage for the backward algorithm
//	Vector<double> rho;				// recombination map (set to size L-1)
	Vector<double> PAC;				// storage for the PAC likelihoods (set to size norders)

	Vector<bool> indel;				// indicates that the indel model is to be used at that site
	int maxCodon;					// 61 if no indels, 62 otherwise
	double indelLambda;				// current value of the indel rate parameter

	int iter;						// current step in the MCMC
	double oldLikelihood;			// likelihood of the last accepted move in the MCMC chain
	oEventStart eventStart;			// records the initial state of the MCMC
	Vector<oEvent> event;			// records the accepted and rejected changes in the MCMC
	double memTime;					// records the time taken for memory allocation (in seconds)
	double compTime;				// records the time taken for the MCMC chain (in minutes)

	/* Output */
	string datafile;				// filename for encoded data output file
	string outfile;					// filename for text output file
	Vector<int> field;				// stores the variables to output to a text file
	int thinning;					// stores the thinning interval for the text file

	/* Misc */
	static const double SQRT2PI;

public:
	omegaMapBase(); // defined in omegaMap.cpp
	/* the following are defined in io.cpp */
	omegaMapBase& headings(ofstream& out);
	omegaMapBase& record(ofstream& out, oEvent& event);
	omegaMapBase& headings(ofstream& out, ostream& omanip(ostream&));
	omegaMapBase& record(ofstream& out, oEvent& event, ostream& omanip(ostream&));
	double thetaS();
};

/* functions for storing and retrieving omegaMapBase classes, defined in io.cpp */
ostream& tab(ostream& out);
ostream& ENDCONTROLWIZARD(ostream& out);
ostream& operator<<(ostream& out, omegaMapBase& oM);
istream& operator>>(istream& in,  omegaMapBase& oM);
ostream& operator<<(ostream& out, oEventBinary& event);
istream& operator>>(istream& in,  oEventBinary& event);
ostream& operator<<(ostream& out, oEvent& event);
istream& operator>>(istream& in,  oEvent& event);
ostream& operator<<(ostream& out, oEventStartBinary& event);
istream& operator>>(istream& in,  oEventStartBinary& event);
ostream& operator<<(ostream& out, oEventStart& event);
istream& operator>>(istream& in,  oEventStart& event);

class omegaMap : public omegaMapBase {
public:	/* Member functions */

	/***************************/
	/* Initialization routines */		/*  omegaMap.cpp  */
	/***************************/

	omegaMap();
	omegaMap& initialize(const char* inifile, Random &r);
	omegaMap& initialize(const int argc, const char* argv[], const char *inifile, Random &r);
	/* For maximum likelihood estimation - simplified initialization */
	omegaMap& mle_initialize(const int argc, const char* argv[], const char *inifile, Random &r);
	~omegaMap();
	omegaMap& dnaToCodons();
	string codonToTriplet(const int a);
	string codonToTriplet61(const int a);
	int tripletToCodon(string &tri);
	inline int removeStopCodons(const int a);
	//omegaMap& setOrders();
	vector<int> setOrders();
	omegaMap& setOrders(vector<int> &orders);

	/***************************/
	/*     Output routines     */		/*     io.cpp     */
	/***************************/

	omegaMap& outputToText(const char* filename);

	/***************************/
	/*   Likelihood routines   */		/* likelihood.cpp */
	/***************************/

	omegaMap& diagonalizeMatrix(oMatrix &oM, const double mu, const double kappa, const double omega);
	omegaMap& buildMatrixA(Matrix<double> &A, const double mu, const double kappa, const double omega);
	double mutProb(int site, int nhaps, int i, int j);		/* it is assumed that i<=j */
	omegaMap& forward(int _5prime, int _3prime);			/* updates alpha from _5prime to _3prime inclusive */
	omegaMap& backward(int _5prime, int _3prime);			/* updates beta from _3prime to _5prime inclusive */
	double likelihood(int pos);
	double star_likelihood();

	/***************************/
	/**     MCMC routines     **/		/*    mcmc.cpp    */
	/***************************/

	omegaMap& go(Random &r, const char* inifile);
	omegaMap& go(const int argc, const char* argv[], Random &r, const char* inifile);
	omegaMap& debug();
	omegaMap& propose();
	omegaMap& redraw();
	
	/***************************/
	/**  omega MCMC routines  **/		/*    mcmc.cpp    */
	/***************************/

	double omegaSplitRate(const int nblo);					/* relative rate at which splits are proposed */
	double omegaMergeRate(const int nblo);					/* relative rate at which merges are proposed */
	omegaMap& propose_change_oBlock();						/* Proposal distribution for changing a block's omega */
	omegaMap& propose_extend_oBlock();						/* Proposal distribution for extending a block 5' or 3' */
	omegaMap& propose_split_oBlock();						/* Proposal distribution for splitting an existing block */
	omegaMap& propose_merge_oBlock();						/* Proposal distribution for merging two adjacent blocks */

	/***************************/
	/**   rho MCMC routines   **/		/*    mcmc.cpp    */
	/***************************/

	double rhoSplitRate(const int nblo);					/* relative rate at which splits are proposed */
	double rhoMergeRate(const int nblo);					/* relative rate at which merges are proposed */
	omegaMap& propose_change_rBlock();						/* Proposal distribution for changing a block's rho */
	omegaMap& propose_extend_rBlock();						/* Proposal distribution for extending a rho block 5' or 3' */
	omegaMap& propose_split_rBlock();						/* Proposal distribution for splitting an existing rho block */
	omegaMap& propose_merge_rBlock();						/* Proposal distribution for merging two adjacent rho blocks */

	/***************************/
	/**  other MCMC routines  **/		/*    mcmc.cpp    */
	/***************************/

	omegaMap& propose_change_mu();							/* Proposal distribution for changing mu */
	omegaMap& propose_change_kappa();						/* Proposal distribution for changing kappa */
	omegaMap& propose_change_indelLambda();					/* Proposal distribution for changing indelLambda */
	omegaMap& propose_change_order();						/* Proposal distribution for changing the orderings */

	/***************************/
	/**    prior functions    **/		/*     inline     */
	/***************************/

	inline double logRatioPriors_change(const double theta, const double thetaPrime, const enum Priors prior, const double priorParam[2]) {
		switch(prior) {
		case uniform:
			return (thetaPrime>=priorParam[0] && thetaPrime<=priorParam[1]) ? 0.0 : -numeric_limits<double>::max();
		case improper_uniform:
			return 0.0;
		case exponential:
			return -priorParam[0] * (thetaPrime - theta);
		case gamma:
			return (priorParam[0]-1.0) * (log(thetaPrime) - log(theta)) - (thetaPrime-theta)/priorParam[1];
		case exponential_ratio:
			return 2.0 * (log(1.0+theta) - log(1.0+thetaPrime));
		case inverse:
			return (thetaPrime>=priorParam[0] && thetaPrime<=priorParam[1]) ? log(theta) - log(thetaPrime) : -numeric_limits<double>::max();
		case improper_inverse:
			return log(theta) - log(thetaPrime);
		case log_normal: {
			const double logTheta = log(theta), logThetaPrime = log(thetaPrime);
			return logTheta - logThetaPrime + 0.5*(pow((logTheta-priorParam[0])/priorParam[1],2.0)-pow((logThetaPrime-priorParam[0])/priorParam[1],2.0));
						 }
		default:
			error("omegaMap::logRatioPriors_change(): distribution not found");
		}
	}
	inline double logRatioPriors_split(const double theta, const double thetaPrime1, const double thetaPrime2, const enum Priors prior, const double priorParam[2], const double blockPriorRatio) {
		switch(prior) {
		case uniform:
			return (thetaPrime1>=priorParam[0] && thetaPrime1<=priorParam[1]
				&&	thetaPrime2>=priorParam[0] && thetaPrime2<=priorParam[1]) 
			  //				? blockPriorRatio - log(priorParam[1]-priorParam[0]) : -numeric_limits<double>::max();
			  ? log(blockPriorRatio/(1.-blockPriorRatio)/(priorParam[1]-priorParam[0])) : -numeric_limits<double>::max();
		case improper_uniform:
			error("omegaMap::logRatioPriors_split(): cannot use improper uniform distribution in RJMCMC");
		case exponential:
		  //			return blockPriorRatio + log(priorParam[0]) - priorParam[0] * (thetaPrime1 + thetaPrime2 - theta);
			return log(blockPriorRatio/(1.-blockPriorRatio) * priorParam[0]) - priorParam[0] * (thetaPrime1 + thetaPrime2 - theta);
		case gamma:
		  //return blockPriorRatio + (priorParam[0]-1.0) * (log(thetaPrime1) + log(thetaPrime2) - log(theta))
		  // - (thetaPrime1+thetaPrime2-theta)/priorParam[1]
		  // - priorParam[0] * log(priorParam[1]) - lnGAMMA(priorParam[0]);
		  return log(blockPriorRatio/(1.-blockPriorRatio)) + (priorParam[0]-1.0) * (log(thetaPrime1) + log(thetaPrime2) - log(theta))
		    - (thetaPrime1+thetaPrime2-theta)/priorParam[1]
		    - priorParam[0] * log(priorParam[1]) - lnGAMMA(priorParam[0]);
		case exponential_ratio:
			return log(blockPriorRatio/(1.-blockPriorRatio)) + 2.0 * (log(1.0+theta) - log(1.0+thetaPrime1) - log(1.0+thetaPrime2));
		case inverse:
			return (thetaPrime1>=priorParam[0] && thetaPrime1<=priorParam[1]
				&&	thetaPrime2>=priorParam[0] && thetaPrime2<=priorParam[1]) 
				? log(blockPriorRatio/(1.-blockPriorRatio)/(log(priorParam[1])-log(priorParam[0])))
				+ log(theta) - log(thetaPrime1) - log(thetaPrime2) : -numeric_limits<double>::max();
		case improper_inverse:
			error("omegaMap::logRatioPriors_split(): cannot use improper inverse distribution in RJMCMC");
		case log_normal:
			return log(blockPriorRatio/(1.-blockPriorRatio)/priorParam[1]/SQRT2PI)
				+ 0.5*(pow((log(theta)-priorParam[0])/priorParam[1],2.0)
				-pow((log(thetaPrime1)-priorParam[0])/priorParam[1],2.0)
				-pow((log(thetaPrime2)-priorParam[0])/priorParam[1],2.0));
		default:
			error("omegaMap::logRatioPriors_split(): distribution not found");
		}
	}
	inline double logRatioPriors_merge(const double theta1, const double theta2, const double thetaPrime, const enum Priors prior, const double priorParam[2], const double blockPriorRatio) {
		switch(prior) {
		case uniform:
		  //return -blockPriorRatio + log(priorParam[1]-priorParam[0]);
		  return log((1.-blockPriorRatio)/blockPriorRatio * (priorParam[1]-priorParam[0]));
		case improper_uniform:
			error("omegaMap::logRatioPriors_merge(): cannot use improper uniform distribution in RJMCMC");
		case exponential:
		  //			return -blockPriorRatio - log(priorParam[0]) + priorParam[0] * (theta1 + theta2 - thetaPrime);
			return log((1.-blockPriorRatio)/blockPriorRatio / priorParam[0]) - priorParam[0] * (thetaPrime - theta1 - theta2);
		case gamma:
			return log((1.-blockPriorRatio)/blockPriorRatio) + (priorParam[0]-1.0) * (log(thetaPrime) - log(theta1) - log(theta2))
				+ (theta1 + theta2 - thetaPrime) / priorParam[1]
				+ priorParam[0] * log(priorParam[1]) + lnGAMMA(priorParam[0]);
		case exponential_ratio:
			return log((1.-blockPriorRatio)/blockPriorRatio) + 2.0 * (log(1.0+theta1) + log(1.0+theta2) - log(1.0+thetaPrime));
		case inverse:
			return log((1.-blockPriorRatio)/blockPriorRatio*(log(priorParam[1])-log(priorParam[0])))
				+ log(theta1) + log(theta2) - log(thetaPrime);
		case improper_inverse:
			error("omegaMap::logRatioPriors_merge(): cannot use improper inverse distribution in RJMCMC");
		case log_normal:
			return log((1.-blockPriorRatio)/blockPriorRatio*priorParam[1]*SQRT2PI)
				+ 0.5*(pow((log(theta1)-priorParam[0])/priorParam[1],2.0)
				+ pow((log(theta2)-priorParam[0])/priorParam[1],2.0)
				- pow((log(thetaPrime)-priorParam[0])/priorParam[1],2.0));
		default:
			error("omegaMap::logRatioPriors_merge(): distribution not found");
		}
	}

	Priors stringToEnumPriors(string &s);
	double lnGAMMA(const double xx);
};


#endif // _OMEGA_MAP_H_
