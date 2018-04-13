/************************************/
/*	omegaMap.cpp 23rd February 2005	*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*	www.danielwilson.me.uk			*/
/************************************/

#include "omegaMap.h"
#include <controlwizard.h>
#include <argumentwizard.h>
#include <time.h>
#include <sstream>

/***************************/
/*  Housekeeping routines  */
/***************************/

const double omegaMapBase::SQRT2PI = 2.5066282746310;

omegaMapBase::omegaMapBase() {
	ran = NULL;
}

omegaMap::omegaMap() {
	ran = NULL;
	omegaVersion = "omegaMap0.5";
	initialized = 0;
	baseToInt['T'] = 1;
	baseToInt['U'] = baseToInt['T'];
	baseToInt['C'] = 2;
	baseToInt['A'] = 3;
	baseToInt['G'] = 4;
	baseToInt['-'] = 5;
	intToBase[0] = 'N';
	intToBase[1] = 'T';
	intToBase[2] = 'C';
	intToBase[3] = 'A';
	intToBase[4] = 'G';
	intToBase[5] = '-';
	efficientLikelihood = true;
	stickyDebugEval = 0;
}

omegaMap& omegaMap::initialize(const char *inifile, Random &r) {
	return initialize(0,0,inifile,r);
}

omegaMap& omegaMap::initialize(const int argc, const char* argv[], const char *inifile, Random &r) {
	if(initialized) error("omegaMap::initialize(): already initialized");

	/***************************/
	/*      Read ini file      */
	/***************************/

	ran = &r;
	ControlWizard con;
	ArgumentWizard arg;

	// choose a bizarre number to indicate that an option has not been set
	// set required fields to 'notSet'
	const double notSet = -28091980.;
	norders = niter = thinning = (int)notSet;
	int seed = (int)notSet;
	muPrior = kaPrior = omPrior = rhPrior = inPrior = none;
	string muPriorString, kaPriorString, omPriorString, rhPriorString, inPriorString;
	muPriorString = kaPriorString = omPriorString = rhPriorString = inPriorString = "";
	vector<double> muPriorParamIn, kaPriorParamIn, omPriorParamIn, rhPriorParamIn, inPriorParamIn;
	double muStart,kaStart,inStart;
	muStart = kaStart = inStart = notSet;
	vector<double> omStart,rhStart;
	oBlockPrior = rBlockPrior = notSet;
	double local_piIndel = notSet;
	string temp_coutput = "true";
	omegaModel = rhoModel = "";

	// set defaults
	vector<double> local_pi;
	vector<int> orders;
	maxCodon = 61;
	writeBinary = false;
	string temp_writeBinary = "";
	string temp_efficient = "";
	string overwrite = "false";
	vector<int> temp_field;

	// set MCMC defaults
	prMOVE = Vector<double>(14,0.0);
	int i;
	for(i=2;i<=12;i++) prMOVE[i] = 0.1;
	prMOVE[4] = prMOVE[10] = 0.2;
	prMOVE[5] = prMOVE[11] = 0.0;
	prMOVE[6] = prMOVE[7] = prMOVE[12] = 0.2/3.0;
	prMOVE[13] = 0.0;

	// read in initialization file
	con.coutput = false;											arg.coutput = false;
	con.add_item("FASTA",TP_STRING,&hapfile);						arg.add_item("FASTA",TP_STRING,&hapfile);
	con.add_item("pi",TP_VEC_DOUBLE,&local_pi);						arg.add_item("pi",TP_VEC_DOUBLE,&local_pi);
	con.add_item("piindel",TP_DOUBLE,&local_piIndel);				arg.add_item("piindel",TP_DOUBLE,&local_piIndel);
	con.add_item("norders",TP_INT,&norders);						arg.add_item("norders",TP_INT,&norders);
	con.add_item("orders",TP_VEC_INT,&orders);						arg.add_item("orders",TP_VEC_INT,&orders);
	con.add_item("seed",TP_INT,&seed);								arg.add_item("seed",TP_INT,&seed);
	
	con.add_item("muPrior",TP_STRING,		&muPriorString);		arg.add_item("muPrior",TP_STRING,		&muPriorString);
	con.add_item("kappaPrior",TP_STRING,	&kaPriorString);		arg.add_item("kappaPrior",TP_STRING,	&kaPriorString);
	con.add_item("omegaPrior",TP_STRING,	&omPriorString);		arg.add_item("omegaPrior",TP_STRING,	&omPriorString);
	con.add_item("rhoPrior",TP_STRING,		&rhPriorString);		arg.add_item("rhoPrior",TP_STRING,		&rhPriorString);
	con.add_item("indelPrior",TP_STRING,	&inPriorString);		arg.add_item("indelPrior",TP_STRING,	&inPriorString);
	con.add_item("muParam",TP_VEC_DOUBLE,	&muPriorParamIn);		arg.add_item("muParam",TP_VEC_DOUBLE,	&muPriorParamIn);
	con.add_item("kappaParam",TP_VEC_DOUBLE,&kaPriorParamIn);		arg.add_item("kappaParam",TP_VEC_DOUBLE,&kaPriorParamIn);
	con.add_item("omegaParam",TP_VEC_DOUBLE,&omPriorParamIn);		arg.add_item("omegaParam",TP_VEC_DOUBLE,&omPriorParamIn);
	con.add_item("rhoParam",TP_VEC_DOUBLE,	&rhPriorParamIn);		arg.add_item("rhoParam",TP_VEC_DOUBLE,	&rhPriorParamIn);
	con.add_item("indelParam",TP_VEC_DOUBLE,&inPriorParamIn);		arg.add_item("indelParam",TP_VEC_DOUBLE,&inPriorParamIn);
	con.add_item("muStart",TP_DOUBLE,		&muStart);				arg.add_item("muStart",TP_DOUBLE,		&muStart);
	con.add_item("kappaStart",TP_DOUBLE,	&kaStart);				arg.add_item("kappaStart",TP_DOUBLE,	&kaStart);
	con.add_item("omegaStart",TP_VEC_DOUBLE,&omStart);				arg.add_item("omegaStart",TP_VEC_DOUBLE,&omStart);
	con.add_item("rhoStart",TP_VEC_DOUBLE,	&rhStart);				arg.add_item("rhoStart",TP_VEC_DOUBLE,	&rhStart);
	con.add_item("indelStart",TP_DOUBLE,	&inStart);				arg.add_item("indelStart",TP_DOUBLE,	&inStart);

/*	con.add_item("mu",TP_DOUBLE,&muPrior);							arg.add_item("mu",TP_DOUBLE,&muPrior);
	con.add_item("kappa",TP_DOUBLE,&kappaPrior);					arg.add_item("kappa",TP_DOUBLE,&kappaPrior);
	con.add_item("omega",TP_DOUBLE,&omegaPrior);					arg.add_item("omega",TP_DOUBLE,&omegaPrior);
	con.add_item("rho",TP_DOUBLE,&rhoPrior);						arg.add_item("rho",TP_DOUBLE,&rhoPrior);
	con.add_item("indel",TP_DOUBLE,&indelLambdaPrior);				arg.add_item("indel",TP_DOUBLE,&indelLambdaPrior);*/

	con.add_item("oblock",TP_DOUBLE,&oBlockPrior);					arg.add_item("oblock",TP_DOUBLE,&oBlockPrior);
	con.add_item("rblock",TP_DOUBLE,&rBlockPrior);					arg.add_item("rblock",TP_DOUBLE,&rBlockPrior);
	con.add_item("niter",TP_INT,&niter);							arg.add_item("niter",TP_INT,&niter);
	con.add_item("OMEGA_MODEL",TP_STRING,&omegaModel);				arg.add_item("OMEGA_MODEL",TP_STRING,&omegaModel);
	con.add_item("RHO_MODEL",TP_STRING,&rhoModel);					arg.add_item("RHO_MODEL",TP_STRING,&rhoModel);
	con.add_item("binary",TP_STRING,&temp_writeBinary);				arg.add_item("binary",TP_STRING,&temp_writeBinary);
	con.add_item("CHANGE_OBLOCK",TP_DOUBLE,&(prMOVE[2]));			arg.add_item("CHANGE_OBLOCK",TP_DOUBLE,&(prMOVE[2]));
	con.add_item("EXTEND_OBLOCK",TP_DOUBLE,&(prMOVE[3]));			arg.add_item("EXTEND_OBLOCK",TP_DOUBLE,&(prMOVE[3]));
	con.add_item("SPLITMERGE_OBLOCK",TP_DOUBLE,&(prMOVE[4]));		arg.add_item("SPLITMERGE_OBLOCK",TP_DOUBLE,&(prMOVE[4]));
	con.add_item("CHANGE_MU",TP_DOUBLE,&(prMOVE[6]));				arg.add_item("CHANGE_MU",TP_DOUBLE,&(prMOVE[6]));
	con.add_item("CHANGE_KAPPA",TP_DOUBLE,&(prMOVE[7]));			arg.add_item("CHANGE_KAPPA",TP_DOUBLE,&(prMOVE[7]));
	con.add_item("CHANGE_RBLOCK",TP_DOUBLE,&(prMOVE[8]));			arg.add_item("CHANGE_RBLOCK",TP_DOUBLE,&(prMOVE[8]));
	con.add_item("EXTEND_RBLOCK",TP_DOUBLE,&(prMOVE[9]));			arg.add_item("EXTEND_RBLOCK",TP_DOUBLE,&(prMOVE[9]));
	con.add_item("SPLITMERGE_RBLOCK",TP_DOUBLE,&(prMOVE[10]));		arg.add_item("SPLITMERGE_RBLOCK",TP_DOUBLE,&(prMOVE[10]));
	con.add_item("CHANGE_INDEL",TP_DOUBLE,&(prMOVE[12]));			arg.add_item("CHANGE_INDEL",TP_DOUBLE,&(prMOVE[12]));
	con.add_item("CHANGE_ORDER",TP_DOUBLE,&(prMOVE[13]));			arg.add_item("CHANGE_ORDER",TP_DOUBLE,&(prMOVE[13]));
	con.add_item("datafile",TP_STRING,&datafile);					arg.add_item("datafile",TP_STRING,&datafile);
	con.add_item("outfile",TP_STRING,&outfile);						arg.add_item("outfile",TP_STRING,&outfile);
	con.add_item("fields",TP_VEC_INT,&temp_field);					arg.add_item("fields",TP_VEC_INT,&temp_field);
	con.add_item("thinning",TP_INT,&thinning);						arg.add_item("thinning",TP_INT,&thinning);
	con.add_item("efficient",TP_STRING,&temp_efficient);			arg.add_item("efficient",TP_STRING,&temp_efficient);
	con.add_item("evaluate",TP_INT,&stickyDebugEval);				arg.add_item("evaluate",TP_INT,&stickyDebugEval);
	con.add_item("coutput",TP_STRING,&temp_coutput);				arg.add_item("coutput",TP_STRING,&temp_coutput);
	con.add_item("overwrite",TP_STRING,&overwrite);					arg.add_item("overwrite",TP_STRING,&overwrite);
	con.read_input(inifile);										arg.read_input(argc,argv);

	// generate a seed for the random number generator if specified
	// (otherwise accept the seed automatically generated from the system clock)
	if(seed!=(int)notSet) {
		cout << "Setting seed for random number generator to " << seed << endl;
		ran->setseed(seed);
	}

	// check whether constant or variable models have been chosen for omega and rho
	for(i=0;i<(int)omegaModel.length();i++) omegaModel[i] = tolower(omegaModel[i]);
	for(i=0;i<(int)rhoModel.length();i++) rhoModel[i] = tolower(rhoModel[i]);
/*	#ifdef __GNUC__
		transform(omegaModel.begin(),omegaModel.end(),omegaModel.begin(),(int(*)(int))std::tolower);
		transform(rhoModel.begin(),rhoModel.end(),rhoModel.begin(),(int(*)(int))std::tolower);
	#else
		transform(omegaModel.begin(),omegaModel.end(),omegaModel.begin(),tolower);
		transform(rhoModel.begin(),rhoModel.end(),rhoModel.begin(),tolower);
	#endif*/
	if(omegaModel=="") error("omegaMap::initialize(): OMEGA_MODEL is a required option");
	if(omegaModel=="constant") {
		prMOVE[3] = prMOVE[4] = prMOVE[5] = 0.0;
	}
	else if(omegaModel=="variable") {
		if(oBlockPrior==notSet) error("omegaMap::initialize(): for variable omega model, must give option oBlock");
	}
	else if(omegaModel=="independent") {
		if(oBlockPrior!=notSet && oBlockPrior!=1.0) warning("omegaMap::initialize(): for independent omega model oBlock should not be set. Ignoring oBlock");
		oBlockPrior = 1.0;
	}
	else error("omegaMap::initialize(): OMEGA_MODEL must equal CONSTANT, VARIABLE or INDEPENDENT");
	if(rhoModel=="") error("omegaMap::initialize(): RHO_MODEL is a required option");
	if(rhoModel=="constant") {
		prMOVE[9] = prMOVE[10] = prMOVE[11] = 0.0;
	}
	else if(rhoModel=="variable") {
		if(rBlockPrior==notSet) error("omegaMap::initialize(): for variable rho model, must give option rBlock");
	}
	else if(rhoModel=="independent") {
		if(rBlockPrior!=notSet && rBlockPrior!=1.0) warning("omegaMap::initialize(); for independent rho model rBlock should not be set. Ignoring rBlock");
		rBlockPrior = 1.0;
	}
	else error("omegaMap::initialize(): RHO_MODEL must equal CONSTANT, VARIABLE or INDEPENDENT");
	oBlockPriorRatio = log(oBlockPrior) - log(1.0-oBlockPrior);
	rBlockPriorRatio = log(rBlockPrior) - log(1.0-rBlockPrior);

	if(local_piIndel==0.0) {
		prMOVE[12] = 0.0;
		maxCodon = 61;
	}

	// check that outfile and datafile can be written to OK
	for(i=0;i<(int)overwrite.size();i++) overwrite[i] = tolower(overwrite[i]);
	if(overwrite!="true" && overwrite!="false") warning("omegaMap::initialize(): overwrite must be TRUE or FALSE. Defaulting to false");
	if(overwrite!="true" && outfile!="") {
		ifstream exists(outfile.c_str());
		if(exists.is_open()) {
			stringstream eText;
			eText << "omegaMap::initialize(): " << outfile << " already exists";
			exists.close();
			error(eText.str().c_str());
		}
		else exists.close();
	}
	if(overwrite!="true" && datafile!="") {
		ifstream exists(datafile.c_str());
		if(exists.is_open()) {
			stringstream eText;
			eText << "omegaMap::initialize(): " << datafile << " already exists";
			exists.close();
			error(eText.str().c_str());
		}
		else exists.close();
	}
	if(outfile!="") {
		ofstream checkWrite(outfile.c_str());
		if(!checkWrite.is_open()) {
			stringstream eText;
			eText << "omegaMap::initialize(): " << outfile << " cannot open for writing";
			checkWrite.close();
			error(eText.str().c_str());
		}
		else checkWrite.close();
	}
	if(datafile!="") {
		ofstream checkWrite(datafile.c_str());
		if(!checkWrite.is_open()) {
			stringstream eText;
			eText << "omegaMap::initialize(): " << datafile << " cannot open for writing";
			checkWrite.close();
			error(eText.str().c_str());
		}
		else checkWrite.close();
	}
	if(outfile=="" && datafile=="") error("omegaMap::initialize(): no outfile and no datafile");
	// read in the haplotypes from the FASTA file
	if(hapfile=="") error("omegaMap::initialize(): must specify a name for FASTA file");
	dna.readFASTA(hapfile.c_str());
	// check for STOP codons, and convert into codon numbers
	dnaToCodons();

	// check the priors have been specified properly
	enum Priors* vecPrior[5] = {&muPrior,&kaPrior,&omPrior,&rhPrior,&inPrior};
	string* vecPriorString[5] = {&muPriorString,&kaPriorString,&omPriorString,&rhPriorString,&inPriorString};
	vector<double>* vecPriorParamIn[5] = {&muPriorParamIn,&kaPriorParamIn,&omPriorParamIn,&rhPriorParamIn,&inPriorParamIn};
	double* vecPriorParam[5] = {muPriorParam,kaPriorParam,omPriorParam,rhPriorParam,inPriorParam};
	const char* vecName[5] = {"mu","kappa","omega","rho","indel"};
	const int reqParams[10] = {0, 1, 2, 0, 1, 2, 0, 2, 0, 2};

	int j;
	for(i=0;i<5;i++) {
		stringstream errorText;
		errorText << "omegaMap::initialize(): ";
		*vecPrior[i] = stringToEnumPriors(*vecPriorString[i]);
		if(vecPriorParamIn[i]->size()!=reqParams[*vecPrior[i]]) {
			if(*vecPrior[i]==fixed && i==2 && vecPriorParamIn[i]->size()==L) {}
			else if(*vecPrior[i]==fixed && i==3 && vecPriorParamIn[i]->size()==(L-1)) {}
			else {
				errorText << "wrong # parameters for prior on " << vecName[i];
				error(errorText.str().c_str());
			}
		}
		if(*vecPrior[i]==fixed)
			for(j=0;j<(int)vecPriorParamIn[i]->size();j++)
				if((*vecPriorParamIn[i])[0]<0.0) {
					errorText << vecName[i] << "[" << j+1 << "] is fixed at a negative number";
					error(errorText.str().c_str());
				}
		if(*vecPrior[i]==uniform && (*vecPriorParamIn[i])[0]<0.0) {
			errorText << vecName[i] << " uniform prior left bound is negative";
			error(errorText.str().c_str());
		}
		if(*vecPrior[i]==uniform && (*vecPriorParamIn[i])[1]<0.0) {
			errorText << vecName[i] << " uniform prior right bound is negative";
			error(errorText.str().c_str());
		}
		if(*vecPrior[i]==uniform && (*vecPriorParamIn[i])[0]>=(*vecPriorParamIn[i])[1]) {
			errorText << vecName[i] << " uniform right bound is less than left";
			error(errorText.str().c_str());
		}
		if(*vecPrior[i]==exponential && (*vecPriorParamIn[i])[0]<=0.0) {
			errorText << vecName[i] << " exponential parameter must be positive";
			error(errorText.str().c_str());
		}
		if(*vecPrior[i]==gamma && (*vecPriorParamIn[i])[0]<=0.0) {
			errorText << vecName[i] << " gamma 1st (shape) parameter must be positive";
			error(errorText.str().c_str());
		}
		if(*vecPrior[i]==gamma && (*vecPriorParamIn[i])[1]<=0.0) {
			errorText << vecName[i] << " gamma 2nd (scale) parameter must be positive";
			error(errorText.str().c_str());
		}
		for(j=0;j<reqParams[*vecPrior[i]];j++) vecPriorParam[i][j] = (*vecPriorParamIn[i])[j];
	}
	if(omegaModel=="variable" && (omPrior==improper_uniform || omPrior==improper_inverse)) error("omegaMap::initialize(): improper prior on omega under variable model");
	if(rhoModel=="variable" && (rhPrior==improper_uniform || rhPrior==improper_inverse)) error("omegaMap::initialize(): improper prior on rho under variable model");
	/* Using improper priors, the initial values of the parameter must be specified */
	if((muPrior==improper_uniform || muPrior==improper_inverse) && muStart==notSet) error("omegaMap::initialize(): improper prior requires that muStart be set");
	if((kaPrior==improper_uniform || kaPrior==improper_inverse) && kaStart==notSet) error("omegaMap::initialize(): improper prior requires that kappaStart be set");
	if((inPrior==improper_uniform || inPrior==improper_inverse) && inStart==notSet) error("omegaMap::initialize(): improper prior requires that indelStart be set");
	if((omPrior==improper_uniform || omPrior==improper_inverse) && omStart.size()==0) error("omegaMap::initialize(): improper prior requires that omegaStart be set");
	if((rhPrior==improper_uniform || rhPrior==improper_inverse) && rhStart.size()==0) error("omegaMap::initialize(): improper prior requires that rhoStart be set");

/*	if(muPrior<=0.0) error("omegaMap::initialize() mu must be greater than zero");
	else muPrior = 1.0/muPrior;
	if(kappaPrior<=0.0) error("omegaMap::initialize() kappa must be greater than zero");
	else kappaPrior = 1.0/kappaPrior;
	if(omegaPrior<=0.0) error("omegaMap::initialize() omega must be greater than zero");
	else omegaPrior = 1.0/omegaPrior;
	if(rhoPrior<=0.0) error("omegaMap::initialize() rho must be greater than zero");
	else rhoPrior = 1.0/rhoPrior;
	if(indelLambdaPrior!=notSet && indelLambdaPrior<=0.0) error("omegaMap::initialize() indel must be greater than zero");*/
	if(omegaModel=="variable") {
		if(oBlockPrior<1.0) error("omegaMap::initialize(): oBlock must be greater than (or equal to) one");
		else if(oBlockPrior==1.0) {
			warning("omegaMap::initialize(): oBlock equals one forces each site to have its own omega block.\nEquivalent to INDEPENDENT model");
			prMOVE[3] = 0.0;
			prMOVE[4] = 0.0;		
		}
		else oBlockPrior = 1.0/oBlockPrior;
	}
	else {
		prMOVE[3] = prMOVE[4] = 0.0;
	}
	if(rhoModel=="variable") {
		if(rBlockPrior<1.0) error("omegaMap::initialize(): rBlock must be greater than (or equal to) one");
		else if(rBlockPrior==1.0) {
			warning("omegaMap::initialize(): rBlock equals one forces each site to have its own rho block.\nEquivalent to INDEPENDENT model");
			prMOVE[9] = 0.0;
			prMOVE[10] = 0.0;		
		}
		else rBlockPrior = 1.0/rBlockPrior;
	}
	else {
		prMOVE[9] = prMOVE[10] = 0.0;
	}
	if(muPrior==fixed) for(j=0;j<(int)muPriorParamIn.size();j++) if(muPriorParamIn[j]<=0.0) error("omegaMap::initialize(): mu cannot be fixed <= 0");
	if(kaPrior==fixed) for(j=0;j<(int)kaPriorParamIn.size();j++) if(kaPriorParamIn[j]<=0.0) error("omegaMap::initialize(): kappa cannot be fixed <= 0");
	if(inPrior==fixed) for(j=0;j<(int)inPriorParamIn.size();j++) if(inPriorParamIn[j]<=0.0) error("omegaMap::initialize(): indel cannot be fixed <= 0");
	if(muPrior==fixed && muStart!=notSet) warning("omegaMap::initialize(): muPrior = fixed will override muStart");
	if(kaPrior==fixed && kaStart!=notSet) warning("omegaMap::initialize(): kappaPrior = fixed will override kappaStart");
	if(omPrior==fixed && omStart.size()!=0) warning("omegaMap::initialize(): omegaPrior = fixed will override omegaStart");
	if(rhPrior==fixed && rhStart.size()!=0) warning("omegaMap::initialize(): rhoPrior = fixed will override rhoStart");
	if(inPrior==fixed && inStart!=notSet) warning("omegaMap::initialize(): indelPrior = fixed will override indelStart");
	if(muPrior==fixed) prMOVE[6] = 0.0;
	if(kaPrior==fixed) prMOVE[7] = 0.0;
	if(omPrior==fixed) prMOVE[2] = prMOVE[3] = prMOVE[4] = 0.0;
	if(rhPrior==fixed) prMOVE[8] = prMOVE[9] = prMOVE[10] = 0.0;
	if(inPrior==fixed) prMOVE[12] = 0.0;
	if(muStart!=notSet && muStart<=0.0) error("omegaMap::initialize(): muStart must be positive");
	if(kaStart!=notSet && kaStart<=0.0) error("omegaMap::initialize(): kappaStart must be positive");
	if(omStart.size()!=0) {
		if(omegaModel=="constant" && omStart.size()!=1) error("omegaMap::initialize(): for constant omegaModel, omegaStart must have length 1");
		else if(omegaModel=="variable") error("omegaMap::initialize(): cannot specify omegaStart for variable omegaModel");
		else if(omegaModel=="independent" && omStart.size()==1) omStart = vector<double>(omStart[0],L);
		else if(omegaModel=="independent" && omStart.size()!=L) error("omegaMap::initialize(): for independent omegaModel, omegaStart must have length L");
	}
	for(i=0;i<(int)omStart.size();i++) if(omStart[i]<0.0) error("omegaMap::initialize(): omegaStart must be non-negative");
	if(rhStart.size()!=0) {
		if(rhoModel=="constant" && rhStart.size()!=1) error("omegaMap::initialize(): for constant rhoModel, rhoStart must have length 1");
		else if(rhoModel=="variable") error("omegaMap::initialize(): cannot specify rhoStart for variable rhoModel");
		else if(rhoModel=="independent" && rhStart.size()==1) rhStart = vector<double>(rhStart[0],L-1);
		else if(rhoModel=="independent" && rhStart.size()!=(L-1)) error("omegaMap::initialize(): for independent rhoModel, rhoStart must have length L-1");
	}
	for(i=0;i<(int)rhStart.size();i++) if(rhStart[i]<0.0) error("omegaMap::initialize(): rhoStart must be non-negative");
	if(inStart!=notSet && inStart<=0.0) error("omegaMap::initialize(): indelStart must be positive");

	// Adjust priors in the case that L==1 or L==2
	if(L==1) {
		oBlockPrior = rBlockPrior = 1.0;
		prMOVE[3] = prMOVE[4] = prMOVE[8] = prMOVE[9] = prMOVE[10] = 0.0;
		warning("omegaMap::initialize(): L=1. There will be one omega block and no recombination.\nAdjusting MCMC proposal probabilities");
	}
	else if(L==2) {
		prMOVE[9] = prMOVE[10] = 0.0;
		warning("omegaMap::initialize(): L=2. There will be one rho block.\nAdjusting MCMC proposal probabilities");
	}

	// check the MCMC move probabilities have been properly specified
	double prTOTAL = 0.0;
	for(i=2;i<prMOVE.size();i++) {
		if(prMOVE[i]<0.0) error("omegaMap::initialize(): negative MCMC move probability specified");
		prTOTAL += prMOVE[i];
	}
	if(prTOTAL==0.0) error("omegaMap::initialize(): MCMC move probabilities sum to zero");
	if(prTOTAL!=1.0)
		for(i=2;i<prMOVE.size();i++) prMOVE[i]/=prTOTAL;
	if(prMOVE[2]==0.0 && omPrior!=fixed) error("omegaMap::initialize(): CHANGE_OBLOCK proposal probability equals zero");
	if(omegaModel=="variable" && omPrior!=fixed) {
		if(prMOVE[3]==0.0) error("omegaMap::initialize(): EXTEND_OBLOCK proposal probability equals zero");
		if(prMOVE[4]==0.0) error("omegaMap::initialize(): SPLITMERGE_OBLOCK proposal probability equals zero");
	}
	if(prMOVE[6]==0.0 && muPrior!=fixed) error("omegaMap::initialize(): CHANGE_MU proposal probability equals zero");
	if(prMOVE[7]==0.0 && kaPrior!=fixed) error("omegaMap::initialize(): CHANGE_KAPPA proposal probability equals zero");
	if(prMOVE[8]==0.0 && rhPrior!=fixed) error("omegaMap::initialize(): CHANGE_RBLOCK proposal probability equals zero");
	if(rhoModel=="variable" && rhPrior!=fixed) {
		if(prMOVE[9]==0.0) error("omegaMap::initialize(): EXTEND_RBLOCK proposal probability equals zero");
		if(prMOVE[10]==0.0) error("omegaMap::initialize(): SPLITMERGE_RBLOCK proposal probability equals zero");
	}

	// Deal with indels
	bool anyIndels = false;
	indel = Vector<bool> (L,false);
	Vector<double> obsPi(62,0.0);
	double obsIndel[2] = {0.0,0.0};
	int dan1,dan2;
	for(dan1=0;dan1<codon.size();dan1++) {
		for(dan2=0;dan2<codon[dan1].size();dan2++) {
			if(codon[dan1][dan2]==61) {
				anyIndels = true;
				indel[dan2] = true;
				obsIndel[0] += 1.0;
			}
			else {
				obsPi[codon[dan1][dan2]] += 1.0;
				obsPi[61] += 1.0;
			}
//			cout << codon[dan1][dan2] << " ";
		}
//		cout << endl;
	}
//	cout << endl;
	for(i=0;i<62;i++) obsPi[i] /= obsPi[61];
	for(i=0;i<L;i++) if(indel[i]) obsIndel[1] += 1.0;
	obsIndel[0] /= (obsIndel[1]*n);
	maxCodon = (anyIndels) ? 62 : 61;
	if(anyIndels && inPrior==none) error("omegaMap::initialize(): indels detected, so indel is a required option");
//	else indelLambdaPrior = 1.0/indelLambdaPrior;

	// check the equilibrium codon and indel frequencies have been correctly specified
	if(local_pi.size()==0) {
		// estimate from data
		local_pi = vector<double>(61,0.0);
		cout << "Estimating codon usage from data:" << endl;
		for(i=0;i<61;i++) {
			local_pi[i] = obsPi[i];
			cout << local_pi[i] << " ";
		}	cout << endl;
	}
	else if(local_pi.size()!=61) error("omegaMap::initialize(): pi must be of length 61");
	if(anyIndels) {
		if(local_piIndel==notSet) {
			// estimate from data
			local_piIndel = obsIndel[0];
			cout << "Estimating indel frequency from the data: " << local_piIndel << endl;
			if(local_piIndel<=0.0) error("omegaMap::initialize(): piindel was estimated to be negative. This is an internal error");
			if(local_piIndel>=1.0) error("omegaMap::initialize(): piindel was estimated to be one or greater. This is an internal error");
		}
		else {
			if(local_piIndel<0.0) error("omegaMap::initialize(): piindel must be non-negative");
			if(local_piIndel==0.0) warning("omegaMap::initialize(): piindel has been set to equal zero");
			if(local_piIndel>=1.0) error("omegaMap::initialize(): piindel must be less than one");
		}
		if(prMOVE[12]==0.0 && inPrior!=fixed) error("omegaMap::initialize(): CHANGE_INDEL proposal probability equals zero, but indels found");
	}
	pi = Vector<double>(maxCodon,0.0);
	double tot = 0.0;
	for(i=0;i<61;i++) {
		if(local_pi[i]==0.0) {
			stringstream wText;
			wText << "omegaMap::initialize(): pi[" << i << "](" << codonToTriplet61(i) << ") = 0.0";
			warning(wText.str().c_str());
		}
		tot += local_pi[i];
	}
	if(fabs(tot-1.0)>1.e-3) {
		stringstream wText;
		wText << "omegaMap::initialize(): pi does not sum exactly to one (sums to " << tot << ")";
		warning(wText.str().c_str());
		for(i=0;i<61;i++) local_pi[i]/=tot;
	}
	for(i=0;i<61;i++) pi[i] = local_pi[i];
	for(i=61;i<maxCodon;i++) pi[i] = local_piIndel;
	sqrtPi = pi;
	for(i=0;i<maxCodon;i++) sqrtPi[i] = sqrt(sqrtPi[i]);

	// check for efficient likelihood calculations
	for(i=0;i<(int)temp_efficient.length();i++) temp_efficient[i] = tolower(temp_efficient[i]);
	if(temp_efficient=="") efficientLikelihood = true;
	else if(temp_efficient=="false") {
		efficientLikelihood = false;
		if(stickyDebugEval<0 || stickyDebugEval>=L) error("omegaMap::initialize(): evaluate must lie in the interval [0,L-1]");
	}
	else if(temp_efficient=="true") efficientLikelihood = true;
	else {
		warning("omegaMap::initialize(): efficient must be true or false. Setting to true");
		efficientLikelihood = true;
	}

	// check the PAC orderings have been properly specified
	if(norders==(int)notSet) error("omegaMap::initialize(): norders is a required option");
	if(norders<=0) error("omegaMap::initialize(): norders must be positive");
	if(orders.size()==0) {
		warning("Orderings not set. The following are used...");
		orders = setOrders();
		for(i=0;i<(int)orders.size();i++) cout << orders[i] << ", ";
		cout << "\b\b \n";
	}
	else {
		if(orders.size()!=norders*n) error("omegaMap::initialize(): number of elements of orders must equal norders*n");
		setOrders(orders);
	}

	for(i=0;i<(int)temp_writeBinary.length();i++) temp_writeBinary[i] = tolower(temp_writeBinary[i]);
	for(i=0;i<(int)temp_coutput.length();i++) temp_coutput[i] = tolower(temp_coutput[i]);
/*	#ifdef __GNUC__
		transform(temp_writeBinary.begin(),temp_writeBinary.end(),temp_writeBinary.begin(),(int(*)(int))std::tolower);
	#else
		transform(temp_writeBinary.begin(),temp_writeBinary.end(),temp_writeBinary.begin(),tolower);
	#endif*/
	if(temp_writeBinary=="true") writeBinary = true;
	else if(temp_writeBinary=="false") writeBinary = false;
	else if(temp_writeBinary=="") writeBinary = false;
	else {
		warning("omegaMap::initialize(): binary should equal true or false. Set to false");
		writeBinary = false;
	}
	if(writeBinary) error("omegaMap::initialize(): cannot write in binary");
	if(temp_coutput=="true") coutput = true;
	else if(temp_coutput=="false") coutput = false;
	else if(temp_coutput=="") coutput = false;
	else coutput = true;

	// check the output file and fields have been correctly specified
	field = Vector<int>((int)temp_field.size());
	for(i=0;i<field.size();i++) field[i] = temp_field[i];

	if(outfile=="" &&
		(thinning!=(int) notSet || field.size()!=0))
		warning("omegaMap::initialize(): outfile not specified. There will be no text output");
	if(outfile!="") {
		if(field.size()==0) {
			//warning("omegaMap::initialize(): no fields specified for outfile. By default all fields will be output");
			temp_field = vector<int>(12);
			for(i=-1;i<11;i++) temp_field[i+1] = i;

			if(omegaModel=="constant") temp_field.push_back(11+L);
			else for(i=11;i<11+2*L;i++) temp_field.push_back(i);

			if(rhoModel=="constant") temp_field.push_back(10+3*L);
			else for(i=11+2*L;i<9+4*L;i++) temp_field.push_back(i);

			field = Vector<int>((int)temp_field.size());
			for(i=0;i<field.size();i++) field[i] = temp_field[i];
		}
		if(thinning==(int) notSet) {
			warning("omegaMap::initialize: thinning interval not specified. There will be no text output");
			outfile = "";
		}
		if(thinning!=(int) notSet && thinning<=0) {
			warning("omegaMap::initialize(): thinning interval must be positive. There will be no text output");
			outfile = "";
		}
	}

	/***************************/
	/*    Memory allocation    */
	/***************************/

	clock_t start = clock();
	size_t memSize = 2*norders*(1+(n-1)*(1+L*(1+(n+2)/2)));
	memSize += 2*L*(1830*n+3785);
	memSize *= sizeof(double);
	cout << "Required memory size is approximately " << (double)(memSize)/1.e6 << " Megabytes" << endl;

	oBlock default_oBlock;
	default_oBlock.start = default_oBlock.end = -1;
	default_oBlock._5prime = default_oBlock._3prime = 0;
	default_oBlock.oMat = 0;
	_block = Vector<oBlock>(L,default_oBlock);
	block = Vector<oBlock*>(L,0);

	rBlock default_rBlock;
	default_rBlock.start = default_rBlock.end = -1;
	default_rBlock._5prime = default_rBlock._3prime = 0;
	default_rBlock.rho = 0.0;
	_rblock = Vector<rBlock>(L-1,default_rBlock);
	rblock = Vector<rBlock*>(L-1,0);

	oMatrix default_oMatrix;
	Vector<double> default_oMatrix_gamma_Vector = Vector<double> (n,-1.0);
	default_oMatrix.gamma = 0;
	default_oMatrix.gamma2 = 0;
	default_oMatrix.R = Matrix<double>(61,61);
	default_oMatrix.lambda = Vector<double>(61);
	default_oMatrix.mu = default_oMatrix.kappa = default_oMatrix.omega = 0.0;
	pamlWork = Matrix<double>(61,61);
	for(i=0;i<L;i++) {
		_block[i].oMat = new oMatrix(default_oMatrix);
		_block[i].oMat->gamma = new LowerTriangularMatrix< Vector<double> >(61,default_oMatrix_gamma_Vector);
		if(maxCodon==62)
			_block[i].oMat->gamma2 = new LowerTriangularMatrix< Vector<double> >(62,default_oMatrix_gamma_Vector);
	}
	oMatTemp = Vector<oMatrix*>(L,0);
	for(i=0;i<L;i++) {
		oMatTemp[i] = new oMatrix(default_oMatrix);
		oMatTemp[i]->gamma = new LowerTriangularMatrix< Vector<double> >(61,default_oMatrix_gamma_Vector);
		if(maxCodon==62) 
			oMatTemp[i]->gamma2 = new LowerTriangularMatrix< Vector<double> >(62,default_oMatrix_gamma_Vector);
	}

	/*************************************/
	/*    Initialize the Markov Chain    */
	/*************************************/

	eventStart.nblocks = 0;
	eventStart.nrblocks = 0;
	if(muPrior==fixed)					eventStart.mu = muPriorParam[0];
	else if(muStart!=notSet)			eventStart.mu = muStart;
	else if(muPrior==uniform)			eventStart.mu = ran->uniform(muPriorParam[0],muPriorParam[1]);
	else if(muPrior==improper_uniform)	eventStart.mu = ran->uniform(0.0,1.0e6);
	else if(muPrior==exponential)		eventStart.mu = ran->exponential(1.0/muPriorParam[0]);
	else if(muPrior==gamma)				eventStart.mu = ran->gamma(muPriorParam[1],muPriorParam[0]);
	else if(muPrior==exponential_ratio) eventStart.mu = ran->exponential_ratio();
	else if(muPrior==inverse)			eventStart.mu = ran->inverse(muPriorParam[0],muPriorParam[1]);
	else if(muPrior==improper_inverse)	eventStart.mu = ran->inverse(1.0e-6,1.0e6);
	else if(muPrior==log_normal)		eventStart.mu = ran->log_normal(muPriorParam[0],muPriorParam[1]);
	else error("omegaMap::initialize(): cannot draw from the prior on mu");

	if(kaPrior==fixed)					eventStart.kappa = kaPriorParam[0];
	else if(kaStart!=notSet)			eventStart.kappa = kaStart;
	else if(kaPrior==uniform)			eventStart.kappa = ran->uniform(kaPriorParam[0],kaPriorParam[1]);
	else if(kaPrior==improper_uniform)	eventStart.kappa = ran->uniform(0.0,1.0e6);
	else if(kaPrior==exponential)		eventStart.kappa = ran->exponential(1.0/kaPriorParam[0]);
	else if(kaPrior==gamma)				eventStart.kappa = ran->gamma(kaPriorParam[1],kaPriorParam[0]);
	else if(kaPrior==exponential_ratio) eventStart.kappa = ran->exponential_ratio();
	else if(kaPrior==inverse)			eventStart.kappa = ran->inverse(kaPriorParam[0],kaPriorParam[1]);
	else if(kaPrior==improper_inverse)	eventStart.kappa = ran->inverse(1.0e-6,1.0e6);
	else if(kaPrior==log_normal)		eventStart.kappa = ran->log_normal(kaPriorParam[0],kaPriorParam[1]);
	else error("omegaMap::initialize(): cannot draw from the prior on kappa");

	if(maxCodon==61) eventStart.indelLambda = 0.0;
	else if(inPrior==fixed)				eventStart.indelLambda = inPriorParam[0];
	else if(inStart!=notSet)			eventStart.indelLambda = inStart;
	else if(inPrior==uniform)			eventStart.indelLambda = ran->uniform(inPriorParam[0],inPriorParam[1]);
	else if(inPrior==improper_uniform)	eventStart.indelLambda = ran->uniform(0.0,1.0e6);
	else if(inPrior==exponential)		eventStart.indelLambda = ran->exponential(1.0/inPriorParam[0]);
	else if(inPrior==gamma)				eventStart.indelLambda = ran->gamma(inPriorParam[1],inPriorParam[0]);
	else if(inPrior==exponential_ratio) eventStart.indelLambda = ran->exponential_ratio();
	else if(inPrior==inverse)			eventStart.indelLambda = ran->inverse(inPriorParam[0],inPriorParam[1]);
	else if(inPrior==improper_inverse)	eventStart.indelLambda = ran->inverse(1.0e-6,1.0e6);
	else if(inPrior==log_normal)		eventStart.indelLambda = ran->log_normal(inPriorParam[0],inPriorParam[1]);
	else error("omegaMap::initialize(): cannot draw from the prior on indel");
	indelLambda = eventStart.indelLambda;

	nblocks = 0;
	oBlock* prevBlock = 0;
	int pos = 0;
	while(pos<=L-1) {
		++nblocks;
		++eventStart.nblocks;
		_block[pos].start = pos;
		_block[pos].end = L-1;
		_block[pos]._5prime = prevBlock;
		if(prevBlock!=0) {
			_block[pos]._5prime->_3prime = &(_block[pos]);
			_block[pos]._5prime->end = pos-1;
		}
		for(i=pos;i<L;i++) block[i] = &(_block[pos]);

		if(omPrior==fixed){
			if(omPriorParamIn.size()==1)	_block[pos].oMat->omega = omPriorParam[0];
			else							_block[pos].oMat->omega = omPriorParamIn[pos];
		}
		else if(omStart.size()!=0)			_block[pos].oMat->omega = omStart[pos];
		else if(omPrior==uniform)			_block[pos].oMat->omega = ran->uniform(omPriorParam[0],omPriorParam[1]);
		else if(omPrior==improper_uniform)	_block[pos].oMat->omega = ran->uniform(0.0,1.0e6);
		else if(omPrior==exponential)		_block[pos].oMat->omega = ran->exponential(1.0/omPriorParam[0]);
		else if(omPrior==gamma)				_block[pos].oMat->omega = ran->gamma(omPriorParam[1],omPriorParam[0]);
		else if(omPrior==exponential_ratio) _block[pos].oMat->omega = ran->exponential_ratio();
		else if(omPrior==inverse)			_block[pos].oMat->omega = ran->inverse(omPriorParam[0],omPriorParam[1]);
		else if(omPrior==improper_inverse)	_block[pos].oMat->omega = ran->inverse(1.0e-6,1.0e6);
		else if(omPrior==log_normal)		_block[pos].oMat->omega = ran->log_normal(omPriorParam[0],omPriorParam[1]);
		else error("omegaMap::initialize(): cannot draw from the prior on omega");

		_block[pos].oMat->mu = eventStart.mu;
		_block[pos].oMat->kappa = eventStart.kappa;
//		_block[pos].oMat->omega = ran->exponential(1.0/omegaPrior);
		diagonalizeMatrix(*(_block[pos].oMat),_block[pos].oMat->mu,_block[pos].oMat->kappa,_block[pos].oMat->omega);

		prevBlock = &(_block[pos]);
		pos += (omegaModel=="constant") ? L : 1 + ran->geometric(oBlockPrior);
		if(oBlockPrior==1.0) ++pos;
	}

	nrblocks = 0;
	rBlock* prevRBlock = 0;
	pos = 0;
	while(pos<=L-2) {
		++nrblocks;
		++eventStart.nrblocks;
		_rblock[pos].start = pos;
		_rblock[pos].end = L-2;
		_rblock[pos]._5prime = prevRBlock;
		if(prevRBlock!=0) {
			_rblock[pos]._5prime->_3prime = &(_rblock[pos]);
			_rblock[pos]._5prime->end = pos-1;
		}
		for(i=pos;i<L-1;i++) rblock[i] = &(_rblock[pos]);

		if(rhPrior==fixed){
			if(rhPriorParamIn.size()==1)	_rblock[pos].rho = rhPriorParam[0];
			else							_rblock[pos].rho = rhPriorParamIn[pos];
		}
		else if(rhStart.size()!=0)			_rblock[pos].rho = rhStart[pos];
		else if(rhPrior==uniform)			_rblock[pos].rho = ran->uniform(rhPriorParam[0],rhPriorParam[1]);
		else if(rhPrior==improper_uniform)	_rblock[pos].rho = ran->uniform(0.0,1.0e6);
		else if(rhPrior==exponential)		_rblock[pos].rho = ran->exponential(1.0/rhPriorParam[0]);
		else if(rhPrior==gamma)				_rblock[pos].rho = ran->gamma(rhPriorParam[1],rhPriorParam[0]);
		else if(rhPrior==exponential_ratio) _rblock[pos].rho = ran->exponential_ratio();
		else if(rhPrior==inverse)			_rblock[pos].rho = ran->inverse(rhPriorParam[0],rhPriorParam[1]);
		else if(rhPrior==improper_inverse)	_rblock[pos].rho = ran->inverse(1.0e-6,1.0e6);
		else if(rhPrior==log_normal)		_rblock[pos].rho = ran->log_normal(rhPriorParam[0],rhPriorParam[1]);
		else error("omegaMap::initialize(): cannot draw from the prior on rho");

//		_rblock[pos].rho = ran->exponential(1.0/rhoPrior);

		prevRBlock = &(_rblock[pos]);
		pos += (rhoModel=="constant") ? L : 1 + ran->geometric(rBlockPrior);
		if(rBlockPrior==1.0) ++pos;
	}

	/***************************/
	/*  Record initial state   */
	/***************************/

	oEvent default_oEvent;
	default_oEvent.type = oEvent::DEFAULT;
	default_oEvent.accepted = false;
	default_oEvent.param[0] = default_oEvent.param[1] = default_oEvent.likelihood = default_oEvent.alpha = 0.0;
	if(datafile!="") event = Vector<oEvent>(niter,default_oEvent);
	else event = Vector<oEvent>(1,default_oEvent);
	iter = 0;

	eventStart.type = oEvent::START;
	eventStart.oBlockStart = Vector<int>(nblocks,0);
	eventStart.rBlockStart = Vector<int>(nrblocks,0);
	eventStart.omega = Vector<double>(nblocks,0.0);
	eventStart.rho = Vector<double>(nrblocks,0.0);
	eventStart.alpha = 0.0;

	oBlock* oB;
	for(i=0, oB=block[0]; oB!=0; i++) {
		eventStart.oBlockStart[i] = oB->start;
		eventStart.omega[i] = oB->oMat->omega;
		oB = oB->_3prime;
	}
	rBlock* rB;
	if(nrblocks>0)
	for(i=0, rB=rblock[0]; rB!=0; i++) {
		eventStart.rBlockStart[i] = rB->start;
		eventStart.rho[i] = rB->rho;
		rB = rB->_3prime;
	}
	eventStart.accepted = true;


	/***************************/
	/*    Memory allocation    */
	/***************************/

	//rho = Vector<double>(L,eventStart.rho[0]);
	PAC = Vector<double>(norders,0.0);

	size_t size_mem_dbl = norders*L*(n*n+n-2) * sizeof(double);
	size_t size_mem_dblp = 2*norders*L*(n-1) * sizeof(double*);
	size_t size_mem_dblpp = 2*norders*(n-1) * sizeof(double**);
	size_t size_mem_dblppp = 2*norders * sizeof(double***);

	/*double* mem3 = new double[norders*(n-1)*L*(n+2)];
	if(!mem3) error("omegaMap::initialize(): failed to allocate mem3");

	double** mem2 = new double*[2*norders*(n-1)*L];
	if(!mem2) error("omegaMap::initialize(): failed to allocate mem2");

	double*** mem1 = new double**[2*norders*(n-1)];
	if(!mem1) error("omegaMap::initialize(): failed to allocate mem1");

	alpha = new double***[norders];
	if(!alpha) error("omegaMap::initialize(): failed to allocate alpha");

	beta = new double***[norders];
	if(!beta) error("omegaMap::initialize(): failed to allocate beta");
	cout << "Memory allocation complete";

	int ord,k,x;
	int mem1ctr,mem2ctr,mem3ctr;
	mem1ctr=mem2ctr=mem3ctr=0;

	for(ord=0;ord<norders;ord++) {
		alpha[ord] = &(mem1[mem1ctr]);
		--alpha[ord];
		for(k=1;k<n;k++,mem1ctr++) {
			alpha[ord][k] = &(mem2[mem2ctr]);
			for(pos=0;pos<L;pos++,mem2ctr++) {
				alpha[ord][k][pos] = &(mem3[mem3ctr]);
				for(x=0;x<=k;x++,mem3ctr++)
					alpha[ord][k][pos][x] = 0.0;
			}
		}
	}
	for(ord=0;ord<norders;ord++) {
		beta[ord] = &(mem1[mem1ctr]);
		--beta[ord];
		for(k=1;k<n;k++,mem1ctr++) {
			beta[ord][k] = &(mem2[mem2ctr]);
			for(pos=0;pos<L;pos++,mem2ctr++) {
				beta[ord][k][pos] = &(mem3[mem3ctr]);
				for(x=0;x<=k;x++,mem3ctr++)
					beta[ord][k][pos][x] = 0.0;
			}
		}
	}*/

	alpha = new double***[norders];
	if(!alpha) error("omegaMap::initialize(): failed to allocate alpha");
	beta = new double***[norders];
	if(!beta) error("omegaMap::initialize(): failed to allocate beta");
	int ord,k,x;
	for(ord=0;ord<norders;ord++) {
		alpha[ord] = new double**[n-1];					// only elements 1..n-1 will be access-	
		if(!alpha[ord]) error("omegaMap::initialize(): failed to allocate alpha[]");
		--alpha[ord];									// ible. This saves memory.				
		beta[ord] = new double**[n-1];
		if(!beta[ord]) error("omegaMap::initialize(): failed to allocate beta[]");
		--beta[ord];
		for(k=1;k<n;k++) {
			alpha[ord][k] = new double*[L];
			if(!alpha[ord][k]) error("omegaMap::initialize(): failed to allocate alpha[][]");
			beta[ord][k] = new double*[L];
			if(!beta[ord][k]) error("omegaMap::initialize(): failed to allocate beta[][]");
			for(pos=0;pos<L;pos++) {
				alpha[ord][k][pos] = new double[k+1];	// the kth element is the sum of the	
				if(!alpha[ord][k][pos]) error("omegaMap::initialize(): failed to allocate alpha[][][]");
				beta[ord][k][pos] = new double[k+1];	// first k-1 elements					
				if(!beta[ord][k][pos]) error("omegaMap::initialize(): failed to allocate beta[][][]");
				for(x=0;x<=k;x++) {
					alpha[ord][k][pos][x] = 1.0;
					beta[ord][k][pos][x] = 1.0;
				}
			}
		}
	}

	memTime = (double)(clock()-start)/CLOCKS_PER_SEC;
	cout << "\rMemory allocation took " << memTime << " seconds." << endl;

	/***************************/
	/* Initialize likelihoods  */
	/***************************/

	alphaMargin = 0;
	betaMargin = L-1;
	pos = (L-1)/2;
	forward(alphaMargin,pos);
	backward(pos,betaMargin);
	oldLikelihood = likelihood(pos);
	eventStart.likelihood = oldLikelihood;
	alphaMargin = pos+1;
	betaMargin = pos-1;

	return *this;
}

omegaMap::~omegaMap() {
	if(initialized) {
		int ord,k,pos;
		for(ord=0;ord<norders;ord++) {
			for(k=1;k<n;k++) {
				for(pos=0;pos<L;pos++) {
					delete [] alpha[ord][k][pos];
					delete [] beta[ord][k][pos];
				}
				delete [] alpha[ord][k];
				delete [] beta[ord][k];
			}
			++alpha[ord];
			delete [] alpha[ord];
			++beta[ord];
			delete [] beta[ord];
		}
		alpha = beta = 0;
		int i;
		for(i=0;i<L;i++) {
			delete _block[i].oMat->gamma;
			if(maxCodon==62) delete _block[i].oMat->gamma2;
			delete _block[i].oMat;
		}
		for(i=0;i<L;i++) {
			delete oMatTemp[i]->gamma;
			if(maxCodon==62) delete oMatTemp[i]->gamma2;
			delete oMatTemp[i];
		}
	}
}

/* Returns 0-60 for non-STOP codons, 61 for indels */
omegaMap& omegaMap::dnaToCodons() {
	n = dna.nseq;
	L = (int)floor((double)dna.lseq/(double)3.0);
	if(n<=0) error("omegaMap::dnaToCodons(): n must be positive");
	if(L<=0) error("omegaMap::dnaToCodons(): L must be positive");
	
	Vector<int> default_codon_sequence(L,-1);
	codon = Vector< Vector<int> >(n,default_codon_sequence);
	int i,j,t2c;
	string triplet,tr1,tr2,tr3;
	for(i=0;i<n;i++)
		for(j=0;j<L;j++) {
			tr1 = string(1,toupper(dna[i][3*j]));
			tr2 = string(1,toupper(dna[i][3*j+1]));
			tr3 = string(1,toupper(dna[i][3*j+2]));
			//triplet = dna[i][3*j]+dna[i][3*j+1]+dna[i][3*j+2];
			triplet = tr1 + tr2 + tr3;
			if(triplet!="---" && (tr1=="-" || tr2=="-" || tr3=="-")) {
				stringstream wrnTxt;
				wrnTxt << "omegaMap::dnaToCodons(): '" << triplet << "' found in sequence " << i+1 << " position " << j+1 << endl;
				wrnTxt << "will be interpreted as an indel '---'";
				warning(wrnTxt.str().c_str());
				triplet = "---";
			}
			t2c = tripletToCodon(triplet);
			if(t2c<0 || t2c>64) {
				stringstream errTxt;
				errTxt << "omegaMap::dnaToCodons(): '" << triplet << "' found in sequence " << i+1 << " position " << j+1 << endl;
				errTxt << "is not understood.";
				error(errTxt.str().c_str());
			}
			codon[i][j] = removeStopCodons(t2c);
		}
	
	return *this;
}

string omegaMap::codonToTriplet(const int a) {
	if(a==64) return string(3,'-');
	string tri(3,'N');
	tri[0] = intToBase[a/16 + 1];
	tri[1] = intToBase[(a%16)/4 + 1];
	tri[2] = intToBase[(a%16)%4 + 1];
	return tri;
}

string omegaMap::codonToTriplet61(const int b) {
	int a = b;
	if(a==61) return string(3,'-');
	if(a>=10) a+=2;
	if(a>=14) a+=1;
	string tri(3,'N');
	tri[0] = intToBase[a/16 + 1];
	tri[1] = intToBase[(a%16)/4 + 1];
	tri[2] = intToBase[(a%16)%4 + 1];
	return tri;
}

/* Returns 0-63 for codons and 64 for indels */
int omegaMap::tripletToCodon(string &tri) {
	int a = baseToInt[tri[0]];
	int b = baseToInt[tri[1]];
	int c = baseToInt[tri[2]];
	bool indel = false;
	if(a==5) indel = true;
	if(b==5) indel = true;
	if(c==5) indel = true;
	if(indel==true) {
		if(a==5 && b==5 && c==5) return 64;
		else return -1;
	}
	/* return a value from 0 to 63 */
	return (a-1)*16 + (b-1)*4 + c - 1;
}

int omegaMap::removeStopCodons(const int a) {
	if(a<0) error("omegaMap::removeStopCodons(): invalid codon number (<0)");
	if(a<10) return a;
	if(a==10) error("omegaMap::removeStopCodons(): found stop codon TAA/UAA");
	if(a==11) error("omegaMap::removeStopCodons(): found stop codon TAG/UAG");
	if(a<14) return a-2;
	if(a==14) error("omegaMap::removeStopCodons(): found stop codon TGA/UGA");
	/* a==64 is for indels. Change this to 61 */
	if(a>64) error("omegaMap::removeStopCodons(): invalid codon number (>64)");
	return a-3;
}

/*omegaMap& omegaMap::setOrders() {
	H = Matrix< Vector<int>* >(norders,n,0);
	int i,j,k,rnum,hap;
	Vector<int> pool(n);
	for(i=0;i<norders;i++) {
		for(j=0;j<n;j++) pool[j] = j;
		for(j=0,k=n-1;j<n;j++,k--) {
			rnum = ran->discrete(0,k);
			hap = pool[rnum];
			pool[rnum] = pool[k];
			H[i][j] = &(codon[hap]);
		}
	}	
	return *this;
}*/

vector<int> omegaMap::setOrders() {
	vector<int> orders(norders*n);
	H = Matrix< Vector<int>* >(norders,n,0);
	int i,j,k,rnum,hap,ctr;
	Vector<int> pool(n);
	for(i=0,ctr=0;i<norders;i++) {
		for(j=0;j<n;j++) pool[j] = j;
		for(j=0,k=n-1;j<n;j++,k--,ctr++) {
			rnum = ran->discrete(0,k);
			hap = pool[rnum];
			pool[rnum] = pool[k];
			orders[ctr] = hap;
			H[i][j] = &(codon[hap]);
		}
	}	
	return orders;
}

omegaMap& omegaMap::setOrders(vector<int> &orders) {
	H = Matrix< Vector<int>* >(norders,n,0);
	int i,j,k;
	for(i=0,k=0;i<norders;i++)
		for(j=0;j<n;j++,k++)
			H[i][j] = &(codon[orders[k]]);
	return *this;
}

/***Debug error-checking*************/
omegaMap& omegaMap::debug() {
	int j;
	if((oldLikelihood==0.0)
		|| (oldLikelihood!=oldLikelihood)) {
		cout << "Likelihood problem\n";
	}
	oBlock* oB = block[0];
	int end;
	while(oB!=0) {
		end = oB->end;
		if(oB->_5prime!=0)
			if(oB->_5prime->_3prime!=oB) {
				warning("omegaMap::debug():5prime oBlock continuity problem");
			}
		if(oB->_3prime!=0)
			if(oB->_3prime->_5prime!=oB) {
				warning("omegaMap::debug():3prime oBlock continuity problem");
			}
		if(oB->_5prime==oB || oB->_3prime==oB) {
			warning("omegaMap::debug():recursive oBlock problem");
		}
		if(oB->start>oB->end) {
			warning("omegaMap::debug():oBlock order consistency problem");
		}
		if(oB->_5prime!=0 && oB->_5prime->end!=oB->start-1) {
			warning("omegaMap::debug():oBlock contiguity problem");
		}
		oB = oB->_3prime;
	}
	if(end!=L-1) {
		warning("omegaMap::debug():oBlock terminal-sequence problem");
	}
	for(j=0;j<L;j++)
		if(block[j]->start>j || block[j]->end<j) {
			warning("omegaMap::debug():oBlock positioning problem");
		}

/*	int i,x,y;
	oB = block[0];
	while(oB!=0) {
		for(i=0;i<maxCodon;i++) for(j=0;j<=i;j++)
			for(x=0;x<=i;x++) for(y=0;y<j;y++) {
				if(&(oB->oMat->gamma[i][j])==&(oB->oMat->gamma[x][y])) {
					warning("omegaMap::debug():memory problem");
				}
			}
		oB = oB->_3prime;
	}*/

	return *this;
}

omegaMap& omegaMap::redraw() {
	double newMu;// = ran->exponential(1.0/muPrior);
	if(muPrior!=fixed) {
		if(muPrior==uniform)				newMu = ran->uniform(muPriorParam[0],muPriorParam[1]);
		else if(muPrior==improper_uniform)	newMu = ran->uniform(0.0,1.0e6);
		else if(muPrior==exponential)		newMu = ran->exponential(1.0/muPriorParam[0]);
		else if(muPrior==gamma)				newMu = ran->gamma(muPriorParam[1],muPriorParam[0]);
		else if(muPrior==exponential_ratio) newMu = ran->exponential_ratio();
		else if(muPrior==inverse)			newMu = ran->inverse(muPriorParam[0],muPriorParam[1]);
		else if(muPrior==improper_inverse)	newMu = ran->inverse(1.0e-6,1.0e6);
		else if(muPrior==log_normal)		newMu = ran->log_normal(muPriorParam[0],muPriorParam[1]);
		else error("omegaMap::redraw(): cannot draw from the prior on mu");
	}
	double newKappa;// = ran->exponential(1.0/kappaPrior);
	if(kaPrior!=fixed) {
		if(kaPrior==uniform)				newKappa = ran->uniform(kaPriorParam[0],kaPriorParam[1]);
		else if(kaPrior==improper_uniform)	newKappa = ran->uniform(0.0,1.0e6);
		else if(kaPrior==exponential)		newKappa = ran->exponential(1.0/kaPriorParam[0]);
		else if(kaPrior==gamma)				newKappa = ran->gamma(kaPriorParam[1],kaPriorParam[0]);
		else if(kaPrior==exponential_ratio) newKappa = ran->exponential_ratio();
		else if(kaPrior==inverse)			newKappa = ran->inverse(kaPriorParam[0],kaPriorParam[1]);
		else if(kaPrior==improper_inverse)	newKappa = ran->inverse(1.0e-6,1.0e6);
		else if(kaPrior==log_normal)		newKappa = ran->log_normal(kaPriorParam[0],kaPriorParam[1]);
		else error("omegaMap::redraw(): cannot draw from the prior on kappa");
	}
	double newIndel;
	if(maxCodon==62 && inPrior!=fixed) {
		if(inPrior==uniform)				newIndel = ran->uniform(inPriorParam[0],inPriorParam[1]);
		else if(inPrior==improper_uniform)	newIndel = ran->uniform(0.0,1.0e6);
		else if(inPrior==exponential)		newIndel = ran->exponential(1.0/inPriorParam[0]);
		else if(inPrior==gamma)				newIndel = ran->gamma(inPriorParam[1],inPriorParam[0]);
		else if(inPrior==exponential_ratio) newIndel = ran->exponential_ratio();
		else if(inPrior==inverse)			newIndel = ran->inverse(inPriorParam[0],inPriorParam[1]);
		else if(inPrior==improper_inverse)	newIndel = ran->inverse(1.0e-6,1.0e6);
		else if(inPrior==log_normal)		newIndel = ran->log_normal(inPriorParam[0],inPriorParam[1]);
		else error("omegaMap::redraw(): cannot draw from the prior on indel");
		indelLambda = newIndel;
	}
	else indelLambda = 0.0;

	int i,pos;
	if(omPrior!=fixed) {
		nblocks = 0;
		oBlock* prevBlock = 0;
		pos = 0;
		while(pos<=L-1) {
			++nblocks;
			_block[pos].start = pos;
			_block[pos].end = L-1;
			_block[pos]._5prime = prevBlock;
			if(prevBlock!=0) {
				_block[pos]._5prime->_3prime = &(_block[pos]);
				_block[pos]._5prime->end = pos-1;
			}
			for(i=pos;i<L;i++) block[i] = &(_block[pos]);

			if(omPrior==uniform)				_block[pos].oMat->omega = ran->uniform(omPriorParam[0],omPriorParam[1]);
			else if(omPrior==improper_uniform)	_block[pos].oMat->omega = ran->uniform(0.0,1.0e6);
			else if(omPrior==exponential)		_block[pos].oMat->omega = ran->exponential(1.0/omPriorParam[0]);
			else if(omPrior==gamma)				_block[pos].oMat->omega = ran->gamma(omPriorParam[1],omPriorParam[0]);
			else if(omPrior==exponential_ratio) _block[pos].oMat->omega = ran->exponential_ratio();
			else if(omPrior==inverse)			_block[pos].oMat->omega = ran->inverse(omPriorParam[0],omPriorParam[1]);
			else if(omPrior==improper_inverse)	_block[pos].oMat->omega = ran->inverse(1.0e-6,1.0e6);
			else if(omPrior==log_normal)		_block[pos].oMat->omega = ran->log_normal(omPriorParam[0],omPriorParam[1]);
			else error("omegaMap::redraw(): cannot draw from the prior on omega");

			_block[pos].oMat->mu = newMu;
			_block[pos].oMat->kappa = newKappa;
	//		_block[pos].oMat->omega = ran->exponential(1.0/omegaPrior);
			diagonalizeMatrix(*(_block[pos].oMat),_block[pos].oMat->mu,_block[pos].oMat->kappa,_block[pos].oMat->omega);

			prevBlock = &(_block[pos]);
			pos += (omegaModel=="constant") ? L : 1 + ran->geometric(oBlockPrior);
			if(oBlockPrior==1.0) ++pos;
		}
	}

	if(rhPrior!=fixed) {
		nrblocks = 0;
		rBlock* prevRBlock = 0;
		pos = 0;
		while(pos<=L-2) {
			++nrblocks;
			_rblock[pos].start = pos;
			_rblock[pos].end = L-2;
			_rblock[pos]._5prime = prevRBlock;
			if(prevRBlock!=0) {
				_rblock[pos]._5prime->_3prime = &(_rblock[pos]);
				_rblock[pos]._5prime->end = pos-1;
			}
			for(i=pos;i<L-1;i++) rblock[i] = &(_rblock[pos]);

			if(rhPrior==uniform)				_rblock[pos].rho = ran->uniform(rhPriorParam[0],rhPriorParam[1]);
			else if(rhPrior==improper_uniform)	_rblock[pos].rho = ran->uniform(0.0,1.0e6);
			else if(rhPrior==exponential)		_rblock[pos].rho = ran->exponential(1.0/rhPriorParam[0]);
			else if(rhPrior==gamma)				_rblock[pos].rho = ran->gamma(rhPriorParam[1],rhPriorParam[0]);
			else if(rhPrior==exponential_ratio) _rblock[pos].rho = ran->exponential_ratio();
			else if(rhPrior==inverse)			_rblock[pos].rho = ran->inverse(rhPriorParam[0],rhPriorParam[1]);
			else if(rhPrior==improper_inverse)	_rblock[pos].rho = ran->inverse(1.0e-6,1.0e6);
			else if(rhPrior==log_normal)		_rblock[pos].rho = ran->log_normal(rhPriorParam[0],rhPriorParam[1]);
			else error("omegaMap::redraw(): cannot draw from the prior on rho");

//			_rblock[pos].rho = ran->exponential(1.0/rhoPrior);

			prevRBlock = &(_rblock[pos]);
			pos += (rhoModel=="constant") ? L : 1 + ran->geometric(rBlockPrior);
			if(rBlockPrior==1.0) ++pos;
		}
	}

	return *this;
}

omegaMap::Priors omegaMap::stringToEnumPriors(string &s) {
	int i;
	for(i=0;i<(int)s.length();i++) s[i] = tolower(s[i]);
	string values[10] = {
		"",
		"fixed",
		"uniform",
		"improper_uniform",
		"exponential",
		"gamma",
		"exponential_ratio",
		"inverse",
		"improper_inverse",
		"log_normal"
	};
	if(s==values[0])		return none;
	else if(s==values[1])	return fixed;
	else if(s==values[2])	return uniform;
	else if(s==values[3])	return improper_uniform;
	else if(s==values[4])	return exponential;
	else if(s==values[5])	return gamma;
	else if(s==values[6])	return exponential_ratio;
	else if(s==values[7])	return inverse;
	else if(s==values[8])	return improper_inverse;
	else if(s==values[9])	return log_normal;
	stringstream eText;
	eText << "omegaMap::initialize(): " << s << " unrecognised prior distribution";
	error(eText.str().c_str());
	return none;
}

/* This function from Numerical Recipes in C++ */
double omegaMap::lnGAMMA(const double xx) {
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}