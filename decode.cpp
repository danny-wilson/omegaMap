/************************************/
/*	decode.cpp 4th January 2006		*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*	www.danielwilson.me.uk			*/
/************************************/

#include <myutils.h>
#include "decode.h"
#include <fstream>
#include <algorithm>
#include <vector>
#include <functional>
#include <sstream>
#include <controlwizard.h>

using namespace myutils;

double mean(vector<double> &x);
void HPD(vector<double> &x, const double alpha, double &lo, double &hi);

class sort_by_vector_double : public std::binary_function<int,int,bool>
{
	const vector<double> &sort_by;
public:
	sort_by_vector_double(const vector<double> &sort_by_in) : sort_by(sort_by_in) {}

	bool operator()(int a, int b) const
	{
		return (sort_by.at(a)<sort_by.at(b));
	}
};

int main(int argc, char* argv[]) {
	if(argc!=4) error("SYNTAX: data-file out-file thinning");
	int thinning = atoi(argv[3]);
	Vector<int> fields(0);

	omegaMapAnalyse oM;
	cout << "Reading " << argv[1] << "..." << endl;
	oM.open(argv[1]);
	oM.initialize(argv[2],fields,thinning);
	cout << "Reconstructing Markov chain..." << endl;
	oM.extract();

	return 0;
}

/* produces the 100(1-alpha)% HPD interval */
/* assumes the posterior is unimodal (Chen and Shao 1999) */
void HPD(vector<double> &x, const double alpha, double &lo, double &hi) {
	lo = hi = 0.0;
	int n = x.size();
	double m = MAX(1.0,ceil(alpha*(double)n));
	//vector<double> y = x;
	sort(x.begin(),x.end());
	vector<double> y((int)m,0.0);
	int i;
	for(i=0;i<(int)m;i++) y[i] = x[n-(int)m+i]-x[i];
	vector<double> z((int)m,0.0);
	for(i=0;i<(int)m;i++) z[i] = (double)i;
	/* sort z by y */
	sort(z.begin(),z.end(),sort_by_vector_double(y));
	int j = (int)z[0];
	lo = x[j];
	hi = x[n-m+j];
}

double mean(vector<double> &x) {
	double res = 0.0;
	int i;
	for(i=0;i<x.size();i++) res+=x[i];
	return res/(double)x.size();
}


omegaMapAnalyse& omegaMapAnalyse::go(const char* datafile, const char* inifile) {
	cout << "Reading " << datafile << "..." << endl;
	open(datafile);
	initialize(inifile);
	cout << "Reconstructing Markov chain..." << endl;
	extract();
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::open(const char* datafile) {
	ifstream in(datafile);
	in >> (*this);
	in.close();
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::initialize(const char* inifile) {
	if(analyseInitialized) recycle();
	int i;

	ControlWizard *con = new ControlWizard;
	con->coutput = false;
	vector<int> temp_field;
	con->add_ITEM("outfile",TP_STRING,&outfile);
	con->add_item("fields",TP_VEC_INT,&temp_field);
	con->add_ITEM("thinning",TP_INT,&thinning);
	con->read_input(inifile);
	if(!con->check_required()) error("omegaMapAnalyse::initialize(): fatal error");

	field = Vector<int>(temp_field.size());
	for(i=0;i<field.size();i++) field[i] = temp_field[i];
	if(field.size()==0) {
		warning("omegaMap::initialize(): no fields specified for outfile. By default all fields will be output");
		field = Vector<int>(10+4*L);
		for(i=0;i<field.size();i++) field[i] = i-1;
	}
	if(thinning<=0) error("omegaMapAnalyse::initialize(): thinning interval must be positive");

	analyseInitialized = true;
	delete con;
	return recycle();
}

omegaMapAnalyse& omegaMapAnalyse::initialize(const char* outfile_in, Vector<int> &fields_in, const int thinning_in) {
	if(analyseInitialized) recycle();
	int i;

	outfile = outfile_in;
	thinning = thinning_in;

	if(field.size()==0) {
		warning("omegaMap::initialize(): no fields specified for outfile. By default all fields will be output");
		field = Vector<int>(10+4*L);
		for(i=0;i<field.size();i++) field[i] = i-1;
	}
	if(thinning<=0) error("omegaMapAnalyse::initialize(): thinning interval must be positive");

	analyseInitialized = true;
	return recycle();
}

omegaMapAnalyse& omegaMapAnalyse::recycle() {
	if(!analyseInitialized) error("omegaMapAnalyse::recycle(): not initialized");
	int i,j;

	oBlock default_oBlock;
	default_oBlock.start = default_oBlock.end = -1;
	default_oBlock._5prime = default_oBlock._3prime = 0;
	default_oBlock.oMat = NULL;
	_block = Vector<oBlock>(L,default_oBlock);
	block = Vector<oBlock*>(L,0);

	rBlock default_rBlock;
	default_rBlock.start = default_rBlock.end = -1;
	default_rBlock._5prime = default_rBlock._3prime = 0;
	default_rBlock.rho = 0.0;
	_rblock = Vector<rBlock>(L-1,default_rBlock);
	rblock = Vector<rBlock*>(L-1,0);

	oMatrix default_oMatrix;
	Vector<double> default_oMatrix_gamma_Vector = Vector<double> (0);
	default_oMatrix.gamma = 0;
	default_oMatrix.gamma2 = 0;
	default_oMatrix.R = Matrix<double>(0,0);
	default_oMatrix.lambda = Vector<double>(0);
	default_oMatrix.mu = default_oMatrix.kappa = default_oMatrix.omega = 0.0;
	pamlWork = Matrix<double>(0,0);
	for(i=0;i<L;i++) {
		if(_block[i].oMat==NULL) _block[i].oMat = new oMatrix(default_oMatrix);
	}
	oMatTemp = Vector<oMatrix*>(0);

	indelLambda = eventStart.indelLambda;

	nblocks = 0;
	oBlock* prevBlock = 0;
	int pos = 0;
	int posmax;
	if(eventStart.nblocks!=eventStart.oBlockStart.size()) error("omegaMapAnalyse::go(): nblocks and oBlockStart.size() do not match");
	for(i=0;i<eventStart.oBlockStart.size();i++) {
		++nblocks;
		pos = eventStart.oBlockStart[i];
		posmax = (i==eventStart.oBlockStart.size()-1) ? L : eventStart.oBlockStart[i+1];
		_block[pos].start = pos;
		_block[pos].end = posmax-1;
		_block[pos]._5prime = prevBlock;
		if(prevBlock!=0) {
			_block[pos]._5prime->_3prime = &(_block[pos]);
			_block[pos]._5prime->end = pos-1;
		}
		for(j=pos;j<posmax;j++) block[j] = &(_block[pos]);

		_block[pos].oMat->mu = eventStart.mu;
		_block[pos].oMat->kappa = eventStart.kappa;
		_block[pos].oMat->omega = eventStart.omega[i];

		prevBlock = &(_block[pos]);
	}

	nrblocks = 0;
	rBlock* prevRBlock = 0;
	pos = 0;
	for(i=0;i<eventStart.rBlockStart.size();i++) {
		++nrblocks;
		pos = eventStart.rBlockStart[i];
		posmax = (i==eventStart.rBlockStart.size()-1) ? L-1 : eventStart.rBlockStart[i+1];
		_rblock[pos].start = pos;
		_rblock[pos].end = posmax-1;
		_rblock[pos]._5prime = prevRBlock;
		if(prevRBlock!=0) {
			_rblock[pos]._5prime->_3prime = &(_rblock[pos]);
			_rblock[pos]._5prime->end = pos-1;
		}
		for(j=pos;j<posmax;j++) rblock[j] = &(_rblock[pos]);

		_rblock[pos].rho = eventStart.rho[i];

		prevRBlock = &(_rblock[pos]);
	}
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::extract() {
//	return extract(NULL,NULL);
//}
//
//omegaMapAnalyse& omegaMapAnalyse::extract(void (*f)(void*),void *farg) {
	ofstream out(outfile.c_str());
	headings(out);
	iter = -1;
	record(out,(oEvent&)eventStart);
	cout << "Done 0 of " << niter << " iterations";
	for(iter=0;iter<niter;iter++) {
		if(event[iter].accepted)
			switch(event[iter].type) {
			case oEvent::DEFAULT: error("omegaMapAnalyse::go(): DEFAULT event type cannot be interpreted");
			case oEvent::START: error("omegaMapAnalyse::go(): START event reserved for first iteration");
			case oEvent::CHANGE_OBLOCK:	change_oBlock();		break;
			case oEvent::EXTEND_OBLOCK:	extend_oBlock();		break;
			case oEvent::SPLIT_OBLOCK:	split_oBlock();			break;
			case oEvent::MERGE_OBLOCK:	merge_oBlock();			break;
			case oEvent::CHANGE_MU:		change_mu();			break;
			case oEvent::CHANGE_KAPPA:	change_kappa();			break;
			case oEvent::CHANGE_RBLOCK:	change_rBlock();		break;
			case oEvent::EXTEND_RBLOCK:	extend_rBlock();		break;
			case oEvent::SPLIT_RBLOCK:	split_rBlock();			break;
			case oEvent::MERGE_RBLOCK:	merge_rBlock();			break;
			case oEvent::CHANGE_INDEL:	change_indelLambda();	break;
			case oEvent::CHANGE_ORDER:	break;
			default: error("omegaMapAnalyse::go(): Event cannot be interpreted");
			}
		if((iter+1)%thinning==0) {
			record(out,event[iter]);
			cout << "\rDone " << iter+1 << " of " << niter << " iterations";
		}
//		if(f!=NULL) f(farg);
	}
	cout << "\rDone " << niter << " of " << niter << " iterations" << endl;
	out.close();
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::iterate(const int iter_in) {
	if(iter_in<0) error("omegaMapAnalyse::iterate(): iter < 0");
	if(iter_in>=niter) error("omegaMapAnalyse::iterate(): iter >= niter");
	iter = iter_in;
	if(event[iter].accepted)
		switch(event[iter].type) {
		case oEvent::DEFAULT: error("omegaMapAnalyse::go(): DEFAULT event type cannot be interpreted");
		case oEvent::START: error("omegaMapAnalyse::go(): START event reserved for first iteration");
		case oEvent::CHANGE_OBLOCK:	change_oBlock();		break;
		case oEvent::EXTEND_OBLOCK:	extend_oBlock();		break;
		case oEvent::SPLIT_OBLOCK:	split_oBlock();			break;
		case oEvent::MERGE_OBLOCK:	merge_oBlock();			break;
		case oEvent::CHANGE_MU:		change_mu();			break;
		case oEvent::CHANGE_KAPPA:	change_kappa();			break;
		case oEvent::CHANGE_RBLOCK:	change_rBlock();		break;
		case oEvent::EXTEND_RBLOCK:	extend_rBlock();		break;
		case oEvent::SPLIT_RBLOCK:	split_rBlock();			break;
		case oEvent::MERGE_RBLOCK:	merge_rBlock();			break;
		case oEvent::CHANGE_INDEL:	change_indelLambda();	break;
		case oEvent::CHANGE_ORDER:	break;
		default: error("omegaMapAnalyse::go(): Event cannot be interpreted");
		}
	return *this;
}

/*~*~*~#=#~*~*~*/
/* MCMC  moves */
/*~*~*~#=#~*~*~*/

omegaMapAnalyse& omegaMapAnalyse::change_oBlock() {
	oBlock* oB = block[event[iter].param[0]];
	oB->oMat->omega = event[iter].param[1];
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::extend_oBlock() {
	oBlock* oB;
	int newpos,oldpos,j;
	if(event[iter].param[0]>event[iter].param[1]) {// move start left
		oB = block[event[iter].param[0]];
		newpos = event[iter].param[1];
		_block[newpos].start = newpos;
		_block[newpos].end = oB->end;
		if(oB->_5prime!=0) oB->_5prime->_3prime = &(_block[newpos]);
		_block[newpos]._5prime = oB->_5prime;
		_block[newpos]._5prime->end = newpos - 1;
		if(oB->_3prime!=0) oB->_3prime->_5prime = &(_block[newpos]);
		_block[newpos]._3prime = oB->_3prime;
		SWAP(_block[newpos].oMat,oB->oMat);
		for(j=newpos;j<=oB->end;j++) block[j] = &(_block[newpos]);
	}
	else if(event[iter].param[0]<event[iter].param[1]) {// move start right
		oB = block[event[iter].param[0]-1];
		newpos = event[iter].param[1]-1;
		oldpos = oB->end+1;
		oB->end = newpos;
		oB->_3prime = &(_block[newpos+1]);
		_block[newpos+1].start = newpos + 1;
		_block[newpos+1].end = _block[oldpos].end;
		_block[newpos+1]._5prime = oB;
		_block[newpos+1]._3prime = _block[oldpos]._3prime;
		if(_block[newpos+1]._3prime!=0) _block[newpos+1]._3prime->_5prime = &(_block[newpos+1]);
		SWAP(_block[newpos+1].oMat,_block[oldpos].oMat);
		for(j=oldpos;j<=newpos;j++) block[j] = oB;
		for(j=newpos+1;j<=_block[newpos+1].end;j++) block[j] = &(_block[newpos+1]);
	}
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::split_oBlock() {
	int pos = event[iter].param[0];
	double U = event[iter].param[1];
	oBlock* oB = block[pos];
	if(pos==oB->start) {
		cout << "Problem" << endl;
	}
	_block[pos].start = pos;
	_block[pos].end = oB->end;
	oB->end = pos-1;
	_block[pos]._5prime = oB;
	_block[pos]._3prime = oB->_3prime;
	oB->_3prime = &(_block[pos]);
	if(_block[pos]._3prime!=0)
		_block[pos]._3prime->_5prime = &(_block[pos]);
	int j;
	for(j=pos;j<=_block[pos].end;j++) block[j] = &(_block[pos]);
	double a = (double)(_block[pos].start - oB->start)/(double)(_block[pos].end - oB->start + 1);
	double oldOmega = oB->oMat->omega;
	double rat = U/(1.-U);
	/*oB->oMat->omega = oldOmega * pow(U/(1.-U),1.-a);
	_block[pos].oMat->omega = oldOmega * pow((1.-U)/U,a);*/
	oB->oMat->omega = oldOmega * pow(rat,1.-a);
	_block[pos].oMat->omega = oldOmega * pow(rat,-a);
	++nblocks;
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::merge_oBlock() {
	int pos = event[iter].param[0];
	oBlock* oB1 = block[pos-1];
	oBlock* oB2 = block[pos];
	oB1->end = oB2->end;
	oB1->_3prime = oB2->_3prime;
	if(oB1->_3prime!=0)
		oB1->_3prime->_5prime = oB1;
	int j;
	for(j=oB2->start;j<=oB2->end;j++) block[j] = oB1;
	oB1->oMat->omega = event[iter].param[1];
	--nblocks;
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::change_rBlock() {
	rBlock* rB = rblock[event[iter].param[0]];
	rB->rho = event[iter].param[1];
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::extend_rBlock() {
	rBlock* rB;
	int newpos,oldpos,j;
	if(event[iter].param[0]>event[iter].param[1]) {// move start left
		rB = rblock[event[iter].param[0]];
		newpos = event[iter].param[1];
		_rblock[newpos].start = newpos;
		_rblock[newpos].end = rB->end;
		if(rB->_5prime!=0) rB->_5prime->_3prime = &(_rblock[newpos]);
		_rblock[newpos]._5prime = rB->_5prime;
		_rblock[newpos]._5prime->end = newpos - 1;
		if(rB->_3prime!=0) rB->_3prime->_5prime = &(_rblock[newpos]);
		_rblock[newpos]._3prime = rB->_3prime;
		SWAP(_rblock[newpos].rho,rB->rho);
		for(j=newpos;j<=rB->end;j++) rblock[j] = &(_rblock[newpos]);
	}
	else if(event[iter].param[0]<event[iter].param[1]) {// move start right
		rB = rblock[event[iter].param[0]-1];
		newpos = event[iter].param[1]-1;
		oldpos = rB->end+1;
		rB->end = newpos;
		rB->_3prime = &(_rblock[newpos+1]);
		_rblock[newpos+1].start = newpos + 1;
		_rblock[newpos+1].end = _rblock[oldpos].end;
		_rblock[newpos+1]._5prime = rB;
		_rblock[newpos+1]._3prime = _rblock[oldpos]._3prime;
		if(_rblock[newpos+1]._3prime!=0) _rblock[newpos+1]._3prime->_5prime = &(_rblock[newpos+1]);
		SWAP(_rblock[newpos+1].rho,_rblock[oldpos].rho);
		for(j=oldpos;j<=newpos;j++) rblock[j] = rB;
		for(j=newpos+1;j<=_rblock[newpos+1].end;j++) rblock[j] = &(_rblock[newpos+1]);
	}
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::split_rBlock() {
	int pos = event[iter].param[0];
	double U = event[iter].param[1];
	rBlock* rB = rblock[pos];
	if(pos==rB->start) {
		cout << "Problem" << endl;
	}
	_rblock[pos].start = pos;
	_rblock[pos].end = rB->end;
	rB->end = pos-1;
	_rblock[pos]._5prime = rB;
	_rblock[pos]._3prime = rB->_3prime;
	rB->_3prime = &(_rblock[pos]);
	if(_rblock[pos]._3prime!=0)
		_rblock[pos]._3prime->_5prime = &(_rblock[pos]);
	int j;
	for(j=pos;j<=_rblock[pos].end;j++) rblock[j] = &(_rblock[pos]);
	double a = (double)(_rblock[pos].start - rB->start)/(double)(_rblock[pos].end - rB->start + 1);
	double oldRho = rB->rho;
	double rat = U/(1.-U);
	/*rB->rho = oldRho * pow(U/(1.-U),1.-a);
	_rblock[pos].rho = oldRho * pow((1.-U)/U,a);*/
	rB->rho = oldRho * pow(rat,1.-a);
	_rblock[pos].rho = oldRho * pow(rat,-a);
	++nrblocks;
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::merge_rBlock() {
	int pos = event[iter].param[0];
	rBlock* rB1 = rblock[pos-1];
	rBlock* rB2 = rblock[pos];
	rB1->end = rB2->end;
	rB1->_3prime = rB2->_3prime;
	if(rB1->_3prime!=0)
		rB1->_3prime->_5prime = rB1;
	int j;
	for(j=rB2->start;j<=rB2->end;j++) rblock[j] = rB1;
	rB1->rho = event[iter].param[1];
	--nrblocks;
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::change_mu() {
	_block[0].oMat->mu = event[iter].param[0];
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::change_kappa() {
	_block[0].oMat->kappa = event[iter].param[0];
	return *this;
}

omegaMapAnalyse& omegaMapAnalyse::change_indelLambda() {
	indelLambda = event[iter].param[0];
	return *this;
}
