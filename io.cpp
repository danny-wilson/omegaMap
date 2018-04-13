/************************************/
/*	io.cpp 5th February 2006		*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*	www.danielwilson.me.uk			*/
/************************************/

#include "omegaMap.h"
#include <controlwizard.h>
#include <fstream>
#include <sstream>
#include <iomanip>

/******************************************/
/*  Input/output routines : omegaMapBase  */
/******************************************/

omegaMapBase& omegaMapBase::headings(ofstream& out) {
	return headings(out,tab);
}

omegaMapBase& omegaMapBase::record(ofstream& out, oEvent& event) {
	return record(out,event,tab);
}

omegaMapBase& omegaMapBase::headings(ofstream& out, ostream& omanip(ostream&)) {
	int i,j,k;
	bool validField;
	out << "iter";
	for(i=0;i<field.size();i++) {
		validField = true;
/*		if(field[i]<9) {
			j = field[i];
		}
		else {
			j = field[i]-9;
			k = j % L;
			j = 9 + j / L;
		}*/
		if(field[i]<11) {
			j = field[i];
		}
		else if(field[i]<11+L) {
			j = 11;
			k = field[i]-j;
		}
		else if(field[i]<11+2*L) {
			j = 12;
			k = field[i]-11-L;
		}
		else if(field[i]<10+3*L) {
			j = 13;
			k = field[i]-11-2*L;
		}
		else if(field[i]<9+4*L) {
			j = 14;
			k = field[i]-10-3*L;
		}
		else {
			validField = false;
		}
		if(validField) 	out << omanip;
		switch(j) {
		case -1:out << "thetaS"; break;
		case 0: out << "mu"; break;
		case 1: out << "kappa"; break;
		case 2: out << "phi"; break;
		case 3: out << "moveType"; break;
		case 4: out << "accepted"; break;
		case 5: out << "likelihood"; break;
		case 6: out << "alpha"; break;
		case 7: out << "nblocks"; break;
		case 8: out << "nrblocks"; break;
		case 9: out << "param0"; break;
		case 10: out << "param1"; break;
		case 11: out << "oBlockStart" << k; break;
		case 12: out << "omega" << k; break;
		case 13: out << "rBlockStart" << k; break;
		case 14: out << "rho" << k; break;
		}
	}
	out << endl;
	return *this;
}

omegaMapBase& omegaMapBase::record(ofstream& out, oEvent& event, ostream& omanip(ostream&)) {
	int i,j,k;
	bool validField;
	out << (iter+1);
	for(i=0;i<field.size();i++) {
		validField = true;
		if(field[i]<11) {
			j = field[i];
		}
		else if(field[i]<11+L) {
			j = 11;
			k = field[i]-j;
		}
		else if(field[i]<11+2*L) {
			j = 12;
			k = field[i]-11-L;
		}
		else if(field[i]<10+3*L) {
			j = 13;
			k = field[i]-11-2*L;
		}
		else if(field[i]<9+4*L) {
			j = 14;
			k = field[i]-10-3*L;
		}
		else {
			validField = false;
		}
		if(validField) 	out << omanip;
		switch(j) {
		case -1:out << thetaS(); break;
		case 0: out << _block[0].oMat->mu; break;
		case 1: out << _block[0].oMat->kappa; break;
		case 2: out << indelLambda; break;
		case 3: out << (int)event.type; break;
		case 4: out << event.accepted; break;
		case 5: out << event.likelihood; break;
		case 6: out << event.alpha; break;
		case 7: out << nblocks; break;
		case 8: out << nrblocks; break;
		case 9: out << event.param[0]; break;
		case 10: out << event.param[1]; break;
		case 11: out << (block[k]->start==k) ? 1 : 0; break;
		case 12: out << block[k]->oMat->omega; break;
		case 13: out << (rblock[k]->start==k) ? 1 : 0; break;
		case 14: out << rblock[k]->rho; break;
		}
	}
	out << endl;
	if(out.fail()) error("omegaMapBase::record(): write fail");
	return *this;
}

/******************************************/
/*    Input/output routines : omegaMap    */
/******************************************/

/*omegaMap& omegaMap::outputToText(const char* filename) {
	ofstream out(filename);
	string omegaVersion = "omegaMap0.01";
	out << omegaVersion << tab << n << tab << L << tab << niter << tab << ran->getseed() << tab;
	int i;
	for(i=0;i<61;i++) out << pi[i] << tab;
	out << norders << tab << muPrior << tab << kappaPrior << tab << omegaPrior << tab << rhoPrior << tab << oBlockPrior << tab << oSplitMergeProbability << tab;
	out << eventStart;
	for(i=0;i<niter;i++) out << event[i];
	out << endl;
	out.close();
	return *this;
}*/

/*omegaMapData& omegaMapData::inputFromText(const char* filename){
	ifstream in(filename);
	in >> omegaVersion;
	if(omegaVersion!="omegaMap0.01") {
		stringstream wText;
		wText << "omegaMapData::inputFromText(): expecting omegaMap0.01 data but received " << omegaVersion;
		warning(wText.str().c_str());
	}
	in >> n >> L >> niter >> seed;
	pi = Vector<double>(61,-1.0);
	event = Vector<oEvent>(niter);
	int i;
	for(i=0;i<61;i++) in >> pi[i];
	in >> norders >> muPrior >> kappaPrior >> omegaPrior >> rhoPrior >> oBlockPrior >> oSplitMergeProbability;
	in >> eventStart;
	for(i=0;i<niter && !in.eof();i++) in >> event[i];
	in.close();
	return *this;
}*/

/******************************************/
/*     Input/output routines : global     */
/******************************************/

inline ostream& tab(ostream& out) {
	return out << "\t";
}

inline ostream& ENDCONTROLWIZARD(ostream& out) {
	return out << '*';
}

ostream& operator<<(ostream& out, omegaMapBase& oM) {
	int i;
	out << setprecision(16);
	out << "omegaVersion = " << oM.omegaVersion << endl;
	out << "initialized = " << (int)oM.initialized << endl;
	out << "coutput = " << (int)oM.coutput << endl;
	out << "hapfile = " << oM.hapfile << endl;
	out << "efficientLikelihood = " << (int)oM.efficientLikelihood << endl;
	out << "stickyDebugEval = " << oM.stickyDebugEval << endl;
/*	out << "muPrior = " << oM.muPrior << endl;
	out << "kappaPrior = " << oM.kaPrior << endl;
	out << "omegaPrior = " << oM.omPrior << endl;
	out << "rhoPrior = " << oM.rhPrior << endl;
	out << "indelLambdaPrior = " << oM.inPrior << endl;*/
	out << "oBlockPrior = " << oM.oBlockPrior << endl;
	out << "oBlockPriorRatio = " << oM.oBlockPriorRatio << endl;
	out << "rBlockPrior = " << oM.rBlockPrior << endl;
	out << "rBlockPriorRatio = " << oM.rBlockPriorRatio << endl;
	out << "writeBinary = " << (int)oM.writeBinary << endl;
	out << "prMOVE = ";
	for(i=0;i<oM.prMOVE.size();i++) {
		if(i>0) out << ", ";
		out << oM.prMOVE[i];
	}	out << endl;
	out << "omegaModel = " << oM.omegaModel << endl;
	out << "rhoModel = " << oM.rhoModel << endl;
	out << "muPrior = " << (int)oM.muPrior << endl;
	out << "kaPrior = " << (int)oM.kaPrior << endl;
	out << "omPrior = " << (int)oM.omPrior << endl;
	out << "rhPrior = " << (int)oM.rhPrior << endl;
	out << "inPrior = " << (int)oM.inPrior << endl;
	out << "muPriorParam = " << oM.muPriorParam[0] << ", " << oM.muPriorParam[1] << endl;
	out << "kaPriorParam = " << oM.kaPriorParam[0] << ", " << oM.kaPriorParam[1] << endl;
	out << "omPriorParam = " << oM.omPriorParam[0] << ", " << oM.omPriorParam[1] << endl;
	out << "rhPriorParam = " << oM.rhPriorParam[0] << ", " << oM.rhPriorParam[1] << endl;
	out << "inPriorParam = " << oM.inPriorParam[0] << ", " << oM.inPriorParam[1] << endl;
	out << "n = " << oM.n << endl;
	out << "L = " << oM.L << endl;
	out << "norders = " << oM.norders << endl;
	out << "pi = ";
	for(i=0;i<oM.pi.size();i++) {
		if(i>0) out << ", ";
		out << oM.pi[i];
	}	out << endl;
	out << "niter = " << oM.niter << endl;
	out << "seed = " << oM.ran->getseed() << endl;
	out << "maxCodon = " << oM.maxCodon << endl;
	out << "memTime = " << oM.memTime << endl;
	out << "compTime = " << oM.compTime << endl;
	out << ENDCONTROLWIZARD;
	if(!oM.writeBinary) {
		out << endl;
		out << oM.eventStart;
		for(i=0;i<oM.niter;i++) {
			out << oM.event[i];
			if(out.fail()) error("ostream& operator<<(ostream&, omegaMapBase&): write failure");
		}
	}
	else {
		#ifdef __GNUC__
			out.setf((_Ios_Fmtflags&)ios::binary);
			out << (oEventStartBinary&)oM.eventStart;
			for(i=0;i<oM.niter;i++) {
				out << (oEventBinary&)oM.event[i];
				if(out.fail()) error("ostream& operator<<(ostream&, omegaMapBase&): write failure");
			}
			out.unsetf((_Ios_Fmtflags&)ios::binary);
			out << endl;
		#else
			out.setf(ios::binary);
			out << (oEventStartBinary&)oM.eventStart;
			for(i=0;i<oM.niter;i++) {
				out << (oEventBinary&)oM.event[i];
				if(out.fail()) error("ostream& operator<<(ostream&, omegaMapBase&): write failure");
			}
			out.unsetf(ios::binary);
			out << endl;
		#endif
	}
	return out;
}

istream& operator>>(istream& in, omegaMapBase& oM) {
	int i;
	ControlWizard con;
	con.coutput = false;
	con.eof_delimiters.push_back('*');
	int temp_initialized,temp_efficientLikelihood,temp_seed,temp_writeBinary,temp_coutput;
	int varPrior[5] = {0,0,0,0,0};
	vector< vector<double> > varPriorParam(5);
	vector<double> temp_prMOVE,temp_pi;
	enum omegaMapBase::Priors* refPrior[5] = {&(oM.muPrior),&(oM.kaPrior),&(oM.omPrior),&(oM.rhPrior),&(oM.inPrior)};
	double* refPriorParam[5] = {oM.muPriorParam,oM.kaPriorParam,oM.omPriorParam,oM.rhPriorParam,oM.inPriorParam};
	con.add_ITEM("omegaVersion",TP_STRING,&(oM.omegaVersion));
	con.add_ITEM("initialized",TP_INT,&temp_initialized);
	con.add_ITEM("coutput",TP_INT,&temp_coutput);
	con.add_ITEM("hapfile",TP_STRING,&(oM.hapfile));
	con.add_ITEM("efficientLikelihood",TP_INT,&temp_efficientLikelihood);
	con.add_ITEM("stickyDebugEval",TP_INT,&(oM.stickyDebugEval));
/*	con.add_ITEM("muPrior",TP_DOUBLE,&(oM.muPrior));
	con.add_ITEM("kappaPrior",TP_DOUBLE,&(oM.kappaPrior));
	con.add_ITEM("omegaPrior",TP_DOUBLE,&(oM.omegaPrior));
	con.add_ITEM("rhoPrior",TP_DOUBLE,&(oM.rhoPrior));
	con.add_ITEM("indelLambdaPrior",TP_DOUBLE,&(oM.indelLambdaPrior));*/
	con.add_ITEM("oBlockPrior",TP_DOUBLE,&(oM.oBlockPrior));
	con.add_ITEM("oBlockPriorRatio",TP_DOUBLE,&(oM.oBlockPriorRatio));
	con.add_ITEM("rBlockPrior",TP_DOUBLE,&(oM.rBlockPrior));
	con.add_ITEM("rBlockPriorRatio",TP_DOUBLE,&(oM.rBlockPriorRatio));
	con.add_ITEM("writeBinary",TP_INT,&temp_writeBinary);
	con.add_ITEM("prMOVE",TP_VEC_DOUBLE,&temp_prMOVE);
	con.add_ITEM("omegaModel",TP_STRING,&(oM.omegaModel));
	con.add_ITEM("rhoModel",TP_STRING,&(oM.rhoModel));
	con.add_ITEM("muPrior",TP_INT,&varPrior[0]);
	con.add_ITEM("kaPrior",TP_INT,&varPrior[1]);
	con.add_ITEM("omPrior",TP_INT,&varPrior[2]);
	con.add_ITEM("rhPrior",TP_INT,&varPrior[3]);
	con.add_ITEM("inPrior",TP_INT,&varPrior[4]);
	con.add_ITEM("muPriorParam",TP_VEC_DOUBLE,&varPriorParam[0]);
	con.add_ITEM("kaPriorParam",TP_VEC_DOUBLE,&varPriorParam[1]);
	con.add_ITEM("omPriorParam",TP_VEC_DOUBLE,&varPriorParam[2]);
	con.add_ITEM("rhPriorParam",TP_VEC_DOUBLE,&varPriorParam[3]);
	con.add_ITEM("inPriorParam",TP_VEC_DOUBLE,&varPriorParam[4]);
	con.add_ITEM("n",TP_INT,&(oM.n));
	con.add_ITEM("L",TP_INT,&(oM.L));
	con.add_ITEM("norders",TP_INT,&(oM.norders));
	con.add_ITEM("pi",TP_VEC_DOUBLE,&temp_pi);
	con.add_ITEM("niter",TP_INT,&(oM.niter));
	con.add_ITEM("seed",TP_INT,&temp_seed);
	con.add_ITEM("maxCodon",TP_INT,&(oM.maxCodon));
	con.add_ITEM("memTime",TP_DOUBLE,&(oM.memTime));
	con.add_ITEM("compTime",TP_DOUBLE,&(oM.compTime));
	con.read_input((std::ifstream&)in);
	if(oM.omegaVersion!="omegaMap0.5") {
		stringstream wText;
		wText << "omegaMapData::inputFromText(): expecting omegaMap0.5 data but received " << oM.omegaVersion;
		warning(wText.str().c_str());
	}
	oM.initialized = !(temp_initialized==0);
	oM.coutput = !(temp_coutput==0);
	oM.efficientLikelihood = !(temp_efficientLikelihood==0);
	oM.writeBinary = !(temp_writeBinary==0);
	oM.prMOVE = Vector<double>((int)temp_prMOVE.size());
	for(i=0;i<oM.prMOVE.size();i++) oM.prMOVE[i] = temp_prMOVE[i];
	oM.pi = Vector<double>((int)temp_pi.size());
	for(i=0;i<oM.pi.size();i++) oM.pi[i] = temp_pi[i];
	if(oM.ran!=NULL) oM.ran->setseed(temp_seed);
	for(i=0;i<5;i++) {
		if(varPrior[i]==0)		*refPrior[i] = omegaMapBase::none;
		else if(varPrior[i]==1) *refPrior[i] = omegaMapBase::fixed;
		else if(varPrior[i]==2) *refPrior[i] = omegaMapBase::uniform;
		else if(varPrior[i]==3) *refPrior[i] = omegaMapBase::improper_uniform;
		else if(varPrior[i]==4) *refPrior[i] = omegaMapBase::exponential;
		else if(varPrior[i]==5) *refPrior[i] = omegaMapBase::gamma;
		else if(varPrior[i]==6) *refPrior[i] = omegaMapBase::exponential_ratio;
		else if(varPrior[i]==7) *refPrior[i] = omegaMapBase::inverse;
		else if(varPrior[i]==8) *refPrior[i] = omegaMapBase::improper_inverse;
		else if(varPrior[i]==9) *refPrior[i] = omegaMapBase::log_normal;
		else error("istream& operator>>(istream&, omegaMapBase&): unrecognised prior");

		if(varPriorParam[i].size()!=2) error("istream& operator>>(istream&, omegaMapBase&): PriorParam not found");
		refPriorParam[i][0] = varPriorParam[i][0];
		refPriorParam[i][1] = varPriorParam[i][1];
	}

	oEvent oEventDefault;
	oEventDefault.accepted = false;
	oEventDefault.alpha = 0.0;
	oEventDefault.likelihood = 0.0;
	oEventDefault.param[0] = oEventDefault.param[1] = 0.0;
	oEventDefault.type = oEvent::DEFAULT;
	oM.event = Vector<oEvent>(oM.niter,oEventDefault);

	if(!oM.writeBinary) {
		in >> oM.eventStart;
		for(i=0;i<oM.niter && !in.eof();i++) {
			in >> oM.event[i];
		}
	}
	else {
		error("Binary representation currently not programmed");
//		int ch = in.get();
//		while(ch!='\n') ch = in.get();
		#ifdef __GNUC__
			in.setf((_Ios_Fmtflags&)ios::binary);
			in >> (oEventStartBinary&)oM.eventStart;
			for(i=0;i<oM.niter && !in.eof();i++) in >> (oEventBinary&)oM.event[i];
			in.unsetf((_Ios_Fmtflags&)ios::binary);
		#else
			in.setf(ios::binary);
			in >> (oEventStartBinary&)oM.eventStart;
			for(i=0;i<oM.niter && !in.eof();i++) in >> (oEventBinary&)oM.event[i];
			in.unsetf(ios::binary);
		#endif
	}
	return in;
}

/* binary versions */
ostream& operator<<(ostream& out, oEventStartBinary& event) {
//	out.write((char*)&event,sizeof(oEventStartBinary));
	out.write((char*)&event.accepted,sizeof(event.accepted));
	out.write((char*)&event.likelihood,sizeof(event.likelihood));
	out.write((char*)&event.alpha,sizeof(event.alpha));
	out.write((char*)&event.nblocks,sizeof(event.nblocks));
	out.write((char*)&event.nrblocks,sizeof(event.nrblocks));
	out.write((char*)&event.mu,sizeof(event.mu));
	out.write((char*)&event.kappa,sizeof(event.kappa));
	out.write((char*)&event.indelLambda,sizeof(event.indelLambda));
	int i;
	for(i=0;i<event.nblocks;i++) out.write((char*)&event.oBlockStart[i],sizeof(event.oBlockStart[i]));
	for(i=0;i<event.nblocks;i++) out.write((char*)&event.omega[i],sizeof(event.omega[i]));
	for(i=0;i<event.nrblocks;i++) out.write((char*)&event.rBlockStart[i],sizeof(event.rBlockStart[i]));
	for(i=0;i<event.nrblocks;i++) out.write((char*)&event.rho[i],sizeof(event.rho[i]));
	return out;
}

istream& operator>>(istream& in, oEventStartBinary& event) {
//	in.read((char*)&event,sizeof(oEventStartBinary));
	in.read((char*)&event.accepted,sizeof(event.accepted));
	in.read((char*)&event.likelihood,sizeof(event.likelihood));
	in.read((char*)&event.alpha,sizeof(event.alpha));
	in.read((char*)&event.nblocks,sizeof(event.nblocks));
	in.read((char*)&event.nrblocks,sizeof(event.nrblocks));
	in.read((char*)&event.mu,sizeof(event.mu));
	in.read((char*)&event.kappa,sizeof(event.kappa));
	in.read((char*)&event.indelLambda,sizeof(event.indelLambda));
	int i;
	for(i=0;i<event.nblocks;i++) in.read((char*)&event.oBlockStart[i],sizeof(event.oBlockStart[i]));
	for(i=0;i<event.nblocks;i++) in.read((char*)&event.omega[i],sizeof(event.omega[i]));
	for(i=0;i<event.nrblocks;i++) in.read((char*)&event.rBlockStart[i],sizeof(event.rBlockStart[i]));
	for(i=0;i<event.nrblocks;i++) in.read((char*)&event.rho[i],sizeof(event.rho[i]));
	return in;
}

ostream& operator<<(ostream& out, oEventBinary& event) {
//	out.write((char*)&event,sizeof(oEventBinary));
	out.write((char*)&event.type,sizeof(event.type));
	out.write((char*)&event.accepted,sizeof(event.accepted));
	out.write((char*)&event.param[0],sizeof(event.param[0]));
	out.write((char*)&event.param[1],sizeof(event.param[1]));
	out.write((char*)&event.likelihood,sizeof(event.likelihood));
	out.write((char*)&event.alpha,sizeof(event.alpha));
	return out;
}

istream& operator>>(istream& in, oEventBinary& event) {
//	in.read((char*)&event,sizeof(oEventBinary));
	in.read((char*)&event.type,sizeof(event.type));
	in.read((char*)&event.accepted,sizeof(event.accepted));
	in.read((char*)&event.param[0],sizeof(event.param[0]));
	in.read((char*)&event.param[1],sizeof(event.param[1]));
	in.read((char*)&event.likelihood,sizeof(event.likelihood));
	in.read((char*)&event.alpha,sizeof(event.alpha));
	return in;
}

/* non-binary (safer) versions */
ostream& operator<<(ostream& out, oEventStart& event) {
	out << (int)event.type << tab << (int)event.accepted << tab << event.likelihood << tab << event.alpha << tab << event.nblocks << tab << event.nrblocks << tab << event.mu << tab << event.kappa << tab << event.indelLambda << tab;
	int i;
	for(i=0;i<event.nblocks;i++) out << event.oBlockStart[i] << tab;
	for(i=0;i<event.nblocks;i++) out << event.omega[i] << tab;
	for(i=0;i<event.nrblocks;i++) out << event.rBlockStart[i] << tab;
	for(i=0;i<event.nrblocks;i++) out << event.rho[i] << tab;
	return out;
}

ostream& operator<<(ostream& out, oEvent& event) {
	out << (int)event.type << tab << (int)event.accepted << tab << event.param[0] << tab << event.param[1] << tab << event.likelihood << tab << event.alpha << tab;
	return out;
}

istream& operator>>(istream& in, oEvent& event) {
	int typeDump;
	in >> typeDump;
	if(in.fail()) error("istream& operator>>(istream&, oEvent&): read failure");
	switch(typeDump) {
	case 0: event.type = oEvent::DEFAULT; break;
	case 1: event.type = oEvent::START; break;
	case 2: event.type = oEvent::CHANGE_OBLOCK; break;
	case 3: event.type = oEvent::EXTEND_OBLOCK; break;
	case 4: event.type = oEvent::SPLIT_OBLOCK; break;
	case 5: event.type = oEvent::MERGE_OBLOCK; break;
	case 6: event.type = oEvent::CHANGE_MU; break;
	case 7: event.type = oEvent::CHANGE_KAPPA; break;
	case 8: event.type = oEvent::CHANGE_RBLOCK; break;
	case 9: event.type = oEvent::EXTEND_RBLOCK; break;
	case 10: event.type = oEvent::SPLIT_RBLOCK; break;
	case 11: event.type = oEvent::MERGE_RBLOCK; break;
	case 12: event.type = oEvent::CHANGE_INDEL; break;
	case 13: event.type = oEvent::CHANGE_ORDER; break;
	default:
		event.type = oEvent::DEFAULT;
		warning("istream& operator>>(istream&, oEvent&): Event type not recognised. Set to default");
	}
	in >> event.accepted >> event.param[0] >> event.param[1] >> event.likelihood >> event.alpha;
	if(in.fail()) {
		//error("istream& operator>>(istream&, oEvent&): read failure");
		in.clear(); // just try to keep going
	}
	return in;
}

istream& operator>>(istream& in, oEventStart& event) {
	int typeDump;
	bool acceptedDump;
	in >> typeDump >> acceptedDump >> event.likelihood >> event.alpha >> event.nblocks >> event.nrblocks >> event.mu >> event.kappa >> event.indelLambda;
	if(in.fail()) error("istream& operator>>(istream&, oEventStart&): read failure");
	if(event.nblocks<=0) error("istream& operator>>(istream&, oEventStart&): nblocks must be positive");
	if(event.nrblocks<=0) error("istream& operator>>(istream&, oEventStart&): nrblocks must be positive");
	event.oBlockStart = Vector<int>(event.nblocks,-1);
	event.rBlockStart = Vector<int>(event.nrblocks,-1);
	event.omega = Vector<double>(event.nblocks,-1.0);
	event.rho = Vector<double>(event.nrblocks,-1.0);
	int i;
	for(i=0;i<event.nblocks;i++) in >> event.oBlockStart[i];
	for(i=0;i<event.nblocks;i++) in >> event.omega[i];
	for(i=0;i<event.nrblocks;i++) in >> event.rBlockStart[i];
	for(i=0;i<event.nrblocks;i++) in >> event.rho[i];
	if(in.fail()) error("istream& operator>>(istream&, oEventStart&): read failure");
	return in;
}
