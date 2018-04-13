/************************************/
/*	mcmc.cpp 17th November 2005		*/
/*	Part of omegaMap v0.5			*/
/*	(c) Danny Wilson.				*/
/*	www.danielwilson.me.uk			*/
/************************************/

#include "omegaMap.h"
#include <math.h>
#include <fstream>
#include <time.h>

/***************************/
/***************************/
/**     MCMC routines     **/
/***************************/
/***************************/

omegaMap& omegaMap::go(Random &r, const char* inifile) {
	return go(0,0,r,inifile);
}

omegaMap& omegaMap::go(const int argc, const char* argv[], Random &r, const char* inifile) {
	initialize(argc,argv,inifile,r);
	bool no_datafile = (datafile=="");
	clock_t start = clock();
	if(outfile!="") {
		ofstream out(outfile.c_str());
		headings(out);
		iter = -1;
		record(out,(oEvent&)eventStart);
		if(coutput) cout << "Completed 0 of " << niter << " iterations";
		for(iter=0;iter<niter;iter++) {
			propose();
#ifdef _DEBUG
			debug();
#endif
			if((iter+1)%thinning==0) record(out,event[iter]);
			if(no_datafile) --event.element;
			if(coutput) cout << "\rCompleted " << iter + 1 << " of " << niter << " iterations";
		}
		out.close();
	}
	else {
		for(iter=0;iter<niter;iter++) {
			propose();
			if(no_datafile) --event.element;
			if(coutput) cout << "\rCompleted " << iter + 1 << " of " << niter << " iterations";
		}
	}
	compTime = (double)(clock()-start)/CLOCKS_PER_SEC/60.0;
	cout << "\rCompleted " << niter << " of " << niter << " in " << compTime << " minutes" << endl;
	if(!no_datafile) {
		ofstream data(datafile.c_str());
		data << (omegaMapBase&)*this;
		data.close();
	}
	if(no_datafile) event.element += niter;
	return *this;
}

omegaMap& omegaMap::propose() {
//	redraw();
	double U = ran->U();
	if(U<=prMOVE[6]) return propose_change_mu();
	U -= prMOVE[6];
	if(U<=prMOVE[7]) return propose_change_kappa();
	U -= prMOVE[7];
	if(U<=prMOVE[12]) return propose_change_indelLambda();
	U -= prMOVE[12];
	if(U<=prMOVE[2]) return propose_change_oBlock();
	U -= prMOVE[2];
	if(U<=prMOVE[3]) return propose_extend_oBlock();
	U -= prMOVE[3];
	if(U<=prMOVE[8]) return propose_change_rBlock();
	U -= prMOVE[8];
	if(U<=prMOVE[9]) return propose_extend_rBlock();
	U -= prMOVE[9];
	if(U<=prMOVE[13]) return propose_change_order();
	U -= prMOVE[13];
	
	double h[6];
	if(oBlockPrior!=1.0) {
		h[0] = omegaSplitRate(nblocks);
		h[1] = omegaMergeRate(nblocks);
		h[2] = h[0] + h[1];
		h[0] *= prMOVE[4]/h[2];
		h[1] *= prMOVE[4]/h[2];
		if(U<=h[0]) return propose_split_oBlock();
		U -= h[0];
		if(U<=h[1]) return propose_merge_oBlock();
		U -= h[1];
	}
	if(rBlockPrior!=1.0) {
		h[3] = rhoSplitRate(nrblocks);
		h[4] = rhoMergeRate(nrblocks);
		h[5] = h[3] + h[4];
		h[3] *= prMOVE[10]/h[5];
		h[4] *= prMOVE[10]/h[5];
		if(U<=h[3]) return propose_split_rBlock();
		U -= h[3];
		if(U<=h[4]) return propose_merge_rBlock();
		U -= h[4];
	}
	warning("omegaMap::propose_move(): no move proposed");
	return propose();
}

/***************************/
/***************************/
/**  omega MCMC routines  **/
/***************************/
/***************************/

/* relative rate at which splits are proposed */
 double omegaMap::omegaSplitRate(const int nblo) {
	return MIN((double) 1.0,(double)(L-nblo) / (double)(nblo)
		* oBlockPrior / (1.0-oBlockPrior));
}

/* relative rate at which merges are proposed */
 double omegaMap::omegaMergeRate(const int nblo) {
	return MIN(1.0,(double)(nblo-1) / (double)(L-nblo+1)
		* (1.0-oBlockPrior) / oBlockPrior);
}

/* Proposal distribution for changing a block's omega */
omegaMap& omegaMap::propose_change_oBlock(){
	event[iter].type = oEvent::CHANGE_OBLOCK;

	int bNum = ran->discrete(0,nblocks-1);
	oBlock *oB = block[0];
	int i;
	for(i=1;i<=bNum;i++) oB = oB->_3prime;

	double U = ran->uniform(-1.0,1.0);
	diagonalizeMatrix(*oMatTemp[0],oB->oMat->mu,oB->oMat->kappa,oB->oMat->omega*exp(-U));
	//double omegaPrime = ran->exponential(1./omegaPrior);
	//diagonalizeMatrix(*oMatTemp[0],oB->oMat->mu,oB->oMat->kappa,omegaPrime);

	SWAP(oB->oMat,oMatTemp[0]);	/* revert this if rejected */

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - oB->start;
		rightMargin = oB->end - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(oB->start<L-oB->end) {
				goForward = true;
				forward(alphaMargin,oB->end+1);
				backward(oB->end+1,betaMargin);
				betaMargin = oB->end;
				newLikelihood = likelihood(oB->end+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,oB->start-1);
				alphaMargin = oB->start;
				backward(oB->start-1,betaMargin);
				newLikelihood = likelihood(oB->start-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>oB->start) alphaMargin = oB->start;
			forward(alphaMargin,oB->end+1);
			newLikelihood = likelihood(oB->end+1);
		}
		else {
			goForward = false;
			if(betaMargin<oB->end) betaMargin = oB->end;
			backward(oB->start-1,betaMargin);
			newLikelihood = likelihood(oB->start-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].param[0] = oB->start;
	event[iter].param[1] = oB->oMat->omega;
	event[iter].likelihood = newLikelihood;

	lnAlpha = newLikelihood - oldLikelihood + log(oB->oMat->omega) - log(oMatTemp[0]->omega) 
		+ logRatioPriors_change(oMatTemp[0]->omega,oB->oMat->omega,omPrior,omPriorParam);
	//lnAlpha = newLikelihood - oldLikelihood +log(oB->oMat->omega) -
	//	log(oMatTemp[0]->omega) - omegaPrior*(oB->oMat->omega - oMatTemp[0]->omega);
	//lnAlpha = newLikelihood - oldLikelihood +log(oB->oMat->omega) - log(oMatTemp[0]->omega); // UNIFORM
	//lnAlpha = 0.0;
	event[iter].alpha = lnAlpha;
	U = log(ran->U());
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? oB->end+1 : oB->start-1;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = oB->start-1;
		else betaMargin = oB->end+1;
		SWAP(oB->oMat,oMatTemp[0]);
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for extending a block 5' or 3' */
omegaMap& omegaMap::propose_extend_oBlock() {
	event[iter].type = oEvent::EXTEND_OBLOCK;

	int bNum = ran->discrete(0,nblocks-1);
	oBlock *oB = block[0];
	int i;
	for(i=1;i<=bNum;i++) oB = oB->_3prime;

	/* true for left, false for right */
	bool goLeft = ran->bernoulliTF(0.5);
	int newpos,leftUpdate,rightUpdate,oldpos;
	if(goLeft) {
		if(oB->_5prime==0) {
			event[iter].param[0] = oB->start;
			event[iter].param[1] = -1;
			event[iter].accepted = false;
			return *this;
		}
		newpos = oB->start-(ran->geometric(0.9)+1);
		if(newpos<=oB->_5prime->start) {
			event[iter].param[0] = oB->start;
			event[iter].param[1] = -1;
			event[iter].accepted = false;
			return *this;
		}

		event[iter].param[0] = oB->start;
		event[iter].param[1] = newpos;

		leftUpdate = newpos;
		rightUpdate = oB->start - 1;

		/* will need to revert the following if rejected */
		_block[newpos].start = newpos;
		_block[newpos].end = oB->end;
		if(oB->_5prime!=0) oB->_5prime->_3prime = &(_block[newpos]);
		_block[newpos]._5prime = oB->_5prime;
		_block[newpos]._5prime->end = newpos - 1;
		if(oB->_3prime!=0) oB->_3prime->_5prime = &(_block[newpos]);
		_block[newpos]._3prime = oB->_3prime;
		SWAP(_block[newpos].oMat,oB->oMat);
		for(i=newpos;i<=oB->end;i++) block[i] = &(_block[newpos]);
	}
	else {
		if(oB->_3prime==0) {
			event[iter].param[0] = L;
			event[iter].param[1] = -2;
			event[iter].accepted = false;
			return *this;
		}
		newpos = oB->end+(ran->geometric(0.9)+1);
		if(newpos>=oB->_3prime->end) {
			event[iter].param[0] = oB->end + 1;
			event[iter].param[1] = -2;
			return *this;
		}

		event[iter].param[0] = oB->end + 1;
		event[iter].param[1] = newpos + 1;

		leftUpdate = oB->end + 1;
		rightUpdate = newpos;
		
		/* will need to revert the following if rejected */
		oldpos = oB->end+1;
		oB->end = newpos;
		oB->_3prime = &(_block[newpos+1]);
		_block[newpos+1].start = newpos + 1;
		_block[newpos+1].end = _block[oldpos].end;
		_block[newpos+1]._5prime = oB;
		_block[newpos+1]._3prime = _block[oldpos]._3prime;
		if(_block[newpos+1]._3prime!=0) _block[newpos+1]._3prime->_5prime = &(_block[newpos+1]);
		SWAP(_block[newpos+1].oMat,_block[oldpos].oMat);
		for(i=oldpos;i<=newpos;i++) block[i] = oB;
		for(i=newpos+1;i<=_block[newpos+1].end;i++) block[i] = &(_block[newpos+1]);
	}

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	/* Hastings ratio = 1 and ratio of priors = 1 */
	double lnAlpha = newLikelihood - oldLikelihood;
	double U = log(ran->U());
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate - 1;
		else betaMargin = rightUpdate + 1;

		if(goLeft) {
			if(oB->_5prime!=0) oB->_5prime->_3prime = oB;
			oB->_5prime->end = oB->start - 1;
			if(oB->_3prime!=0) oB->_3prime->_5prime = oB;
			SWAP(_block[newpos].oMat,oB->oMat);
			for(i=newpos;i<oB->start;i++) block[i] = oB->_5prime;
			for(i=oB->start;i<=oB->end;i++) block[i] = oB;
		}
		else {
			oB->end = oldpos-1;
			oB->_3prime = &(_block[oldpos]);
			if(_block[newpos+1]._3prime!=0) _block[newpos+1]._3prime->_5prime = &(_block[oldpos]);
			SWAP(_block[newpos+1].oMat,_block[oldpos].oMat);
			for(i=oldpos;i<=_block[oldpos].end;i++) block[i] = &(_block[oldpos]);
		}
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for splitting an existing block */
omegaMap& omegaMap::propose_split_oBlock() {
	/* Choose a position uniformly on [0,L-2] excluding current splits */
	/* The split occurs to the right of the chosen numbered site */
	event[iter].type = oEvent::SPLIT_OBLOCK;

	const int K = (L - 1) - (nblocks - 1);
	int pos;
	oBlock *oB;
	/* The constant 0.2 determines the efficiency of the search routine */
	if((double)K/(double)(L-1)>=0.2) {
		pos = ran->discrete(0,L-2);
		while(true) {
			if(pos != block[pos]->end) break;
			pos = ran->discrete(0,L-2);
		}
	}
	else {
		pos = ran->discrete(0,K-1);
		oB = block[0];
		while(oB!=0) {
			if(pos>=oB->end) ++pos;
			else break;
			oB = oB->_3prime;
		}
	}
	/* now the split occurs to the left of pos */
	++pos;
	oB = block[pos];
	if(oB==&(_block[pos])){
		cout << "Major problem" << endl;
	}
	int i;

	event[iter].param[0] = pos;
	double U = ran->U();
	event[iter].param[1] = U;
	if(U==0.0 || U==1.0) {
		event[iter].accepted = false;			
		return *this;
	}

	_block[pos].start = pos;
	_block[pos].end = oB->end;
	oB->end = pos - 1;

	_block[pos]._5prime = oB;
	_block[pos]._3prime = oB->_3prime;
	oB->_3prime = &(_block[pos]);
	if(_block[pos]._3prime!=0)
		_block[pos]._3prime->_5prime = &(_block[pos]);
	for(i=pos;i<=_block[pos].end;i++) block[i] = &(_block[pos]);

	/* create the new oMatrices */
	double a = (double)(_block[pos].start - oB->start)/(double)(_block[pos].end - oB->start + 1);
	/*const double h[3] = {
		oB->oMat->omega,
		oB->oMat->omega * pow(U/(1.-U),1.-a),
		oB->oMat->omega * pow((1.-U)/U,a)
	};*/
	const double rat = U/(1.-U);
	const double h[3] = {
		oB->oMat->omega,
		oB->oMat->omega * pow(rat,1.-a),
		oB->oMat->omega * pow(rat,-a)
	};
	diagonalizeMatrix(*oMatTemp[0],oB->oMat->mu,oB->oMat->kappa,h[1]);
	diagonalizeMatrix(*oMatTemp[1],oB->oMat->mu,oB->oMat->kappa,h[2]);

	/* may need to revert these two lines of code */
	SWAP(oB->oMat,oMatTemp[0]);
	SWAP(_block[pos].oMat,oMatTemp[1]);
	
	int leftUpdate = oB->start;
	int rightUpdate = _block[pos].end;

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	const double lnRatioPriors = logRatioPriors_split(h[0],h[1],h[2],omPrior,omPriorParam,oBlockPrior);
	//const double lnRatioPriors = log(oBlockPrior/(1.-oBlockPrior) * omegaPrior)
	//	- omegaPrior * (h[1] + h[2] - h[0]);
	//const double lnRatioPriors = log(oBlockPrior/(1.-oBlockPrior)); // QUASI-UNIFORM
	//const double RatioPriors = exp(lnRatioPriors);
	const double eta[4] = {
		omegaSplitRate(nblocks+1),
		omegaMergeRate(nblocks+1),
		omegaSplitRate(nblocks),
		omegaMergeRate(nblocks)
	};
	const double Hastings = (double)K / (double)nblocks * eta[1] * (eta[2]+eta[3]) / eta[2] / (eta[0]+eta[1]);
	const double Jacobian = pow(h[1]+h[2],2.)/h[0];
	const double lnAlpha = newLikelihood - oldLikelihood
		+ lnRatioPriors + log(Hastings * Jacobian);
	U = log(ran->U());
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		++nblocks;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate - 1;
		else betaMargin = rightUpdate + 1;
		oB->end = _block[pos].end;
		oB->_3prime = _block[pos]._3prime;
		if(oB->_3prime!=0)
			oB->_3prime->_5prime = oB;
		for(i=pos;i<=oB->end;i++) block[i] = oB;
		SWAP(oB->oMat,oMatTemp[0]);
		SWAP(_block[pos].oMat,oMatTemp[1]);
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for merging two adjacent blocks */
omegaMap& omegaMap::propose_merge_oBlock() {
	event[iter].type = oEvent::MERGE_OBLOCK;
	const int K = (L - 1) - (nblocks - 2);
	int bNum = ran->discrete(0,nblocks-2); /* merge to right so only allow first nblocks-1 to do so */
	oBlock *oB1 = block[0];
	int i;
	for(i=1;i<=bNum;i++) oB1 = oB1->_3prime;
	oBlock *oB2 = oB1->_3prime;
	event[iter].param[0] = oB2->start;

	/* need to specify how to get new omega */
	double U;

	oB1->end = oB2->end;
	oB1->_3prime = oB2->_3prime;
	if(oB1->_3prime!=0)
		oB1->_3prime->_5prime = oB1;
	for(i=oB2->start;i<=oB2->end;i++) block[i] = oB1;

	double a = (double)(oB2->start - oB1->start)/(double)(oB2->end - oB1->start + 1);
	const double h[3] = {
		//pow(oB1->oMat->omega,1.-a)*pow(oB2->oMat->omega,a),
		pow(oB1->oMat->omega,a)*pow(oB2->oMat->omega,1.-a),
		oB1->oMat->omega,
		oB2->oMat->omega
	};
	event[iter].param[1] = h[0];
	diagonalizeMatrix(*oMatTemp[0],oB1->oMat->mu,oB1->oMat->kappa,h[0]);

	/* revert this if rejected */
	SWAP(oB1->oMat,oMatTemp[0]);

	int leftUpdate = oB1->start;
	int rightUpdate = oB2->end;

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	const double lnRatioPriors = logRatioPriors_merge(h[1],h[2],h[0],omPrior,omPriorParam,oBlockPrior);
	//const double lnRatioPriors = log((1.-oBlockPrior)/oBlockPrior / omegaPrior)
	//	- omegaPrior * (h[0] - h[1] - h[2]);
	//const double lnRatioPriors = log((1.-oBlockPrior)/oBlockPrior); // QUASI-UNIFORM. IN FACT NOT BAYESIAN (LIKELIHOOD ONLY)
	const double eta[4] = {
		omegaSplitRate(nblocks-1),
		omegaMergeRate(nblocks-1),
		omegaSplitRate(nblocks),
		omegaMergeRate(nblocks)
	};
	const double Hastings = (double)(nblocks-1) / (double)K * eta[0] * (eta[2]+eta[3]) / eta[3] / (eta[0]+eta[1]);
	const double Jacobian = h[0] / pow(h[1]+h[2],2.);
	const double lnAlpha = newLikelihood - oldLikelihood
		+ lnRatioPriors + log(Hastings * Jacobian);
	event[iter].alpha = lnAlpha;
	U = log(ran->U());
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		--nblocks;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate - 1;
		else betaMargin = rightUpdate + 1;

		oB1->end = oB2->start - 1;
		oB1->_3prime = oB2;
		if(oB2->_3prime!=0)
			oB2->_3prime->_5prime = oB2;
		for(i=oB2->start;i<=oB2->end;i++) block[i] = oB2;
		SWAP(oB1->oMat,oMatTemp[0]);
		event[iter].accepted = false;
	}

	return *this;
}

/***************************/
/***************************/
/**   rho MCMC routines   **/
/***************************/
/***************************/

/* relative rate at which splits are proposed */
 double omegaMap::rhoSplitRate(const int nblo) {
	return MIN(1.0,(double)(L-nblo-1) / (double)(nblo)
		* rBlockPrior / (1.0-rBlockPrior));
}

/* relative rate at which merges are proposed */
 double omegaMap::rhoMergeRate(const int nblo) {
	return MIN(1.0,(double)(nblo-1) / (double)(L-nblo)
		* (1.0-rBlockPrior) / rBlockPrior);
}

/* Proposal distribution for changing a block's rho */
omegaMap& omegaMap::propose_change_rBlock(){
	event[iter].type = oEvent::CHANGE_RBLOCK;

	int bNum = ran->discrete(0,nrblocks-1);
	rBlock *rB = rblock[0];
	int i;
	for(i=1;i<=bNum;i++) rB = rB->_3prime;

	double U = ran->uniform(-1.0,1.0);
	double oldRho = rB->rho;
	rB->rho = oldRho*exp(-U); /* revert this if rejected */

	double leftUpdate = rB->start;
	double rightUpdate = rB->end;
	++leftUpdate;	// convert from rho units to omega units (rightUpdate is unchanged)

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-1-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].param[0] = rB->start;
	event[iter].param[1] = rB->rho;
	event[iter].likelihood = newLikelihood;

	lnAlpha = newLikelihood - oldLikelihood + log(rB->rho) - log(oldRho)
		+ logRatioPriors_change(oldRho,rB->rho,rhPrior,rhPriorParam);
	//lnAlpha = newLikelihood - oldLikelihood +log(rB->rho) -
	//	log(oldRho) - rhoPrior*(rB->rho - oldRho);
	//lnAlpha = newLikelihood - oldLikelihood +log(rB->rho) - log(oldRho); // UNIFORM
	//lnAlpha = 0.0;
	event[iter].alpha = lnAlpha;
	U = log(ran->U());
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate-1;
		else betaMargin = rightUpdate+1;
		rB->rho = oldRho;
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for extending a rho block 5' or 3' */
omegaMap& omegaMap::propose_extend_rBlock() {
	event[iter].type = oEvent::EXTEND_RBLOCK;

	int bNum = ran->discrete(0,nrblocks-1);
	rBlock *rB = rblock[0];
	int i;
	for(i=1;i<=bNum;i++) rB = rB->_3prime;

	/* true for left, false for right */
	bool goLeft = ran->bernoulliTF(0.5);
	int newpos,leftUpdate,rightUpdate,oldpos;
	if(goLeft) {
		if(rB->_5prime==0) {
			event[iter].param[0] = rB->start;
			event[iter].param[1] = -1;
			event[iter].accepted = false;
			return *this;
		}
		newpos = rB->start-(ran->geometric(0.9)+1);
		if(newpos<=rB->_5prime->start) {
			event[iter].param[0] = rB->start;
			event[iter].param[1] = -1;
			event[iter].accepted = false;
			return *this;
		}

		event[iter].param[0] = rB->start;
		event[iter].param[1] = newpos;

		leftUpdate = newpos;
		rightUpdate = rB->start - 1;

		/* will need to revert the following if rejected */
		_rblock[newpos].start = newpos;
		_rblock[newpos].end = rB->end;
		if(rB->_5prime!=0) rB->_5prime->_3prime = &(_rblock[newpos]);
		_rblock[newpos]._5prime = rB->_5prime;
		_rblock[newpos]._5prime->end = newpos - 1;
		if(rB->_3prime!=0) rB->_3prime->_5prime = &(_rblock[newpos]);
		_rblock[newpos]._3prime = rB->_3prime;
		SWAP(_rblock[newpos].rho,rB->rho);
		for(i=newpos;i<=rB->end;i++) rblock[i] = &(_rblock[newpos]);
	}
	else {
		if(rB->_3prime==0) {
			event[iter].param[0] = L;
			event[iter].param[1] = -2;
			event[iter].accepted = false;
			return *this;
		}
		newpos = rB->end+(ran->geometric(0.9)+1);
		if(newpos>=rB->_3prime->end) {
			event[iter].param[0] = rB->end + 1;
			event[iter].param[1] = -2;
			return *this;
		}

		event[iter].param[0] = rB->end + 1;
		event[iter].param[1] = newpos + 1;

		leftUpdate = rB->end + 1;
		rightUpdate = newpos;
		
		/* will need to revert the following if rejected */
		oldpos = rB->end+1;
		rB->end = newpos;
		rB->_3prime = &(_rblock[newpos+1]);
		_rblock[newpos+1].start = newpos + 1;
		_rblock[newpos+1].end = _rblock[oldpos].end;
		_rblock[newpos+1]._5prime = rB;
		_rblock[newpos+1]._3prime = _rblock[oldpos]._3prime;
		if(_rblock[newpos+1]._3prime!=0) _rblock[newpos+1]._3prime->_5prime = &(_rblock[newpos+1]);
		SWAP(_rblock[newpos+1].rho,_rblock[oldpos].rho);
		for(i=oldpos;i<=newpos;i++) rblock[i] = rB;
		for(i=newpos+1;i<=_rblock[newpos+1].end;i++) rblock[i] = &(_rblock[newpos+1]);
	}

	++leftUpdate;	// convert from rho units to omega units (rightUpdate is unchanged)

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-1-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	/* Hastings ratio = 1 and ratio of priors = 1 */
	double lnAlpha = newLikelihood - oldLikelihood;
	double U = log(ran->U());
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate - 1;
		else betaMargin = rightUpdate + 1;

		if(goLeft) {
			if(rB->_5prime!=0) rB->_5prime->_3prime = rB;
			rB->_5prime->end = rB->start - 1;
			if(rB->_3prime!=0) rB->_3prime->_5prime = rB;
			SWAP(_rblock[newpos].rho,rB->rho);
			for(i=newpos;i<rB->start;i++) rblock[i] = rB->_5prime;
			for(i=rB->start;i<=rB->end;i++) rblock[i] = rB;
		}
		else {
			rB->end = oldpos-1;
			rB->_3prime = &(_rblock[oldpos]);
			if(_rblock[newpos+1]._3prime!=0) _rblock[newpos+1]._3prime->_5prime = &(_rblock[oldpos]);
			SWAP(_rblock[newpos+1].rho,_rblock[oldpos].rho);
			for(i=oldpos;i<=_rblock[oldpos].end;i++) rblock[i] = &(_rblock[oldpos]);
		}
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for splitting an existing rho block */
omegaMap& omegaMap::propose_split_rBlock() {
	/* Choose a position uniformly on [0,L-3] excluding current splits */
	/* The split occurs to the right of the chosen numbered site */
	event[iter].type = oEvent::SPLIT_RBLOCK;

	const int K = (L - 2) - (nrblocks - 1);
	int pos;
	rBlock *rB;
	/* The constant 0.2 determines the efficiency of the search routine */
	if((double)K/(double)(L-2)>=0.2) {
		pos = ran->discrete(0,L-3);
		while(true) {
			if(pos != rblock[pos]->end) break;
			pos = ran->discrete(0,L-3);
		}
	}
	else {
		pos = ran->discrete(0,K-1);
		rB = rblock[0];
		while(rB!=0) {
			if(pos>=rB->end) ++pos;
			else break;
			rB = rB->_3prime;
		}
	}
	/* now the split occurs to the left of pos */
	++pos;
	rB = rblock[pos];
	if(rB==&(_rblock[pos])){
		cout << "Major problem" << endl;
	}
	int i;

	event[iter].param[0] = pos;
	double U = ran->U();
	event[iter].param[1] = U;
	if(U==0.0 || U==1.0) {
		event[iter].accepted = false;			
		return *this;
	}

	_rblock[pos].start = pos;
	_rblock[pos].end = rB->end;
	rB->end = pos - 1;

	_rblock[pos]._5prime = rB;
	_rblock[pos]._3prime = rB->_3prime;
	rB->_3prime = &(_rblock[pos]);
	if(_rblock[pos]._3prime!=0)
		_rblock[pos]._3prime->_5prime = &(_rblock[pos]);
	for(i=pos;i<=_rblock[pos].end;i++) rblock[i] = &(_rblock[pos]);

	/* create the new oMatrices */
	double a = (double)(_rblock[pos].start - rB->start)/(double)(_rblock[pos].end - rB->start + 1);
/*	const double h[3] = {
		rB->rho,
		rB->rho * pow(U/(1.-U),1.-a),
		rB->rho * pow((1.-U)/U,a)
	};*/
	const double rat = U/(1.-U);
	const double h[3] = {
		rB->rho,
		rB->rho * pow(rat,1.-a),
		rB->rho * pow(rat,-a)
	};

	double oldRho = rB->rho;
	/* may need to revert these two lines of code */
	rB->rho = h[1];
	_rblock[pos].rho = h[2];
	
	int leftUpdate = rB->start;
	int rightUpdate = _rblock[pos].end;
	++leftUpdate;	// convert from rho units to omega units (rightUpdate is unchanged)

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	const double lnRatioPriors = logRatioPriors_split(h[0],h[1],h[2],rhPrior,rhPriorParam,rBlockPrior);
	//const double lnRatioPriors = log(rBlockPrior/(1.-rBlockPrior) * rhoPrior)
	//	- rhoPrior * (h[1] + h[2] - h[0]);
	//const double lnRatioPriors = log(rBlockPrior/(1.-rBlockPrior)); // QUASI-UNIFORM
	//const double RatioPriors = exp(lnRatioPriors);
	const double eta[4] = {
		rhoSplitRate(nrblocks+1),
		rhoMergeRate(nrblocks+1),
		rhoSplitRate(nrblocks),
		rhoMergeRate(nrblocks)
	};
	const double Hastings = (double)K / (double)nrblocks * eta[1] * (eta[2]+eta[3]) / eta[2] / (eta[0]+eta[1]);
	const double Jacobian = pow(h[1]+h[2],2.)/h[0];
	const double lnAlpha = newLikelihood - oldLikelihood
		+ lnRatioPriors + log(Hastings * Jacobian);
	U = log(ran->U());
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		++nrblocks;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate - 1;
		else betaMargin = rightUpdate + 1;
		rB->end = _rblock[pos].end;
		rB->_3prime = _rblock[pos]._3prime;
		if(rB->_3prime!=0)
			rB->_3prime->_5prime = rB;
		for(i=pos;i<=rB->end;i++) rblock[i] = rB;
		rB->rho = oldRho;
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for merging two adjacent rho blocks */
omegaMap& omegaMap::propose_merge_rBlock() {
	event[iter].type = oEvent::MERGE_RBLOCK;
	const int K = (L - 2) - (nrblocks - 2);
	int bNum = ran->discrete(0,nrblocks-2); /* merge to right so only allow first nrblocks-1 to do so */
	rBlock *rB1 = rblock[0];
	int i;
	for(i=1;i<=bNum;i++) rB1 = rB1->_3prime;
	rBlock *rB2 = rB1->_3prime;
	event[iter].param[0] = rB2->start;

	/* need to specify how to get new omega */
	double U;

	rB1->end = rB2->end;
	rB1->_3prime = rB2->_3prime;
	if(rB1->_3prime!=0)
		rB1->_3prime->_5prime = rB1;
	for(i=rB2->start;i<=rB2->end;i++) rblock[i] = rB1;

	double a = (double)(rB2->start - rB1->start)/(double)(rB2->end - rB1->start + 1);
	const double h[3] = {
	//	pow(rB1->rho,1.-a)*pow(rB2->rho,a),
		pow(rB1->rho,a)*pow(rB2->rho,1.-a),
		rB1->rho,
		rB2->rho
	};
	event[iter].param[1] = h[0];

	double oldRho = rB1->rho;
	/* revert this if rejected */
	rB1->rho = h[0];

	int leftUpdate = rB1->start;
	int rightUpdate = rB2->end;
	++leftUpdate;	// convert from rho units to omega units (rightUpdate is unchanged)

	bool goForward;
	int leftMargin,rightMargin;
	double newLikelihood;
	if(efficientLikelihood) {
		leftMargin = alphaMargin - leftUpdate;
		rightMargin = rightUpdate - betaMargin;
		if(leftMargin<=0 && rightMargin<=0) {
			if(leftUpdate<L-rightUpdate) {
				goForward = true;
				forward(alphaMargin,rightUpdate+1);
				backward(rightUpdate+1,betaMargin);
				betaMargin = rightUpdate;
				newLikelihood = likelihood(rightUpdate+1);
			}
			else {
				goForward = false;
				forward(alphaMargin,leftUpdate-1);
				alphaMargin = leftUpdate;
				backward(leftUpdate-1,betaMargin);
				newLikelihood = likelihood(leftUpdate-1);
			}
		}
		else if(leftMargin<=rightMargin) {
			goForward = true;
			if(alphaMargin>leftUpdate) alphaMargin = leftUpdate;
			forward(alphaMargin,rightUpdate+1);
			newLikelihood = likelihood(rightUpdate+1);
		}
		else {
			goForward = false;
			if(betaMargin<rightUpdate) betaMargin = rightUpdate;
			backward(leftUpdate-1,betaMargin);
			newLikelihood = likelihood(leftUpdate-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	const double lnRatioPriors = logRatioPriors_merge(h[1],h[2],h[0],rhPrior,rhPriorParam,rBlockPrior);
	//const double lnRatioPriors = log((1.-rBlockPrior)/rBlockPrior / rhoPrior)
	//	- rhoPrior * (h[0] - h[1] - h[2]);
	//const double lnRatioPriors = log((1.-rBlockPrior)/rBlockPrior); // QUASI-UNIFORM.
	const double eta[4] = {
		rhoSplitRate(nrblocks-1),
		rhoMergeRate(nrblocks-1),
		rhoSplitRate(nrblocks),
		rhoMergeRate(nrblocks)
	};
	const double Hastings = (double)(nrblocks-1) / (double)K * eta[0] * (eta[2]+eta[3]) / eta[3] / (eta[0]+eta[1]);
	const double Jacobian = h[0] / pow(h[1]+h[2],2.);
	const double lnAlpha = newLikelihood - oldLikelihood
		+ lnRatioPriors + log(Hastings * Jacobian);
	event[iter].alpha = lnAlpha;
	U = log(ran->U());
	if(lnAlpha>=0.0 || lnAlpha>=U) /* then accept */ {
		alphaMargin = betaMargin = (goForward) ? rightUpdate+1 : leftUpdate-1;
		oldLikelihood = newLikelihood;
		--nrblocks;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = leftUpdate - 1;
		else betaMargin = rightUpdate + 1;

		rB1->end = rB2->start - 1;
		rB1->_3prime = rB2;
		if(rB2->_3prime!=0)
			rB2->_3prime->_5prime = rB2;
		for(i=rB2->start;i<=rB2->end;i++) rblock[i] = rB2;
		rB1->rho = oldRho;
		event[iter].accepted = false;
	}

	return *this;
}

/***************************/
/***************************/
/**  other MCMC routines  **/
/***************************/
/***************************/

/* Proposal distribution for changing mu */
omegaMap& omegaMap::propose_change_mu() {
	event[iter].type = oEvent::CHANGE_MU;
	oBlock *oB = &(_block[0]);
	double muPrime = oB->oMat->mu*exp(-ran->uniform(-1.0,1.0));
	while(oB!=0) {
		diagonalizeMatrix(*oMatTemp[oB->start],muPrime,oB->oMat->kappa,oB->oMat->omega);
		SWAP(oB->oMat,oMatTemp[oB->start]);
		oB = oB->_3prime;
	}
	event[iter].param[0] = muPrime;

	bool goForward;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		if(alphaMargin<L-betaMargin) {
			goForward = true;
			alphaMargin = 0;
			forward(0,L);
			newLikelihood = likelihood(L);
		}
		else {
			goForward = false;
			betaMargin = L-1;
			backward(-1,L-1);
			newLikelihood = likelihood(-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	lnAlpha = newLikelihood - oldLikelihood + log(muPrime) - log(oMatTemp[0]->mu)
		+ logRatioPriors_change(oMatTemp[0]->mu,muPrime,muPrior,muPriorParam);
	//lnAlpha = newLikelihood - oldLikelihood + log(muPrime) -
	//	log(oMatTemp[0]->mu) - muPrior*(muPrime - oMatTemp[0]->mu);
	//lnAlpha = newLikelihood - oldLikelihood + log(muPrime) - log(oMatTemp[0]->mu); // UNIFORM
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=log(ran->U())) /* then accept */ {
		alphaMargin = betaMargin = (goForward)? L-1 : 0;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = 0;
		else betaMargin = L-1;
		oB = &(_block[0]);
		while(oB!=0) {
			SWAP(oB->oMat,oMatTemp[oB->start]);
			oB = oB->_3prime;
		}
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for changing kappa */
omegaMap& omegaMap::propose_change_kappa() {
	event[iter].type = oEvent::CHANGE_KAPPA;
	oBlock *oB = &(_block[0]);
	double kappaPrime = oB->oMat->kappa*exp(-ran->uniform(-1.0,1.0));
	while(oB!=0) {
		diagonalizeMatrix(*oMatTemp[oB->start],oB->oMat->mu,kappaPrime,oB->oMat->omega);
		SWAP(oB->oMat,oMatTemp[oB->start]);
		oB = oB->_3prime;
	}
	event[iter].param[0] = kappaPrime;

	bool goForward;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		if(alphaMargin<L-betaMargin) {
			goForward = true;
			alphaMargin = 0;
			forward(0,L);
			newLikelihood = likelihood(L);
		}
		else {
			goForward = false;
			betaMargin = L-1;
			backward(-1,L-1);
			newLikelihood = likelihood(-1);
		}
	}
	else {
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	lnAlpha = newLikelihood - oldLikelihood + log(kappaPrime) - log(oMatTemp[0]->kappa)
		+ logRatioPriors_change(oMatTemp[0]->kappa,kappaPrime,kaPrior,kaPriorParam);
	//lnAlpha = newLikelihood - oldLikelihood + log(kappaPrime) -
	//	log(oMatTemp[0]->kappa) - kappaPrior*(kappaPrime - oMatTemp[0]->kappa);
	//lnAlpha = newLikelihood - oldLikelihood + log(kappaPrime) - log(oMatTemp[0]->kappa); // UNIFORM
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=log(ran->U())) /* then accept */ {
		alphaMargin = betaMargin = (goForward)? L-1 : 0;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = 0;
		else betaMargin = L-1;
		oB = &(_block[0]);
		while(oB!=0) {
			SWAP(oB->oMat,oMatTemp[oB->start]);
			oB = oB->_3prime;
		}
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for changing rho */
/* assumes that rho is the same along the sequence */
/* omegaMap& propose_change_rho() {
	event[iter].type = oEvent::CHANGE_RHO;
	oBlock *oB = &(_block[0]);
	double oldRho = rho[0];
	double rhoPrime = oldRho*exp(-ran->uniform(-1.0,1.0));
	rho = Vector<double>(L,rhoPrime);
	event[iter].param[0] = rhoPrime;

	bool goForward;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		if(alphaMargin<L-betaMargin) {
			goForward = true;
			alphaMargin = 0;
			forward(0,L);
			newLikelihood = likelihood(L);
		}
		else {
			goForward = false;
			betaMargin = L-1;
			backward(-1,L-1);
			newLikelihood = likelihood(-1);
		}
	}
	else {
		//************ likelihood debug ************
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		//******************************************
	}

	event[iter].likelihood = newLikelihood;
	lnAlpha = newLikelihood - oldLikelihood + log(rhoPrime) -
		log(oldRho) - rhoPrior*(rhoPrime - oldRho);
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=log(ran->U())) { // then accept
		alphaMargin = betaMargin = (goForward)? L-1 : 0;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else { // else reject
		if(goForward) alphaMargin = 0;
		else betaMargin = L-1;
		rho = Vector<double>(L,oldRho);
		event[iter].accepted = false;
	}

	return *this;
}*/

/* Proposal distribution for changing indelLambda */
omegaMap& omegaMap::propose_change_indelLambda() {
	if(maxCodon!=62) return *this;
	event[iter].type = oEvent::CHANGE_INDEL;
	double oldIndelLambda = indelLambda;
	double indelLambdaPrime = indelLambda*exp(-ran->uniform(-1.0,1.0));
	indelLambda = indelLambdaPrime;
	oBlock* oB = &(_block[0]);
	int i,j,k;
	while(oB!=0) {
		//copy matrix to temp, then swap it, then set gamma[][61] and gamma[61][] to -1.0
/*			oMatTemp[oB->start]->gamma = oB->oMat->gamma;
		for(i=0;i<61;i++)
			for(j=0;j<=i;j++)
				for(k=0;k<n;k++)
					(*oMatTemp[oB->start]->gamma)[i][j][k] = (*oB->oMat->gamma)[i][j][k];*/
		for(i=0;i<62;i++)
			for(j=0;j<=i;j++)
				for(k=0;k<n;k++)
					(*oMatTemp[oB->start]->gamma2)[i][j][k] = -1.0;

/*			oMatTemp[oB->start]->R = oB->oMat->R;
		oMatTemp[oB->start]->lambda = oB->oMat->lambda;
		oMatTemp[oB->start]->mu = oB->oMat->mu;
		oMatTemp[oB->start]->kappa = oB->oMat->kappa;
		oMatTemp[oB->start]->omega = oB->oMat->omega;
		SWAP(oB->oMat,oMatTemp[oB->start]);*/
		SWAP(oB->oMat->gamma2,oMatTemp[oB->start]->gamma2);
		oB = oB->_3prime;
	}
	event[iter].param[0] = indelLambdaPrime;

	bool goForward;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		if(alphaMargin<L-betaMargin) {
			goForward = true;
			alphaMargin = 0;
			forward(0,L);
			newLikelihood = likelihood(L);
		}
		else {
			goForward = false;
			betaMargin = L-1;
			backward(-1,L-1);
			newLikelihood = likelihood(-1);
		}
	}
	else {
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	lnAlpha = newLikelihood - oldLikelihood + log(indelLambdaPrime) - log(oldIndelLambda)
		+ logRatioPriors_change(oldIndelLambda,indelLambdaPrime,inPrior,inPriorParam);
	//lnAlpha = newLikelihood - oldLikelihood + log(indelLambdaPrime) -
	//	log(oldIndelLambda) - indelLambdaPrior*(indelLambdaPrime - oldIndelLambda);
	//lnAlpha = newLikelihood - oldLikelihood + log(indelLambdaPrime) - log(oldIndelLambda); // UNIFORM
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=log(ran->U())) /* then accept */ {
		alphaMargin = betaMargin = (goForward)? L-1 : 0;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = 0;
		else betaMargin = L-1;
		indelLambda = oldIndelLambda;
		oB = &(_block[0]);
		while(oB!=0) {
			//SWAP(oB->oMat,oMatTemp[oB->start]);
			SWAP(oB->oMat->gamma2,oMatTemp[oB->start]->gamma2);
			oB = oB->_3prime;
		}
		event[iter].accepted = false;
	}

	return *this;
}

/* Proposal distribution for changing the orderings */
omegaMap& omegaMap::propose_change_order(){
	event[iter].type = oEvent::CHANGE_ORDER;

	Matrix< Vector<int>* > Hold = H;
	//event[iter].param[0] = (double)ran->getidum();	/* the ordering can be reconstructed from this */
	//setOrders();

	int swap0 = ran->discrete(0,n-1);
	int swap1 = ran->discrete(0,n-2);
	if(swap1==swap0) swap1 = n-1;
	event[iter].param[0] = (double)swap0;
	event[iter].param[1] = (double)swap1;

	int ord;
	for(ord=0;ord<norders;ord++) swap(H[ord][swap0],H[ord][swap1]);

	bool goForward;
	double newLikelihood,lnAlpha;
	if(efficientLikelihood) {
		if(alphaMargin<L-betaMargin) {
			goForward = true;
			alphaMargin = 0;
			forward(0,L);
			newLikelihood = likelihood(L);
		}
		else {
			goForward = false;
			betaMargin = L-1;
			backward(-1,L-1);
			newLikelihood = likelihood(-1);
		}
	}
	else{
		/************ likelihood debug ************/
		forward(0,stickyDebugEval);
		backward(stickyDebugEval,L-1);
		newLikelihood = likelihood(stickyDebugEval);
		/******************************************/
	}

	event[iter].likelihood = newLikelihood;
	lnAlpha = newLikelihood - oldLikelihood;
	event[iter].alpha = lnAlpha;
	if(lnAlpha>=0.0 || lnAlpha>=log(ran->U())) /* then accept */ {
		alphaMargin = betaMargin = (goForward)? L-1 : 0;
		oldLikelihood = newLikelihood;
		event[iter].accepted = true;
	}
	else /* else reject */ {
		if(goForward) alphaMargin = 0;
		else betaMargin = L-1;
		//H = Hold;
		for(ord=0;ord<norders;ord++) swap(H[ord][swap0],H[ord][swap1]);
		event[iter].accepted = false;
	}

	return *this;
}

