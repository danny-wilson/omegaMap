/********************************************/
/*	DNA.h 23rd February 2005				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _DNA_H_
#define _DNA_H_

#pragma warning(disable: 4786)

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <myutils.h>
#include <map>

using namespace std;
using namespace myutils;

class DNA {
public:
	vector<string> label;
	vector<string> sequence;
	int nseq;
	int lseq;
	vector<double> ntimes;
	bool coutput;

protected:
	vector<int> _uniqueHaps;
	vector<int> _sites;
	LowerTriangularMatrix<int> _B_B;
	vector<int> _M;
	vector<double> _F;
	vector<double> _four;
	LowerTriangularMatrix< vector<double> > _G;
	LowerTriangularMatrix<double> _A;
	LowerTriangularMatrix<double> B_B;
	LowerTriangularMatrix<double> _CC;
	Matrix<double> _D;

public:
	DNA() {
		coutput = false;
	}
	DNA(const char* filename) {
		coutput = false;
		readFASTA(filename);
	}
	DNA& readFASTA(const char* filename) {
		ifstream* in = new ifstream(filename);
		if(in->is_open()==false) {
			string errmsg = "DNA::readFASTA(): File ";
			errmsg += string(filename);
			errmsg += " not found";
			error(errmsg.c_str());
		}

		int str;
		nseq = 0;
		while(!in->eof()) {
			str = in->get();
			if((char)str=='>') {
				++nseq;
			}
		}
		in->close();
		delete in;
		if(coutput) cout << "Read in " << nseq << " sequence" << endl;
		if(nseq==0) {
			lseq = 0;
			return *this;
		}

		in = new ifstream(filename);
		if(in->is_open()==false)error("File not found");

		lseq = 0;
		string junk;
		while(!in->eof()) {
			str = in->get();
			if((char)str=='>') {
				getline(*in,junk);
				while(!in->eof()) {
					str = in->get();
					if((char)str=='>') break;
					if((char)str!='\n' && (char)str!='\r')
						++lseq;
				}
				if(coutput) cout << "Sequences are " << lseq << " long" << endl;
				break;
			}
		}
		in->close();
		delete in;

		string blank(lseq,' ');
		sequence.resize(nseq,blank);
		label.resize(nseq);
		ntimes.resize(nseq,0.0);

		in = new ifstream(filename);
		int NSEQ = 0; int LSEQ = 0;
		while(true) {
			str = in->get();
			if(in->eof()) error("Cannot find sequences!");
			if((char)str=='>') {
				getline(*in,label[NSEQ]);
				break;
			}
		}
		while(true) {
			str = in->get();
			if(in->eof()) break;
			if(LSEQ<lseq)
				sequence[NSEQ][LSEQ] = (char)str;
			if((char)str!='\n' && (char)str!='\r')
				++LSEQ;
			if((char)str=='>') {
				++NSEQ;
				getline(*in,label[NSEQ]);
				LSEQ=0;
			}
		}
		in->close();

		if(coutput) for(NSEQ=0;NSEQ<nseq;NSEQ++) {
			cout << label[NSEQ] << endl;
			cout << sequence[NSEQ] << endl;
		}
		
		return *this;
	}
	DNA& writeFASTA(const char* filename) {
		ofstream fout(filename);
		int n,pos;
		for(n=0;n<nseq;n++)
		{
			fout << ">" << label[n] << endl;
			for(pos=0;pos<lseq;pos++)
				fout << sequence[n][pos];
			fout << endl;
		}
		fout.close();
		return *this;
	}
	/*Subscript operator*/
	inline string& operator[](int pos) {
		return sequence[pos];
	}
	/*Resize wipes current entries*/
	DNA& resize(const int NSEQ, const int LSEQ) {
		if(NSEQ<0) error("DNA: NSEQ must be non-negative");
		if(LSEQ<0) error("DNA: LSEQ must be non-negative");
		nseq = NSEQ;
		lseq = LSEQ;

		string blank(lseq,' ');
		string empty = "";
		
		sequence.resize(nseq,blank);
		label.resize(nseq,empty);
		ntimes.resize(nseq);
		int i;
		for(i=0;i<nseq;i++) {
			sequence[i] = blank;
			label[i] = empty;
			ntimes[i] = 0.0;
		}
		return *this;
	}
	DNA& resize(const int NSEQ, const string& str) {
		if(NSEQ<0) error("DNA: NSEQ must be non-negative");
		nseq = NSEQ;
		sequence.resize(nseq,str);
		lseq = (int)sequence[0].size();
		string empty = "";
		label.resize(nseq,empty);
		ntimes.resize(nseq,0.0);
		return *this;
	}
	DNA& clear() {
		string blank(lseq,' ');
		sequence.resize(nseq,blank);
		string empty = "";
		label.resize(nseq,empty);
		ntimes.resize(nseq,0.0);
		return *this;
	}
public:
	/* Number of segregating sites */
	double S() {
		double result = 0.0;
		if(nseq==0) return 0.0;

		int i,j;
		for(j=0;j<lseq;j++) {
			char hap = sequence[0][j];
			for(i=1;i<nseq;i++)
				if(sequence[i][j]!=hap) {
					++result;
					break;
				}
		}

		return result;
	}

	/* Number of unique haplotypes */
	double H() {
		int result = 1;

		if(nseq==0) return 0.0;
		vector<int> uniqueHaps(nseq,-1);
		uniqueHaps[0] = 0;

		int i,ii,j;
		bool unique;
		for(i=1;i<nseq;i++) {
			unique = true;
			for(ii=0;ii<result;ii++) {
				for(j=0;j<lseq;j++)
					if(sequence[i][j]!=sequence[uniqueHaps[ii]][j]) break;
				if(j==lseq) unique = false;
			}
			if(unique==true) {
				uniqueHaps[result] = i;
				++result;
			}
		}

		return (double)result;
	}

	/* Average number of pairwise differences */
	double pi() {
		double result = 0.0;

		int i,j,k;
		for(i=0;i<nseq;i++)
			for(j=0;j<i;j++)
				for(k=0;k<lseq;k++)
					result += (sequence[i][k]==sequence[j][k]) ? 0.0 : 1.0;
		result *= 2.0/(double)(nseq)/(double)(nseq-1);

		return result;
	}

	/* Variance in number of pairwise differences */
	double Varpi() {
		double E,EE,pi;
		int i,j,k;

		E = EE = 0.0;
		for(i=0;i<nseq;i++)
			for(j=0;j<i;j++) {
				pi = 0.0;
				for(k=0;k<lseq;k++)
					pi += (sequence[i][k]==sequence[j][k]) ? 0.0 : 1.0;
				E += pi;
				EE += pi*pi;
			}
		E *= 2.0/(double)(nseq)/(double)(nseq-1);
		EE *= 2.0/(double)(nseq)/(double)(nseq-1);

		double result = EE - E*E;
		return result;
	}

	double Tajima() {
		double D = 0.0;
		int i,j,k,n,L;
		n = nseq;
		L = lseq;
		double a1,a2,b1,b2,c1,c2,e1,e2,khat,S;
		bool segregating;
		khat = S = 0.0;
		for(k=0;k<L;k++) {
			segregating = false;
			for(i=0;i<n;i++)
				for(j=0;j<i;j++)
					if(sequence[i][k]!=sequence[j][k]) {
						++khat;
						segregating = true;
					}
			if(segregating) ++S;
		}
		if(S==0) return 0.0;
		khat /= (double)(n*(n-1)/2);
		a1 = a2 = 0.0;
		for(i=1;i<=n-1;i++) {
			a1 += 1./(double)i;
			a2 += 1./(double)(i*i);
		}
		b1 = (double)(n+1)/(double)(3*(n-1));
		b2 = (double)(2*(n*n+n+3))/(double)(9*n*(n-1));
		c1 = b1 - 1./a1;
		c2 = b2 - (double)(n+2)/a1/(double)(n) + a2/a1/a1;
		e1 = c1/a1;
		e2 = c2/(a1*a1+a2);
		D = (khat - S/a1)/sqrt(e1*S+e2*S*(S-1.));
		return D;
	}
	/* This function counts the average number of pairwise differences, where the matrix
	   diff defines those differences using 0 = different or 1 = identical.
	   If diff is the identity matrix then this function is equivalent to double pi(). */
	double pi(Matrix<int> &diff, map<char,int> &chmap) {
		double result = 0.0;

		int i,j,k;
		for(i=0;i<nseq;i++)
			for(j=0;j<i;j++)
				for(k=0;k<lseq;k++)
					result += (diff[chmap[(char)sequence[i][k]]][chmap[(char)sequence[j][k]]]==0);
		result *= 2.0/(double)(nseq)/(double)(nseq);

		return result;
	}
	/* Hudson and Kaplan's Rm, the minimum # recombinations. See Myers and Griffiths(2003)*/
	double Rm() {
		if(nseq==0) return 0.0;
		if(lseq==0) return 0.0;

		/* Determine which sites are biallelic segregating */
		_sites = vector<int>(lseq,0);
		int i,j,k;
		int S = 0;
		char hap0,hap1;
		bool segregating;
		for(j=0;j<lseq;j++) {
			segregating = false;
			hap0 = sequence[0][j];
			for(i=1;i<nseq;i++) {
				if(!segregating && sequence[i][j]!=hap0) {
					segregating = true;
					hap1 = sequence[i][j];
				}
				else if(segregating && sequence[i][j]!=hap0 && sequence[i][j]!=hap1) {
					segregating = false;	// define segregating only for biallelic sites
					break;
				}
			}
			if(segregating) {
				_sites[S] = j;
				++S;
			}
		}
		if(S<2) return 0.0;

		/* Calculate the compatibility matrix */
		_B_B = LowerTriangularMatrix<int>(S,0);	// so j>=k always
		// _B_B[j][k] = 0 for compatible, 1 for incompatible
		bool comb[3];
		for(j=0;j<S;j++)
			for(k=0;k<j;k++)
			{
				hap0 = sequence[0][_sites[j]];
				hap1 = sequence[0][_sites[k]];
				comb[0] = false;				// hap0  hap1'
				comb[1] = false;				// hap0' hap1
				comb[2] = false;				// hap0' hap1'
				for(i=1;i<nseq;i++) {
					if(sequence[i][_sites[j]]==hap0 && sequence[i][_sites[k]]!=hap1) comb[0] = true;
					if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]==hap1) comb[1] = true;
					if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]!=hap1) comb[2] = true;
					if(comb[0] && comb[1] && comb[2]) break;
				}
				_B_B[j][k] = (comb[0] && comb[1] && comb[2]) ? 1 : 0;			
			}

		/* Calculate the dynamic programming partition matrix */
		_M = vector<int>(S,0);
		int maxM = 0;
		_M[S-1] = 0;
		_M[S-2] = _B_B[S-1][S-2];
		for(i=S-3;i>=0;i--) {
			_M[i] = _B_B[i+1][i] + _M[i+1];
			for(k=i+2;k<S;k++) if(_B_B[k][i]+_M[k]>_M[i]) _M[i] = _B_B[k][i]+_M[k];
		}

		return (double)_M[0];
	}
	void RecCorrelations(double* result) {
		RecCorrelations(result,true);
	}
	void RecCovariances(double* result) {
		RecCorrelations(result,false);
	}
	void RecCorrelations(double* result, bool normalize) {
		result[0] = result[1] = result[2] = 0.0;

		if(nseq==0) return;
		if(lseq==0) return;

		/* Determine which sites are biallelic segregating */
		_sites = vector<int>(lseq,0);
		int i,j,k;
		int S = 0;
		char hap0,hap1;
		bool segregating;
		for(j=0;j<lseq;j++) {
			segregating = false;
			hap0 = sequence[0][j];
			for(i=1;i<nseq;i++) {
				if(!segregating && sequence[i][j]!=hap0) {
					segregating = true;
					hap1 = sequence[i][j];
				}
				else if(segregating && sequence[i][j]!=hap0 && sequence[i][j]!=hap1) {
					segregating = false;	// define segregating only for biallelic sites
					break;
				}
			}
			if(segregating) {
				_sites[S] = j;
				++S;
			}
		}
		if(S<3) return;
		
		/* Calculate frequency statistics */
		_F = vector<double>(S,1.0);							/* _F is the marginal frequency of hap0 at site j */
		for(j=0;j<S;j++) {
			hap0 = sequence[0][_sites[j]];
			for(i=1;i<nseq;i++)
				if(sequence[i][_sites[j]]==hap0) _F[j]++;
			_F[j] /= (double)nseq;
		}

		_four = vector<double>(4,0.0);							/* _G[j][k] is the frequency of AB (_G[j][k][0]),	*/
		_G = LowerTriangularMatrix< vector<double> >(S,_four);	/* Ab (1), aB (2), ab (3) for sites j and k		*/
		for(j=0;j<S;j++)
		  for(k=0;k<j;k++) {
			  hap0 = sequence[0][_sites[j]];
			  hap1 = sequence[0][_sites[k]];
			  for(i=0;i<nseq;i++) {
				  if(sequence[i][_sites[j]]==hap0 && sequence[i][_sites[k]]==hap1) ++_G[j][k][0];
				  else if(sequence[i][_sites[j]]==hap0 && sequence[i][_sites[k]]!=hap1) ++_G[j][k][1];
				  else if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]==hap1) ++_G[j][k][2];
				  else if(sequence[i][_sites[j]]!=hap0 && sequence[i][_sites[k]]!=hap1) ++_G[j][k][3];
				  else warning("Unexpected choice");
			  }
			  for(i=0;i<4;i++) _G[j][k][i] /= (double)nseq;
		  }
  
		/* Calculate LD statistics for pairs of sites */
		_A = LowerTriangularMatrix<double>(S,0.0);			//	rsq
		B_B = LowerTriangularMatrix<double>(S,0.0);			//	Dprime
		_CC = LowerTriangularMatrix<double>(S,0.0);			//	G4
		_D = Matrix<double>(S,S,0.0);

		double temp;
		for(i=0;i<S;i++) {
			for(j=0;j<i;j++) {
				temp = _G[i][j][0] - _F[i]*_F[j];
				_A[i][j] = pow(temp,2.0)/(_F[i]*(1.-_F[i])*_F[j]*(1.-_F[j]));
				B_B[i][j] = (temp < 0.0) ? -temp/MIN(_F[i]*_F[j],(1.-_F[i])*(1.-_F[j])) : temp/MIN(_F[i]*(1.-_F[j]),(1.-_F[i])*_F[j]);
				_CC[i][j] = (_G[i][j][0]>0.0 && _G[i][j][1]>0.0 && _G[i][j][2]>0.0 && _G[i][j][3]>0.0) ? 1.0 : 0.0;
				_D[i][j] = _D[j][i] = _sites[i] - _sites[j];
			}
		}

		double  E[4] = {0.0,0.0,0.0,0.0};
		double EE[4] = {0.0,0.0,0.0,0.0};
		double ED[3] = {0.0,0.0,0.0};
		int ctr;
		for(i=0,ctr=0;i<S;i++)
			for(j=0;j<i;j++,ctr++) {
				E[0] += _A[i][j]; E[1] += B_B[i][j]; E[2] += _CC[i][j]; E[3] += _D[i][j];
				EE[0] += _A[i][j]*_A[i][j]; EE[1] += B_B[i][j]*B_B[i][j]; EE[2] += _CC[i][j]*_CC[i][j]; EE[3] += _D[i][j]*_D[i][j];
				ED[0] += _A[i][j]*_D[i][j]; ED[1] += B_B[i][j]*_D[i][j]; ED[2] += _CC[i][j]*_D[i][j];
			}

		if(normalize)	// Calculate correlation
			for(k=0;k<3;k++)
				result[k] = (ED[k]-E[k]*E[3]/(double)(ctr))/sqrt((EE[k]-E[k]*E[k]/(double)(ctr)))/sqrt((EE[3]-E[3]*E[3]/(double)(ctr)));
		else			// Calculate covariance
			for(k=0;k<3;k++)
				result[k] = (ED[k]-E[k]*E[3]/(double)(ctr))/(double)(ctr);

		return;
	}
};

#endif // _DNA_H_
