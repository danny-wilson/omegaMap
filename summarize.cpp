/************************************************
	summarize takes output from one or more
	omegaMap MCMC runs and performs simple
	analyses of the (merged) data
************************************************/

#include <myutils.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <functional>
#include <sstream>
#include <controlwizard.h>
//#include <float.h>

using namespace myutils;
using namespace std;

double mean(vector<double> &x);
double geomean(vector<double> &x);
void HPD(vector<double> &x, const double alpha, double &lo, double &hi);
void geoHPD(vector<double> &x, const double alpha, double &lo, double &hi);
void summary(double &mean, double &lo, double &hi, double &count, const double alpha, vector<double> &x, vector<double> &y, vector<double> &z);

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

int main(const int argc, const char* argv[]) {
	if(argc<3) error("SYNTAX: burn-in file1 [file2 [file3 [... ");
	const int burnin = atoi(argv[1]);
	const int nfiles = argc-2;
	vector<ifstream*> file(nfiles);
	int i;
	for(i=0;i<nfiles;i++) {
		file[i] = new ifstream(argv[2+i]);
		if(!file[i]->is_open()) {
			stringstream errText;
			errText << "Cannot find " << argv[2+i] << endl;
			error(errText.str().c_str());
		}
	}

	vector<int> ncol(nfiles,0);
	vector<string> fieldname(0);
	map<string, int> fieldnum;
	string this_fieldname = "";
	int character = 0;
	i = 0;
	while(!file[i]->eof()) {
		character = tolower(file[i]->get());
		if(character=='\t'||character=='\r'||character=='\n') {
			fieldnum[this_fieldname]=fieldname.size();
			fieldname.push_back(this_fieldname);
			this_fieldname = "";
			if(character=='\r'||character=='\n') {
				character = file[i]->peek();
				if(character=='\r'||character=='\n') file[i]->get();
				break;
			}
		}
		else this_fieldname += (char)character;
	}
	for(i=0;i<fieldname.size();i++) if(fieldnum[fieldname[i]]!=i) error("Duplicate headings found");
	for(i=1;i<nfiles;i++) {
		int which_field = 0;
		while(!file[i]->eof()) {
			character = tolower(file[i]->get());
			if(character=='\t'||character=='\r'||character=='\n') {
				if(fieldnum[this_fieldname]!=which_field || fieldname[which_field]!=this_fieldname) {
					stringstream errText;
					errText << "Inconsistency between " << argv[2] << " and " << argv[2+i] << ".\n";
					errText << "Field number " << which_field << " corresponds, in the two files, to\n";
					errText << fieldname[which_field] << " and " << this_fieldname << " respectively.\n";
					errText << "Whereas field name " << this_fieldname << " corresponds, in the two files, to\n";
					errText << fieldnum[this_fieldname] << " and " << which_field << " respectively.\n";
					/*int fNUM = fieldnum[this_fieldname];
					string fNAME = fieldname[which_field];
					bool firST = fieldnum[this_fieldname]!=which_field;
					bool secOND = fieldname[which_field]!=this_fieldname;
					for(int j=0;j<fieldname.size();j++) cout << fieldnum[fieldname[j]] << " ";
					cout << endl;*/
					error(errText.str().c_str());
				}
				this_fieldname = "";
				if(character=='\r'||character=='\n') {
					character = file[i]->peek();
					if(character=='\r'||character=='\n') file[i]->get();
					break;
				}
				++which_field;
			}
			else this_fieldname += (char)character;
		}
	}
	
	/* Details of memory-allocation */
	const int blocksize = 10000;
	const int nfields = fieldname.size();
	unsigned long nan[2]={0xffffffff, 0x7fffffff};
	double NOTANUMBER = *( double* )nan;
	const int fiter = fieldnum["iter"];
	vector<double> defaultVec(blocksize,0.0);
	vector< vector<double> > A(nfields,defaultVec);
	int row = 0, col;
	for(i=0;i<nfiles;i++) {
		int frow = 0;
		col = 0;
		string word = "";
		while(!file[i]->eof()) {
			character = file[i]->get();
			if(character=='\t'||character=='\r'||character=='\n') {
				if(col==nfields) {
					stringstream errText;
					errText << "In file " << argv[2+i] << " row " << frow << endl;
					errText << "Number of rows exceeded expected number of " << nfields << endl;
					error(errText.str().c_str());
				}
				if(word=="1.#INF" || word=="-1.#INF" || word=="1.#QNAN" || word=="-1.#QNAN" || word=="1.#IND" || word=="-1.#IND")
					A[col][row] = NOTANUMBER;
				else A[col][row] = atof(word.c_str());
				word="";
				++col;
				if(character=='\r'||character=='\n') {
					character = file[i]->peek();
					if(character=='\r'||character=='\n') file[i]->get();
					++frow;
					if(A[fiter][row]>=burnin)++row;
					col = 0;
					if(row==A[0].size()) {
						int j;
						for(j=0;j<nfields;j++) A[j].resize(A[j].size()+blocksize);
					}
				}
			}
			else word += (char)character;
		}
		/*if(frow<burnin) {
			stringstream errText;
			errText << "In file " << argv[2+i] << " there were " << frow << " rows," << endl;
			errText << "which is less than the burn-in of " << burnin << endl;
			error(errText.str().c_str());
		}*/
	}
	const int nrows = row;
	for(i=0;i<nfields;i++) A[i].resize(nrows);
	vector<double> synonymous_theta(nrows);
	const int mu = fieldnum["mu"];
	const int kappa = fieldnum["kappa"];
	const int indel = fieldnum["phi"];
	for(i=0;i<nrows;i++) synonymous_theta[i] = A[mu][i]/155.*(6.+5.*A[kappa][i]);
	const double alpha = 0.05;
	vector<double> y(nrows),z(nrows);
	const int oBeg = fieldnum["omega0"];
	const int oEnd = fieldnum["rblockstart0"]-1;
	const int rBeg = fieldnum["rho0"];
	const int rEnd = rBeg + (oEnd-oBeg-1);
    
	printf("\t\tomega\t\tPosterior\t\trho\t\t\t\tSynonymous theta\t\t\tkappa\t\t\tphi\t\n");
	printf("Site\tLower\tPoint\tHigher\tprob of +ve\tLower\tPoint\tHigher\t\tLower\tPoint\tHigher\tLower\tPoint\tHigher\tLower\tPoint\tHigher\n");
	printf("\t95%c HPD\testimate\t95%c HPD\tselection\t95%c HPD\testimate\t95%c HPD\t\t95%c HPD\testimate\t95%c HPD\t95%c HPD\testimate\t95%c HPD\t95%c HPD\testimate\t95%c HPD\n",'%','%','%','%','%','%','%','%','%','%');
	i=oBeg; int j=rBeg;
	double olo,oav,ohi,pos,rlo,rav,rhi,rpos;
	double slo,sav,shi,spos,klo,kav,khi,kpos,plo,pav,phi,ppos;
	summary(oav,olo,ohi,pos,alpha,A[i],y,z);
	summary(rav,rlo,rhi,rpos,alpha,A[j],y,z);
	summary(sav,slo,shi,spos,alpha,synonymous_theta,y,z);
	summary(kav,klo,khi,kpos,alpha,A[kappa],y,z);
	summary(pav,plo,phi,ppos,alpha,A[indel],y,z);
	printf("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",0,olo,oav,ohi,pos,rlo,rav,rhi,slo,sav,shi,klo,kav,khi,plo,pav,phi);
	i++; j++;
	for(;i<oEnd;i++,j++) {
		summary(oav,olo,ohi,pos,alpha,A[i],y,z);
		summary(rav,rlo,rhi,rpos,alpha,A[j],y,z);
		printf("%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",i-oBeg,olo,oav,ohi,pos,rlo,rav,rhi);
	}
	for(;i<=oEnd;i++) {
		summary(oav,olo,ohi,pos,alpha,A[i],y,z);
		printf("%d\t%g\t%g\t%g\t%g\n",i-oBeg,olo,oav,ohi,pos);
	}
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

void geoHPD(vector<double> &x, const double alpha, double &lo, double &hi) {
	lo = hi = 0.0;
	int n = x.size();
	double m = MAX(1.0,ceil(alpha*(double)n));
	//vector<double> y = x;
	sort(x.begin(),x.end());
	vector<double> y((int)m,0.0);
	int i;
	for(i=0;i<(int)m;i++) y[i] = log(x[n-(int)m+i])-log(x[i]);
	vector<double> z((int)m,0.0);
	for(i=0;i<(int)m;i++) z[i] = (double)i;
	/* sort z by y */
	sort(z.begin(),z.end(),sort_by_vector_double(y));
	int j = (int)z[0];
	lo = x[j];
	hi = x[n-m+j];
}

double mean(vector<double> &x) {
	double num = 0.0, denom = 0.0;
	vector<double>::iterator X;
	for(X=x.begin();X!=x.end();X++)
		if((*X)==(*X)) {
			num += *X;
			denom += 1.0;
		}
	return num/denom;
}

double geomean(vector<double> &x) {
	double num = 0.0, denom = 0.0;
	vector<double>::iterator X;
	double val;
	for(X=x.begin();X!=x.end();X++) {
		val = log(*X);
		if(val==val) {
			num += val;
			denom += 1.0;
		}
	}
	return exp(num/denom);
}

/*	Results are stored in mean, lo, hi and count
	The column to summarize is stored in x
	alpha is for the 100(1-alpha)% HPD intervals
	Vectors y and z are for temporary storage	*/
void summary(double &mean, double &lo, double &hi, double &count, const double alpha, vector<double> &x, vector<double> &y, vector<double> &z) {
	mean = lo = hi = count = 0.0;
	double val;
	int n = x.size();
	double m = MAX(1.0,ceil(alpha*(double)n));
	sort(x.begin(),x.end());
	y.resize((int)m);
	z.resize((int)m);
	int i;
	for(i=0;i<(int)m;i++) {
		val = log(x[i]);
		mean += val;
		if(val>=0.0) ++count;
		y[i] = log(x[n-(int)m+i])-val;
		z[i] = (double)i;
	}
	for(;i<n;i++) {
		val = log(x[i]);
		mean += val;
		if(val>=0.0) ++count;
	}
	mean = exp(mean/(double)n);
	count /= (double)n;
	/* sort z by y */
	sort(z.begin(),z.end(),sort_by_vector_double(y));
	int j = (int)z[0];
	lo = x[j];
	hi = x[n-m+j];
}
