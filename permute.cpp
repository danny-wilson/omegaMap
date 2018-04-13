/*	permute.exe	*/

/*	permute reads in the FASTA file specified in the command line
	and performs a permutation test for recombination. The null
	hypothesis of exchangeability between sites only holds for zero
	recombination or complete linkage equilibrium.

	There are 3 test statistics. Each measures the correlation co-
	efficient between physical distance and a measure of linkage
	disequilibrium:

	Rsq		R-squared
	Dpr		D-prime
	G4		The four-gamete test compatibility

	A p value of less than 0.05 indicates that the chance of observing
	a correlation coefficient as or more extreme under the null
	hypothesis is less than 5%, and so the null hypothesis can be
	rejected.

	*/

#include <myutils.h>
#include <fstream>
#include <sstream>

using namespace myutils;
using namespace std;

class PermutationTest {
public:
	int niter;
	double v[3];
	double p[3];
	DNA *dna;
	DNA sdna;
	bool coutput;
public:
	PermutationTest() {
		coutput = false;
		niter = 999;
		v[0] = v[1] = v[2] = p[0] = p[1] = p[2] = 0.0;
		dna = 0;
	}
	PermutationTest& permute(const char* filename, const int niter, Random &ran);
	PermutationTest& permute(DNA *dna_in, const int niter, Random &ran);
	PermutationTest& printResults();
};

int main(const int argc, const char* argv[]) {
  if(argc<2 || argc>3) error("SYNTAX: file [no.permutations]");
  int niter = (argc==3) ? atoi(argv[2]) : 999;

  Random ran;
  PermutationTest P;
  //P.coutput = true;

  int ctr,i,j;

  P.permute(argv[1],niter,ran);
  P.printResults();
  //myutils::pause();
  return 0;

  /*ofstream out("__out.txt");
  char tab = '\t';
  out << "rsq" << tab << "Dpr" << tab << "G4" << endl;
  for(ctr=0,i=0;i<100;i++)
	  for(j=0;j<1;j++,ctr++) {
		  stringstream fname;
		  //fname << argv[1] << "SIM" << i << "." << j << ".txt";
		  fname << argv[1] << "SIM" << i << ".txt";
		  P.permute(fname.str().c_str(),niter,ran);
		  out << P.p[0] << tab << P.p[1] << tab << P.p[2] << endl;
		  cout << "\rDone " << ctr << " of 1200 runs" << flush;
	  }
  cout << endl;
  out.close();

  myutils::pause();
  return 0;*/
}

PermutationTest& PermutationTest::printResults() {
  cout << "Read " << dna->nseq << " seqs " << dna->lseq << " bp long with " << sdna.lseq << " segregating sites" << endl;
  cout << endl;
  cout << "Tests for recombination. H0 is no/free recombination." << endl;
  cout << "Correlation between distance and r-squared p = " << p[0] << " (" << v[0] << ")" << endl;
  cout << "Correlation between distance and |D'|      p = " << p[1] << " (" << v[1] << ")" << endl;
  cout << "Meunier and Eyre-Walker's (2001) G4        p = " << p[2] << " (" << -v[2] << ")" << endl;
  cout << endl;
  cout << "   When the correlation (indicated in brakets) is negative, linkage" << endl;
  cout << "   disequilibrium decreases with physical distance, a pattern consistent" << endl;
  cout << "   with intragenic recombination." << endl;
  cout << endl;
  if(p[0]<=0.05 || p[1]<=0.05 || p[2]<=0.05) {
	  if((p[0]<=0.05 && v[0]>0.0) ||
		 (p[1]<=0.05 && v[1]>0.0) ||
		 (p[2]<=0.05 && -v[2]>0.0)) {
			 cout << "!! A significant positive correlation was detected between LD and distance." << endl;
			 cout << "!! Such a pattern cannot be explained by recombination, and is biologically" << endl;
			 cout << "!! unlikely. This may indicate a problem in the permutation test." << endl;
		 }
	  else {
		  cout << "** Detectable levels of recombination were found in your sequences," << endl;
		  cout << "** indicated by a significant negative correlation between one or more" << endl;
		  cout << "** LD statistics and physical distance." << endl;
	  }
  }
  return *this;
}


PermutationTest& PermutationTest::permute(const char* filename, const int niter, Random &ran) {
  if(dna==0) dna = new DNA;
  dna->readFASTA(filename);
  return permute(dna,niter,ran);
}

PermutationTest& PermutationTest::permute(DNA *dna_in, const int niter, Random &ran) {
  dna = dna_in;
  int i,j,k,pos;
  
  /* Find out which sites are bi-allelic segregating */
  /* (segregating for short) */
  char state1,state2;
  vector<bool> segregating(dna->lseq);
  for(i=0;i<dna->lseq;i++) {
    segregating[i] = false;
    state1 = (*dna)[0][i];
    for(j=1;j<dna->nseq;j++)
      if((*dna)[j][i]!=state1) {
	state2 = (*dna)[j][i];
	segregating[i] = true;
	j++;
	for(;j<dna->nseq;j++)
	  if((*dna)[j][i]!=state1 && (*dna)[j][i]!=state2) {
	    /* Tri-allelic segregating */
	    segregating[i] = false;
	    break;
	  }
	break;
      }
  }

  /* Count number of segregating sites */
  int nseg = 0;
  for(i=0;i<dna->lseq;i++)
    if(segregating[i]) nseg++;
  if(nseg==0) error("No segregating sites found");

  /* Copy segregating sites to a new DNA object */
  sdna.resize(dna->nseq,nseg);
  vector<double> d(nseg,0.0);
  for(i=0,pos=0;i<dna->lseq;i++)
    if(segregating[i]) {
      for(j=0;j<dna->nseq;j++)
	sdna[j][pos] = (*dna)[j][i];
	  d[pos] = (double)i;
      ++pos;
    }

  /* Calculate frequency statistics */
  vector<double> F(sdna.lseq,1.0);								/* F is the marginal frequency of the dna[0][i] (A) allele at site i */
  vector<double> four(4,0.0);									/* G[i][j] is the frequency of AB (G[i][j][0]), Ab (1), aB (2), ab (3) */
  LowerTriangularMatrix< vector<double> > G(sdna.lseq,four);
  for(i=0;i<sdna.lseq;i++) {
    state1 = sdna[0][i];
    for(j=1;j<sdna.nseq;j++)
      if(sdna[j][i]==state1) F[i]++;
    F[i] /= (double)sdna.nseq;
  }
  for(i=0;i<sdna.lseq;i++)
	  for(j=0;j<i;j++) {
		  state1 = sdna[0][i];
		  state2 = sdna[0][j];
		  for(k=0;k<sdna.nseq;k++) {
			  if(sdna[k][i]==state1 && sdna[k][j]==state2) ++G[i][j][0];
			  else if(sdna[k][i]==state1 && sdna[k][j]!=state2) ++G[i][j][1];
			  else if(sdna[k][i]!=state1 && sdna[k][j]==state2) ++G[i][j][2];
			  else if(sdna[k][i]!=state1 && sdna[k][j]!=state2) ++G[i][j][3];
			  else warning("Unexpected choice");
		  }
		  for(k=0;k<4;k++) G[i][j][k] /= (double)sdna.nseq;
	  }
  
  /* Calculate LD statistics for pairs of sites */
  LowerTriangularMatrix<double> A(nseg,0.0);
  LowerTriangularMatrix<double> B(nseg,0.0);
  LowerTriangularMatrix<double> C(nseg,0.0);
  Matrix<double> D(nseg,nseg,0.0);
  double temp;
//  ofstream out("__out.txt");
//  char tab = '\t';
//  out << "locusA" << tab << "locusB" << tab << "rsq" << tab << "Dprime" << tab << "G4" << "dist" << endl;
  for(i=0;i<nseg;i++) {
	  for(j=0;j<i;j++) {
		  temp = G[i][j][0] - F[i]*F[j];
		  A[i][j] = pow(temp,2.0)/(F[i]*(1.-F[i])*F[j]*(1.-F[j]));
		  B[i][j] = (temp < 0.0) ? -temp/MIN(F[i]*F[j],(1.-F[i])*(1.-F[j])) : temp/MIN(F[i]*(1.-F[j]),(1.-F[i])*F[j]);
		  C[i][j] = (G[i][j][0]>0.0 && G[i][j][1]>0.0 && G[i][j][2]>0.0 && G[i][j][3]>0.0) ? 1.0 : 0.0;
		  D[i][j] = D[j][i] = d[i] - d[j];
//		  out << i << tab << j << tab << A[i][j] << tab << B[i][j] << tab << C[i][j] << tab << D[i][j] << endl;
	  }
  }
//  out.close();

  /* Calculate remaining statistics for correlation coefficients */
  double Abar,Bbar,Cbar,Dbar;
  double Adev,Bdev,Cdev,Ddev;
  Abar=Bbar=Cbar=Dbar=Adev=Bdev=Cdev=Ddev=0.0;
  double ctr = 0.0;
  for(i=0;i<nseg;i++)
	for(j=0;j<i;j++,ctr++) {
	  Abar += A[i][j]; Bbar += B[i][j]; Cbar += C[i][j]; Dbar += D[i][j];
	}
  Abar /= ctr; Bbar /= ctr; Cbar /= ctr; Dbar /= ctr;
  for(i=0;i<nseg;i++)
	  for(j=0;j<i;j++) {
		  Adev += pow(A[i][j]-Abar,2.0);
		  Bdev += pow(B[i][j]-Bbar,2.0);
		  Cdev += pow(C[i][j]-Cbar,2.0);
		  Ddev += pow(D[i][j]-Dbar,2.0);
	  }


  /* Perform the permutation tests */
  double a,b,c,pa,pb,pc;
  double PA,PB,PC,dist;

	a = b = c = 0.0;
	for(j=0;j<nseg;j++)
		for(k=0;k<j;k++) {
			a += A[j][k] * D[j][k];
			b += B[j][k] * D[j][k];
			c += C[j][k] * D[j][k];
	}
	double Aadd = ctr*Abar*Dbar; double Adiv = sqrt(Adev*Ddev);
	double Badd = ctr*Bbar*Dbar; double Bdiv = sqrt(Bdev*Ddev);
	double Cadd = ctr*Cbar*Dbar; double Cdiv = sqrt(Cdev*Ddev);
	a -= Aadd; a /= Adiv;
	b -= Badd; b /= Bdiv;
	c -= Cadd; c /= Cdiv;

  PA = PB = PC = 0.0;
  vector<int> order(sdna.lseq);
  vector<int> pool(sdna.lseq);
  for(i=0;i<niter;i++) {
    for(j=0;j<nseg;j++) pool[j] = j;
    for(j=nseg-1;j>=0;j--) {
      pos = ran.discrete(0,j);
      order[j] = pool[pos];
      pool[pos] = pool[j];
    }
	pa = pb = pc = 0.0;
	for(j=0;j<nseg;j++)
		for(k=0;k<j;k++) {
			dist = D[order[j]][order[k]];
			pa += A[j][k] * dist;
			pb += B[j][k] * dist;
			pc += C[j][k] * dist;
		}
	pa -= Aadd; pa /= Adiv; pb -= Badd; pb /= Bdiv;
	pc -= Cadd; pc /= Cdiv;
	if(pa<=a) ++PA;
	if(pb<=b) ++PB;
	if(pc>=c) ++PC;
  }

  p[0] = (PA+1.)/(double)(niter+1);	v[0] = a;
  p[1] = (PB+1.)/(double)(niter+1);	v[1] = b;
  p[2] = (PC+1.)/(double)(niter+1);	v[2] = c;
  if(coutput) printResults();
  return *this;
}
