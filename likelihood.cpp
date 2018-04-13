/****************************************/
/*	likelihood.cpp 20th December 2005	*/
/*	Part of omegaMap v0.5				*/
/*	(c) Danny Wilson.					*/
/*	www.danielwilson.me.uk				*/
/****************************************/

#include "omegaMap.h"
#include <math.h>
#include <limits>
extern "C" {
	#include "paml.h"
}

/***************************/
/*  omegaMapBase routines  */
/***************************/

double omegaMapBase::thetaS() {
	double A = 0.0;
	double mu = _block[0].oMat->mu;
	double kappa = _block[0].oMat->kappa;
	A += pi[0]*pi[1]*kappa;
	A += pi[2]*pi[3]*kappa;
	A += pi[2]*pi[15]*kappa;
	A += pi[3]*pi[16]*kappa;
	A += pi[4]*pi[5]*kappa;
	A += pi[4]*pi[6];
	A += pi[4]*pi[7];
	A += pi[5]*pi[6];
	A += pi[5]*pi[7];
	A += pi[6]*pi[7]*kappa;
	A += pi[8]*pi[9]*kappa;
	A += pi[10]*pi[11]*kappa;
	A += pi[13]*pi[14]*kappa;
	A += pi[13]*pi[15];
	A += pi[13]*pi[16];
	A += pi[14]*pi[15];
	A += pi[14]*pi[16];
	A += pi[15]*pi[16]*kappa;
	A += pi[17]*pi[18]*kappa;
	A += pi[17]*pi[19];
	A += pi[17]*pi[20];
	A += pi[18]*pi[19];
	A += pi[18]*pi[20];
	A += pi[19]*pi[20]*kappa;
	A += pi[21]*pi[22]*kappa;
	A += pi[23]*pi[24]*kappa;
	A += pi[25]*pi[26]*kappa;
	A += pi[25]*pi[27];
	A += pi[25]*pi[28];
	A += pi[26]*pi[27];
	A += pi[26]*pi[28];
	A += pi[27]*pi[28]*kappa;
	A += pi[27]*pi[43];
	A += pi[28]*pi[44];
	A += pi[29]*pi[30]*kappa;
	A += pi[29]*pi[31];
	A += pi[30]*pi[31];
	A += pi[33]*pi[34]*kappa;
	A += pi[33]*pi[35];
	A += pi[33]*pi[36];
	A += pi[34]*pi[35];
	A += pi[34]*pi[36];
	A += pi[35]*pi[36]*kappa;
	A += pi[37]*pi[38]*kappa;
	A += pi[39]*pi[40]*kappa;
	A += pi[41]*pi[42]*kappa;
	A += pi[43]*pi[44]*kappa;
	A += pi[45]*pi[46]*kappa;
	A += pi[45]*pi[47];
	A += pi[45]*pi[48];
	A += pi[46]*pi[47];
	A += pi[46]*pi[48];
	A += pi[47]*pi[48]*kappa;
	A += pi[49]*pi[50]*kappa;
	A += pi[49]*pi[51];
	A += pi[49]*pi[52];
	A += pi[50]*pi[51];
	A += pi[50]*pi[52];
	A += pi[51]*pi[52]*kappa;
	A += pi[53]*pi[54]*kappa;
	A += pi[55]*pi[56]*kappa;
	A += pi[57]*pi[58]*kappa;
	A += pi[57]*pi[59];
	A += pi[57]*pi[60];
	A += pi[58]*pi[59];
	A += pi[58]*pi[60];
	A += pi[59]*pi[60]*kappa;
	A *= 4.0*mu;
	return A;
}

/***************************/
/*   Likelihood routines   */
/***************************/
omegaMap& omegaMap::diagonalizeMatrix(oMatrix &oM, const double mu, const double kappa, const double omega) {
#ifdef _TESTPRIOR
	oM.mu = mu;
	oM.kappa = kappa;
	oM.omega = omega;
	return *this;
#endif
	buildMatrixA(oM.R,mu,kappa,omega);
	eigenRealSym(oM.R.array,61,oM.lambda.element,pamlWork.array);
	oM.mu = mu;
	oM.kappa = kappa;
	oM.omega = omega;
	int i,j,k;
	for(i=0;i<61;i++)
		for(j=0;j<=i;j++)
			for(k=0;k<n;k++)
				(*oM.gamma)[i][j][k] = -1.0;
	if(maxCodon==62)
		for(i=0;i<62;i++)
			for(j=0;j<=i;j++)
				for(k=0;k<n;k++)
					(*oM.gamma2)[i][j][k] = -1.0;

	return *this;
}

omegaMap& omegaMap::buildMatrixA(Matrix<double> &A, const double mu, const double kappa, const double omega) {
	int i,j;
	/*Initialize to zero*/
	for(i=0;i<61;i++){for(j=0;j<61;j++)A[i][j]=0.0;}
	
	A[0][1] = sqrtPi[0]*sqrtPi[1]*kappa*mu;
	A[0][2] = sqrtPi[0]*sqrtPi[2]*omega*mu;
	A[0][3] = sqrtPi[0]*sqrtPi[3]*omega*mu;
	A[0][4] = sqrtPi[0]*sqrtPi[4]*kappa*omega*mu;
	A[0][8] = sqrtPi[0]*sqrtPi[8]*omega*mu;
	A[0][10] = sqrtPi[0]*sqrtPi[10]*omega*mu;
	A[0][13] = sqrtPi[0]*sqrtPi[13]*kappa*omega*mu;
	A[0][29] = sqrtPi[0]*sqrtPi[29]*omega*mu;
	A[0][45] = sqrtPi[0]*sqrtPi[45]*omega*mu;
	A[1][2] = sqrtPi[1]*sqrtPi[2]*omega*mu;
	A[1][3] = sqrtPi[1]*sqrtPi[3]*omega*mu;
	A[1][5] = sqrtPi[1]*sqrtPi[5]*kappa*omega*mu;
	A[1][9] = sqrtPi[1]*sqrtPi[9]*omega*mu;
	A[1][11] = sqrtPi[1]*sqrtPi[11]*omega*mu;
	A[1][14] = sqrtPi[1]*sqrtPi[14]*kappa*omega*mu;
	A[1][30] = sqrtPi[1]*sqrtPi[30]*omega*mu;
	A[1][46] = sqrtPi[1]*sqrtPi[46]*omega*mu;
	A[2][3] = sqrtPi[2]*sqrtPi[3]*kappa*mu;
	A[2][6] = sqrtPi[2]*sqrtPi[6]*kappa*omega*mu;
	A[2][15] = sqrtPi[2]*sqrtPi[15]*kappa*mu;
	A[2][31] = sqrtPi[2]*sqrtPi[31]*omega*mu;
	A[2][47] = sqrtPi[2]*sqrtPi[47]*omega*mu;
	A[3][7] = sqrtPi[3]*sqrtPi[7]*kappa*omega*mu;
	A[3][12] = sqrtPi[3]*sqrtPi[12]*omega*mu;
	A[3][16] = sqrtPi[3]*sqrtPi[16]*kappa*mu;
	A[3][32] = sqrtPi[3]*sqrtPi[32]*omega*mu;
	A[3][48] = sqrtPi[3]*sqrtPi[48]*omega*mu;
	A[4][5] = sqrtPi[4]*sqrtPi[5]*kappa*mu;
	A[4][6] = sqrtPi[4]*sqrtPi[6]*mu;
	A[4][7] = sqrtPi[4]*sqrtPi[7]*mu;
	A[4][8] = sqrtPi[4]*sqrtPi[8]*omega*mu;
	A[4][10] = sqrtPi[4]*sqrtPi[10]*omega*mu;
	A[4][17] = sqrtPi[4]*sqrtPi[17]*kappa*omega*mu;
	A[4][33] = sqrtPi[4]*sqrtPi[33]*omega*mu;
	A[4][49] = sqrtPi[4]*sqrtPi[49]*omega*mu;
	A[5][6] = sqrtPi[5]*sqrtPi[6]*mu;
	A[5][7] = sqrtPi[5]*sqrtPi[7]*mu;
	A[5][9] = sqrtPi[5]*sqrtPi[9]*omega*mu;
	A[5][11] = sqrtPi[5]*sqrtPi[11]*omega*mu;
	A[5][18] = sqrtPi[5]*sqrtPi[18]*kappa*omega*mu;
	A[5][34] = sqrtPi[5]*sqrtPi[34]*omega*mu;
	A[5][50] = sqrtPi[5]*sqrtPi[50]*omega*mu;
	A[6][7] = sqrtPi[6]*sqrtPi[7]*kappa*mu;
	A[6][19] = sqrtPi[6]*sqrtPi[19]*kappa*omega*mu;
	A[6][35] = sqrtPi[6]*sqrtPi[35]*omega*mu;
	A[6][51] = sqrtPi[6]*sqrtPi[51]*omega*mu;
	A[7][12] = sqrtPi[7]*sqrtPi[12]*omega*mu;
	A[7][20] = sqrtPi[7]*sqrtPi[20]*kappa*omega*mu;
	A[7][36] = sqrtPi[7]*sqrtPi[36]*omega*mu;
	A[7][52] = sqrtPi[7]*sqrtPi[52]*omega*mu;
	A[8][9] = sqrtPi[8]*sqrtPi[9]*kappa*mu;
	A[8][10] = sqrtPi[8]*sqrtPi[10]*kappa*omega*mu;
	A[8][21] = sqrtPi[8]*sqrtPi[21]*kappa*omega*mu;
	A[8][37] = sqrtPi[8]*sqrtPi[37]*omega*mu;
	A[8][53] = sqrtPi[8]*sqrtPi[53]*omega*mu;
	A[9][11] = sqrtPi[9]*sqrtPi[11]*kappa*omega*mu;
	A[9][22] = sqrtPi[9]*sqrtPi[22]*kappa*omega*mu;
	A[9][38] = sqrtPi[9]*sqrtPi[38]*omega*mu;
	A[9][54] = sqrtPi[9]*sqrtPi[54]*omega*mu;
	A[10][11] = sqrtPi[10]*sqrtPi[11]*kappa*mu;
	A[10][12] = sqrtPi[10]*sqrtPi[12]*omega*mu;
	A[10][25] = sqrtPi[10]*sqrtPi[25]*kappa*omega*mu;
	A[10][41] = sqrtPi[10]*sqrtPi[41]*omega*mu;
	A[10][57] = sqrtPi[10]*sqrtPi[57]*omega*mu;
	A[11][12] = sqrtPi[11]*sqrtPi[12]*omega*mu;
	A[11][26] = sqrtPi[11]*sqrtPi[26]*kappa*omega*mu;
	A[11][42] = sqrtPi[11]*sqrtPi[42]*omega*mu;
	A[11][58] = sqrtPi[11]*sqrtPi[58]*omega*mu;
	A[12][28] = sqrtPi[12]*sqrtPi[28]*kappa*omega*mu;
	A[12][44] = sqrtPi[12]*sqrtPi[44]*omega*mu;
	A[12][60] = sqrtPi[12]*sqrtPi[60]*omega*mu;
	A[13][14] = sqrtPi[13]*sqrtPi[14]*kappa*mu;
	A[13][15] = sqrtPi[13]*sqrtPi[15]*mu;
	A[13][16] = sqrtPi[13]*sqrtPi[16]*mu;
	A[13][17] = sqrtPi[13]*sqrtPi[17]*kappa*omega*mu;
	A[13][21] = sqrtPi[13]*sqrtPi[21]*omega*mu;
	A[13][25] = sqrtPi[13]*sqrtPi[25]*omega*mu;
	A[13][29] = sqrtPi[13]*sqrtPi[29]*omega*mu;
	A[13][45] = sqrtPi[13]*sqrtPi[45]*omega*mu;
	A[14][15] = sqrtPi[14]*sqrtPi[15]*mu;
	A[14][16] = sqrtPi[14]*sqrtPi[16]*mu;
	A[14][18] = sqrtPi[14]*sqrtPi[18]*kappa*omega*mu;
	A[14][22] = sqrtPi[14]*sqrtPi[22]*omega*mu;
	A[14][26] = sqrtPi[14]*sqrtPi[26]*omega*mu;
	A[14][30] = sqrtPi[14]*sqrtPi[30]*omega*mu;
	A[14][46] = sqrtPi[14]*sqrtPi[46]*omega*mu;
	A[15][16] = sqrtPi[15]*sqrtPi[16]*kappa*mu;
	A[15][19] = sqrtPi[15]*sqrtPi[19]*kappa*omega*mu;
	A[15][23] = sqrtPi[15]*sqrtPi[23]*omega*mu;
	A[15][27] = sqrtPi[15]*sqrtPi[27]*omega*mu;
	A[15][31] = sqrtPi[15]*sqrtPi[31]*omega*mu;
	A[15][47] = sqrtPi[15]*sqrtPi[47]*omega*mu;
	A[16][20] = sqrtPi[16]*sqrtPi[20]*kappa*omega*mu;
	A[16][24] = sqrtPi[16]*sqrtPi[24]*omega*mu;
	A[16][28] = sqrtPi[16]*sqrtPi[28]*omega*mu;
	A[16][32] = sqrtPi[16]*sqrtPi[32]*omega*mu;
	A[16][48] = sqrtPi[16]*sqrtPi[48]*omega*mu;
	A[17][18] = sqrtPi[17]*sqrtPi[18]*kappa*mu;
	A[17][19] = sqrtPi[17]*sqrtPi[19]*mu;
	A[17][20] = sqrtPi[17]*sqrtPi[20]*mu;
	A[17][21] = sqrtPi[17]*sqrtPi[21]*omega*mu;
	A[17][25] = sqrtPi[17]*sqrtPi[25]*omega*mu;
	A[17][33] = sqrtPi[17]*sqrtPi[33]*omega*mu;
	A[17][49] = sqrtPi[17]*sqrtPi[49]*omega*mu;
	A[18][19] = sqrtPi[18]*sqrtPi[19]*mu;
	A[18][20] = sqrtPi[18]*sqrtPi[20]*mu;
	A[18][22] = sqrtPi[18]*sqrtPi[22]*omega*mu;
	A[18][26] = sqrtPi[18]*sqrtPi[26]*omega*mu;
	A[18][34] = sqrtPi[18]*sqrtPi[34]*omega*mu;
	A[18][50] = sqrtPi[18]*sqrtPi[50]*omega*mu;
	A[19][20] = sqrtPi[19]*sqrtPi[20]*kappa*mu;
	A[19][23] = sqrtPi[19]*sqrtPi[23]*omega*mu;
	A[19][27] = sqrtPi[19]*sqrtPi[27]*omega*mu;
	A[19][35] = sqrtPi[19]*sqrtPi[35]*omega*mu;
	A[19][51] = sqrtPi[19]*sqrtPi[51]*omega*mu;
	A[20][24] = sqrtPi[20]*sqrtPi[24]*omega*mu;
	A[20][28] = sqrtPi[20]*sqrtPi[28]*omega*mu;
	A[20][36] = sqrtPi[20]*sqrtPi[36]*omega*mu;
	A[20][52] = sqrtPi[20]*sqrtPi[52]*omega*mu;
	A[21][22] = sqrtPi[21]*sqrtPi[22]*kappa*mu;
	A[21][23] = sqrtPi[21]*sqrtPi[23]*omega*mu;
	A[21][24] = sqrtPi[21]*sqrtPi[24]*omega*mu;
	A[21][25] = sqrtPi[21]*sqrtPi[25]*kappa*omega*mu;
	A[21][37] = sqrtPi[21]*sqrtPi[37]*omega*mu;
	A[21][53] = sqrtPi[21]*sqrtPi[53]*omega*mu;
	A[22][23] = sqrtPi[22]*sqrtPi[23]*omega*mu;
	A[22][24] = sqrtPi[22]*sqrtPi[24]*omega*mu;
	A[22][26] = sqrtPi[22]*sqrtPi[26]*kappa*omega*mu;
	A[22][38] = sqrtPi[22]*sqrtPi[38]*omega*mu;
	A[22][54] = sqrtPi[22]*sqrtPi[54]*omega*mu;
	A[23][24] = sqrtPi[23]*sqrtPi[24]*kappa*mu;
	A[23][27] = sqrtPi[23]*sqrtPi[27]*kappa*omega*mu;
	A[23][39] = sqrtPi[23]*sqrtPi[39]*omega*mu;
	A[23][55] = sqrtPi[23]*sqrtPi[55]*omega*mu;
	A[24][28] = sqrtPi[24]*sqrtPi[28]*kappa*omega*mu;
	A[24][40] = sqrtPi[24]*sqrtPi[40]*omega*mu;
	A[24][56] = sqrtPi[24]*sqrtPi[56]*omega*mu;
	A[25][26] = sqrtPi[25]*sqrtPi[26]*kappa*mu;
	A[25][27] = sqrtPi[25]*sqrtPi[27]*mu;
	A[25][28] = sqrtPi[25]*sqrtPi[28]*mu;
	A[25][41] = sqrtPi[25]*sqrtPi[41]*omega*mu;
	A[25][57] = sqrtPi[25]*sqrtPi[57]*omega*mu;
	A[26][27] = sqrtPi[26]*sqrtPi[27]*mu;
	A[26][28] = sqrtPi[26]*sqrtPi[28]*mu;
	A[26][42] = sqrtPi[26]*sqrtPi[42]*omega*mu;
	A[26][58] = sqrtPi[26]*sqrtPi[58]*omega*mu;
	A[27][28] = sqrtPi[27]*sqrtPi[28]*kappa*mu;
	A[27][43] = sqrtPi[27]*sqrtPi[43]*mu;
	A[27][59] = sqrtPi[27]*sqrtPi[59]*omega*mu;
	A[28][44] = sqrtPi[28]*sqrtPi[44]*mu;
	A[28][60] = sqrtPi[28]*sqrtPi[60]*omega*mu;
	A[29][30] = sqrtPi[29]*sqrtPi[30]*kappa*mu;
	A[29][31] = sqrtPi[29]*sqrtPi[31]*mu;
	A[29][32] = sqrtPi[29]*sqrtPi[32]*omega*mu;
	A[29][33] = sqrtPi[29]*sqrtPi[33]*kappa*omega*mu;
	A[29][37] = sqrtPi[29]*sqrtPi[37]*omega*mu;
	A[29][41] = sqrtPi[29]*sqrtPi[41]*omega*mu;
	A[29][45] = sqrtPi[29]*sqrtPi[45]*kappa*omega*mu;
	A[30][31] = sqrtPi[30]*sqrtPi[31]*mu;
	A[30][32] = sqrtPi[30]*sqrtPi[32]*omega*mu;
	A[30][34] = sqrtPi[30]*sqrtPi[34]*kappa*omega*mu;
	A[30][38] = sqrtPi[30]*sqrtPi[38]*omega*mu;
	A[30][42] = sqrtPi[30]*sqrtPi[42]*omega*mu;
	A[30][46] = sqrtPi[30]*sqrtPi[46]*kappa*omega*mu;
	A[31][32] = sqrtPi[31]*sqrtPi[32]*kappa*omega*mu;
	A[31][35] = sqrtPi[31]*sqrtPi[35]*kappa*omega*mu;
	A[31][39] = sqrtPi[31]*sqrtPi[39]*omega*mu;
	A[31][43] = sqrtPi[31]*sqrtPi[43]*omega*mu;
	A[31][47] = sqrtPi[31]*sqrtPi[47]*kappa*omega*mu;
	A[32][36] = sqrtPi[32]*sqrtPi[36]*kappa*omega*mu;
	A[32][40] = sqrtPi[32]*sqrtPi[40]*omega*mu;
	A[32][44] = sqrtPi[32]*sqrtPi[44]*omega*mu;
	A[32][48] = sqrtPi[32]*sqrtPi[48]*kappa*omega*mu;
	A[33][34] = sqrtPi[33]*sqrtPi[34]*kappa*mu;
	A[33][35] = sqrtPi[33]*sqrtPi[35]*mu;
	A[33][36] = sqrtPi[33]*sqrtPi[36]*mu;
	A[33][37] = sqrtPi[33]*sqrtPi[37]*omega*mu;
	A[33][41] = sqrtPi[33]*sqrtPi[41]*omega*mu;
	A[33][49] = sqrtPi[33]*sqrtPi[49]*kappa*omega*mu;
	A[34][35] = sqrtPi[34]*sqrtPi[35]*mu;
	A[34][36] = sqrtPi[34]*sqrtPi[36]*mu;
	A[34][38] = sqrtPi[34]*sqrtPi[38]*omega*mu;
	A[34][42] = sqrtPi[34]*sqrtPi[42]*omega*mu;
	A[34][50] = sqrtPi[34]*sqrtPi[50]*kappa*omega*mu;
	A[35][36] = sqrtPi[35]*sqrtPi[36]*kappa*mu;
	A[35][39] = sqrtPi[35]*sqrtPi[39]*omega*mu;
	A[35][43] = sqrtPi[35]*sqrtPi[43]*omega*mu;
	A[35][51] = sqrtPi[35]*sqrtPi[51]*kappa*omega*mu;
	A[36][40] = sqrtPi[36]*sqrtPi[40]*omega*mu;
	A[36][44] = sqrtPi[36]*sqrtPi[44]*omega*mu;
	A[36][52] = sqrtPi[36]*sqrtPi[52]*kappa*omega*mu;
	A[37][38] = sqrtPi[37]*sqrtPi[38]*kappa*mu;
	A[37][39] = sqrtPi[37]*sqrtPi[39]*omega*mu;
	A[37][40] = sqrtPi[37]*sqrtPi[40]*omega*mu;
	A[37][41] = sqrtPi[37]*sqrtPi[41]*kappa*omega*mu;
	A[37][53] = sqrtPi[37]*sqrtPi[53]*kappa*omega*mu;
	A[38][39] = sqrtPi[38]*sqrtPi[39]*omega*mu;
	A[38][40] = sqrtPi[38]*sqrtPi[40]*omega*mu;
	A[38][42] = sqrtPi[38]*sqrtPi[42]*kappa*omega*mu;
	A[38][54] = sqrtPi[38]*sqrtPi[54]*kappa*omega*mu;
	A[39][40] = sqrtPi[39]*sqrtPi[40]*kappa*mu;
	A[39][43] = sqrtPi[39]*sqrtPi[43]*kappa*omega*mu;
	A[39][55] = sqrtPi[39]*sqrtPi[55]*kappa*omega*mu;
	A[40][44] = sqrtPi[40]*sqrtPi[44]*kappa*omega*mu;
	A[40][56] = sqrtPi[40]*sqrtPi[56]*kappa*omega*mu;
	A[41][42] = sqrtPi[41]*sqrtPi[42]*kappa*mu;
	A[41][43] = sqrtPi[41]*sqrtPi[43]*omega*mu;
	A[41][44] = sqrtPi[41]*sqrtPi[44]*omega*mu;
	A[41][57] = sqrtPi[41]*sqrtPi[57]*kappa*omega*mu;
	A[42][43] = sqrtPi[42]*sqrtPi[43]*omega*mu;
	A[42][44] = sqrtPi[42]*sqrtPi[44]*omega*mu;
	A[42][58] = sqrtPi[42]*sqrtPi[58]*kappa*omega*mu;
	A[43][44] = sqrtPi[43]*sqrtPi[44]*kappa*mu;
	A[43][59] = sqrtPi[43]*sqrtPi[59]*kappa*omega*mu;
	A[44][60] = sqrtPi[44]*sqrtPi[60]*kappa*omega*mu;
	A[45][46] = sqrtPi[45]*sqrtPi[46]*kappa*mu;
	A[45][47] = sqrtPi[45]*sqrtPi[47]*mu;
	A[45][48] = sqrtPi[45]*sqrtPi[48]*mu;
	A[45][49] = sqrtPi[45]*sqrtPi[49]*kappa*omega*mu;
	A[45][53] = sqrtPi[45]*sqrtPi[53]*omega*mu;
	A[45][57] = sqrtPi[45]*sqrtPi[57]*omega*mu;
	A[46][47] = sqrtPi[46]*sqrtPi[47]*mu;
	A[46][48] = sqrtPi[46]*sqrtPi[48]*mu;
	A[46][50] = sqrtPi[46]*sqrtPi[50]*kappa*omega*mu;
	A[46][54] = sqrtPi[46]*sqrtPi[54]*omega*mu;
	A[46][58] = sqrtPi[46]*sqrtPi[58]*omega*mu;
	A[47][48] = sqrtPi[47]*sqrtPi[48]*kappa*mu;
	A[47][51] = sqrtPi[47]*sqrtPi[51]*kappa*omega*mu;
	A[47][55] = sqrtPi[47]*sqrtPi[55]*omega*mu;
	A[47][59] = sqrtPi[47]*sqrtPi[59]*omega*mu;
	A[48][52] = sqrtPi[48]*sqrtPi[52]*kappa*omega*mu;
	A[48][56] = sqrtPi[48]*sqrtPi[56]*omega*mu;
	A[48][60] = sqrtPi[48]*sqrtPi[60]*omega*mu;
	A[49][50] = sqrtPi[49]*sqrtPi[50]*kappa*mu;
	A[49][51] = sqrtPi[49]*sqrtPi[51]*mu;
	A[49][52] = sqrtPi[49]*sqrtPi[52]*mu;
	A[49][53] = sqrtPi[49]*sqrtPi[53]*omega*mu;
	A[49][57] = sqrtPi[49]*sqrtPi[57]*omega*mu;
	A[50][51] = sqrtPi[50]*sqrtPi[51]*mu;
	A[50][52] = sqrtPi[50]*sqrtPi[52]*mu;
	A[50][54] = sqrtPi[50]*sqrtPi[54]*omega*mu;
	A[50][58] = sqrtPi[50]*sqrtPi[58]*omega*mu;
	A[51][52] = sqrtPi[51]*sqrtPi[52]*kappa*mu;
	A[51][55] = sqrtPi[51]*sqrtPi[55]*omega*mu;
	A[51][59] = sqrtPi[51]*sqrtPi[59]*omega*mu;
	A[52][56] = sqrtPi[52]*sqrtPi[56]*omega*mu;
	A[52][60] = sqrtPi[52]*sqrtPi[60]*omega*mu;
	A[53][54] = sqrtPi[53]*sqrtPi[54]*kappa*mu;
	A[53][55] = sqrtPi[53]*sqrtPi[55]*omega*mu;
	A[53][56] = sqrtPi[53]*sqrtPi[56]*omega*mu;
	A[53][57] = sqrtPi[53]*sqrtPi[57]*kappa*omega*mu;
	A[54][55] = sqrtPi[54]*sqrtPi[55]*omega*mu;
	A[54][56] = sqrtPi[54]*sqrtPi[56]*omega*mu;
	A[54][58] = sqrtPi[54]*sqrtPi[58]*kappa*omega*mu;
	A[55][56] = sqrtPi[55]*sqrtPi[56]*kappa*mu;
	A[55][59] = sqrtPi[55]*sqrtPi[59]*kappa*omega*mu;
	A[56][60] = sqrtPi[56]*sqrtPi[60]*kappa*omega*mu;
	A[57][58] = sqrtPi[57]*sqrtPi[58]*kappa*mu;
	A[57][59] = sqrtPi[57]*sqrtPi[59]*mu;
	A[57][60] = sqrtPi[57]*sqrtPi[60]*mu;
	A[58][59] = sqrtPi[58]*sqrtPi[59]*mu;
	A[58][60] = sqrtPi[58]*sqrtPi[60]*mu;
	A[59][60] = sqrtPi[59]*sqrtPi[60]*kappa*mu;

	/*Fill in the lower triangle*/
	for(i=0;i<61;i++){
		for(j=i+1;j<61;j++)A[j][i]=A[i][j];}

	/*Compute the diagonal*/
	for(i=0;i<61;i++)
	{
		double rowsum=0.0;
		for(j=0;j<61;j++) rowsum+=A[i][j];
		A[i][i]=-rowsum;
	}

	return *this;
}

/* it is assumed that i<=j */
double omegaMap::mutProb(int site, int nhaps, int i, int j) {
#ifdef _DEBUG
	if(j>i) warning("omegaMap::mutProb(): j must be less than or equal to i");
	if(j<0) warning("omegaMap::mutProb(): j must be greater than or equal to 0");
	if(i>=maxCodon) warning("omegaMap::mutProb(): i must be less than maxCodon");
#endif
	oMatrix &oM = *(block[site]->oMat);
	bool changed = false;
	if(!indel[site]) {
		double &pr = (*oM.gamma)[i][j][nhaps];
		if(pr==-1.0) {
			changed = true;
			pr = 0.0;
			int k;
			for(k=0;k<61;k++) pr += oM.R[i][k] * oM.R[j][k] * sqrtPi[j] / sqrtPi[i] * (double)nhaps / ((double)nhaps - 2.0 * oM.lambda[k]);
			pr *= (i==j) ? pi[i] : 2.0*pi[i];
			if(pr<=0.0) pr = numeric_limits<double>::epsilon();
		}
		return pr;
	}
	else {
		double &pr = (*oM.gamma2)[i][j][nhaps];
		if(pr==-1.0) {
			changed = true;
			pr = 0.0;
			int k;
			if(i<61) { /* both are codons */
				pr = (1.-pi[61]) + (double)nhaps*pi[61]/((double)nhaps+2.*oM.omega*indelLambda) 
					- (double)nhaps/((double)nhaps+2.*pi[61]*oM.omega*indelLambda);
				pr *= pi[j];
				for(k=0;k<61;k++) pr += oM.R[i][k] * oM.R[j][k] * sqrtPi[j] / sqrtPi[i] 
					* (double)nhaps / ((double)nhaps + 2.*pi[61]*oM.omega*indelLambda - 2.0 * oM.lambda[k]);
				pr *= (i==j) ? pi[i]*(1.-pi[61]) : 2.*pi[i]*(1.-pi[61]);
			}
			else if(j<61) { /* one indel, one codon */
				pr = 2.*pi[i]*pi[j]*(1.-pi[61])*(1. - (double)nhaps/((double)nhaps+2.*oM.omega*indelLambda));
			}
			else { /* both are indels */
				pr = pi[61]*(pi[61] + (double)nhaps*(1.-pi[61])/((double)nhaps+2.*oM.omega*indelLambda));
			}
			if(pr<=0.0) pr = numeric_limits<double>::epsilon();
		}
		return pr;
	}
}

/* updates alpha from _5prime to _3prime inclusive */
omegaMap& omegaMap::forward(int _5prime, int _3prime) {
#ifdef _TESTPRIOR
	return *this;
#endif
#ifdef _STARLIKELIHOOD
	return *this;
#endif
	int ord,k,pos,x,i,j;
	double lnk,diff,P,lnP,OneMinusP,lnOneMinusP;

	int _5primeMin = MAX(0,_5prime);
	int _3primeMax = MIN(L-1,_3prime);
	for(ord=0;ord<norders;ord++) {
		for(k=1;k<n;k++) {
			lnk = log((double)k);
			for(pos=_5primeMin;pos<=_3primeMax;pos++) {
				if(pos>0) {
					lnP = -rblock[pos-1]->rho/(double)k;
					if(lnP>=0.0) error("omegaMap::forward(): P=1.0");
					P = exp(lnP);
					OneMinusP = 1.0 - P;
					lnOneMinusP = log(OneMinusP);
				}
				alpha[ord][k][pos][k] = 0.0; /* this is the sum of elements [0]..[k-1] */
				if(pos==0) {
					for(x=0;x<k;x++) {
						if((*(H[ord][x]))[pos]<=(*(H[ord][k]))[pos]) {
							i = (*(H[ord][k]))[pos];
							j = (*(H[ord][x]))[pos];
						}
						else {
							i = (*(H[ord][x]))[pos];
							j = (*(H[ord][k]))[pos];
						}
						alpha[ord][k][pos][x] = log(mutProb(pos,k,i,j)) - lnk;
						diff = alpha[ord][k][pos][x] - alpha[ord][k][pos][k];
						if(x==0) alpha[ord][k][pos][k] = alpha[ord][k][pos][x];
						else if(diff<0.0) alpha[ord][k][pos][k] += log(1.0 + exp(diff));
						else alpha[ord][k][pos][k] = alpha[ord][k][pos][x] + log(1.0 + exp(-diff));
					}
				}
				else {
					for(x=0;x<k;x++) {
						if((*(H[ord][x]))[pos]<=(*(H[ord][k]))[pos]) {
							i = (*(H[ord][k]))[pos];
							j = (*(H[ord][x]))[pos];
						}
						else {
							i = (*(H[ord][x]))[pos];
							j = (*(H[ord][k]))[pos];
						}
						if(P>=0.5) {
								alpha[ord][k][pos][x] = log(mutProb(pos,k,i,j)) + lnP + alpha[ord][k][pos-1][x]
									+ log(1.0 + exp(lnOneMinusP-lnP+alpha[ord][k][pos-1][k]-alpha[ord][k][pos-1][x]-lnk));
								diff = alpha[ord][k][pos][x]-alpha[ord][k][pos][k];
								if(x==0) alpha[ord][k][pos][k] = alpha[ord][k][pos][x];
								else if(diff<0.0) alpha[ord][k][pos][k] += log(1.0 + exp(diff));
								else alpha[ord][k][pos][k] = alpha[ord][k][pos][x] + log(1.0 + exp(-diff));
						}
						else {
								alpha[ord][k][pos][x] = log(mutProb(pos,k,i,j)) + lnOneMinusP - lnk + alpha[ord][k][pos-1][k]
									+ log(1.0 + exp(lnP-lnOneMinusP+alpha[ord][k][pos-1][x]+lnk-alpha[ord][k][pos-1][k]));
								diff = alpha[ord][k][pos][x]-alpha[ord][k][pos][k];
								if(x==0) alpha[ord][k][pos][k] = alpha[ord][k][pos][x];
								else if(diff<0.0) alpha[ord][k][pos][k] += log(1.0 + exp(diff));
								else alpha[ord][k][pos][k] = alpha[ord][k][pos][x] + log(1.0 + exp(-diff));
						}
					}
				}
			}
		}
	}

	return *this;
}

/* updates beta from _3prime to _5prime inclusive */
omegaMap& omegaMap::backward(int _5prime, int _3prime) {
#ifdef _TESTPRIOR
	return *this;
#endif
#ifdef _STARLIKELIHOOD
	return *this;
#endif
	int ord,k,pos,x,i,j;
	double lnk,diff,P,lnP,OneMinusP,lnOneMinusP;

	int _5primeMin = MAX(0,_5prime);
	int _3primeMax = MIN(L-1,_3prime);
	for(ord=0;ord<norders;ord++) {
		for(k=1;k<n;k++) {
			lnk = log((double)k);
			for(pos=_3primeMax;pos>=_5primeMin;pos--) {
				if(pos<L-1) {
					lnP = -rblock[pos]->rho/(double)k;
					if(lnP>=0.0) error("omegaMap::backward(): P=1.0");
					P = exp(lnP);
					OneMinusP = 1.0 - P;
					lnOneMinusP = log(OneMinusP);
				}
				beta[ord][k][pos][k] = 0.0; /* this is the sum of elements [0]..[k-1] */
				if(pos==L-1) {
					for(x=0;x<k;x++) {
						if((*(H[ord][x]))[pos]<=(*(H[ord][k]))[pos]) {
							i = (*(H[ord][k]))[pos];
							j = (*(H[ord][x]))[pos];
						}
						else {
							i = (*(H[ord][x]))[pos];
							j = (*(H[ord][k]))[pos];
						}
						beta[ord][k][pos][x] = log(mutProb(pos,k,i,j)) - lnk;
						diff = beta[ord][k][pos][x] - beta[ord][k][pos][k];
						if(x==0) beta[ord][k][pos][k] = beta[ord][k][pos][x];
						else if(diff<0.0) beta[ord][k][pos][k] += log(1.0 + exp(diff));
						else beta[ord][k][pos][k] = beta[ord][k][pos][x] + log(1.0 + exp(-diff));
					}
				}
				else {
					for(x=0;x<k;x++) {
						if((*(H[ord][x]))[pos]<=(*(H[ord][k]))[pos]) {
							i = (*(H[ord][k]))[pos];
							j = (*(H[ord][x]))[pos];
						}
						else {
							i = (*(H[ord][x]))[pos];
							j = (*(H[ord][k]))[pos];
						}
						if(P>=0.5) {
								beta[ord][k][pos][x] = log(mutProb(pos,k,i,j)) + lnP + beta[ord][k][pos+1][x]
									+ log(1.0 + exp(lnOneMinusP-lnP+beta[ord][k][pos+1][k]-beta[ord][k][pos+1][x]-lnk));
								diff = beta[ord][k][pos][x]-beta[ord][k][pos][k];
								if(x==0) beta[ord][k][pos][k] = beta[ord][k][pos][x];
								else if(diff<0.0) beta[ord][k][pos][k] += log(1.0 + exp(diff));
								else beta[ord][k][pos][k] = beta[ord][k][pos][x] + log(1.0 + exp(-diff));
						}
						else {
								beta[ord][k][pos][x] = log(mutProb(pos,k,i,j)) + lnOneMinusP - lnk + beta[ord][k][pos+1][k]
									+ log(1.0 + exp(lnP-lnOneMinusP+beta[ord][k][pos+1][x]+lnk-beta[ord][k][pos+1][k]));
								diff = beta[ord][k][pos][x]-beta[ord][k][pos][k];
								if(x==0) beta[ord][k][pos][k] = beta[ord][k][pos][x];
								else if(diff<0.0) beta[ord][k][pos][k] += log(1.0 + exp(diff));
								else beta[ord][k][pos][k] = beta[ord][k][pos][x] + log(1.0 + exp(-diff));
						}
					}
				}
			}
		}
	}

	return *this;
}

double omegaMap::likelihood(int pos) {
#ifdef _TESTPRIOR
	return 1.0;
#endif
#ifdef _STARLIKELIHOOD
	return star_likelihood();
#endif
	double Lik = 0.0;
	int ord,k,x,i,j;
	double sum,summand,diff,lnk;

	for(ord=0;ord<norders;ord++) {
		PAC[ord] = (double)L * log(0.5);
		for(k=1;k<n;k++) {
			if(pos<=0) sum = beta[ord][k][0][k];
			else if(pos>=L-1) sum = alpha[ord][k][L-1][k];
			else {
				lnk = log((double)k);
				sum = 0.0;
				for(x=0;x<k;x++) {
					if((*(H[ord][x]))[pos]<=(*(H[ord][k]))[pos]) {
						i = (*(H[ord][k]))[pos];
						j = (*(H[ord][x]))[pos];
					}
					else {
						i = (*(H[ord][x]))[pos];
						j = (*(H[ord][k]))[pos];
					}
					summand = alpha[ord][k][pos][x] + beta[ord][k][pos][x]
						- log(mutProb(pos,k,i,j)) + lnk;
					/* the lnk comes in so that the 1/k isn't counted twice */
					diff = summand - sum;
					if(x==0) sum = summand;
					else if(diff<0.0) sum += log(1.0 + exp(diff));
					else sum = summand + log(1.0 + exp(-diff));
				}
			}
			PAC[ord] += sum;
		}
		diff = PAC[ord] - Lik;
		if(ord==0) Lik = PAC[ord];
		else if(diff<0.0) Lik += log(1.0 + exp(diff));
		else {
			Lik = log(1.0 + exp(-diff));
			Lik += PAC[ord];
		}
	}
	Lik -= log((double)norders);
	return Lik;
}

double omegaMap::star_likelihood() {
	double Lik = 0.0;
	int pos,hap,anc,i,j,k;
	double margLik,diff;
	for(pos=0;pos<L;pos++) {
		oMatrix &oM = *(block[pos]->oMat);
		if(!indel[pos]) {
			for(anc=0;anc<61;anc++) {
				margLik = log(pi[anc]);
				for(hap=0;hap<n;hap++) {
					i = anc;
					j = (*H[0][hap])[pos];
					if(j>i) SWAP(i,j);
					double &pr = (*oM.gamma)[i][j][0];
					if(pr==-1.0) {
						pr = 0.0;
						for(k=0;k<61;k++) pr += oM.R[i][k] * oM.R[j][k] * sqrtPi[j] / sqrtPi[i] * exp(oM.lambda[k]*2.0);
						if(pr<=0.0) pr = numeric_limits<double>::epsilon();
					}
					margLik += log(pr);
				}
				diff = margLik = Lik;
				if(anc==0) Lik = margLik;
				else if(diff<0.0) Lik += log(1.0 + exp(diff));
				else Lik = margLik + log(1.0 + exp(-diff));
			}
		}
		else {
			for(anc=0;anc<61;anc++) {
				margLik = log(pi[anc]);
				for(hap=0;hap<n;hap++) {
					i = anc;
					j = (*H[0][hap])[pos];
					if(j>i) SWAP(i,j);
					double &pr = (*oM.gamma2)[i][j][0];
					if(pr==-1.0) {
						pr = 0.0;
						if(i<61) { /* both are codons */
							pr = (1.-pi[61]) + pi[61]*exp(-2.*oM.omega*indelLambda) 
								- exp(-pi[61]*2.*oM.omega*indelLambda);
							pr *= pi[j];
							for(k=0;k<61;k++) pr += oM.R[i][k] * oM.R[j][k] * sqrtPi[j] / sqrtPi[i] 
								* exp(2.*(oM.lambda[k] - pi[61]*oM.omega*indelLambda));
						}
						else if(j<61) { /* i indel, j codon */
							pr = pi[j]*(1.-pi[61])*(1. - exp(-2.*oM.omega*indelLambda));
						}
						else { /* both are indels */
							pr = pi[61] + (1.-pi[61])*exp(-2.*oM.omega*indelLambda);
						}
						if(pr<=0.0) pr = numeric_limits<double>::epsilon();
					}
					margLik += log(pr);
				}
				diff = margLik = Lik;
				if(anc==0) Lik = margLik;
				else if(diff<0.0) Lik += log(1.0 + exp(diff));
				else Lik = margLik + log(1.0 + exp(-diff));
			}
		}
	}
	return Lik;
}