#include <myerror.h>
#include <random.h>
#include <vector.h>
#include <utils.h>
#include <iostream>

using namespace myutils;
using namespace std;

int main(const int argc, const char* argv[]) {
	if(argc!=3) error("SYNTAX: sample-size number-of-orderings");
	const int n = atoi(argv[1]);
	const int norder = atoi(argv[2]);

	Random ran;
	Vector<int> pool(n);
	int o,i,rnum;
	for(o=0;o<norder;o++) {
		for(i=0;i<n;i++) pool[i] = i;
		for(i=n-1;i>=0;i--) {
            rnum = ran.discrete(0,i);
			cout << pool[rnum] << ", ";
			SWAP(pool[i],pool[rnum]);
		}
	}
	cout << "\b\b " << endl;
	return 0;
}