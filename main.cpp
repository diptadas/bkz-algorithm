#include <bits/stdc++.h>
#include <NTL/vec_xdouble.h>
#include <NTL/vec_double.h>
#include <NTL/mat_RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/new.h>

using namespace std;
using namespace NTL;

#include "FPLLL.h"
#include "MyTools.h"
#include "NTL_LLL_XD.h"
#include "NTL_BKZ_XD.h"
#include "MyEnum.h"
#include "MyBKZ.h"

#define getTime(time_1) double(clock() - time_1) / CLOCKS_PER_SEC


int main()
{
	int n, k, seed, option;

	cout << "1. BKZ Pruned Enumeration\n2. RBKZ\n3. RBKZ DP\n4. BKZ Randomize Input\n";
	cout << "Option: ";	
	cin >> option;

	cout << "Dim(n): ";
	cin >> n;
	cout << "Blocksize(k): ";
	cin >> k;
	cout << "Seed: ";
	cin >> seed;

	mat_ZZ B;
	generateInput(B, n, seed);

	cout << endl << endl;

	if(option == 1) MyBKZ(B, n, k);
	else if(option == 2) MyBKZ_Rec(B, n, k, true);
	else if(option == 3) MyBKZ_Dp(B, n, k, true);
	else if(option == 4) MyBKZ_Rand(B, n, k, n/3, true, true);
	else cout << "No reduction";

	cout << endl << endl;

	cout << "Shortest norm after BKZ: " << getEuclidNorm(B(1), n) << endl;
	cout << "Shortest vector:" << endl << B(1) << endl;
	//cout << "Time: " << getTime(time_1) << endl << endl;

	return 0;
}


