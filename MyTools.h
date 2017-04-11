RR getEuclidNorm(vec_ZZ v, int n)
{
	ZZ euclid_norm = to_ZZ(0);
	for (int i = 1; i <= n; i++)
		euclid_norm += v(i) * v(i);
	return sqrt(to_RR(euclid_norm));
}

void generate_random_HNF(vec_ZZ& out, long n, long bit, ZZ seed)
{
	SetSeed(seed);
	out.SetLength(n);
	clear(out);
	ZZ p;
	GenPrime(p, bit * n);
	out(1) = p;
	for (int i = 2; i <= n; i++)
	{
		RandomBnd(out(i), p);
	}
}

//*******generate SVP chalange input********//
void generateInput(mat_ZZ& B, int n, int seed)
{
	B.SetDims(n, n);
	clear(B);

	vec_ZZ v;
	generate_random_HNF(v, n, 10, to_ZZ(seed));

	B(1, 1) = v(1);
	for (int i = 2; i <= n; i++)
	{
		B(i, 1) = v(i);
		B(i, i) = 1;
	}
}

//********randomize basis**********//
void randomizeBasis(mat_ZZ& BR, mat_ZZ& B, mat_ZZ& U_rand, long d, long seed)
{
	srand(seed);

	mat_ZZ U1;
	mat_ZZ U2;
	mat_ZZ U_1;
	mat_ZZ U_2;

	U1.SetDims(d, d);
	clear(U1);
	U2.SetDims(d, d);
	clear(U2);
	U_1.SetDims(d, d);
	clear(U_1);
	U_2.SetDims(d, d);
	clear(U_2);

	//generating unimodular matrix//
	ident(U_1, d);
	ident(U_2, d);

	for (int i = 1; i <= d; i++)
	{
		for (int j = 1; j <= d; j++)
		{
			U1(i, j) = (rand() % 3) - 1;
			U1(j, j) = 1;
			U2(i, j) = (rand() % 3) - 1;
			U2(j, j) = 1;
		}

		U_1 = U_1 * U1;
		U_2 = U_2 * U2;

		clear(U1);
		clear(U2);
	}

	U_rand = U_1 * U_2;

	//randomize the input basis
	BR = U_rand * B;
}


//********generating block GSO matrix***********//
void genBlockGSO(mat_RR& B_k, vec_RR c, mat_RR mu, long d, long jj, long kk)
{
	mat_RR ccc;
	ccc.SetDims(d, d);
	clear(ccc);

	mat_RR muu;
	muu.SetDims(d, d);
	clear(muu);

	for (int i = 1, j = jj; i <= d, j <= kk; i++, j++)
		conv(ccc(i, i), sqrt(c(j)));

	for (int i = 1, i_i = jj; i <= d, i_i <= kk; i++, i_i++)
	{
		for (int j = 1, j_j = jj; j <= i, j_j <= i_i; j++, j_j++)
			conv(muu(i, j), mu(i_i, j_j));

		muu(i, i) = 1;
	}

	B_k = ccc * muu; //block GSO matrix
}


//**********convert mat_RR to mat_ZZ for LLL*********//
void conv_RR_to_ZZ(mat_ZZ& BB, mat_RR& B_k, long d)
{
	//******K = max(C*d^2/|B_k^*_d|, 1) as scaling factor for conversion***********//

	RR norm = to_RR(0);

	for (int i = 1; i <= d; i++) // norm squared of B_k^*_d
		norm = norm + B_k(d, i) * B_k(d, i);

	RR b_d = sqrt(norm);

	long C = 50;

	ZZ K;
	conv(K, (C * to_RR(d)*to_RR(d)) / b_d);

	if (K < 1) K = 1; // K = max(K,1);

	for (int i = 1; i <= d; i++)
		for (int j = 1; j <= d; j++)
			conv(BB(i, j), to_RR(K)*B_k(i, j)); //convert RR type to ZZ type
}


//******number of tours calculation********//
xdouble tour(long n, long beta)
{
	xdouble tr, tr1, tr2;
	xdouble T;

	conv(tr1, (double(n) / double(beta)));
	tr2 = log((n * n) / 50); //log2 or ln ?

	T = tr1 * tr1 * tr2;

	return T;

}


//******************* log vol. of ball calculation*********************//
double func_vol_ball(long k, double R)
{
	const double PI = 3.1416;
	double  r, kk, x, y, z;
	int i;

	x = 0;
	kk = k;

	R = sqrt(R); // sqrt of radius is taken as R comes here as squared, not be in general case!!!

	//gamma function
	if ((k % 2) != 0)
	{
		for (r = 0.5; r <= kk / 2.0; r++)
			x = x + log(r);
		x = x + 0.5 * log(PI);
	}
	else
	{
		for (i = 1; i <= k / 2; i++)
			x = x + log(i);
	}

	//volume of ball of radius 1
	y = (double(k) / 2.0) * log(PI) + (double(k)) * log(R);
	z = (y - x);


	return (z);
}


// ************Enum. radius R = min(\sqrt(1.1)*{GH(L)}, |b*_j|***************)
xdouble enum_R(vec_RR& c, long k)
{
	xdouble lgdet;

	xdouble *x;
	x = new xdouble[k + 1];

	int i;

	for (i = 1; i <= k; i++)
		x[i] = log(sqrt(to_xdouble(c(i)))); //log|b*_i|

	// lgdet = sum(x, k);
	lgdet = 0;
	for (i = 1; i <= k; i++)
		lgdet = lgdet + x[i];

	xdouble rad, R;

	rad = log(1.04) - (1.0 / double(k) * func_vol_ball(k, 1.0)) + ( 1.0 / double(k) * lgdet);

	if ( exp(to_double(rad)) > exp(to_double(x[1])) ) //store the minimum
		R = exp(to_double(x[1]));
	else
		R = exp(to_double(rad));

	delete [] x;

	return R;
}


//*********log factorial of n numbers*************//
double lgfact(long l)
{

	double x;
	x = 0;
	long i;

	for (i = 1; i <= l; i++)
		x = x + log(i);

	return x;
}


//************Successive Integration Calculation********************//
RR SuccInt(long d, long r, mat_RR& t)
{

	RR I, cons;

	vec_RR p;
	p.SetLength(1000);

	long j, i;

	for (i = 1; i <= d; i++)
		p[i] = 0.0;

	p[1] = 1.0;

	for (j = d; j >= 1; j--)
	{
		for (i = (d - j + 1); i >= 1; i--)
		{
			p[i + 1] = -(1.0 / to_RR(i)) * (p[i]);
		}

		p[1] = 0.0;
		cons = -p[d - j + 2];

		for (i = (d - j + 1); i >= 1; i--)
			cons = cons * t[r][j] - p[i];

		p[1] = cons;
	}

	I = p[1];

	return I;

}


//********Calculating exact value of Psucc******************//
void Init_ExtPrune_P(long k, mat_RR& t,  vec_RR& lgP)
{
	mat_RR  td;
	td.SetDims(k + 1, k + 1);

	RR  I;

	long l, i, r;

	for (r = 2; r <= k; r++)
	{
		l = r / 2 - 1;

		for (i = 1; i <= l; i++)
			td[r][i] = t[r][2 * i];

		I = SuccInt(l, r, td);

		lgP[r] = lgfact(l) + log(I);

	}
}
