//********Insert new vector and remove linear dependency via LLL, update Gram coefficient (steps 5-8 in CN11)*********//
static
long BKZ_rest(mat_ZZ& B, vec_RR& cc, mat_RR& muu, mat_ZZ& U, vec_ZZ& vc, long k, long jj, long kk, long chk)
{
	long i, j;
	long h, s, new_m;
	long count = 0;

	xdouble t1;

	ZZ MU;

	long m = B.NumRows();
	long n = B.NumCols();
	long m_orig = m;

	B.SetDims(m + 1, n);
	U.SetDims(m + 1, n);

	LLLCheckFct check = 0;
	long quit = 0;
	xdouble delta = to_xdouble(0.999);

	//Steps 5-8 of BKZ (Algorithm 1 in [CN11]) reduction

	xdouble **B1;
	typedef xdouble *xdoubleptr;
	B1 = NTL_NEW_OP xdoubleptr[m + 2];
	for (i = 1; i <= m + 1; i++)
		B1[i] = NTL_NEW_OP xdouble[n + 1];

	xdouble *b; // squared lengths of basis vectors
	b = NTL_NEW_OP xdouble[m + 2];

	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			conv(B1[i][j], B(i, j));

	for (i = 1; i <= m; i++)
		b[i] = InnerProduct(B1[i], B1[i], n);

	xdouble **mu;
	mu = NTL_NEW_OP xdoubleptr[m + 2];

	for (i = 1; i <= m + 1; i++)
	{
		mu[i] = NTL_NEW_OP xdouble[m + 2];
	}

	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			conv(mu[i][j], muu(i, j));

	xdouble *c; // squared lengths of Gramm-Schmidt basis vectors
	c = NTL_NEW_OP xdouble[m + 1];

	for (i = 1; i <= m; i++)
		conv(c[i], cc(i));

	long *v;
	v = NTL_NEW_OP long[m + 1];

	for (i = 1; i <= m; i++)
		conv(v[i], vc(i));

	for (i = jj + 1; i <= kk; i++)
	{
		if (v[i] != 0)
		{
			count++;
			break;
		}
	}

	h = min(kk + 1, m);

	if ((chk == 1) && ( count != 0) ) //checking condition v != 1,0,......,0 and a solution from enumeration
	{

		// we treat the case that the new vector is b_s (jj < s <= kk)
		// as a special case that appears to occur most of the time.

		s = 0;

		for (i = jj + 1; i <= kk; i++)
		{
			if (v[i] != 0)
			{
				if (s == 0) s = i;
				else s = -1;
			}
		}

		if (s == 0) Error("BKZ2.0: internal error 1");

		if (s > 0)
		{
			// special case

			for (i = s; i > jj; i--)
			{
				// swap i, i-1
				swap(B(i - 1), B(i));
				swap(U(i - 1), U(i));

				swap(B1[i - 1], B1[i]);
				swap(b[i - 1], b[i]);
			}

			new_m = ll_LLL_XD(B, &U, delta, 0, check, B1, mu, b, c, h, jj, quit);

			if (new_m != h) Error("BKZ2.0: internal error 2");
			//if (quit) break;
		}
		else
		{
			// the general case

			for (i = 1; i <= n; i++)
				conv(B(m + 1, i), 0);

			for (i = 1; i <= m_orig; i++)
				conv(U(m + 1, i), 0);

			for (i = jj; i <= kk; i++)
			{
				if (v[i] == 0) continue;
				conv(MU, v[i]);
				RowTransform2(B(m + 1), B(i), MU);
				RowTransform2(U(m + 1), U(i), MU);
			}

			for (i = m + 1; i >= jj + 1; i--)
			{
				// swap i, i-1
				swap(B(i - 1), B(i));
				swap(U(i - 1), U(i));

				swap(B1[i - 1], B1[i]);
				swap(b[i - 1], b[i]);
			}

			for (i = 1; i <= n; i++)
				conv(B1[jj][i], B(jj, i));

			b[jj] = InnerProduct(B1[jj], B1[jj], n);

			if (b[jj] == 0) Error("BKZ2.0: internal error 3");

			// remove linear dependencies
			new_m = ll_LLL_XD(B, &U, delta, 0, 0, B1, mu, b, c, kk + 1, jj, quit);

			if (new_m != kk) Error("BKZ2.0: internal error 4");

			// remove zero vector

			for (i = kk + 2; i <= m + 1; i++)
			{
				// swap i, i-1
				swap(B(i - 1), B(i));
				swap(U(i - 1), U(i));

				swap(B1[i - 1], B1[i]);
				swap(b[i - 1], b[i]);
			}

			//if (quit) break;

			if (h > kk)
			{
				// extend reduced basis

				new_m = ll_LLL_XD(B, &U, delta, 0, check, B1, mu, b, c, h, h, quit);

				if (new_m != h) Error("BKZ2.0: internal error 5");
				//if (quit) break;
			}
		}
	}
	else
	{
		new_m = ll_LLL_XD(B, &U, delta, 0, check, B1, mu, b, c, h, h, quit);
		if (new_m != h) Error("BKZ2.0: internal error 6");

	}// if (v[1] == 1........) ends

	if (m_orig > m)
	{
		// for consistency, we move zero vectors to the front

		for (i = m + 1; i <= m_orig; i++)
		{
			swap(B(i), B(i + 1));
			swap(U(i), U(i + 1));
		}

		for (i = 0; i < m; i++)
		{
			swap(B(m_orig - i), B(m - i));
			swap(U(m_orig - i), U(m - i));
		}
	}

	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			conv(muu(i, j), mu[i][j]);

	for (i = 1; i <= m; i++)
		conv(cc(i), c[i]);

	for (i = 1; i <= m; i++)
		conv(vc(i), v[i]);

	B.SetDims(m_orig, n);
	U.SetDims(m_orig, n);

	for (i = 1; i <= m + 1; i++)
		delete [] mu[i];
	delete [] mu;

	for (i = 1; i <= m + 1; i++)
		delete [] B1[i];
	delete [] B1;

	delete [] b;
	delete [] c;
	delete [] v;

	return k;
}

bool dp[100][100];

//**********BKZ with CN11 reduction***********//
static
void MyBKZ(mat_ZZ& B, long n, long k, vec_RR& c, mat_RR& mu, mat_ZZ U, long level, bool isRec, bool isDP, long pointer, bool verb)
{
	if(isDP && dp[pointer][n]) return; //already processed

	long i, j;

	long d_min = 10; //minimum block size for recursive processing
	long step = 10;

	double delta = 0.999; //LLL delta
	double alpha = 0.3; //Pruning parameter(0<alpha<0.5)
	RR eta = to_RR(1.0 / exp(1));

	long TR = to_int(tour(n, k)); //no. of tours count for aborted BKZ preprocessing

	long jj = 1;
	long kk;
	long tr = 1;

	if (verb && level == 0) cout << "Number of Tours: " << TR << endl;

	while (tr <= TR) //no. of tours for aborted BKZ
	{
		if (verb && level == 0 && jj == 1) cout << "Tour no.: " << tr << " ";

		kk = min(jj + k - 1, n);

		long d = kk - jj + 1; //block size

		//branch_factor = to_int(1.0/exp(lgP[d])); //branching factor (1/p_succ)
		long branch_factor = 1;

		long chk = 0;

		mat_RR B_k; // Block GSO
		B_k.SetDims(d, d);
		clear(B_k);
		genBlockGSO(B_k, c, mu, d, jj, kk);

		mat_ZZ BB; //convert RR type to ZZ type
		BB.SetDims(d, d);
		clear(BB);
		conv_RR_to_ZZ(BB, B_k, d);

		RR sum_min = to_RR(10e100);

		vec_ZZ v; //solution vector
		v.SetLength(n + 1);
		clear(v);
		v(jj) = 1;

		//repeat branching factor 1/p_succ times

		for (long seed = 0; seed < branch_factor; seed++)
		{
			// randomize Basis

			mat_ZZ B_RR;
			B_RR.SetDims(d, d);
			clear(B_RR);

			mat_ZZ U_rand;
			U_rand.SetDims(d, d);
			clear(U_rand);
			randomizeBasis(B_RR, BB, U_rand, d, seed);

			// call LLL

			vec_RR c_c;
			c_c.SetLength(d);
			clear(c_c);

			mat_RR m_u;
			m_u.SetDims(d, d);
			clear(m_u);

			mat_ZZ U_p;
			U_p.SetDims(d, d);
			clear(U_p);
			ident(U_p, d);

			LLL_XD(B_RR, c_c, m_u, U_p, delta);


			// apply recursive preprocessing

			mat_ZZ U_prep;
			U_prep.SetDims(d, d);
			clear(U_prep);
			ident(U_prep, d);

			if (isRec && d - step >= d_min)	MyBKZ(B_RR, d, d - step, c_c, m_u, U_prep, level + 1, isRec, isDP, pointer + jj, verb);

			// pruned eumueration with optimal bounding profile R

			vec_ZZ vc;
			vc.SetLength(d);
			clear(vc);

			bool ret = pruned_enum(c_c, m_u, d, vc, alpha, eta);

			// store the best integer combinations that minimizes the short vector

			if (ret == 1)  //solution found
			{
				chk = 1;

				vec_ZZ vk;
				vk.SetLength(d);
				clear(vk);

				vk = ((vc * U_prep) * U_p) * U_rand;

				vec_RR v_RR;
				v_RR.SetLength(d);
				clear(v_RR);

				for (i = 1; i <= d; i++) //vec_ZZ to vec_RR
					conv(v_RR(i), vk(i));

				vec_RR y;
				y.SetLength(d);
				clear(y);

				y = v_RR * B_k;

				RR sum = to_RR(0);

				for (i = 1; i <= d; i++) // inner product
					sum = sum + y(i) * y(i);

				//cout<<"sqr_sum:"<<sqrt(sum)<<endl;//length of the shorest vector

				if (sqrt(sum) < sum_min) //update solution vector w.r.t the basis B_k
				{
					for (i = 1, j = jj; i <= d, j <= kk; i++, j++) //Integer coefficients
					{
						v(j) = vk(i);
					}

					sum_min = sqrt(sum);
				}

			} //if ret ends
		}

		//Steps 5-8 of BKZ (Algorithm 1 in CN11) reduction
		BKZ_rest(B, c, mu, U, v, k, jj, kk, chk);

		jj++;

		if (jj == n)
		{
			if (verb && level == 0) cout << "Norm: " << sqrt(c(1)) << endl;

			jj = 1;
			kk = k;
			tr++;
		}

	}//while tr condition ends

	if(isDP) dp[pointer][n] = true;
}

static
void MyBKZ(mat_ZZ& B, long n, long k, bool rec, bool dp, bool verb)
{
	mat_ZZ U;
	U.SetDims(n + 1, n);
	clear(U);
	ident(U, n);

	vec_RR c;
	c.SetLength(n + 1);
	clear(c);

	mat_RR mu;
	mu.SetDims(n + 2, n + 1);
	clear(mu);

	LLL_XD(B, c, mu, U); //first LLL
	if (verb) cout << "First LLL: " << sqrt(c(1)) << endl;

	clear(U);

	long level = 0; //recursive level
	long pointer = 0;

	MyBKZ(B, n, k, c, mu, U, level, rec, dp, pointer, verb);
}

static
void MyBKZ(mat_ZZ& B, long n, long k, bool verb = false)
{
	MyBKZ(B, n, k, false, false, verb);
}

static
void MyBKZ_Rec(mat_ZZ& B, long n, long k, bool verb = false)
{
	MyBKZ(B, n, k, true, false, verb);
}

static
void MyBKZ_Dp(mat_ZZ& B, long n, long k, bool verb = false)
{
	MyBKZ(B, n, k, true, true, verb);
}

static
void MyBKZ_Rand(mat_ZZ& B, long n, long k, long iteration, bool abort, bool verb = false)
{
	FpBKZ(B, n, k);

	RR norm_base = getEuclidNorm(B(1), n);
	RR norm_min = norm_base;

	mat_ZZ BB = B;

	mat_ZZ BR;
	BR.SetDims(n, n);
	mat_ZZ U_rand;
	U_rand.SetDims(n, n);

	if(verb) cout << "Base: " << norm_base << endl;

	for (int i = 0; i < iteration; i++)
	{
		clear(BR);
		clear(U_rand);

		randomizeBasis(BR, BB, U_rand, n, i);

		FpBKZ(BR, n, k);

		RR norm_cur = getEuclidNorm(BR(1), n);

		if(verb) cout << "Loop: " << i+1 << " Norm: " << norm_cur << endl;

		if(abort && abs(norm_cur - norm_min) < 1) break; //repeated

		if(norm_cur < norm_min)
		{
			norm_min = norm_cur;
			B = BR;
		}

		BB = BR;

	}
}
