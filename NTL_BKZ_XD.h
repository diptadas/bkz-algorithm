static vec_xdouble BKZConstant;

static
void ComputeBKZConstant(long beta, long p)
{
    const double c_PI = 3.14159265358979323846264338328;
    const double LogPI = 1.14472988584940017414342735135;

    BKZConstant.SetLength(beta-1);

    vec_double Log;
    Log.SetLength(beta);


    long i, j, k;
    double x, y;

    for (j = 1; j <= beta; j++)
        Log(j) = log(double(j));

    for (i = 1; i <= beta-1; i++)
    {
        // First, we compute x = gamma(i/2)^{2/i}

        k = i/2;

        if ((i & 1) == 0)   // i even
        {
            x = 0;
            for (j = 1; j <= k; j++)
                x = x + Log(j);

            x = x * (1/double(k));

            x = exp(x);
        }
        else   // i odd
        {
            x = 0;
            for (j = k + 2; j <= 2*k + 2; j++)
                x = x + Log(j);

            x = 0.5*LogPI + x - 2*(k+1)*Log(2);

            x = x * (2.0/double(i));

            x = exp(x);
        }

        // Second, we compute y = 2^{2*p/i}

        y = -(2*p/double(i))*Log(2);
        y = exp(y);

        BKZConstant(i) = x*y/c_PI;
    }
}

static vec_xdouble BKZThresh;

static
void ComputeBKZThresh(xdouble *c, long beta)
{
    BKZThresh.SetLength(beta-1);

    long i;
    double x;

    x = 0;

    for (i = 1; i <= beta-1; i++)
    {
        x += log(c[i-1]);
        BKZThresh(i) = xexp(x/double(i))*BKZConstant(i);
    }
}

static
long NTL_BKZ_XD(mat_ZZ& BB, long beta)
{
    mat_ZZ* UU = 0;
    xdouble delta = to_xdouble(0.999);
    long prune = 0;
    LLLCheckFct check = 0;


    long m = BB.NumRows();
    long n = BB.NumCols();
    long m_orig = m;

    long i, j;
    ZZ MU;

    xdouble t1;
    ZZ T1;
    xdouble *tp;

    init_red_fudge();

    mat_ZZ B;
    B = BB;

    B.SetDims(m+1, n);


    xdouble **B1;  // approximates B

    typedef xdouble *xdoubleptr;

    B1 = NTL_NEW_OP xdoubleptr[m+2];
    if (!B1) Error("BKZ_XD: out of memory");

    for (i = 1; i <= m+1; i++)
    {
        B1[i] = NTL_NEW_OP xdouble[n+1];
        if (!B1[i]) Error("BKZ_XD: out of memory");
    }

    xdouble **mu;
    mu = NTL_NEW_OP xdoubleptr[m+2];
    if (!mu) Error("BKZ_XD: out of memory");

    for (i = 1; i <= m+1; i++)
    {
        mu[i] = NTL_NEW_OP xdouble[m+1];
        if (!mu[i]) Error("BKZ_XD: out of memory");
    }

    xdouble *c; // squared lengths of Gramm-Schmidt basis vectors

    c = NTL_NEW_OP xdouble[m+2];
    if (!c) Error("BKZ_XD: out of memory");

    xdouble *b; // squared lengths of basis vectors

    b = NTL_NEW_OP xdouble[m+2];
    if (!b) Error("BKZ_XD: out of memory");

    xdouble cbar;

    xdouble *ctilda;
    ctilda = NTL_NEW_OP xdouble[m+2];
    if (!ctilda) Error("BKZ_XD: out of memory");

    xdouble *vvec;
    vvec = NTL_NEW_OP xdouble[m+2];
    if (!vvec) Error("BKZ_XD: out of memory");

    xdouble *yvec;
    yvec = NTL_NEW_OP xdouble[m+2];
    if (!yvec) Error("BKZ_XD: out of memory");

    xdouble *uvec;
    uvec = NTL_NEW_OP xdouble[m+2];
    if (!uvec) Error("BKZ_XD: out of memory");

    xdouble *utildavec;
    utildavec = NTL_NEW_OP xdouble[m+2];
    if (!utildavec) Error("BKZ_XD: out of memory");


    long *Deltavec;
    Deltavec = NTL_NEW_OP long[m+2];
    if (!Deltavec) Error("BKZ_XD: out of memory");

    long *deltavec;
    deltavec = NTL_NEW_OP long[m+2];
    if (!deltavec) Error("BKZ_XD: out of memory");

    mat_ZZ Ulocal;
    mat_ZZ *U;

    if (UU)
    {
        Ulocal.SetDims(m+1, m);
        for (i = 1; i <= m; i++)
            conv(Ulocal(i, i), 1);
        U = &Ulocal;
    }
    else
        U = 0;

    long quit;
    long new_m;
    long z, jj, kk;
    long s, t;
    long h;
    xdouble eta;


    for (i = 1; i <=m; i++)
        for (j = 1; j <= n; j++)
            conv(B1[i][j], B(i, j));


    for (i = 1; i <= m; i++)
    {
        b[i] = InnerProduct(B1[i], B1[i], n);
    }

    // cerr << "\n";
    // cerr << "first LLL\n";

    m = ll_LLL_XD(B, U, delta, 0, check, B1, mu, b, c, m, 1, quit);

    unsigned long NumIterations = 0;
    unsigned long NumTrivial = 0;
    unsigned long NumNonTrivial = 0;
    unsigned long NumNoOps = 0;


    if (m < m_orig)
    {
        for (i = m_orig+1; i >= m+2; i--)
        {
            // swap i, i-1

            swap(B(i), B(i-1));
            if (U) swap((*U)(i), (*U)(i-1));
        }
    }

    long clean = 1;

    if (!quit && m > 1)
    {
        // cerr << "continuing\n";
        if (beta > m) beta = m;

        if (prune > 0)
            ComputeBKZConstant(beta, prune);

        z = 0;
        jj = 0;

        while (z < m-1)
        {
            jj++;
            kk = min(jj+beta-1, m);

            if (jj == m)
            {
                jj = 1;
                kk = beta;
                clean = 1;
            }

            // ENUM

            if (prune > 0)
                ComputeBKZThresh(&c[jj], kk-jj+1);

            cbar = c[jj];
            utildavec[jj] = uvec[jj] = 1;

            yvec[jj] = vvec[jj] = 0;
            Deltavec[jj] = 0;


            s = t = jj;
            deltavec[jj] = 1;

            for (i = jj+1; i <= kk+1; i++)
            {
                ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
                Deltavec[i] = 0;
                vvec[i] = 0;
                deltavec[i] = 1;
            }

            while (t <= kk)
            {
                ctilda[t] = ctilda[t+1] +
                            (yvec[t]+utildavec[t])*(yvec[t]+utildavec[t])*c[t];

                if (prune > 0 && t > jj)
                {
                    eta = BKZThresh(t-jj);
                }
                else
                    eta = 0;

                if (ctilda[t] < cbar - eta)
                {
                    if (t > jj)
                    {
                        t--;
                        t1 = 0;
                        for (i = t+1; i <= s; i++)
                        {
                            t1 += utildavec[i]*mu[i][t];
                        }


                        yvec[t] = t1;
                        t1 = -t1;
                        if (t1 >= 0)
                            t1 = ceil(t1-0.5);
                        else
                            t1 = floor(t1+0.5);

                        utildavec[t] = vvec[t] = t1;
                        Deltavec[t] = 0;
                        if (utildavec[t] > -yvec[t])
                            deltavec[t] = -1;
                        else
                            deltavec[t] = 1;
                    }
                    else
                    {
                        cbar = ctilda[jj];
                        for (i = jj; i <= kk; i++)
                        {
                            uvec[i] = utildavec[i];
                        }
                    }
                }
                else
                {
                    t++;
                    s = max(s, t);
                    if (t < s) Deltavec[t] = -Deltavec[t];
                    if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];
                    utildavec[t] = vvec[t] + Deltavec[t];
                }
            }

            NumIterations++;

            h = min(kk+1, m);

            if ((delta-8*red_fudge)*c[jj] > cbar)
            {

                clean = 0;

                // we treat the case that the new vector is b_s (jj < s <= kk)
                // as a special case that appears to occur most of the time.

                s = 0;
                for (i = jj+1; i <= kk; i++)
                {
                    if (uvec[i] != 0)
                    {
                        if (s == 0)
                            s = i;
                        else
                            s = -1;
                    }
                }

                if (s == 0) Error("BKZ_XD: internal error");

                if (s > 0)
                {
                    // special case

                    NumTrivial++;

                    for (i = s; i > jj; i--)
                    {
                        // swap i, i-1
                        swap(B(i-1), B(i));
                        if (U) swap((*U)(i-1), (*U)(i));
                        tp = B1[i-1];
                        B1[i-1] = B1[i];
                        B1[i] = tp;
                        t1 = b[i-1];
                        b[i-1] = b[i];
                        b[i] = t1;
                    }

                    // cerr << "special case\n";
                    new_m = ll_LLL_XD(B, U, delta, 0, check,
                                      B1, mu, b, c, h, jj, quit);
                    if (new_m != h) Error("BKZ_XD: internal error");
                    if (quit) break;
                }
                else
                {
                    // the general case

                    NumNonTrivial++;

                    for (i = 1; i <= n; i++) conv(B(m+1, i), 0);

                    if (U)
                    {
                        for (i = 1; i <= m_orig; i++)
                            conv((*U)(m+1, i), 0);
                    }

                    for (i = jj; i <= kk; i++)
                    {
                        if (uvec[i] == 0) continue;
                        conv(MU, uvec[i]);
                        RowTransform2(B(m+1), B(i), MU);
                        if (U) RowTransform2((*U)(m+1), (*U)(i), MU);
                    }

                    for (i = m+1; i >= jj+1; i--)
                    {
                        // swap i, i-1
                        swap(B(i-1), B(i));
                        if (U) swap((*U)(i-1), (*U)(i));
                        tp = B1[i-1];
                        B1[i-1] = B1[i];
                        B1[i] = tp;
                        t1 = b[i-1];
                        b[i-1] = b[i];
                        b[i] = t1;
                    }

                    for (i = 1; i <= n; i++)
                        conv(B1[jj][i], B(jj, i));

                    b[jj] = InnerProduct(B1[jj], B1[jj], n);

                    if (b[jj] == 0) Error("BKZ_XD: internal error");

                    // remove linear dependencies

                    // cerr << "general case\n";
                    new_m = ll_LLL_XD(B, U, delta, 0, 0, B1, mu, b, c, kk+1, jj, quit);

                    if (new_m != kk) Error("BKZ_XD: internal error");

                    // remove zero vector

                    for (i = kk+2; i <= m+1; i++)
                    {
                        // swap i, i-1
                        swap(B(i-1), B(i));
                        if (U) swap((*U)(i-1), (*U)(i));
                        tp = B1[i-1];
                        B1[i-1] = B1[i];
                        B1[i] = tp;
                        t1 = b[i-1];
                        b[i-1] = b[i];
                        b[i] = t1;
                    }

                    quit = 0;
                    if (check)
                    {
                        for (i = 1; i <= kk; i++)
                            if ((*check)(B(i)))
                            {
                                quit = 1;
                                break;
                            }
                    }

                    if (quit) break;

                    if (h > kk)
                    {
                        // extend reduced basis

                        new_m = ll_LLL_XD(B, U, delta, 0, check,
                                          B1, mu, b, c, h, h, quit);

                        if (new_m != h) Error("BKZ_XD: internal error");
                        if (quit) break;
                    }
                }

                z = 0;
            }
            else
            {
                // LLL_XD
                // cerr << "progress\n";

                NumNoOps++;

                if (!clean)
                {
                    new_m =
                        ll_LLL_XD(B, U, delta, 0, check, B1, mu, b, c, h, h, quit);
                    if (new_m != h) Error("BKZ_XD: internal error");
                    if (quit) break;
                }

                z++;
            }
        }
    }

    // clean up

    if (m_orig > m)
    {
        // for consistency, we move zero vectors to the front

        for (i = m+1; i <= m_orig; i++)
        {
            swap(B(i), B(i+1));
            if (U) swap((*U)(i), (*U)(i+1));
        }

        for (i = 0; i < m; i++)
        {
            swap(B(m_orig-i), B(m-i));
            if (U) swap((*U)(m_orig-i), (*U)(m-i));
        }
    }

    B.SetDims(m_orig, n);
    BB = B;

    if (U)
    {
        U->SetDims(m_orig, m_orig);
        *UU = *U;
    }

    for (i = 1; i <= m_orig+1; i++)
    {
        delete [] B1[i];
    }

    delete [] B1;

    for (i = 1; i <= m_orig+1; i++)
    {
        delete [] mu[i];
    }

    delete [] mu;


    delete [] c;
    delete [] b;
    delete [] ctilda;
    delete [] vvec;
    delete [] yvec;
    delete [] uvec;
    delete [] utildavec;
    delete [] Deltavec;
    delete [] deltavec;

    return m;
}

