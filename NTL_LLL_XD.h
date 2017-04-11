
typedef long (*LLLCheckFct)(const vec_ZZ&);

extern double LLLStatusInterval;
extern char *LLLDumpFile;

static xdouble InnerProduct(xdouble *a, xdouble *b, long n)
{
    xdouble s;
    long i;

    s = 0;
    for (i = 1; i <= n; i++)
        MulAdd(s, s, a[i], b[i]);

    return s;
}


static void RowTransform(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1)
// x = x - y*MU
{
    static ZZ T, MU;
    long k;

    long n = A.length();
    long i;

    MU = MU1;

    if (MU == 1)
    {
        for (i = 1; i <= n; i++)
            sub(A(i), A(i), B(i));

        return;
    }

    if (MU == -1)
    {
        for (i = 1; i <= n; i++)
            add(A(i), A(i), B(i));

        return;
    }

    if (MU == 0) return;

    if (NumTwos(MU) >= NTL_ZZ_NBITS)
        k = MakeOdd(MU);
    else
        k = 0;


    if (MU.WideSinglePrecision())
    {
        long mu1;
        conv(mu1, MU);

        if (k > 0)
        {

            for (i = 1; i <= n; i++)
            {
                mul(T, B(i), mu1);
                LeftShift(T, T, k);
                sub(A(i), A(i), T);
            }

        }
        else
        {

            for (i = 1; i <= n; i++)
            {
                MulSubFrom(A(i), B(i), mu1);
            }

        }
    }
    else
    {
        for (i = 1; i <= n; i++)
        {
            mul(T, B(i), MU);
            if (k > 0) LeftShift(T, T, k);
            sub(A(i), A(i), T);
        }
    }
}

static void RowTransform2(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1)
// x = x + y*MU
{
    static ZZ T, MU;
    long k;

    long n = A.length();
    long i;

    MU = MU1;

    if (MU == 1)
    {
        for (i = 1; i <= n; i++)
            add(A(i), A(i), B(i));

        return;
    }

    if (MU == -1)
    {
        for (i = 1; i <= n; i++)
            sub(A(i), A(i), B(i));

        return;
    }

    if (MU == 0) return;

    if (NumTwos(MU) >= NTL_ZZ_NBITS)
        k = MakeOdd(MU);
    else
        k = 0;

    if (MU.WideSinglePrecision())
    {
        long mu1;
        conv(mu1, MU);

        for (i = 1; i <= n; i++)
        {
            mul(T, B(i), mu1);
            if (k > 0) LeftShift(T, T, k);
            add(A(i), A(i), T);
        }
    }
    else
    {
        for (i = 1; i <= n; i++)
        {
            mul(T, B(i), MU);
            if (k > 0) LeftShift(T, T, k);
            add(A(i), A(i), T);
        }
    }
}

static
void ComputeGS(mat_ZZ& B, xdouble **B1, xdouble **mu, xdouble *b,
               xdouble *c, long k, xdouble bound, long st, xdouble *buf)
{
    long n = B.NumCols();
    long i, j;
    xdouble s, t1, y, t;
    ZZ T1;

    xdouble *mu_k = mu[k];

    if (st < k)
    {
        for (i = 1; i < st; i++)
            buf[i] = mu_k[i]*c[i];
    }

    for (j = st; j <= k-1; j++)
    {
        if (b[k]*b[j] < NTL_FDOUBLE_PRECISION*NTL_FDOUBLE_PRECISION)
        {
            double z = 0;
            xdouble *B1_k = B1[k];
            xdouble *B1_j = B1[j];

            for (i = 1; i <= n; i++)
                z += B1_k[i].x * B1_j[i].x;

            s = z;
        }
        else
        {
            s = InnerProduct(B1[k], B1[j], n);

            if (s*s <= b[k]*b[j]/bound)
            {
                InnerProduct(T1, B(k), B(j));
                conv(s, T1);
            }
        }

        xdouble *mu_j = mu[j];

        t1 = 0;
        for (i = 1; i <= j-1; i++)
            MulAdd(t1, t1, mu_j[i], buf[i]);

        mu_k[j] = (buf[j] = (s - t1))/c[j];
    }

    s = 0;
    for (j = 1; j <= k-1; j++)
        MulAdd(s, s, mu_k[j], buf[j]);

    c[k] = b[k] - s;
}

static xdouble red_fudge = to_xdouble(0);
static long log_red = 0;

static void init_red_fudge()
{
    long i;

    log_red = long(0.50*NTL_DOUBLE_PRECISION);
    red_fudge = 1;

    for (i = log_red; i > 0; i--)
        red_fudge = red_fudge*0.5;
}

static void inc_red_fudge()
{

    red_fudge = red_fudge * 2;
    log_red--;

    cerr << "LLL_XD: warning--relaxing reduction (" << log_red << ")\n";

    if (log_red < 4)
        Error("LLL_XD: can not continue...sorry");
}


static
long ll_LLL_XD(mat_ZZ& B, mat_ZZ* U, xdouble delta, long deep,
               LLLCheckFct check, xdouble **B1, xdouble **mu,
               xdouble *b, xdouble *c,
               long m, long init_k, long &quit)
{
    //delta = 0.999;

    long n = B.NumCols();


    long i, j, k, Fc1;
    ZZ MU;
    xdouble mu1;

    xdouble t1;
    ZZ T1;
    xdouble *tp;

    static xdouble bound = to_xdouble(0);

    if (bound == 0)
    {
        // we tolerate a 15% loss of precision in computing
        // inner products in ComputeGS.

        bound = 1;
        for (i = 2*long(0.15*NTL_DOUBLE_PRECISION); i > 0; i--)
        {
            bound = bound * 2;
        }
    }


    xdouble half = to_xdouble(0.5);
    xdouble half_plus_fudge = 0.5 + red_fudge;

    quit = 0;
    k = init_k;

    vec_long st_mem;
    st_mem.SetLength(m+2);
    long *st = st_mem.elts();

    for (i = 1; i < k; i++)
        st[i] = i;

    for (i = k; i <= m+1; i++)
        st[i] = 1;

    xdouble *buf;
    buf = NTL_NEW_OP xdouble [m+1];
    if (!buf) Error("out of memory in lll_LLL_XD");

    long rst;
    long counter;

    long trigger_index;
    long small_trigger;
    long cnt;

    long max_k = 0;


    while (k <= m)
    {
        if (k > max_k)
        {
            max_k = k;
        }

        if (st[k] == k)
            rst = 1;
        else
            rst = k;

        if (st[k] < st[k+1]) st[k+1] = st[k];
        ComputeGS(B, B1, mu, b, c, k, bound, st[k], buf);
        st[k] = k;

        counter = 0;
        trigger_index = k;
        small_trigger = 0;
        cnt = 0;

        do
        {
            // size reduction

            counter++;
            if (counter > 10000)
            {
                cerr << "LLL_XD: warning--possible infinite loop\n";
                counter = 0;
            }

            Fc1 = 0;

            for (j = rst-1; j >= 1; j--)
            {
                t1 = fabs(mu[k][j]);
                if (t1 > half_plus_fudge)
                {
                    if (!Fc1)
                    {
                        if (j > trigger_index ||
                                (j == trigger_index && small_trigger))
                        {
                            cnt++;

                            if (cnt > 10)
                            {
                                inc_red_fudge();
                                half_plus_fudge = 0.5 + red_fudge;
                                cnt = 0;
                            }
                        }

                        trigger_index = j;
                        small_trigger = (t1 < 4);
                    }

                    Fc1 = 1;

                    mu1 = mu[k][j];
                    if (mu1 >= 0)
                        mu1 = ceil(mu1-half);
                    else
                        mu1 = floor(mu1+half);


                    xdouble *mu_k = mu[k];
                    xdouble *mu_j = mu[j];

                    if (mu1 == 1)
                    {
                        for (i = 1; i <= j-1; i++)
                            mu_k[i] -= mu_j[i];
                    }
                    else if (mu1 == -1)
                    {
                        for (i = 1; i <= j-1; i++)
                            mu_k[i] += mu_j[i];
                    }
                    else
                    {
                        for (i = 1; i <= j-1; i++)
                            MulSub(mu_k[i], mu_k[i], mu1, mu_j[i]);
                    }

                    mu_k[j] -= mu1;

                    conv(MU, mu1);

                    // cout << j << " " << MU << "\n";

                    RowTransform(B(k), B(j), MU);
                    if (U) RowTransform((*U)(k), (*U)(j), MU);
                }
            }

            if (Fc1)
            {
                for (i = 1; i <= n; i++)
                    conv(B1[k][i], B(k, i));

                b[k] = InnerProduct(B1[k], B1[k], n);
                ComputeGS(B, B1, mu, b, c, k, bound, 1, buf);
            }
        }
        while (Fc1);

        if (check && (*check)(B(k)))
            quit = 1;

        if (b[k] == 0)
        {
            for (i = k; i < m; i++)
            {
                // swap i, i+1
                swap(B(i), B(i+1));
                tp = B1[i];
                B1[i] = B1[i+1];
                B1[i+1] = tp;
                t1 = b[i];
                b[i] = b[i+1];
                b[i+1] = t1;
                if (U) swap((*U)(i), (*U)(i+1));
            }

            for (i = k; i <= m+1; i++) st[i] = 1;

            m--;
            if (quit) break;
            continue;
        }

        if (quit) break;

        if (deep > 0)
        {
            // deep insertions

            xdouble cc = b[k];
            long l = 1;
            while (l <= k-1 && delta*c[l] <= cc)
            {
                cc = cc - mu[k][l]*mu[k][l]*c[l];
                l++;
            }

            if (l <= k-1 && (l <= deep || k-l <= deep))
            {
                // deep insertion at position l

                for (i = k; i > l; i--)
                {
                    // swap rows i, i-1
                    swap(B(i), B(i-1));
                    tp = B1[i];
                    B1[i] = B1[i-1];
                    B1[i-1] = tp;
                    tp = mu[i];
                    mu[i] = mu[i-1];
                    mu[i-1] = tp;
                    t1 = b[i];
                    b[i] = b[i-1];
                    b[i-1] = t1;
                    if (U) swap((*U)(i), (*U)(i-1));
                }

                k = l;
                continue;
            }
        } // end deep insertions

        // test LLL reduction condition

        if (k > 1 && delta*c[k-1] > c[k] + mu[k][k-1]*mu[k][k-1]*c[k-1])
        {
            // swap rows k, k-1
            swap(B(k), B(k-1));
            tp = B1[k];
            B1[k] = B1[k-1];
            B1[k-1] = tp;
            tp = mu[k];
            mu[k] = mu[k-1];
            mu[k-1] = tp;
            t1 = b[k];
            b[k] = b[k-1];
            b[k-1] = t1;
            if (U) swap((*U)(k), (*U)(k-1));

            k--;

            // cout << "- " << k << "\n";
        }
        else
        {
            k++;
            // cout << "+ " << k << "\n";
        }
    }

    delete [] buf;

    return m;
}



long LLL_XD(mat_ZZ& B, vec_RR& cc, mat_RR& muu, mat_ZZ& U, double delta = 0.999, long deep = 0,
            LLLCheckFct check = 0)
{
    long m = B.NumRows();
    long n = B.NumCols();

    long i, j;
    long new_m, dep, quit;
    xdouble s;
    ZZ MU;
    xdouble mu1;

    xdouble t1;
    ZZ T1;

    init_red_fudge();

    // if (U) ident(*U, m);

    xdouble **B1;  // approximates B

    typedef xdouble *xdoubleptr;

    B1 = NTL_NEW_OP xdoubleptr[m+1];
    if (!B1) Error("LLL_XD: out of memory");

    for (i = 1; i <= m; i++)
    {
        B1[i] = NTL_NEW_OP xdouble[n+1];
        if (!B1[i]) Error("LLL_XD: out of memory");
    }

    xdouble **mu;
    mu = NTL_NEW_OP xdoubleptr[m+1];
    if (!mu) Error("LLL_XD: out of memory");

    for (i = 1; i <= m; i++)
    {
        mu[i] = NTL_NEW_OP xdouble[m+1];
        if (!mu[i]) Error("LLL_XD: out of memory");
    }

    xdouble *c; // squared lengths of Gramm-Schmidt basis vectors

    c = NTL_NEW_OP xdouble[m+1];
    if (!c) Error("LLL_XD: out of memory");

    xdouble *b; // squared lengths of basis vectors

    b = NTL_NEW_OP xdouble[m+1];
    if (!b) Error("LLL_XD: out of memory");

    for (i = 1; i <=m; i++)
        for (j = 1; j <= n; j++)
            conv(B1[i][j], B(i, j));

    for (i = 1; i <= m; i++)
    {
        b[i] = InnerProduct(B1[i], B1[i], n);
    }

    new_m = ll_LLL_XD(B, &U, to_xdouble(delta), deep, check, B1, mu, b, c, m, 1, quit);
    dep = m - new_m;
    m = new_m;

    if (dep > 0)
    {
        // for consistency, we move all of the zero rows to the front

        for (i = 0; i < m; i++)
        {
            swap(B(m+dep-i), B(m-i));
            swap(U(m+dep-i), U(m-i));
        }
    }

    for (i = 1; i <=m; i++)
        for (j = 1; j <= n; j++)
            conv(muu(i, j), mu[i][j]);

    for (i = 1; i <=m; i++)
        conv(cc(i), c[i]);


    // clean-up

    for (i = 1; i <= m+dep; i++)
    {
        delete [] B1[i];
    }

    delete [] B1;

    for (i = 1; i <= m+dep; i++)
    {
        delete [] mu[i];
    }

    delete [] mu;

    delete [] c;

    delete [] b;

    return m;
}
