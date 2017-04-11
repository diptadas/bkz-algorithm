//********* CN11 reduction *************//
bool pruned_enum(vec_RR c, mat_RR mu, long long n, vec_ZZ& v, double alpha, RR eta, long bound_type = 2)
{
    long cnt = 0;

    xdouble R[n + 1];

    if (bound_type == 1) //piecewise linear bounding function
    {
        for (int i = 1; i <= n; i++)
        {
            if (i <= n / 2) conv(R[i], sqrt(2 * i * alpha * c(1) / n));
            else conv(R[i], ((2 * alpha - 1 + 2 * i * (1 - alpha) / n)*c(1)));
        }
    }
    else if (bound_type == 2) //modified PLBF
    {
        xdouble R_opt;
        R_opt =  enum_R(c, n) * enum_R(c, n);
        for (int i = 1; i <= n; i++)
        {
            if (i <= double(n)*to_double(eta))
                R[i] = ((double(i) * alpha) / (to_double(eta) * double(n))) * to_double(R_opt);
            else
                R[i] = ( alpha +  (1.0 - alpha) / (double(n) - to_double(eta) * double(n)) * (double(i) - to_double(eta) * double(n)) ) * to_double(R_opt);
        }
    }

    xdouble sig[n + 2][n + 1];

    for (int i = 1; i <= n + 1; i++)
    {
        for (int j = 1; j <= n; j++) conv(sig[i][j], 0);
    }

    xdouble rho[n + 1]; //partial norm
    xdouble cen[n + 1]; //centers
    long long r[n + 1];
    long long w[n + 1]; //jumps

    for (int i = 1; i <= n; i++)
    {
        r[i] = i;
        rho[i] = 0;
        cen[i] = 0;
        v(i) = 0;
        w[i] = 0;
    }

    r[0] = 0;
    rho[n + 1] = 0;
    v(1) = 1;

    long long k = 1;
    long long last_nonzero = 1;

    while (true)
    {
        cnt++;

        //compute norm-squared of current node
        rho[k] = rho[k + 1] + ((to_double(v(k)) - cen[k]) * (to_double(v(k)) - cen[k])) * to_double(c(k));

        if (rho[k] <= R[n + 1 - k]) //we are below the bound
        {
            if (k == 1) return true; //solution found

            else
            {
                k--; //going down the tree

                r[k - 1] = max(r[k - 1], r[k]); //to maintain the invariant for j<k

                for (int i = r[k]; i >= (k + 1); i--)
                {
                    sig[i][k] = sig[i + 1][k] + to_double(v(i)) * to_double(mu(i, k));
                }

                cen[k] = -(sig[k + 1][k]);

                xdouble t = cen[k];

                if (t >= 0)
                    t = ceil(t - 0.5);
                else
                    t = floor(t + 0.5);

                //v(k) = to_ZZ(t); //nearest integer of cen[k]
                v(k) = to_ZZ(cen[k] + 0.5); //nearest integer of cen[k]
                w[k] = 1;
            }
        }
        else
        {
            k++; //going up the tree

            if (k == (n + 1)) return false; //no solution

            r[k - 1] = k;

            //update v[k]
            if (k >= last_nonzero)
            {
                last_nonzero = k;
                v(k) = v(k) + 1;
            }
            else
            {
                if (to_double(v(k)) > cen[k]) v(k) = v(k) - w[k];
                else v(k) = v(k) + w[k];

                w[k] = w[k] + 1;
            }
        }
    }
}



