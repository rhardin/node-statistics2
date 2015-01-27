/* This is a work in progress. Sharing now to solicit feedback and pull requests */
(function () {
    var SQ2PI = Math.sqrt(2 * Math.PI),
        LOG_2PI = Math.log(2 * Math.PI),
        N = 8,
        B0  = 1.0,
        B1  = -1.0 / 2.0,
        B2  = 1.0 / 6.0,
        B4  = -1.0 / 30.0,
        B6  =  1.0 / 42.0,
        B8  = -1.0 / 30.0,
        B10 =  5.0 / 66.0,
        B12 = -691.0 / 2730.0,
        B14 =  7.0 / 6.0,
        B16 = -3617.0 / 510.0,

        // Returns the integral of normal distribution over (-Infty, x].
        normaldist = function (z) {
            return p_nor(z);
        },

        // Returns the P-value of normaldist(x).
        pnormaldist = function (y) {
            return pnorm(y);
        },

        // Returns the integral of Chi-squared distribution with n degrees of freedom over [0, x].
        chi2dist = function (n, x) { return 1.0 - q_chi2(n, x); },

        // Returns the P-value of chi2dist().
        pchi2dist = function (n, y) { return pchi2(n, 1.0 - y); },

        // Returns the integral of t-distribution with n degrees of freedom over (-Infty, x].
        tdist = function (n, t) { return p_t(n, t); },

        // Returns the P-value of tdist().
        ptdist = function (n, y) {
            if (y > 0.5) { return pt(2.0 - y*2.0, n); }
            return -pt(y*2.0, n);
        },

        // inverse of F-distribution ([2])
        pfsub = function (x, y, z) { return (Math.sqrt(z) - y) / x / 2.0; },

        // Returns the integral of F-distribution with n1 and n2 degrees of freedom over [0, x].
        fdist = function (n1, n2, f) { return 1.0 - q_f(n1, n2, f); },

        // Returns the P-value of fdist().
        pfdist = function (n1, n2, y) { return pf(1.0 - y, n1, n2); },

        /*
         * discrete distributions
         */

        perm = function (n, x) {
            var r;
            x = typeof x !== 'undefined' ? x : n;
            if (n < 0 || x < 0) { console.error("RangeError: perm"); }

            r = 1;
            while (x >= 1) {
                r = r * n;
                n = n - 1;
                x = x - 1;
            }
            return r;
        },

        combi = function (n, x) {
            if (n < 0 || x < 0) { console.error("RangeError: combi"); }
            if (z*2 > n) { x = n - x; }
            return perm(n, x) / perm(x, x);
        },

        bindens = function (n, p, x) {
            q = 1.0 - p;
            return combi(n, x) * Math.pow(p, x) * Math.pow(q, (n - x));
        },

        bindist = function (n, p, x) {
            var i = 0, acc = 0.0;
            for(i; i < x; i++) { acc = acc + bindens(n, p, i); }
            return acc;
        },

        poissondens = function (m, x) {
            if (x < 0) { return 0.0; }
            return Math.pow(m, x) * Math.exp((-m) / perm(x));
        },

        poissondist = function (m, x) {
            var i = 0, acc = 0.0;
            for (i; i < x; i++) { acc = acc + poissondens(m, i); }
            return acc;
        },

        /*
         * normal-distribution
         */

        // Returns the integral of normal distribution over (-Infty, x].
        normalxXX_ = function (z) { return normaldist(z); },

        // Returns the integral of normal distribution over [0, x].
        normal__X_ = function (z) { return normaldist(z) - 0.5; },

        // Returns the integral of normal distribution over [x, Infty).
        normal___x = function (z) { return 1.0 - normaldist(z); },

        // Returns the integral of normal distribution over (-Infty, -x] + [x, Infty).
        normalx__x = function (z) { return 2.0 - normaldist(z) * 2.0; },

        /*
         * inverse of normal-distribution
         */

        // Return the P-value of the corresponding integral.
        pnormalxXX_ = function (z) { return pnormaldist(z); },

        // Return the P-value of the corresponding integral.
        pnormal__X_ = function (y) { return pnormalxXX_(y + 0.5); },

        // Return the P-value of the corresponding integral.
        pnormal___x = function (y) { return pnormalxXX_(1.0 - y); },

        // Return the P-value of the corresponding integral.
        pnormalx__x = function (y) { return pnormalxXX_(1.0 - y/2.0); },

        /*
         * chi2-distribution
         */

        // Returns the integral of Chi-squared distribution with n degrees of freedom over [x, Infty).
        chi2_x = function (n, x) { return 1.0 - chi2dist(n, x); },

        // Returns the integral of Chi-squared distribution with n degrees of freedom over [0, x].
        chi2X_ = function (n, x) { return chi2dist(n, x); },

        /*
         * inverse of chi2-distribution
         */

        // Return the P-value of the corresponding integral.
        pchi2_x = function (n, y) { return pchi2dist(n, 1.0 - y); },

        // Return the P-value of the corresponding integral.
        pchi2X_ = function (n, y) { return pchi2dist(n, y); },

        /*
         * t-distribution
         */

        // Returns the integral of normal distribution with n degrees of freedom over (-Infty, -x] + [x, Infty).
        tx__x = function (n, x) { return 2.0 - tdist(n, x) * 2.0; },

        // Returns the integral of t-distribution with n degrees of freedom over (-Infty, x].
        txXX_ = function (n, x) { return tdist(n, x); },

        // Returns the integral of t-distribution with n degrees of freedom over [0, x].
        t__X_ = function (n, x) { return tdist(n, x) - 0.5; },

        // Returns the integral of t-distribution with n degrees of freedom over [x, Infty).
        t___x = function (n, x) { return 1.0 - tdist(n, x); },

        /*
         * inverse of t-distribution
         */

        // Return the P-value of the corresponding integral.
        ptx__x = function (n, y) { return ptdist(n, 1.0 - y / 2.0); },

        // Return the P-value of the corresponding integral.
        ptxXX_ = function (n, y) { return ptdist(n, y); },

        // Return the P-value of the corresponding integral.
        pt__X_ = function (n, y) { return ptdist(n, 0.5 + y); },

        // Return the P-value of the corresponding integral.
        pt___x = function (n, y) { return ptdist(n, 1.0 - y); },

        /*
         * F-distribution
         */

        // Returns the integral of F-distribution with n1 and n2 degrees of freedom over [x, Infty).
        f_x = function (n1, n2, x) { return 1.0 - fdist(n1, n2, x); },

        // Returns the integral of F-distribution with n1 and n2 degrees of freedom over [0, x].
        fX_ = function (n1, n2, x) { return fdist(n1, n2, x); },


        /*
         * inverse of F-distribution
         */

        // Return the P-value of the corresponding integral.
        pf_x = function (n1, n2, x) { return pfdist(n1, n2, 1.0 - x); },

        // Return the P-value of the corresponding integral.
        pfX_ = function (n1, n2, x) { return pfdist(n1, n2, x); },

        // discrete distributions
        binX_ = function (n, p, x) { return bindist(n, p, x); },
        bin_x = function (n, p, x) { return bindist(n, 1.0 - p, n - x); },

        poissonX_ = function (m, x) { return poissondist(m, x); },
        poisson_x = function (m, x) { return 1.0 - poissondist(m, x-1); };


    //////////////

    /*
     * API
     */
    exports.normaldist = normaldist;
    exports.pnormaldist = pnormaldist;
    exports.chi2dist = chi2dist;
    exports.pchi2dist = pchi2dist;
    exports.tdist = tdist;
    exports.ptdist = ptdist;
    exports.pfsub = pfsub;
    exports.fdist = fdist;
    exports.pfdist = pfdist;
    exports.perm = perm;
    exports.combi = combi;
    exports.bindens = bindens;
    exports.bindist = bindist;
    exports.poissondens = poissondens;
    exports.poissondist = poissondist;
    exports.normalxXX_ = normalxXX_;
    exports.normal__X_ = normal__X_;
    exports.normal___x = normal___x;
    exports.normalx__x = normalx__x;
    exports.pnormalxXX_ = pnormalxXX_;
    exports.pnormal__X_ = pnormal__X_;
    exports.pnormal___x = pnormal___x;
    exports.pnormalx__x = pnormalx__x;
    exports.chi2_x = chi2_x;
    exports.chi2X_ = chi2X_;
    exports.pchi2_x = pchi2_x;
    exports.pchi2X_ = pchi2X_;
    exports.tx__x = tx__x;
    exports.txXX_ = txXX_;
    exports.t__X_ = t__X_;
    exports.t___x = t___x;
    exports.ptx__x = ptx__x;
    exports.ptxXX_ = ptxXX_;
    exports.pt__X_ = pt__X_;
    exports.pt___x = pt___x;
    exports.f_x = f_x;
    exports.fX_ = fX_;
    exports.pf_x = pf_x;
    exports.pfX_ = pfX_;
    exports.binX_ = binX_;
    exports.bin_x = bin_x;
    exports.poissonX_ = poissonX_;
    exports.poisson_x = poisson_x;


    /////////////

    // Gamma function
    function loggamma(x) {
        var v = 1.0,
            w = 1.0 / (x * x),
            ret;
        while (x < N) {
            v = v * x;
            x = x + 1;
        }
        ret = B16 / (16 * 15);
        ret = ret * w + B14 / (14 * 13);
        ret = ret * w + B12 / (12 * 11);
        ret = ret * w + B10 / (10 *  9);
        ret = ret * w + B8  / ( 8 *  7);
        ret = ret * w + B6  / ( 6 *  5);
        ret = ret * w + B4  / ( 4 *  3);
        ret = ret * w + B2  / ( 2 *  1);
        ret = ret / x + 0.5 * LOG_2PI - Math.log(v) - x + (x - 0.5) * Math.log(x);
        return ret;
    }

    function gamma(x) {
        if (x < 0.0) {
            return Math.PI / (Math.sin(Math.PI * x) * Math.exp(loggamma(1 - x)));
        }
        return Math.exp(loggamma(x));
    }

    //normal-distribution
    // (-\infty, z]
    function p_nor(z) {
        var e, z2, prev = 0.0, t, p, i = 3;

        if (z < -12) { return 0.0; }
        if (z > 12) { return 1.0; }
        if (z === 0.0) { return 0.5; }
        if (z > 0) {
            e = true;
        } else {
            e = false;
            z = -z;
        }

        z2 = z * z;
        p = z * Math.exp(-0.5 * z2) / SQ2PI;
        t = p;

        for (i; i < 200; i = i + 2) {
            prev = p;
            t = t * (z2 / i);
            p = p + t;
            if (p <= prev) { return (e ? 0.5 + p : 0.5 - p); }
        }

        return (e ? 1.0 : 0.0);
    }

    // inverse of normal distribution ([2])
    // Pr( (-\infty, x] ) = qn -> x
    function pnorm(qn) {
        var b = [1.570796288, 0.03706987906, -0.8364353589e-3,
                -0.2250947176e-3, 0.6841218299e-5, 0.5824238515e-5,
                -0.104527497e-5, 0.8360937017e-7, -0.3231081277e-8,
                0.3657763036e-10, 0.6936233982e-12],
            w1 = qn,
            w3 = -Math.log(4.0 * w1 * (1.0 - w1)),
            i = 1;

        if (qn < 0.0 || 1.0 < qn) { return 0.0; }
        if (qn === 0.5) { return 0.0; }
        if (qn > 0.5) { w1 = 1.0 - w1; }

        w1 = b[0];
        for (i; i < 11; i++) {
            w1 += b[i] * Math.pow(w3, i);
        }

        if (qn > 0.5) { return Math.sqrt(w1 * w3); }

        return -Math.sqrt(w1 * w3);
    }

    // chi-square distribution ([1])
    // [x, \infty)
    function q_chi2(df, chi2) {
        var chi, t = 0, k = 0, s = 0;
        if ((df & 1) !== 0) {
            chi = Math.sqrt(chi2);

            if (df === 1) { return 2 * normal__x(chi); }

            t = chi * Math.exp(-0.5 * chi2) / SQ2PI;
            s = t;
            k = 3;
            while (k < df) {
                t = t * (chi2 / k);
                s = s + t;
                k = k + 2;
            }

            return 2 * (normal__x(chi) + s);
        }

        t = Math.exp(-0.5 * chi2);
        s = t;
        k = 2;
        while (k < df) {
            t = t * (chi2 / k);
            s = s + t;
            k = k + 2;
        }

        return s;
    }

    function chi2dens(n, x) {
        var n2;

        if (n === 1) {
            return 1.0 / Math.sqrt(2 * Math.PI * x) * Math.exp(-x / 2.0);
        } else if (n === 2) {
            return 0.5 * Math.exp(-x / 2.0);
        }

        n2 = n / 2;

        return 1.0 / Math.pow(2, n2) / gamma(n2) * Math.pow(x, (n2 - 1.0)) * Math.exp(-x / 2.0);
    }

    // [x, \infty)
    // Pr([x, \infty)) = y -> x
    function pchi2(n, y) {
        var w, eps, v, s, qe;

        if (n === 1) {
            w = pnorm(1 - y / 2);
            return w * w;
        } else if (n === 2) {
            return -2.0 * Math.log(y);
        }

        eps = 1.0e-5;
        v = 0.0;
        s = 10.0;
        while (true) {
            v = v + 2;

            if (s <= eps) { break; }
            if ((qe = q_chi2(n, v) - y) === 0.0) { break; }
            if (qe < 0.0) {
                v = v - s;
                s = s / 10.0;
            }
        }
        return v;
    }

    // t-distribution ([1])
    // (-\infty, x]
    function p_t(df, t) {
        var c2, s, p, i;

        c2 = df / (df + t * t);
        s = Math.sqrt(1.0 - c2);

        if (t < 0.0) { s = -s; }

        p = 0.0;
        i = df % 2 + 2;
        while (i <= df) {
            p = p + s;
            s = s * ((i - 1) * c2 / i);
            i = i + 2;
        }

        if ((df & 1) !== 0) {
            return 0.5 + (p * Math.sqrt(c2) + Math.atan(t / Math.sqrt(df))) / Math.PI;
        }

        return (1.0 + p) / 2.0;
    }

    // inverse of t-distribution ([2])
    // (-\infty, -q/2] + [q/2, \infty)
    function ptsub(q, n) {
        var eps;

        if (n == 1 && q > 0.001 && q < 0.01) {
            eps = 1.0e-4;
        } else if (n == 2 && q < 0.0001) {
            eps = 1.0e-4;
        } else if (n == 1 && q < 0.001) {
            eps = 1.0e-2;
        } else {
            eps = 1.0e-5;
        }

        s = 10000.0;
        w = 0.0;
        while (true) {
            w = w + s;

            if (s <= eps) { return w; }
            if ((qe = 2.0 - p_t(n, w) * 2.0 - q) === 0.0) { return w; }
            if (qe < 0.0) {
                w = w - s;
                s = s / 10.0;
            }
        }
    }

    function pt(q, n) {
        var f1 = f2 = f3 = f4 = f5 = u = w0 = w1 = w2 = w3 = w4 = w = 0; //leak?

        if (q < 1.0e-5 || q > 1.0 || n < 1) {
            console.error('Error : Illigal parameter in pt()!');
            return 0.0;
        }

        if (n <= 5) { return ptsub(q, n); }
        if (q <= 5.0e-3 && n <= 13) { return ptsub(q, n); }

        f1 = 4.0 * (f = n);
        f5 = (f4 = (f3 = (f2 = f * f) * f) * f) * f;
        f2 = f2 * 96.0;
        f3 = f3 * 384.0;
        f4 = f4 * 92160.0;
        f5 = f5 * 368640.0;
        u = pnormaldist(1.0 - q / 2.0);

        w0 = (u2 = u * u) * u;
        w1 = w0 * u2;
        w2 = w1 * u2;
        w3 = w2 * u2;
        w4 = w3 * u2;
        w = (w0 + u) / f1;
        w = w + ((5.0 * w1 + 16.0 * w0 + 3.0 * u) / f2);
        w = w + ((3.0 * w2 + 19.0 * w1 + 17.0 * w0 - 15.0 * u) / f3);
        w = w + ((79.0 * w3 + 776.0 * w2 + 1482.0 * w1 - 1920.0 * w0 - 9450.0 * u) / f4);
        w = w + ((27.0 * w4 + 339.0 * w3 + 930.0 * w2 - 1782.0 * w1 - 765.0 * w0 + 17955.0 * u) / f5);
        return u + w;
    }

    // F-distribution ([1])
    // [x, \infty)
    function q_f(df1, df2, f) {
        var cos2, sin2, prob, temp, i;

        if (f <= 0.0) { return 1.0; }
        if (df1 % 2 !== 0 && df2 % 2 === 0) { return 1.0 - q_f(df2, df1, 1.0 / f); }

        cos2 = 1.0 / (1.0 + df1 * f / df2);
        sin2 = 1.0 - cos2;

        if (df1 % 2 === 0) {
            prob = Math.pow(cos2, (df2 / 2.0));
            temp = prob;

            i = 2;
            while (i < df1) {
                temp = temp * ((df2 + i - 2) * sin2 / i);
                prob = prob + temp;
                i = i + 2;
            }

            return prob;
        }

        prob = Math.atan(Math.sqrt(df2 / (df1 * f)));
        temp = Math.sqrt(sin2 * cos2);
        i = 3;
        while (i <= df1) {
            prob = prob + temp;
            temp = temp * ((i - 1) * sin2 / i);
            i = i + 2.0;
        }

        temp = temp * df1;
        i = 3;
        while (i <= df2) {
            prob = prob - temp;
            temp = temp * ((df1 + i - 2) * cos2 / i);
            i = i + 2;
        }

        return prob * 2.0 / Math.PI;
    }

    // [x, \infty)
    function pf(q, n1, n2) {
        var eps, fw, s, qe, w1, w2, w3, w4, u, u2, a, b, c, d;

        if (q < 0.0 || q > 1.0 || n1 < 1 || n2 < 1) {
            console.error('Error : Illegal parameter in pf()!');
            return 0.0;
        }

        if (n1 <= 240 || n2 <= 240) { eps = 1.0e-5; }
        if (n2 === 1) { eps = 1.0e-4; }

        fw = 0.0;
        s = 1000.0;
        while (true) {
            fw = fw + s;
            if (s <= eps) { return fw; }
            if ((qe = q_f(n1, n2, fw) - q) === 0.0) { return fw; }
            if (qe < 0.0) {
                fw = fw - s;
                s = s / 10.0;
            }
        }

        eps = 1.0e-6;
        qn = q;

        if (q < 0.5) { qn = 1.0 - q; }

        u = pnorm(qn);
        w1 = 2.0 / n1 / 9.0;
        w2 = 2.0 / n2 / 9.0;
        w3 = 1.0 - w1;
        w4 = 1.0 - w2;
        u2 = u * u;
        a = w4 * w4 - u2 * w2;
        b = -2.0 * w3 * w4;
        c = w3 * w3 - u2 * w1;
        d = b * b - 4 * a * c;

        if (d < 0.0) {
            fw = pfsub(a, b, 0.0);
        } else {
            if (a.abs > eps) {
                fw = pfsub(a, b, d);
            } else {
                if (abs(b) > eps) { return -c / b; }
                fw = pfsub(a, b, 0.0);
            }
        }

        return fw * fw * fw;
    }
})();
