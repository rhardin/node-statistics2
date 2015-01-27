(function () {
    'use strict';

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
        B16 = -3617.0 / 510.0;

        // Returns the integral of normal distribution over (-Infty, x].
        var normaldist = function (z) {
            return p_nor(z);
        };

        // Returns the P-value of normaldist(x).
        var pnormaldist = function (y) {
            return pnorm(y);
        };

        // Returns the integral of Chi-squared distribution with n degrees of freedom over [0, x].
        var chi2dist = function (n, x) { return 1.0 - q_chi2(n, x); };

        // Returns the P-value of chi2dist().
        var pchi2dist = function (n, y) { return pchi2(n, 1.0 - y); };

        // Returns the integral of t-distribution with n degrees of freedom over (-Infty, x].
        var tdist = function (n, t) { return p_t(n, t); }

        // Returns the P-value of tdist().
        var ptdist = function (n, y) {
            if (y > 0.5) { return pt(2.0 - y*2.0, n); }
            return -pt(y*2.0, n);
        }


    //////////


    // Newton approximation
    function newton_a(y, ini, epsilon, limit) {
        var x = ini, i = 0, prev, f, df;
        epsilon = typeof epsilon !== 'undefined' ? epsilon : 1.0e-6;
        limit = typeof limit !== 'undefined' ? limit : 30;
        for (i; i < x; i++) {
            prev = x;
            f, df = yield prev;
            x = (y - f) / df + prev;
            if (abs(x - prev) < epsilon) { return x; }
        }
        console.error("Warning(newton approximation): over limit");
        return x;
    }

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
        var e, z2, prev, t, q, i;
        if (z < -12) { return 0.0; }
        if (z > 12) { return 1.0; }
        if (z === 0.0) { return 0.5; }

        if (z > 0.0) {
            e = true;
        } else {
            e = false;
            z = -z;
        }
        z2 = z * z;
        t = (q = z * Math.exp(-0.5 * z2) / SQ2PI);

        for (i; i < 200; i = i + 2) {
            prev = q;
            t *= z2 / i;
            q += t;
            if (q <= prev) { return (e ? 0.5 + q : 0.5 - q); }
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
            s = (t = chi * Math.exp(-0.5 * chi2) / SQ2PI);
            k = 3;
            while (k < df) {
                t = t * (chi2 / k);
                s = s + t;
                k = k + 2;
            }
            return 2 * (normal__x(chi) + s);
        }
        s = (t = Math.exp(-0.5 * chi2));
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
        return 1.0 / Math.pow(2, n2) / gamma(n2) * Math.pow(x, (n2 - 1.0)) * Math.exp(-x/2.0);
    }

    // [x, \infty)
    // Pr([x, \infty)) = y -> x
    function pchi2(n, y) {
        var w, eps, v, s, qe;
        if (n === 1) {
            w = pnorm(1 - y/2);
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
        if (df & 1 !== 0) {
            return 0.5 + (p * Math.sqrt(c2) + Math.atan(t / Math.sqrt(df))) / Math.PI;
        }
        return (1.0 + p) / 2.0;
    }

    // inverse of t-distribution ([2])
    // (-\infty, -q/2] + [q/2, \infty)
    function ptsub(q, n) {
        var eps;

        if (n == 1 && 0.001 < q && q < 0.01) {
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
            if ((qe = 2.0 - p_t(n, w)*2.0 - q) === 0.0) { return w; }
            if (qe < 0.0) {
                w = w - s;
                s = s / 10.0;
            }
        }
    }

    function pt(q, n) {
        var f1 = f2 = f3 = f4 = f5 = u = w0 = w1 = w2 = w3 = w4 = w = 0;
        if(q < 1.0e-5 || q > 1.0 || n < 1) {
            console.error("Error : Illigal parameter in pt()!");
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
})();
