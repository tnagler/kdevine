#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector kern_gauss(const NumericVector x) {
    NumericVector out(x.size());
    for (int i = 0; i < x.size(); ++i) {
        if ((fabs(x[i]) >= 5)){
            out[i] = 0;
        } else {
            // normalize by 0.9999994267 because of truncation
            out[i] = exp(- 0.5 * pow(x[i], 2)) / (sqrt(2 * M_PI)) / 0.9999994267;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector ikern_gauss(const NumericVector x) {
    NumericVector out(x.size());
    NumericVector xx(1);
    for (int i = 0; i < x.size(); ++i) {
        if (x[i] <= -5.0) {
            out[i] = 0.0;
        } else if (x[i] > 5.0) {
            out[i] = 1.0;
        } else {
            xx[0] = x[i];
            // pnorm(-5.0) = 2.866516e-07
            out[i] = (pnorm(xx)[0] - 2.866516e-07) / 0.9999994267;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector eval_kde1d(const NumericVector xsort,
                         const NumericVector xev,
                         const double xmin,
                         const double xmax,
                         const double bw)
{
    // The vector xsort has to be in ascending order !!!
    NumericVector tmp(xev.size());
    NumericVector out(xev.size());
    NumericVector xeveff(xev.size());
    NumericVector xevdst(xev.size());
    double n = xsort.size();

    for (int i = 0; i < xev.size(); ++i) {
        if ((xev[i] < xmin) | (xev[i] > xmax)) {
            out[i] = 0;
        } else {
            if (xev[i] < xsort[0] - 0.99 * bw) {
                xevdst[i] = xsort[0] - xev[i];
                xeveff[i] = xsort[0] - 0.99 * bw;
            } else if (xev[i] > xsort[n - 1] + 0.99 * bw) {
                xevdst[i] = xev[i] - xsort[n - 1];
                xeveff[i] = xsort[n - 1] + 0.99 * bw;
            } else {
                xeveff[i] = xev[i];
                xevdst[i] = 0;
            }
            if (pow(xevdst[i], 2) > 5 * bw) {
                out[i] = 0;
            } else {
                tmp = kern_gauss((xsort - xeveff[i]) / bw);
                if(!(xmin != xmin))
                    tmp = tmp + kern_gauss((2 * xmin - xsort - xeveff[i]) / bw);
                if(!(xmax != xmax))
                    tmp = tmp + kern_gauss((2 * xmax - xsort - xeveff[i]) / bw);
                out[i] = sum(tmp) / (n * bw) * exp(- 0.5 *  pow(xevdst[i], 2) / bw) ;
            }
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector eval_pkde1d(const NumericVector x,
                          const NumericVector xev,
                          const double xmin,
                          const double xmax,
                          const double bw)
{
    NumericVector tmp(xev.size());
    NumericVector out(xev.size());
    double n = x.size();

    for (int i = 0; i < xev.size(); ++i) {
        if (xev[i] <= xmin) {
            out[i] = 0;
        } else if (xev[i] >= xmax) {
            out[i] = 1;
        } else {
            // integrate kernels on regular data
            tmp = ikern_gauss((xev[i] - x) / bw);
            if (!(xmin != xmin)) {
                // substract integral up to xmin
                tmp += - ikern_gauss((xmin - x) / bw);
                // add integrals for reflected data below xmin
                tmp += ikern_gauss((2 * xmin - x - xmin) / bw);
                // substract integral up to xmin
                tmp += - ikern_gauss((2 * xmin - x - xev[i]) / bw);
            }
            if (!(xmax != xmax)) {
                // add integrals for reflected data above xmax
                tmp += ikern_gauss(-(2 * xmax - x - xev[i]) / bw);
            }
            out[i] = sum(tmp) / n;
        }
    }
    return out;
}

// [[Rcpp::export]]
NumericVector eval_qkde1d(const NumericVector x,
                          const NumericVector qev,
                          const double xmin,
                          const double xmax,
                          const double bw)
{
    NumericVector out(qev.size()), x0, x1, ans, val;
    ans = 0.0, val = 0.0;
    double tol = ::fmax(1e-10 * (x1[0] - x0[0]), 1e-10);

    for (int i = 0; i < qev.size(); ++i) {
        if (::fabs(qev[i]) < 1e-30)
            out[i] = (xmin != xmin) ? R_NegInf : xmin;
        if (::fabs(qev[i] - 1) < 1e-30)
            out[i] = (xmax != xmax) ? R_PosInf : xmin;;

        int br = 0;
        x0 = min(x) - 5 * bw;
        x1 = max(x) + 5 * bw;
        NumericVector ql = eval_pkde1d(x, x0, xmin, xmax, bw);
        NumericVector qh = eval_pkde1d(x, x1, xmin, xmax, bw);
        ql = ql - qev[i];
        qh = qh - qev[i];
        if (::fabs(ql[0]) <= tol) {
            ans = x0;
            br = 1;
        }
        if (::fabs(qh[0]) <= tol) {
            ans = x1;
            br = 1;
        }

        int maxit = 20;
        for (int it = 0; it < maxit; ++it) {
            ans[0] = (x0[0] + x1[0]) / 2.0;
            NumericVector val = eval_pkde1d(x, ans, xmin, xmax, bw);
            val[0] = val[0] - qev[i];
            //stop if values become too close (avoid infinite loop)
            if (::fabs(val[0]) <= tol)
                br = 1;
            if (::fabs(x0[0] - x1[0]) <= tol)
                br = 1;

            if (val[0] > 0.0) {
                x1[0] = ans[0];
                qh[0] = val[0];
            } else {
                x0[0] = ans[0];
                ql[0] = val[0];
            }

            if (br == 1)
                break;
        }

        out[i] = ans[0];
    }

    return out;
}


