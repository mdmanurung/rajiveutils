// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// RobRSVD1_cpp
//
// Exact C++ mirror of R/RobustSVD.R: RobRSVD1()
//
// Computes a single robust rank-1 SVD component via iterative M-estimation
// (Huber weighting).  Every operation, matrix construction, and convergence
// criterion is a direct translation of the R source.  No algebraic
// simplifications are applied.
//
// Argument order matches the R version: data, sinit, uinit, vinit, then the
// optional tuning parameters huberk / niter / tol.
//
// R dimensions:
//   data  : m x n
//   uinit : length m  (first left  singular vector of data)
//   vinit : length n  (first right singular vector of data)
//   sinit : scalar    (first singular value of data)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List RobRSVD1_cpp(
    const arma::mat& data,
    double            sinit,
    const arma::vec&  uinit,
    const arma::vec&  vinit,
    double huberk = 1.345,
    int    niter  = 1000,
    double tol    = 1e-5
) {
    // R: size_data = c(dim(data)); m = size_data[1]; n = size_data[2]
    int m = (int)data.n_rows;
    int n = (int)data.n_cols;

    // R: sold = sinit;  vold = vinit;  uold = sold * uinit
    double   sold = sinit;
    arma::vec vold = vinit;
    arma::vec uold = sold * uinit;

    // R: Appold = uold %*% t(vold)   [m x n]
    arma::mat Appold = uold * vold.t();

    // R: Rmat = data - Appold
    arma::mat Rmat = data - Appold;

    // Initial robust scale estimate; updated every iteration below.
    arma::vec Rvec    = arma::vectorise(Rmat);
    double    mysigma = arma::median(arma::abs(Rvec)) / 0.675;
    if (!std::isfinite(mysigma) || mysigma <= 0.0) {
        mysigma = std::sqrt(std::numeric_limits<double>::epsilon());
    }

    // Hoist identity matrices — their dimensions are fixed throughout the loop
    arma::mat eye_m = arma::eye<arma::mat>(m, m);
    arma::mat eye_n = arma::eye<arma::mat>(n, n);

    int    iter      = 1;
    double localdiff = 9999.0;

    while (localdiff > tol && iter < niter) {

        // W-M7: recompute robust scale each iteration from the current
        // residual to avoid freezing Huber weights at initialization.
        Rvec = arma::vectorise(Rmat);
        mysigma = arma::median(arma::abs(Rvec)) / 0.675;
        if (!std::isfinite(mysigma) || mysigma <= 0.0) {
            mysigma = std::sqrt(std::numeric_limits<double>::epsilon());
        }

        // R: Wmat = huberk / abs(Rmat / mysigma)
        arma::mat Wmat = huberk / arma::abs(Rmat / mysigma);

        // R: Wmat[Wmat > 1] = 1
        Wmat.elem(arma::find(Wmat > 1.0)).fill(1.0);

        // ---- uterm1 ----
        // R: diag(c(vold^2)) %*% t(Wmat)
        //    diag(c(vold^2)) is n x n diagonal;  t(Wmat) is n x m  →  result n x m
        //    colSums(n x m)  →  length-m row-vector  →  m x m diagonal matrix
        arma::rowvec cs_u = arma::sum(
            arma::diagmat(arma::pow(vold, 2)) * Wmat.t(),
            0   // column sums
        );

        // R: (2 * mysigma^2) * (c(t(vold) %*% (diag(n)) %*% vold) * (diag(m))
        //                       - diag(sum(vold^2), m))
        //    t(vold) %*% diag(n) %*% vold  =  scalar (sum of vold^2)
        //    c(...) * diag(m)               =  scalar * m x m identity
        //    diag(sum(vold^2), m)           =  scalar * m x m identity
        //    difference                     =  0  (preserved faithfully)
        double    sum_vold2    = arma::dot(vold, vold);
        arma::mat uterm1_second =
            (2.0 * mysigma * mysigma) *
            (arma::as_scalar(vold.t() * eye_n * vold) * eye_m -
             arma::diagmat(arma::ones<arma::vec>(m) * sum_vold2));

        arma::mat uterm1 = arma::diagmat(cs_u) + uterm1_second;

        // R: uterm2 = (Wmat * data) %*% vold    (* is element-wise in R)
        arma::vec uterm2 = (Wmat % data) * vold;

        // R: unew = solve(uterm1) %*% uterm2
        arma::vec unew = arma::solve(uterm1, uterm2);

        // ---- vterm1 ----
        // R: diag(c(unew^2)) %*% Wmat
        //    diag(c(unew^2)) is m x m diagonal;  Wmat is m x n  →  result m x n
        //    colSums(m x n)  →  length-n row-vector  →  n x n diagonal matrix
        arma::rowvec cs_v = arma::sum(
            arma::diagmat(arma::pow(unew, 2)) * Wmat,
            0   // column sums
        );

        // R: (2 * mysigma^2) * (c(t(unew) %*% (diag(m)) %*% unew) * (diag(n))
        //                       - diag(sum(unew^2), n))
        double    sum_unew2    = arma::dot(unew, unew);
        arma::mat vterm1_second =
            (2.0 * mysigma * mysigma) *
            (arma::as_scalar(unew.t() * eye_m * unew) * eye_n -
             arma::diagmat(arma::ones<arma::vec>(n) * sum_unew2));

        arma::mat vterm1 = arma::diagmat(cs_v) + vterm1_second;

        // R: vterm2 = t(Wmat * data) %*% unew
        arma::vec vterm2 = (Wmat % data).t() * unew;

        // R: vnew = solve(vterm1) %*% vterm2
        arma::vec vnew = arma::solve(vterm1, vterm2);

        // R: Appnew = unew %*% t(vnew)
        arma::mat Appnew = unew * vnew.t();

        // R: Rmat = data - Appnew
        Rmat = data - Appnew;

        // R: localdiff = max(abs(Appnew - Appold))
        localdiff = arma::abs(Appnew - Appold).max();

        // R: Appold = Appnew
        Appold = Appnew;

        // R: uold = sqrt(sum(vnew^2)) * unew;  vold = vnew / sqrt(sum(vnew^2))
        double norm_vnew = std::sqrt(arma::dot(vnew, vnew));
        uold = norm_vnew * unew;
        vold = vnew / norm_vnew;

        iter++;
    }

    // R: v = vold;  s = sqrt(sum(uold^2));  u = uold / sqrt(sum(uold^2))
    arma::vec v   = vold;
    double    s   = std::sqrt(arma::dot(uold, uold));
    arma::vec u   = uold / s;

    // Return plain R numeric vectors (not 1-column matrices).
    // std::vector<double> wraps to a dim-free R vector via Rcpp::wrap.
    return List::create(
        Named("s") = s,
        Named("u") = std::vector<double>(u.begin(), u.end()),
        Named("v") = std::vector<double>(v.begin(), v.end())
    );
}


// ---------------------------------------------------------------------------
// RobRSVD_all_cpp
//
// Exact C++ mirror of R/RobustSVD.R: RobRSVD.all()
//
// Sequentially extracts nrank robust rank-1 components by deflation.
// For each deflated residual, initialize RobRSVD1 from the leading classical
// SVD triplet of that residual (warm-start in the current subproblem).
//
// Returns a named list: list(d = numeric vector,
//                            u = m x nrank matrix,
//                            v = n x nrank matrix)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List RobRSVD_all_cpp(
    const arma::mat& data,
    int              nrank,
    double           sinit1,
    const arma::vec& uinit1,
    const arma::vec& vinit1,
    double huberk = 1.345,
    int    niter  = 1000,
    double tol    = 1e-5
) {
    // R: Rm = min(min(dim(data)), nrank)
    int Rm = std::min({(int)data.n_rows, (int)data.n_cols, nrank});

    // --- First component (rank 1) ---
    // R: data.svd1 <- RobRSVD1(data, sinit=svdinit$d[1],
    //                           uinit=svdinit$u[,1], vinit=svdinit$v[,1])
    List res0 = RobRSVD1_cpp(data, sinit1, uinit1, vinit1, huberk, niter, tol);

    // Accumulate d (vector), U (m x r), V (n x r)
    arma::vec d(Rm);
    arma::mat U(data.n_rows, Rm);
    arma::mat V(data.n_cols, Rm);

    d(0) = Rcpp::as<double>(res0["s"]);
    U.col(0) = Rcpp::as<arma::vec>(res0["u"]);
    V.col(0) = Rcpp::as<arma::vec>(res0["v"]);

    // R: Red = d * u %*% t(v)   (scalar * outer product)
    arma::mat Red = d(0) * U.col(0) * V.col(0).t();

    // --- Ranks 2 .. Rm ---
    // R: for(i in 1:(Rm-1)) { data.svd1 <- RobRSVD1(data - Red, ...) }
    for (int i = 1; i < Rm; i++) {
        arma::mat residual = data - Red;

        arma::mat U0;
        arma::vec s0;
        arma::mat V0;
        bool ok = arma::svd_econ(U0, s0, V0, residual);

        if (!ok || s0.n_elem == 0) {
            break;
        }

        List resi = RobRSVD1_cpp(residual,
                                 s0(0),
                                 U0.col(0),
                                 V0.col(0),
                                 huberk, niter, tol);

        d(i)      = Rcpp::as<double>(resi["s"]);
        U.col(i)  = Rcpp::as<arma::vec>(resi["u"]);
        V.col(i)  = Rcpp::as<arma::vec>(resi["v"]);

        // R: Red <- (u %*% diag(d) %*% t(v))
        Red = U.cols(0, i) * arma::diagmat(d.head(i + 1)) * V.cols(0, i).t();
    }

    // Return d as a plain R numeric vector (arma::vec wraps to 1-column matrix).
    return List::create(
        Named("d") = std::vector<double>(d.begin(), d.end()),
        Named("u") = U,
        Named("v") = V
    );
}
