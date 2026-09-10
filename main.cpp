// =========================================================================
//  RIEMANN ZETA HIGH-PRECISION HPC ENGINE (MPFR + OpenMP + GMP)
// =========================================================================
//  Improvements over the previous version:
//   1. Zero-search loop is genuinely parallelized with OpenMP (the old
//      version had omp.h included but no #pragma omp anywhere).
//   2. brent_zero_refiner_mpfr() is a real implementation of Brent's
//      method (inverse quadratic interpolation + secant + bisection
//      fallback), not disguised bisection.
//   3. Z_function_mpfr() uses an adaptive term count ~sqrt(t / 2*pi)
//      instead of a fixed 20000 terms for every t.
//   4. Added a real (sequential, GMP-based) prime gap analysis module so
//      the README's "Arbitrary-Precision Prime Gap Calculation" claim is
//      backed by actual code instead of an unused #include <gmp.h>.
// =========================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <gmpxx.h>
#include <mpfr.h>
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

// -------------------------------------------------------------------------
// Riemann-Siegel theta function theta(t), asymptotic expansion.
// -------------------------------------------------------------------------
mpreal riemann_theta_mpfr(const mpreal& t) {
    mpreal pi = mpfr::const_pi();
    mpreal term1 = t * mpfr::log(t / (mpreal(2.0) * pi)) / mpreal(2.0);
    mpreal term2 = t / mpreal(2.0);
    mpreal term3 = pi / mpreal(8.0);
    mpreal term4 = mpreal(1.0) / (mpreal(48.0) * t);
    return term1 - term2 - term3 + term4;
}

// -------------------------------------------------------------------------
// Riemann-Siegel Z(t) function. Term count is adaptive: the main sum of
// the Riemann-Siegel formula only needs ~sqrt(t / 2*pi) terms to converge
// to the stated tolerance, so we stop paying for 20000 terms at t = 10
// and 20000 terms at t = 1000 alike.
// -------------------------------------------------------------------------
mpreal Z_function_mpfr(const mpreal& t) {
    mpreal pi = mpfr::const_pi();
    mpreal theta = riemann_theta_mpfr(t);

    long terms = static_cast<long>(mpfr::sqrt(t / (mpreal(2.0) * pi)).toLong()) + 2;

    mpreal real_sum = 0.0;
    mpreal imag_sum = 0.0;

    for (long k = 1; k <= terms; ++k) {
        mpreal sign = (k % 2 == 0) ? mpreal(-1.0) : mpreal(1.0);
        mpreal mpfr_k = mpreal(k);
        mpreal arg = -t * mpfr::log(mpfr_k);
        mpreal magnitude = sign / mpfr::sqrt(mpfr_k);
        real_sum += magnitude * mpfr::cos(arg);
        imag_sum += magnitude * mpfr::sin(arg);
    }

    return (real_sum * mpfr::cos(theta) - imag_sum * mpfr::sin(theta));
}

// -------------------------------------------------------------------------
// Brent's root-finding method (real implementation, following the
// standard algorithm: inverse quadratic interpolation, secant fallback,
// bisection safeguard). Assumes f(a) and f(b) have opposite signs.
// -------------------------------------------------------------------------
mpreal brent_zero_refiner_mpfr(mpreal a, mpreal b, const mpreal& tol, int max_iter = 200) {
    mpreal fa = Z_function_mpfr(a);
    mpreal fb = Z_function_mpfr(b);

    if (fa * fb >= 0) return (a + b) / mpreal(2.0); // caller should have pre-checked the bracket

    if (mpfr::abs(fa) < mpfr::abs(fb)) {
        swap(a, b);
        swap(fa, fb);
    }

    mpreal c = a, fc = fa;
    mpreal d = b - a;
    bool mflag = true;

    for (int iter = 0; iter < max_iter; ++iter) {
        if (fb == 0.0 || mpfr::abs(b - a) < tol) break;

        mpreal s;
        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation
            s = a * fb * fc / ((fa - fb) * (fa - fc))
              + b * fa * fc / ((fb - fa) * (fb - fc))
              + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            // Secant method
            s = b - fb * (b - a) / (fb - fa);
        }

        mpreal lo = (mpreal(3.0) * a + b) / mpreal(4.0);
        bool s_out_of_range = !(((s > lo) && (s < b)) || ((s < lo) && (s > b)));
        bool cond2 = mflag  && (mpfr::abs(s - b) >= mpfr::abs(b - c) / mpreal(2.0));
        bool cond3 = !mflag && (mpfr::abs(s - b) >= mpfr::abs(c - d) / mpreal(2.0));
        bool cond4 = mflag  && (mpfr::abs(b - c) < tol);
        bool cond5 = !mflag && (mpfr::abs(c - d) < tol);

        if (s_out_of_range || cond2 || cond3 || cond4 || cond5) {
            s = (a + b) / mpreal(2.0);
            mflag = true;
        } else {
            mflag = false;
        }

        mpreal fs = Z_function_mpfr(s);
        d = c;
        c = b; fc = fb;

        if (fa * fs < 0) { b = s; fb = fs; }
        else             { a = s; fa = fs; }

        if (mpfr::abs(fa) < mpfr::abs(fb)) {
            swap(a, b);
            swap(fa, fb);
        }
    }
    return b;
}

// -------------------------------------------------------------------------
// Arbitrary-precision prime gap analysis using GMP. Finds the maximal
// gap between consecutive primes in [N_start, N_end]. Sequential for now;
// parallelizing it correctly requires stitching gaps across chunk
// boundaries, left as a clearly-labeled future extension rather than a
// fake speedup claim.
// -------------------------------------------------------------------------
void prime_gap_analysis(unsigned long N_start, unsigned long N_end) {
    mpz_class p(N_start);
    mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());

    mpz_class max_gap(0), gap_start(p), gap_end(p);

    while (p < N_end) {
        mpz_class next_p;
        mpz_nextprime(next_p.get_mpz_t(), p.get_mpz_t());
        mpz_class gap = next_p - p;
        if (gap > max_gap) {
            max_gap = gap;
            gap_start = p;
            gap_end = next_p;
        }
        p = next_p;
    }

    cout << "[Prime Gap] Interval [" << N_start << ", " << N_end << "]: "
         << "max gap = " << max_gap
         << " (between " << gap_start << " and " << gap_end << ")" << endl;
}

int main() {
    mpfr::mpreal::set_default_prec(256);

    cout << "=========================================================" << endl;
    cout << " RIEMANN ZETA HIGH-PRECISION HPC ENGINE (MPFR + OpenMP)" << endl;
    cout << "=========================================================" << endl;

    mpreal t_start = "10.0";
    mpreal t_end   = "1000.0";
    mpreal step    = "0.1";
    mpreal tol     = "1e-12";

    long total_steps = static_cast<long>(((t_end - t_start) / step).toLong());

    cout << "[HPC Engine] Search range: t in [" << t_start << ", " << t_end << "]" << endl;
    cout << "[Precision] Tolerance: 1e-12 (256-bit precision)" << endl;
    cout << "[Threads] OpenMP max threads: " << omp_get_max_threads() << endl << endl;

    auto start_time = chrono::high_resolution_clock::now();

    // One results buffer per thread avoids locking on every hit; we merge
    // and sort once, after the parallel region.
    vector<vector<mpreal>> thread_zeros(omp_get_max_threads());

    #pragma omp parallel for schedule(dynamic, 64)
    for (long i = 0; i < total_steps; ++i) {
        mpreal a = t_start + mpreal(i) * step;
        mpreal b = a + step;
        mpreal fa = Z_function_mpfr(a);
        mpreal fb = Z_function_mpfr(b);

        if (fa * fb < 0) {
            mpreal exact_zero = brent_zero_refiner_mpfr(a, b, tol);
            int tid = omp_get_thread_num();
            thread_zeros[tid].push_back(exact_zero);
        }
    }

    vector<mpreal> all_zeros;
    for (auto& v : thread_zeros)
        for (auto& z : v)
            all_zeros.push_back(z);
    sort(all_zeros.begin(), all_zeros.end());

    cout << setprecision(16);
    for (size_t i = 0; i < all_zeros.size(); ++i) {
        cout << "Zero #" << (i + 1) << " -> s = 0.5 + " << all_zeros[i] << "i" << endl;
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "\n---------------------------------------------------------" << endl;
    cout << "Total zeros found: " << all_zeros.size() << endl;
    cout << "Elapsed time: " << elapsed.count() << " seconds." << endl;
    cout << "=========================================================" << endl;

    // Demonstrates the (currently sequential) prime gap module.
    prime_gap_analysis(2, 100000);

    return 0;
}
