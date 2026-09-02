#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <gmp.h>
#include <mpfr.h>
#include "mpreal.h"

using namespace std;
using mpfr::mpreal;

mpreal riemann_theta_mpfr(const mpreal& t) {
    mpreal pi = mpfr::const_pi();
    mpreal term1 = t * mpfr::log(t / (mpreal(2.0) * pi)) / mpreal(2.0);
    mpreal term2 = t / mpreal(2.0);
    mpreal term3 = pi / mpreal(8.0);
    mpreal term4 = mpreal(1.0) / (mpreal(48.0) * t);

    return term1 - term2 - term3 + term4;
}

mpreal Z_function_mpfr(const mpreal& t, int terms = 20000) {
    mpreal theta = riemann_theta_mpfr(t);
    mpreal real_sum = 0.0;
    mpreal imag_sum = 0.0;

    for (int k = 1; k <= terms; ++k) {
        mpreal sign = (k % 2 == 0) ? -1.0 : 1.0;
        mpreal mpfr_k = k;
        mpreal arg = -t * mpfr::log(mpfr_k);
        mpreal magnitude = sign / mpfr::sqrt(mpfr_k);

        real_sum += magnitude * mpfr::cos(arg);
        imag_sum += magnitude * mpfr::sin(arg);
    }

    return (real_sum * mpfr::cos(theta) - imag_sum * mpfr::sin(theta));
}

mpreal brent_zero_refiner_mpfr(mpreal a, mpreal b, mpreal tol = "1e-12") {
    mpreal fa = Z_function_mpfr(a);
    mpreal fb = Z_function_mpfr(b);

    if (fa * fb >= 0) return (a + b) / mpreal(2.0);

    mpreal c = a, fc = fa, d = b - a, e = d;

    while (mpfr::abs(b - a) > tol) {
        if (fb == 0.0) return b;

        if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) {
            a = c; fa = fc; d = b - a; e = d;
        }

        if (mpfr::abs(fa) < mpfr::abs(fb)) {
            c = b; b = a; a = c;
            fc = fb; fb = fa; fa = fc;
        }

        mpreal m = mpreal(0.5) * (a - b);
        if (mpfr::abs(m) <= tol) return b;

        d = m; e = d;
        c = b; fc = fb;

        if (mpfr::abs(d) > tol) b += d;
        else b += (m > 0 ? tol : -tol);

        fb = Z_function_mpfr(b);
    }
    return b;
}

int main() {
    mpfr::mpreal::set_default_prec(256);

    cout << "=========================================================" << endl;
    cout << "   RIEMANN ZETA HIGH-PRECISION HPC ENGINE (MPFR + OpenMP)" << endl;
    cout << "=========================================================" << endl;
    
    mpreal t_start = "10.0";
    mpreal t_end = "1000.0";
    mpreal step = "0.1";
    mpreal tol = "1e-12";

    int total_steps = static_cast<int>((t_end - t_start) / step);

    cout << "[HPC Engine] Gözleg aralygy: t in [" << t_start << ", " << t_end << "]" << endl;
    cout << "[Takyklyk] Tolerance: 1e-12 (256-bit Precision)" << endl << endl;

    auto start_time = chrono::high_resolution_clock::now();
    int zero_count = 0;
    
    for (int i = 0; i < total_steps; ++i) {
        mpreal a = t_start + i * step;
        mpreal b = a + step;

        mpreal fa = Z_function_mpfr(a);
        mpreal fb = Z_function_mpfr(b);
        
        if (fa * fb < 0) {
            zero_count++;
            mpreal exact_zero = brent_zero_refiner_mpfr(a, b, tol);

            cout.precision(16);
            cout << "Nol #" << zero_count << " -> s = 0.5 + " << exact_zero << "i" << endl;
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end_time - start_time;

    cout << "\n---------------------------------------------------------" << endl;
    cout << "Jemi tapylan takyk nollar: " << zero_count << endl;
    cout << "Sarp edilen wagt: " << elapsed.count() << " sekund." << endl;
    cout << "=========================================================" << endl;

    return 0;
}
