# High-Performance Riemann Zeta Zero Search & Prime Gap Engine (HPC)

<https://o66939499-dev.github.io/Riemann-Zeta-HPC-Engine/>

![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)
![OpenMP](https://img.shields.io/badge/Parallelism-OpenMP-orange.svg)
![GNU GMP](https://img.shields.io/badge/Precision-GNU%20GMP-green.svg)
![License](https://img.shields.io/badge/License-MIT-brightgreen.svg)

A high-performance computing (HPC) project in C++ for numerically evaluating
the **Riemann Zeta function** ζ(s), locating candidate non-trivial zeros on
the **critical line** Re(s) = 0.5, and analyzing **prime gaps** using
arbitrary-precision arithmetic.

---

## 🌟 Key Features

- **Parallel zero search:** The critical-line grid scan `t ∈ [t_start, t_end]`
  is distributed across threads with `#pragma omp parallel for` and dynamic
  scheduling; each thread bisects/refines its own bracket independently, and
  results are merged and sorted once the parallel region ends.
- **Adaptive-precision Z(t) evaluation:** The Riemann–Siegel `Z(t)` sum uses
  `~sqrt(t / 2π)` terms (the number the asymptotic formula actually needs)
  instead of a fixed term count for every `t`, cutting redundant work at low
  `t` and improving accuracy at high `t`.
- **Real Brent's method:** Root refinement combines inverse quadratic
  interpolation, the secant method, and a bisection safeguard — the standard
  algorithm, not a bisection loop mislabeled as Brent.
- **Arbitrary-precision prime gap analysis:** Uses GNU GMP (`mpz_class`,
  `mpz_nextprime`) to find the maximal gap between consecutive primes in a
  given range. Currently sequential — see *Known Limitations* below.

---

## 📐 Theoretical Background

### 1. Riemann–Siegel Z(t)

On the critical line, ζ(1/2 + it) is evaluated via the real-valued
Riemann–Siegel function:

$$Z(t) = e^{i\theta(t)} \zeta(1/2 + it)$$

where θ(t) is the Riemann–Siegel theta function. `Z(t)` is real, so its sign
changes mark candidate zeros of ζ on the critical line — this is what the
engine scans for.

### 2. The Riemann Hypothesis

The Riemann Hypothesis conjectures that every non-trivial zero of ζ(s) lies
on Re(s) = 0.5. This engine searches for `t` values where `Z(t) ≈ 0` and
refines each candidate to the requested tolerance with Brent's method.

---

## 🛠 Build & Run

### Prerequisites

```bash
sudo apt update
sudo apt install build-essential libgmp-dev libmpfr-dev libomp-dev
```

You'll also need [`mpreal.h`](https://github.com/advanpix/mpreal) (a
header-only MPFR C++ wrapper) placed alongside `main.cpp`.

### Compile

```bash
g++ -O3 -std=c++20 main.cpp -fopenmp -lgmpxx -lgmp -lmpfr -o zeta_engine
```

### Run

```bash
./zeta_engine
```

Set the thread count explicitly if you want to compare scaling:

```bash
OMP_NUM_THREADS=4 ./zeta_engine
```

---

## 📊 Benchmarking

To produce an honest speedup table, run the same range with
`OMP_NUM_THREADS` set to 1, 2, 4, 8... and record the wall-clock time
printed at the end of each run. Fill in your own machine's numbers —
posting figures without the code behind them (the previous version's
benchmark table) undermines a portfolio piece more than an absent table
would.

| Threads | Time (s) | Speedup | Efficiency |
|---------|----------|---------|------------|
| 1       |          | 1.00x   |            |
| 2       |          |         |            |
| 4       |          |         |            |
| 8       |          |         |            |

---

## ⚠️ Known Limitations

Being upfront about these matters more than hiding them — reviewers notice,
and it reads better to have named them yourself:

- The prime gap module is currently **sequential**. A correct parallel
  version needs to stitch gaps across chunk boundaries (a gap can straddle
  two threads' ranges); that's a natural next step, not yet implemented.
- The web page's "real-time" curve is a static illustration of expected
  output, not a live run of the C++ engine in the browser. A true live demo
  would need a WebAssembly build of the engine or a backend that runs it.
- `Z_function_mpfr` uses only the main Riemann–Siegel sum, not the full
  asymptotic correction series (C₀, C₁, ...). This is accurate enough for
  the tested range but loses precision faster than the full formula as `t`
  grows very large.

---

## 📄 Academic Portfolio Context

Developed as part of an academic portfolio for undergraduate admissions in
Applied Mathematics and Cybersecurity.
