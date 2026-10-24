// Regression check of the coupled (mixed singlet/triplet) solver against the
// reference tables in data/rashba/tc (Vs=0.5, Vt=2, nl=1.5, t'=0).
// Each sampled reference row (gamma0, Tc, mu, delta-ratio) seeds the bracket
// [0.5x, 2x] exactly as coupled_tabulate1D does with the previous point.
// Usage: regress_rashba <file.dat> <s|t> [stride=5] [Vs=0.5] [Vt=2] [nl=1.5]
//   NCS_BRACKET=f  solve in [x/f, x*f] instead of coupled_get_res's [x/2, 2x]
//   NCS_VERBOSE=1  keep the solver's own messages (e.g. "zbr err: Same signs")
#define NCS_IMPLEMENTATION
#include "../../NCS.h"

#include <cstdio>
#include <cstdlib>
#include <array>

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: %s file.dat s|t [stride] [Vs] [Vt] [nl]\n", argv[0]); return 2; }
    const bool singlet = argv[2][0] == 's';
    const int stride = argc > 3 ? std::atoi(argv[3]) : 5;
    const double Vs = argc > 4 ? std::atof(argv[4]) : 0.5, Vt = argc > 5 ? std::atof(argv[5]) : 2.0;
    const double nl = argc > 6 ? std::atof(argv[6]) : 1.5, tp = 0.0;
    const double bracket = std::getenv("NCS_BRACKET") ? std::atof(std::getenv("NCS_BRACKET")) : 0.0;

    std::ifstream in(argv[1]);
    std::vector<std::array<double, 4>> rows;
    for (std::array<double, 4> r; in >> r[0] >> r[1] >> r[2] >> r[3];) rows.push_back(r);
    std::printf("# %zu reference rows, branch %s, every %d-th\n", rows.size(), singlet ? "singlet" : "triplet", stride);
    std::printf("# %-10s %-14s %-14s %-9s %-13s %-13s %-12s %-12s\n",
                "gamma0", "Tc_ref", "Tc_now", "rel_dTc", "mu_ref", "mu_now", "ratio_ref", "ratio_now");
    double worst = 0.0;
    for (size_t i = 0; i < rows.size(); i += stride) {
        auto& r = rows[i];
        Result_Pair prev = { true, r[1], r[2] };
        std::streambuf* old = std::getenv("NCS_VERBOSE") ? std::cout.rdbuf() : std::cout.rdbuf(nullptr);
        Result_Triple now;
        if (bracket > 0.0) {  // isolate this branch with a tighter bracket [x/f, x*f]
            Result_Pair p = coupled_solve1D_zbr(coupled_gap, r[1] / bracket, r[1] * bracket, 1e-15, r[0], Vs, Vt,
                                                r[2] / bracket, r[2] * bracket, nl, tp);
            now.T = p.some ? p.T : 0.0; now.mi = p.some ? p.mi : 0.0;
            now.delta = singlet ? coupled_deltas(now.T, now.mi, r[0], Vs, tp) : coupled_deltat(now.T, now.mi, r[0], Vt, tp);
        } else {
            now = coupled_get_res(Vs, Vt, nl, tp, r[0], 0, 0, 0, 0, singlet, prev);  // brackets [x/2, 2x]
        }
        std::cout.rdbuf(old);
        double rel = std::fabs(now.T - r[1]) / r[1];
        worst = std::max(worst, rel);
        std::printf("  %-10.8g %-14.9g %-14.9g %-9.1e %-13.9g %-13.9g %-12.6g %-12.6g\n",
                    r[0], r[1], now.T, rel, r[2], now.mi, r[3], now.delta);
        std::fflush(stdout);
    }
    std::printf("# worst relative Tc deviation: %.2e\n", worst);
    return 0;
}
