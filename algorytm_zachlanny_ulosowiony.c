#include "integral.h"
 
int main(int argc, char *argv[]) {
    if (argc<3) {
        fprintf(stderr, "Uzycie: %s n k [seed=-1] [limit_s=60] [tol=1e-6]\n", argv[0]);
        return 1;
    }
    int n=atoi(argv[1]);
    int k=atoi(argv[2]);
    long long seed=(argc>3) ? atoll(argv[3]) : -1;
    double limit=(argc>4) ? atof(argv[4]) : 60.0;
    double tol = (argc > 5) ? atof(argv[5]) : 1e-6;
    if (n < 2 || n > MAXN || k < n - 1 || k > n * (n - 1) / 2 || limit <= 0) {
        fprintf(stderr, "Blad: niepoprawne parametry (n=%d, k=%d, limit=%g)\n", n, k, limit);
        return 1;
    }
    seed=resolve_seed(seed);
    Rng rng={(uint64_t)seed};
    fprintf(stderr, "AKU: n=%d, k=%d, seed=%lld, limit=%.1fs\n", n, k, seed, limit);
 
    double t0=now_s(), deadline = t0 + limit;
    static adj cur, best;
    double best_cost=COST_INF;
    int restarts=0;
    SearchStats st={0, 0, 0};
 
    do {
        generate_connected(&rng, n, k, cur);
        restarts++;
        double c = local_search(n, cur, &rng, deadline, tol, &st);
        if (c < best_cost) {
            best_cost = c;
            memcpy(best, cur, sizeof(adj));
            fprintf(stderr, "[restart %d, t=%.1fs] najlepszy koszt: %.3e\n",restarts, now_s() - t0, best_cost);
        }
    } while (best_cost >= tol && now_s() < deadline);
 
    Verdict v;
    verify_graph(n, best, tol, &v);
    char g6[G6_BUF];
    g6_encode(n, best, g6);
    fprintf(stderr, "Koniec: restartow=%d, koszt=%.3e, calkowity=%d\n",restarts, best_cost, v.integral_graph);
    print_solution_json("AKU", n, k, seed, restarts, &st, best_cost, now_s() - t0, g6, &v);
    return 0;
}