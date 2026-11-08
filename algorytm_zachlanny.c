#include "integral.h"
 
int main(int argc, char *argv[]) {
    if (argc<3) {
        fprintf(stderr, "Uzycie: %s n k [seed=-1] [limit_s=60] [tol=1e-6]\n", argv[0]);
        return 1;
    }
    int n=atoi(argv[1]);
    int k=atoi(argv[2]);
    long long seed=(argc>3) ? atoll(argv[3]):-1;
    double limit=(argc>4) ? atof(argv[4]) :60.0;
    double tol=(argc>5)?atof(argv[5]) : 1e-6;
    if (n<2 || n>MAXN || k<n-1 || k>n*(n-1)/2 || limit<= 0) {
        fprintf(stderr, "Blad: niepoprawne parametry (n=%d, k=%d, limit=%g)\n", n, k, limit);
        return 1;
    }
    seed=resolve_seed(seed);
    Rng rng={(uint64_t)seed};
    fprintf(stderr, "AKZ: n=%d, k=%d, seed=%lld, limit=%.1fs\n", n, k, seed, limit);
 
    double t0=now_s();
    static adj adj;
    char g6[G6_BUF];
    generate_connected(&rng, n, k, adj);
    g6_encode(n, adj, g6);
    fprintf(stderr, "Graf startowy: %s\n", g6);
    SearchStats st={0, 0, 0};
    double cost=local_search(n, adj, NULL, t0 + limit, tol, &st);
    Verdict v;
    verify_graph(n, adj, tol, &v);
    g6_encode(n, adj, g6);
    fprintf(stderr, "Koniec: koszt=%.3e, ruchow=%lld, ocen=%lld, calkowity=%d\n",cost, st.steps, st.evals, v.integral_graph);
    print_solution_json("AKZ", n, k, seed, 1, &st, cost, now_s() - t0, g6, &v);
    return 0;
}