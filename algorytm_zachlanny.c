/**
 * @file algorytm_zachlanny.c
 * @brief AKZ: algorytm zachlanny (wspinaczka, najlepszy ruch z calego sasiedztwa)
 *        od jednego losowego grafu startowego. Deterministyczny dla danego ziarna.
 * Uzycie: ./algorytm_zachlanny n k [seed=-1] [limit_s=60]
 */
#include "integral.h"

int main(int argc, char *argv[]){
    if (argc < 3){ fprintf(stderr, "Uzycie: %s n k [seed=-1] [limit_s=60]\n", argv[0]); return 1; }
    int n = atoi(argv[1]), k = atoi(argv[2]);
    long long seed = (argc > 3) ? atoll(argv[3]) : -1;
    double limit = (argc > 4) ? atof(argv[4]) : 60.0;
    if (n < 2 || n > MAXN_EXACT || k < n-1 || k > n*(n-1)/2 || limit <= 0){
        fprintf(stderr, "Blad: niepoprawne parametry (n=%d (max %d), k=%d, limit=%g)\n", n, MAXN_EXACT, k, limit);
        return 1;
    }
    seed = resolve_seed(seed);
    Rng rng = { (uint64_t)seed };
    fprintf(stderr, "AKZ: n=%d, k=%d, seed=%lld, limit=%.1fs\n", n, k, seed, limit);

    double t0 = now_s();
    static adj a;
    static SearchStats st;
    char g6[G6_BUF];
    generate_connected(&rng, n, k, a);
    g6_encode(n, a, g6);
    fprintf(stderr, "Graf startowy: %s\n", g6);
    double cost = local_search(n, a, NULL, t0 + limit, &st);
    st.restarts = 1;
    if (cost < COST_ZERO) record_found(&st, n, a, NULL);
    Verdict v;
    verify_graph(n, a, &v);
    g6_encode(n, a, g6);
    fprintf(stderr, "Koniec: koszt=%.0f, ruchow=%lld, przejrzano=%lld, znaleziono=%d\n", cost, st.steps, st.evals, st.nfound);
    print_solution_json("AKZ", n, k, seed, 1, &st, cost, now_s() - t0, g6, &v);
    return 0;
}