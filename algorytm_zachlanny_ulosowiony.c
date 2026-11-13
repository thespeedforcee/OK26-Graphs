/**
 * @file algorytm_zachlanny_ulosowiony.c
 * @brief AKU: algorytm zachlanny ulosowiony = wielokrotne restarty z losowego grafu,
 *        wspinaczka z losowa kolejnoscia sasiadow (pierwszy ruch poprawiajacy).
 *        Dziala do limitu czasu (lub do pierwszego grafu, gdy stop_first=1).
 *        Na koncu raportuje: ile grafow przejrzano, ile znaleziono, jakie.
 * Uzycie: ./algorytm_zachlanny_ulosowiony n k [seed=-1] [limit_s=300] [plik_wynikow] [stop_first=0]
 *   plik_wynikow - znalezione grafy (graph6, jedna linia = jeden graf; dopisywane)
 */
#include "integral.h"

int main(int argc, char *argv[]){
    if (argc < 3){ fprintf(stderr, "Uzycie: %s n k [seed=-1] [limit_s=300] [plik_wynikow] [stop_first=0]\n", argv[0]); return 1; }
    int n = atoi(argv[1]), k = atoi(argv[2]);
    long long seed = (argc > 3) ? atoll(argv[3]) : -1;
    double limit = (argc > 4) ? atof(argv[4]) : 300.0;
    const char *path = (argc > 5 && strcmp(argv[5], "-") != 0) ? argv[5] : NULL;
    int stop_first = (argc > 6) ? atoi(argv[6]) : 0;
    if (n < 2 || n > MAXN_EXACT || k < n-1 || k > n*(n-1)/2 || limit <= 0){
        fprintf(stderr, "Blad: niepoprawne parametry (n=%d (max %d), k=%d, limit=%g)\n", n, MAXN_EXACT, k, limit);
        return 1;
    }
    FILE *out = path ? fopen(path, "a") : NULL;
    if (path && !out){ fprintf(stderr, "Blad: nie mozna otworzyc %s\n", path); return 1; }
    seed = resolve_seed(seed);
    Rng rng = { (uint64_t)seed };
    fprintf(stderr, "AKU: n=%d, k=%d, seed=%lld, limit=%.1fs\n", n, k, seed, limit);

    double t0 = now_s(), deadline = t0 + limit, last = t0;
    static adj cur, best;
    static SearchStats st;
    st.best_cost = COST_INF;
    do {
        generate_connected(&rng, n, k, cur);
        st.restarts++;
        double c = local_search(n, cur, &rng, deadline, &st);
        if (c < st.best_cost){
            st.best_cost = c;
            memcpy(best, cur, sizeof(adj));
        }
        if (c < COST_ZERO){
            record_found(&st, n, cur, out);
            if (stop_first) break;
        }
        if (now_s() - last >= 10.0){
            last = now_s();
            fprintf(stderr, "[t=%.0fs] przejrzano=%lld, restartow=%lld, znaleziono=%d, najlepszy koszt=%.0f\n",
                    last - t0, st.evals, st.restarts, st.nfound, st.best_cost);
        }
    } while (now_s() < deadline);
    st.timed_out = now_s() >= deadline;

    Verdict v;
    verify_graph(n, best, &v);
    char g6[G6_BUF];
    g6_encode(n, best, g6);
    fprintf(stderr, "Koniec: przejrzano=%lld grafow, restartow=%lld, ZNALEZIONO=%d, najlepszy koszt=%.0f\n",
            st.evals, st.restarts, st.nfound, st.best_cost);
    print_solution_json("AKU", n, k, seed, (int)st.restarts, &st, st.best_cost, now_s() - t0, g6, &v);
    if (out) fclose(out);
    return 0;
}