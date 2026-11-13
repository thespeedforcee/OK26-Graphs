/**
 * @file metaheurystyka_wyzarzanie.c
 * @brief SA: symulowane wyzarzanie (metaheurystyka). Stan = spojny graf G(n,k);
 *        ruch = zamiana jednej krawedzi na nie-krawedz (liczba krawedzi sie nie zmienia);
 *        ruchy rozlaczajace graf sa odrzucane. Koszt = n - liczba calkowitych wartosci
 *        wlasnych (sito, arytmetyka calkowita). Ruch o rownym lub mniejszym koszcie jest
 *        przyjmowany zawsze, gorszy z prawdopodobienstwem exp(-delta/T) (liczby
 *        zmiennoprzecinkowe: temperatura i losowanie). Temperatura maleje geometrycznie
 *        od T0 do Tmin w kazdym restarcie.
 *        Dziala do limitu czasu; raportuje ile przejrzano, ile znaleziono i jakie.
 * Uzycie: ./metaheurystyka_wyzarzanie n k [seed=-1] [limit_s=300] [plik_wynikow] [stop_first=0]
 *                                     [T0=1.0] [Tmin=0.15] [iter_na_restart=20000]
 */
#include "integral.h"

int main(int argc, char *argv[]){
    if (argc < 3){ fprintf(stderr, "Uzycie: %s n k [seed] [limit_s] [plik|-] [stop_first] [T0] [Tmin] [iter]\n", argv[0]); return 1; }
    int n = atoi(argv[1]), k = atoi(argv[2]);
    long long seed = (argc > 3) ? atoll(argv[3]) : -1;
    double limit = (argc > 4) ? atof(argv[4]) : 300.0;
    const char *path = (argc > 5 && strcmp(argv[5], "-") != 0) ? argv[5] : NULL;
    int stop_first = (argc > 6) ? atoi(argv[6]) : 0;
    double T0 = (argc > 7) ? atof(argv[7]) : 1.0;
    double Tmin = (argc > 8) ? atof(argv[8]) : 0.15;
    long iters = (argc > 9) ? atol(argv[9]) : 20000;
    if (n < 2 || n > MAXN_EXACT || k < n-1 || k > n*(n-1)/2 || limit <= 0 || T0 <= Tmin || Tmin <= 0 || iters < 1){
        fprintf(stderr, "Blad: niepoprawne parametry\n"); return 1;
    }
    FILE *out = path ? fopen(path, "a") : NULL;
    if (path && !out){ fprintf(stderr, "Blad: nie mozna otworzyc %s\n", path); return 1; }
    seed = resolve_seed(seed);
    Rng rng = { (uint64_t)seed };
    fprintf(stderr, "SA: n=%d, k=%d, seed=%lld, limit=%.1fs, T0=%.2f, Tmin=%.2f, iter/restart=%ld\n", n, k, seed, limit, T0, Tmin, iters);

    double t0 = now_s(), deadline = t0 + limit, last = t0;
    double alpha = pow(Tmin/T0, 1.0/(double)iters);       /* wspolczynnik chlodzenia */
    static adj a, best;
    static Pair edges[MAX_PAIRS], non_edges[MAX_PAIRS];
    static SearchStats st;
    st.best_cost = COST_INF;
    int stop = 0;

    while (!stop && now_s() < deadline){
        generate_connected(&rng, n, k, a);
        st.restarts++;
        int ne = 0, nn = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++){
            if (a[i][j]) edges[ne++] = (Pair){i, j}; else non_edges[nn++] = (Pair){i, j};
        }
        double cur = graph_cost(n, a), T = T0;
        st.evals++;
        for (long it = 0; it < iters && !stop; it++, T *= alpha){
            if (cur < st.best_cost){ st.best_cost = cur; memcpy(best, a, sizeof(adj)); }
            if (cur < COST_ZERO){
                record_found(&st, n, a, out);
                if (stop_first) stop = 1;
                break;                                   /* restart po znalezieniu */
            }
            if ((it & 255) == 0 && now_s() >= deadline){ stop = 1; break; }
            int i = rng_int(&rng, ne), j = rng_int(&rng, nn);
            swap_edges(a, edges[i], non_edges[j]);
            double c = graph_cost(n, a);
            st.evals++;
            double delta = c - cur;
            int accept = (c < COST_INF) && (delta <= 0 || exp(-delta/T) > (double)(rng_next(&rng) >> 11) / 9007199254740992.0);
            if (accept){
                Pair t = edges[i]; edges[i] = non_edges[j]; non_edges[j] = t;   /* role par zamienione */
                cur = c; st.steps++;
            } else {
                swap_edges(a, non_edges[j], edges[i]);                          /* cofniecie */
            }
        }
        if (now_s() - last >= 10.0){
            last = now_s();
            fprintf(stderr, "[t=%.0fs] przejrzano=%lld, restartow=%lld, znaleziono=%d, najlepszy koszt=%.0f\n",
                    last - t0, st.evals, st.restarts, st.nfound, st.best_cost);
        }
    }
    st.timed_out = now_s() >= deadline;

    Verdict v;
    verify_graph(n, best, &v);
    char g6[G6_BUF];
    g6_encode(n, best, g6);
    fprintf(stderr, "Koniec: przejrzano=%lld grafow, restartow=%lld, ZNALEZIONO=%d, najlepszy koszt=%.0f\n",
            st.evals, st.restarts, st.nfound, st.best_cost);
    print_solution_json("SA", n, k, seed, (int)st.restarts, &st, st.best_cost, now_s() - t0, g6, &v);
    if (out) fclose(out);
    return 0;
}