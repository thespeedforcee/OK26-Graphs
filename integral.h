/**
 * @file integral.h
 * @brief Wspolna biblioteka: generowanie grafow, kodowanie graph6, dokladne sito
 *        grafow calkowitych (arytmetyka calkowita), przeszukiwanie lokalne.
 *
 * Graf jest calkowity (integralny), gdy wszystkie wartosci wlasne jego macierzy
 * sasiedztwa sa liczbami calkowitymi.
 *
 * SITO (zastepuje dawna metode Jacobiego):
 *  1. wielomian charakterystyczny det(xI - A) liczony DOKLADNIE (Faddeev-Leverrier,
 *     __int128, bez zaokraglen),
 *  2. z wielomianu wyrzucane sa pierwiastki calkowite z [-Delta, Delta]
 *     (|lambda| <= max stopien) dzieleniem Hornera, z krotnosciami.
 *  Koszt grafu = n - (liczba pierwiastkow calkowitych). Koszt 0 <=> graf calkowity.
 *
 * Ograniczenie: n <= MAXN_EXACT (zakres bezpieczny dla __int128).
 *
 */
#ifndef INTEGRAL_H
#define INTEGRAL_H
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 62               /**< maks. rzad grafu w formacie graph6 (wersja 1-bajtowa) */
#define MAXN_EXACT 25         /**< maks. rzad, dla ktorego sito jest uruchamiane */
#define MAX_PAIRS (MAXN*(MAXN-1)/2) /**< maks. liczba par wierzcholkow */
#define G6_BUF 330            /**< rozmiar bufora na lancuch graph6 */
#define COST_INF 1e300        /**< koszt grafu niespojnego */
#define COST_EPS 1e-12        /**< minimalna poprawa kosztu uznawana za zysk */
#define COST_ZERO 0.5         /**< koszt ponizej tej wartosci = graf calkowity (koszt jest calkowity) */

typedef __int128 i128;        /**< 128-bitowa liczba calkowita (GCC) */
typedef int adj[MAXN][MAXN];  /**< macierz sasiedztwa */

/** Para wierzcholkow (krawedz lub nie-krawedz). */
typedef struct {
    int u, v;
} Pair;

/** @brief Czas monotoniczny w sekundach. */
static inline double now_s(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9*(double)ts.tv_nsec;
}

/** Stan generatora pseudolosowego (splitmix64), przekazywany jawnie. */
typedef struct {
    uint64_t s;
} Rng;

/** @brief Nastepna liczba pseudolosowa 64-bitowa. */
static inline uint64_t rng_next(Rng *r){
    uint64_t z = (r->s += 0x9E3779B97F4A7C15ULL);
    z = (z^(z>>30))*0xBF58476D1CE4E5B9ULL;
    z = (z^(z>>27))*0x94D049BB133111EBULL;
    return z^(z>>31);
}

/** @brief Liczba losowa z zakresu [0, m). */
static inline int rng_int(Rng *r, int m){
    return (int)(rng_next(r) % (uint64_t)m);
}

/** @brief Ziarno: wartosc ujemna oznacza ziarno z zegara (ustalane RAZ w main). */
static inline long long resolve_seed(long long seed){
    if (seed < 0) seed = (long long)time(NULL);
    return seed;
}

/** @brief Losowa permutacja tablicy par (Fisher-Yates). */
static inline void shuffle_pairs(Rng *r, Pair *pairs, int num_pairs){
    for (int i = num_pairs-1; i > 0; i--){
        int j = rng_int(r, i+1);
        Pair t = pairs[i]; pairs[i] = pairs[j]; pairs[j] = t;
    }
}

/** @brief Liczba krawedzi grafu. */
static inline int count_edges(int n, const adj a){
    int k = 0;
    for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) k += a[i][j];
    return k;
}

/** @brief Maksymalny stopien wierzcholka. */
static inline int max_degree(int n, const adj a){
    int maxd = 0;
    for (int i = 0; i < n; i++){
        int d = 0;
        for (int j = 0; j < n; j++) d += a[i][j];
        if (d > maxd) maxd = d;
    }
    return maxd;
}

/** @brief Test spojnosci (przeszukiwanie wszerz). @return 1 gdy graf spojny. */
static int is_connected(int n, adj adj){
    int visited[MAXN] = {0}, queue[MAXN];
    int head = 0, tail = 0, cnt = 1;
    visited[0] = 1; queue[tail++] = 0;
    while (head < tail){
        int u = queue[head++];
        for (int v = 0; v < n; v++){
            if (adj[u][v] && !visited[v]){ visited[v] = 1; queue[tail++] = v; cnt++; }
        }
    }
    return cnt == n;
}

/**
 * @brief Dekoduje graph6 do macierzy sasiedztwa.
 * @return rzad grafu n, lub -1 gdy lancuch jest niepoprawny.
 */
static inline int g6_decode(const char *g6, adj adj){
    int n = (unsigned char)g6[0] - 63;
    if (n < 1 || n > MAXN) return -1;
    int need = (n*(n-1)/2 + 5)/6;
    if ((int)strlen(g6) < 1+need) return -1;
    for (int i = 1; i <= need; i++){
        unsigned char c = (unsigned char)g6[i];
        if (c < 63 || c > 126) return -1;
    }
    memset(adj, 0, sizeof(int)*MAXN*MAXN);
    int pos = 1, bit = 5;
    for (int j = 1; j < n; j++){
        for (int i = 0; i < j; i++){
            int v = (((unsigned char)g6[pos]-63) >> bit) & 1;
            adj[i][j] = adj[j][i] = v;
            if (--bit < 0){ bit = 5; pos++; }
        }
    }
    return n;
}

/** @brief Koduje graf do graph6 (bufor min. G6_BUF znakow). */
static inline void g6_encode(int n, adj adj, char *buf){
    int pos = 0, bits = 0, val = 0;
    buf[pos++] = (char)(n+63);
    for (int j = 1; j < n; j++){
        for (int i = 0; i < j; i++){
            val = (val<<1) | adj[i][j];
            if (++bits == 6){ buf[pos++] = (char)(val+63); val = 0; bits = 0; }
        }
    }
    if (bits > 0) buf[pos++] = (char)((val << (6-bits)) + 63);
    buf[pos] = '\0';
}

/**
 * @brief Kopiuje lancuch do bufora z ucieczka znakow \ i " (graph6 moze je zawierac).
 *        Wymagane do poprawnego JSON. Bufor out: min. 2*strlen(s)+1.
 */
static inline void json_escape(const char *s, char *out){
    for (; *s; s++){ if (*s == '\\' || *s == '"') *out++ = '\\'; *out++ = *s; }
    *out = '\0';
}

/**
 * @brief Generator: losowy spojny graf G(n,k). Losowe drzewo + losowanie bez zwracania
 *        brakujacych krawedzi. Zlozonosc O(n^2).
 */
static inline void generate_connected(Rng *r, int n, int k, adj adj){
    memset(adj, 0, sizeof(int)*MAXN*MAXN);
    int perm[MAXN];
    for (int i = 0; i < n; i++) perm[i] = i;
    for (int i = n-1; i > 0; i--){
        int j = rng_int(r, i+1);
        int t = perm[i]; perm[i] = perm[j]; perm[j] = t;
    }
    for (int i = 1; i < n; i++){
        int p = perm[rng_int(r, i)];
        adj[p][perm[i]] = adj[perm[i]][p] = 1;
    }
    static Pair miss[MAX_PAIRS];
    int m = 0;
    for (int u = 0; u < n; u++) for (int v = u+1; v < n; v++) if (!adj[u][v]) miss[m++] = (Pair){u, v};
    for (int i = 0; i < k-(n-1); i++){
        int j = i + rng_int(r, m-i);
        Pair e = miss[j]; miss[j] = miss[i];
        adj[e.u][e.v] = adj[e.v][e.u] = 1;
    }
}

/*  SITO (arytmetyka calkowita)  */

/**
 * @brief Wspolczynniki wielomianu charakterystycznego det(xI - A) (Faddeev-Leverrier).
 *        c[i] to wspolczynnik przy x^i, c[n] = 1. Dzielenie przez k jest dokladne.
 *        Zlozonosc O(n^4).
 * @param n rzad grafu (n <= MAXN_EXACT)
 * @param c (wyjscie) tablica n+1 wspolczynnikow
 */
static inline void char_poly(int n, adj adj, i128 *c){
    int nb[MAXN_EXACT][MAXN_EXACT], deg[MAXN_EXACT];
    i128 m[MAXN_EXACT][MAXN_EXACT], t[MAXN_EXACT][MAXN_EXACT];
    for (int i = 0; i < n; i++){                       /* listy sasiadow: A*M kosztuje 2k*n zamiast n^3 */
        deg[i] = 0;
        for (int j = 0; j < n; j++) if (adj[i][j]) nb[i][deg[i]++] = j;
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) m[i][j] = (i == j);   /* M_1 = I */
    c[n] = 1;
    for (int k = 1; k <= n; k++){
        i128 tr = 0;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                i128 sum = 0;
                for (int q = 0; q < deg[i]; q++) sum += m[nb[i][q]][j];
                t[i][j] = sum;
            }
            tr += t[i][i];
        }
        c[n-k] = -tr/k;                                /* dzielenie dokladne */
        if (k < n)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) m[i][j] = t[i][j] + (i == j ? c[n-k] : 0);
    }
}

/**
 * @brief SITO pierwiastkow calkowitych: dzieli wielomian przez (x - r) dla r z
 *        [-maxdeg, maxdeg] (schemat Hornera), z uwzglednieniem krotnosci.
 * @param spec (wyjscie) znalezione pierwiastki calkowite (rosnaco, z krotnosciami)
 * @return liczba znalezionych pierwiastkow; wynik == n <=> graf calkowity
 */
static inline int integer_roots(int n, const i128 *c, int maxdeg, int *spec){
    i128 p[MAXN_EXACT+1];
    for (int i = 0; i <= n; i++) p[i] = c[i];
    int deg = n, found = 0;
    for (int r = -maxdeg; r <= maxdeg; r++){
        while (deg > 0){
            i128 q[MAXN_EXACT+1], rem = p[deg];
            for (int i = deg-1; i >= 0; i--){ q[i] = rem; rem = p[i] + rem*r; }
            if (rem != 0) break;
            for (int i = 0; i < deg; i++) p[i] = q[i];
            deg--;
            spec[found++] = r;
        }
    }
    return found;
}

/**
 * @brief Funkcja celu minimalizowana przez heurystyki: n - liczba calkowitych
 *        wartosci wlasnych; COST_INF dla grafu niespojnego. Czysto calkowitoliczbowa.
 */
static inline double graph_cost(int n, adj adj){
    if (!is_connected(n, adj)) return COST_INF;
    i128 c[MAXN_EXACT+1];
    int spec[MAXN_EXACT+1];
    char_poly(n, adj, c);
    return (double)(n - integer_roots(n, c, max_degree(n, adj), spec));
}

/*  Weryfikacja wyniku  */

/** Wynik pelnej weryfikacji grafu. */
typedef struct {
    int n, k;                /**< rzad i liczba krawedzi */
    int connected;           /**< czy spojny */
    int cost;                /**< n - liczba pierwiastkow calkowitych */
    int integral_graph;      /**< spojny i calkowity */
    int spec[MAXN_EXACT];    /**< pierwiastki calkowite (cale widmo, gdy graf calkowity) */
    int nspec;               /**< ile pierwiastkow calkowitych znaleziono */
} Verdict;

/** @brief Pelna weryfikacja grafu (n <= MAXN_EXACT). */
static inline void verify_graph(int n, adj adj, Verdict *v){
    i128 c[MAXN_EXACT+1];
    v->n = n;
    v->k = count_edges(n, adj);
    v->connected = is_connected(n, adj);
    char_poly(n, adj, c);
    v->nspec = integer_roots(n, c, max_degree(n, adj), v->spec);
    v->cost = n - v->nspec;
    v->integral_graph = v->connected && v->cost == 0;
}

/** @brief Wypisuje widmo jako tablice JSON (lub null, gdy graf nie jest calkowity). */
static inline void print_spectrum_json(FILE *f, const Verdict *v){
    if (v->cost != 0){ fputs("null", f); return; }
    fputc('[', f);
    for (int i = 0; i < v->n; i++) fprintf(f, "%s%d", i ? "," : "", v->spec[i]);
    fputc(']', f);
}

/*  Przeszukiwanie lokalne  */

#define MAX_FOUND 1000   /**< ile znalezionych grafow pamietamy na liscie (reszta tylko liczona) */

/** Statystyki i wyniki przeszukiwania. */
typedef struct {
    long long evals;       /**< liczba obliczen funkcji celu = liczba "przejrzanych" grafow */
    long long steps;       /**< liczba przyjetych ruchow */
    long long restarts;    /**< liczba startow od losowego grafu */
    int timed_out;         /**< czy przerwano limitem czasu */
    int nfound;            /**< liczba roznych (jako lancuchy graph6) znalezionych grafow calkowitych */
    double best_cost;      /**< najlepszy koszt w calym przebiegu */
    char best_g6[G6_BUF];  /**< graf o najlepszym koszcie */
    char found[MAX_FOUND][G6_BUF]; /**< znalezione grafy calkowite (graph6) */
} SearchStats;

/**
 * @brief Zapisuje znaleziony graf calkowity: lista w pamieci, stderr i opcjonalnie plik
 *        (jeden graph6 na linie, gotowe do `labelg`). Pomija identyczne lancuchy graph6.
 * @return 1 gdy graf jest nowy, 0 gdy juz byl
 */
static inline int record_found(SearchStats *st, int n, adj a, FILE *out){
    char g6[G6_BUF];
    g6_encode(n, a, g6);
    int lim = st->nfound < MAX_FOUND ? st->nfound : MAX_FOUND;
    for (int i = 0; i < lim; i++) if (strcmp(st->found[i], g6) == 0) return 0;
    if (st->nfound < MAX_FOUND) strcpy(st->found[st->nfound], g6);
    st->nfound++;
    fprintf(stderr, "ZNALEZIONO #%d: %s\n", st->nfound, g6);
    if (out){ fprintf(out, "%s\n", g6); fflush(out); }
    return 1;
}

/** @brief Ruch: usuwa krawedz rem i dodaje krawedz add. */
static inline void swap_edges(adj adj, Pair rem, Pair add){
    adj[rem.u][rem.v] = adj[rem.v][rem.u] = 0;
    adj[add.u][add.v] = adj[add.v][add.u] = 1;
}

/**
 * @brief Wspinaczka: sasiedztwo = zamiana jednej krawedzi na nie-krawedz.
 * @param rng  NULL => wersja deterministyczna AKZ (najlepszy ruch z calego sasiedztwa);
 *             niezerowy => AKU (losowa kolejnosc, pierwszy ruch poprawiajacy)
 * @param deadline czas (now_s) po ktorym przerywa
 * @return koszt koncowy; graf zostaje w adj
 */
static inline double local_search(int n, adj adj, Rng *rng, double deadline, SearchStats *st){
    static Pair edges[MAX_PAIRS], non_edges[MAX_PAIRS];
    double cur = graph_cost(n, adj);
    st->evals++;
    while (cur >= COST_ZERO){
        int ne = 0, nn = 0;
        for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++){
            if (adj[i][j]) edges[ne++] = (Pair){i, j}; else non_edges[nn++] = (Pair){i, j};
        }
        if (rng){ shuffle_pairs(rng, edges, ne); shuffle_pairs(rng, non_edges, nn); }
        double best = cur;
        int bi = -1, bj = -1, improved = 0, timeout = 0;
        for (int a = 0; a < ne && !improved && !timeout; a++){
            for (int b = 0; b < nn; b++){
                if (now_s() >= deadline){ timeout = 1; break; }
                swap_edges(adj, edges[a], non_edges[b]);
                double c = graph_cost(n, adj);
                st->evals++;
                if (c < best - COST_EPS){
                    best = c; bi = a; bj = b;
                    if (rng){ improved = 1; break; }   /* zostawiamy ruch */
                }
                swap_edges(adj, non_edges[b], edges[a]);   /* cofniecie */
            }
        }
        if (rng){
            if (improved){ cur = best; st->steps++; }
        } else if (bi >= 0){
            swap_edges(adj, edges[bi], non_edges[bj]);
            cur = best; improved = 1; st->steps++;
        }
        if (timeout){ st->timed_out = 1; break; }
        if (!improved) break;
    }
    return cur;
}

/** @brief Wypisuje wynik dzialania algorytmu jako JSON na stdout. */
static inline void print_solution_json(const char *alg, int n, int k, long long seed, int restarts,
                                       const SearchStats *st, double cost, double time_s,
                                       const char *g6, const Verdict *v){
    char esc[2*G6_BUF];
    json_escape(g6, esc);
    printf("{\n");
    printf("  \"algorithm\": \"%s\",\n", alg);
    printf("  \"n\": %d,\n  \"k\": %d,\n  \"seed\": %lld,\n", n, k, seed);
    printf("  \"restarts\": %d,\n", restarts);
    printf("  \"evaluations\": %lld,\n  \"steps\": %lld,\n", st->evals, st->steps);
    printf("  \"timed_out\": %d,\n", st->timed_out);
    printf("  \"time_s\": %.3f,\n", time_s);
    printf("  \"cost\": %.0f,\n", cost);
    printf("  \"found_count\": %d,\n  \"found\": [", st->nfound);
    for (int i = 0; i < st->nfound && i < MAX_FOUND; i++){
        char e2[2*G6_BUF]; json_escape(st->found[i], e2);
        printf("%s\"%s\"", i ? "," : "", e2);
    }
    printf("],\n");
    printf("  \"graph6\": \"%s\",\n", esc);
    printf("  \"connected\": %d,\n", v->connected);
    printf("  \"integral_graph\": %d,\n", v->integral_graph);
    printf("  \"spectrum\": ");
    print_spectrum_json(stdout, v);
    printf("\n}\n");
}

#endif