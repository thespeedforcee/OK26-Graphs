/**
 * @file rozszerz.c
 * @brief Konstrukcja grafow z mniejszych/innych grafow: dopisywanie i przepinanie wierzcholka.
 *
 * Tryb "dopisz": z grafu rzedu m robi rzad m+1 dodajac jeden wierzcholek polaczony
 *   z d wybranymi wierzcholkami, gdzie d = k_docelowe - k(G). Sprawdza WSZYSTKIE
 *   podzbiory d-elementowe (kod Gospera). Np. z grafow n=17 do n=18.
 * Tryb "przepnij": dla kazdego wierzcholka v usuwa go i dolacza z powrotem na wszystkie
 *   mozliwe sposoby tak, by liczba krawedzi = k_docelowe (nowe, duze sasiedztwo
 *   w porownaniu z zamiana jednej krawedzi).
 *
 * Uzycie: ... | ./rozszerz dopisz|przepnij k_docelowe [limit_s=0 (bez limitu)]
 * Wejscie: graph6, jeden graf na linie. Wyjscie: JSON (jeden graf calkowity na linie).
 *
 */
#include "integral.h"

static long long evals = 0;

/** @brief Probuje wszystkich podzbiorow rozmiaru d zbioru cand[0..m-1] jako sasiadow wierzcholka w. */
static int try_attach(int n, adj a, int w, const int *cand, int m, int d, double deadline, const char *src){
    if (d < 1 || d > m) return 0;
    int found = 0;
    uint32_t s = (1u << d) - 1, limit = 1u << m;
    while (s < limit){
        for (int i = 0; i < m; i++) if (s >> i & 1) a[w][cand[i]] = a[cand[i]][w] = 1;
        evals++;
        if (graph_cost(n, a) < COST_ZERO){
            Verdict v; char g6[G6_BUF], esc[2*G6_BUF];
            verify_graph(n, a, &v);
            g6_encode(n, a, g6); json_escape(g6, esc);
            printf("{\"graph6\":\"%s\",\"n\":%d,\"k\":%d,\"zrodlo\":\"%s\",\"spectrum\":", esc, n, v.k, src);
            print_spectrum_json(stdout, &v);
            printf("}\n"); fflush(stdout);
            found++;
        }
        for (int i = 0; i < m; i++) a[w][cand[i]] = a[cand[i]][w] = 0;
        uint32_t c = s & -s, r = s + c;                /* nastepny podzbior o tej samej liczbie bitow */
        s = (((r ^ s) >> 2) / c) | r;
        if (deadline > 0 && now_s() > deadline) break;
    }
    return found;
}

int main(int argc, char *argv[]){
    if (argc < 3){ fprintf(stderr, "Uzycie: %s dopisz|przepnij k_docelowe [limit_s=0]\n", argv[0]); return 1; }
    int dopisz = strcmp(argv[1], "dopisz") == 0;
    int kt = atoi(argv[2]);
    double limit = (argc > 3) ? atof(argv[3]) : 0;
    double deadline = limit > 0 ? now_s() + limit : 0;
    static adj a;
    char line[G6_BUF];
    long graphs = 0, found = 0;
    while (fgets(line, sizeof line, stdin)){
        line[strcspn(line, "\r\n")] = 0;
        if (!*line || *line == '>') continue;
        int n = g6_decode(line, a);
        if (n < 1) { fprintf(stderr, "pomijam: %s\n", line); continue; }
        graphs++;
        int k = count_edges(n, a);
        if (dopisz){
            int nn = n + 1;
            if (nn > MAXN_EXACT){ fprintf(stderr, "n+1 > %d\n", MAXN_EXACT); continue; }
            int cand[MAXN];
            for (int i = 0; i < n; i++) cand[i] = i;
            for (int i = 0; i < nn; i++) a[n][i] = a[i][n] = 0;
            found += try_attach(nn, a, n, cand, n, kt - k, deadline, line);
        } else {
            for (int v = 0; v < n; v++){
                int dv = 0, cand[MAXN], m = 0;
                int saved[MAXN];
                for (int i = 0; i < n; i++){ saved[i] = a[v][i]; dv += a[v][i]; a[v][i] = a[i][v] = 0; }
                for (int i = 0; i < n; i++) if (i != v) cand[m++] = i;
                found += try_attach(n, a, v, cand, m, kt - (k - dv), deadline, line);
                for (int i = 0; i < n; i++) a[v][i] = a[i][v] = saved[i];
            }
        }
        if (deadline > 0 && now_s() > deadline) break;
    }
    fprintf(stderr, "Grafow wejsciowych: %ld, ocen: %lld, znalezionych calkowitych (z powtorzeniami): %ld\n", graphs, evals, found);
    return 0;
}