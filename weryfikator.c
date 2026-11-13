/**
 * @file weryfikator.c
 * @brief Filtr (sito) grafow calkowitych. Czyta graph6 ze stdin (jeden graf na linie,
 *        np. wyjscie `geng`, lub dowolny JSON z polami "graph6") i dla kazdego
 *        wypisuje JSON z wynikiem dokladnej weryfikacji.
 * Uzycie:  geng -c 12 20:20 | ./weryfikator [-c|-g]
 *   -c  wypisuj tylko grafy calkowite (JSON)
 *   -g  wypisuj tylko grafy calkowite, jako czysty graph6 (jedna linia = graf; dla labelg)
 */
#include "integral.h"

static long total = 0, integral_cnt = 0;

/** @brief Weryfikuje jeden graf i (opcjonalnie) wypisuje wynik. */
static void process(const char *g6, int only_integral, int plain){
    static adj a;
    int n = g6_decode(g6, a);
    if (n < 0 || n > MAXN_EXACT){ fprintf(stderr, "pomijam (zly format lub n>%d): %s\n", MAXN_EXACT, g6); return; }
    Verdict v;
    verify_graph(n, a, &v);
    total++; integral_cnt += v.integral_graph;
    if (only_integral && !v.integral_graph) return;
    if (plain){ printf("%s\n", g6); fflush(stdout); return; }
    char esc[2*G6_BUF];
    json_escape(g6, esc);
    printf("{\"graph6\":\"%s\",\"n\":%d,\"k\":%d,\"connected\":%d,\"cost\":%d,\"integral_graph\":%d,\"spectrum\":",
           esc, v.n, v.k, v.connected, v.cost, v.integral_graph);
    print_spectrum_json(stdout, &v);
    printf("}\n");
}

/** @brief Wczytuje cale wejscie do pamieci. */
static char *slurp(FILE *f){
    size_t cap = 1<<16, len = 0, r;
    char *b = malloc(cap);
    if (!b) return NULL;
    while ((r = fread(b+len, 1, cap-len-1, f)) > 0){
        len += r;
        if (len+1 >= cap){ cap *= 2; char *t = realloc(b, cap); if (!t){ free(b); return NULL; } b = t; }
    }
    b[len] = '\0';
    return b;
}

int main(int argc, char *argv[]){
    int plain = (argc > 1 && strcmp(argv[1], "-g") == 0);
    int only_integral = plain || (argc > 1 && strcmp(argv[1], "-c") == 0);
    char *data = slurp(stdin);
    if (!data){ fprintf(stderr, "brak pamieci\n"); return 1; }
    const char *key = "\"graph6\":\"";
    if (strstr(data, key)){                       /* wejscie JSON (graph6 z ucieczkami \\ i \") */
        char *p = data;
        while ((p = strstr(p, key)) != NULL){
            p += strlen(key);
            char buf[G6_BUF]; int len = 0;
            while (*p && *p != '"' && len < G6_BUF-1){
                if (*p == '\\' && p[1]) p++;       /* usun ucieczke */
                buf[len++] = *p++;
            }
            buf[len] = '\0';
            process(buf, only_integral, plain);
        }
    } else {                                      /* wejscie: jeden graph6 na linie */
        char *line = data;
        while (*line){
            char *nl = strpbrk(line, "\r\n"), *next;
            if (nl){ *nl = '\0'; next = nl+1; } else next = line + strlen(line);
            if (*line && *line != '>') process(line, only_integral, plain);
            line = next;
        }
    }
    fprintf(stderr, "Zweryfikowano %ld grafow, w tym calkowitych: %ld\n", total, integral_cnt);
    free(data);
    return 0;
}