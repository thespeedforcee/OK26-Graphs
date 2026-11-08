#include "integral.h"

static long total=0, integral_cnt=0;

//przetwarza jeden graf i wypisuje json
static void process(const char *g6, double tol){
    static adj adj;
    int n=g6_decode(g6, adj);
    if (n<0){
        fprintf(stderr, "niepoprawny graf: %s\n", g6);
        return;
    }
    Verdict v;
    verify_graph(n, adj, tol, &v);
    total++;
    integral_cnt += v.integral_graph;
    printf("{\"graph6\":\"%s\",\"n\":%d,\"k\":%d,\"connected\":%d,\"cost\":%.3e," "\"jacobi_sweeps\":%d,\"float_integral\":%d,\"exact_integral\":%d," "\"integral_graph\":%d,\"spectrum\":", g6, v.n, v.k, v.connected, v.cost, v.sweeps, v.float_integral, v.exact_integral, v.integral_graph);
    print_spectrum_json(stdout, &v);
    printf("}\n");
    if (v.exact_integral >= 0 && v.exact_integral != v.float_integral){
        fprintf(stderr, "uwaga: rozbieznosc float vs exact: %s (float=%d, exact=%d)\n", g6, v.float_integral, v.exact_integral);
    }
}

static char *slurp(FILE *f) {
    size_t cap=1<<16, len = 0;
    char *b=malloc(cap);
    if (!b){
        return NULL;
    }
    size_t r;
    while ((r = fread(b+len, 1, cap-len-1, f))>0) {
        len+=r;
        if (len+1>=cap) {
            cap*=2;
            char *t=realloc(b, cap);
            if (!t) { 
                free(b); return NULL; 
            }
            b=t;
        }
    }
    b[len]='\0';
    return b;
}

int main(int argc, char *argv[]) {
    double tol=(argc>1)?atof(argv[1]):1e-6;
    char *data=slurp(stdin);
    if (!data) {
        fprintf(stderr, "brak pamieci\n");
        return 1;
    }
    const char *key="\"graph6\":\"";
    if (strstr(data, key)) {
        char *p=data;
        while ((p=strstr(p, key))!=NULL) {
            p+=strlen(key);
            char *e=strchr(p, '"');
            if (!e){
                break;
            }
            *e='\0';
            process(p, tol);
            p=e+1;
        }
    } else {
        char *line=data;
        while (*line) {
            char *nl=strpbrk(line, "\r\n");
            char *next;
            if (nl) { 
                *nl='\0'; 
                next=nl+1; 
            } else {
                next=line+strlen(line);
            }
            if (*line && *line!='>') {
                process(line, tol);
            }
            line=next;
        }
    }
    fprintf(stderr, "Zweryfikowano %ld grafow, w tym calkowitych: %ld\n", total, integral_cnt);
    free(data);
    return 0;
}