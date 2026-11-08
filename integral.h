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

#define MAXN 62 //max rzad grafu
#define MAXN_EXACT 25 //do tego dokl liczy
#define MAX_PAIRS (MAXN*(MAXN-1)/2)
#define G6_BUF 330 //bufor na grapf6
#define COST_INF 1e300 //koszt grafu niespojnego
#define COST_EPS 1e-12 //dokładność obliczeń zmiennoprzecinkowych

typedef __int128 i128;
typedef int adj[MAXN][MAXN]; //macierz sasiedztwa
typedef struct{ 
    int u,v;
} Pair;

//czas w sek, monotoniczny
static inline double now_s(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec+1e-9*(double)ts.tv_nsec;
}
//stan generatora liczb pseudolosowych
typedef struct{
    uint64_t s;
} Rng;
//zwraca nastepna liczbe pseudolosowa
static inline uint64_t rng_next(Rng *r) {
    uint64_t z=(r->s+=0x9E3779B97F4A7C15ULL);
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;
    z=(z^(z>>27))*0x94D049BB133111EBULL;
    return z^(z>>31);
}
//liczba z zakr. [0,m]
static inline int rng_int(Rng *r, int m){
    return (int)(rng_next(r)%(uint64_t)m);
}
//ziarno, jesli na minusie to bierze z zegar
static inline long long resolve_seed(long long seed){
    if (seed<0){
        seed=(long long)time(NULL);
    }
    return seed;
}
//losowa permutacja tab par
static inline void shuffle_pairs(Rng *r, Pair *pairs, int num_pairs){
    for (int i=num_pairs-1; i>0; i--){
        int j=rng_int(r, i+1);
        Pair t=pairs[i];
        pairs[i]=pairs[j];
        pairs[j]=t;
    }
}
//lb krawedzi
static inline int count_edges(int n, const adj a){
    int k=0;
    for (int i=0; i<n; i++){
        for (int j=i+1; j<n; j++){
            k+=a[i][j];
        }
    }
    return k;
}
//max st. wierzcholka
static inline int max_degree(int n, const adj a){
    int maxd=0;
    for (int i=0; i<n; i++){
        int d=0;
        for (int j=0; j<n; j++){
            d+=a[i][j];
        }
        if (d>maxd){
            maxd=d;
        }
    }
    return maxd;
}

//test spojnosci bfs
static int is_connected(int n, adj adj){
    int visited[MAXN]={0};
    int stack[MAXN];
    int top=0;
    int height=0;
    int cnt=1;
    visited[0]=1;
    stack[top++]=0;
    while (height<top){
        int u=stack[height++];
        for (int v=0; v<n; v++){
            if (adj[u][v] && !visited[v]){
                visited[v]=1;
                stack[top++]=v;
                cnt++;
            }
        }
    }
    return cnt==n;
}

//grapf6 -> macierz sasiedztwa
static inline int g6_decode(const char *g6, adj adj){
    int n=(unsigned char)g6[0]-63;
    if (n<1 || n>MAXN){
        return -1;
    }
    int need=(n*(n-1)/2+5)/6;
    if ((int)strlen(g6)<1+need){
        return -1;
    }
    for (int i=1; i<=need; i++){
        unsigned char c=(unsigned char)g6[i];
        if (c<63 || c>126){
            return -1;
        }
    }
    memset(adj,0,sizeof(int)*MAXN*MAXN);
    int pos=1;
    int bit=5;
    for (int j=1; j<n; j++){
        for (int i=0; i<j; i++){
            int v=((unsigned char)g6[pos]-63)>>bit&1;
            adj[i][j]=adj[j][i]=v;
            if (--bit<0){
                bit=5;
                pos++;
            }
        }
    }
    return n;
}

//koduje graf do grpaf6
static inline void g6_encode(int n, adj adj, char *buf){
    int pos=0;
    int bits=0;
    int val=0;
    buf[pos++]=(char)(n+63);
    for (int j=1; j<n; j++){
        for (int i=0; i<j; i++){
            val=(val<<1)|adj[i][j];
            if (++bits==6){
                buf[pos++]=(char)(val+63);
                val=0;
                bits=0;
            }
        }
    }
    if (bits>0){
        buf[pos++]=(char)((val<<(6-bits))+63);
    }
    buf[pos]='\0';
}

//generator: losowy graf spojny G(n,k), losowe drzewo, losowanie bez zwracania
static inline void generate_connected(Rng *r, int n, int k, adj adj){
    memset(adj, 0, sizeof(int)*MAXN*MAXN);
    int perm[MAXN];
    for (int i=0; i<n; i++){
        perm[i]=i;
    }
    for (int i=n-1; i>0; i--){
        int j=rng_int(r, i+1);
        int t=perm[i];
        perm[i]=perm[j];
        perm[j]=t;
    }
    for (int i=1; i<n; i++){
        int p=perm[rng_int(r,i)];
        adj[p][perm[i]]=adj[perm[i]][p]=1;
    }
    Pair miss[MAX_PAIRS];
    int m=0;
    for (int u=0; u<n; u++){
        for (int v=u+1; v<n; v++){
            if (!adj[u][v]){
                miss[m++]=(Pair){u,v};
            }
        }
    }
    for (int i=0; i<k-(n-1); i++){
        int j=i+rng_int(r,m-i);
        Pair e =miss[j];
        miss[j]=miss[i];
        adj[e.u][e.v]=adj[e.v][e.u]=1;
    }
}

//obl zmiennoprzecinkowe, wart wlasne macierzy sasiedztwa cykliczna metoda jacobiego
static inline int jacobi(int n, adj adj, double *lam){
    double a[MAXN][MAXN];
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            a[i][j]=adj[i][j];
        }
    }
    int sweep;
    for (sweep=0; sweep<60; sweep++){
        double off=0.0;
        for (int p=0; p<n; p++){
            for (int q=p+1; q<n; q++){
                off+=a[p][q]*a[p][q];
            }
        }
        if (off<1e-24){
            break;
        }
        for (int p=0; p<n-1; p++){
            for (int q=p+1; q<n; q++){
                if (fabs(a[p][q])<1e-300){
                    continue;
                }
                double theta=(a[q][q]-a[p][p])/(2.0*a[p][q]);
                double t=(theta>=0?1.0:-1.0)/(fabs(theta)+sqrt(theta*theta+1.0));
                double c=1.0/sqrt(t*t+1.0);
                double s=t*c;
                for (int k=0; k<n; k++){ //kolumny
                    double akp=a[k][p];
                    double akq=a[k][q];
                    a[k][p]=c*akp-s*akq;
                    a[k][q]=s*akp+c*akq;
                }
                for (int k=0; k<n; k++){ //wiersze
                    double apk=a[p][k];
                    double aqk=a[q][k];
                    a[p][k]=c*apk-s*aqk;
                    a[q][k]=s*apk+c*aqk;
                }
            }
        }
    } 
    for (int i=0; i<n; i++){
        lam[i]=a[i][i];
    }
    return sweep;      
}

//koszt calkowitosci
static inline double integrality_cost(int n, const double *lam){
    double c=0.0;
    for (int i=0; i<n; i++){
        c+=fabs(lam[i]-round(lam[i]));
    }
    return c;
}

//funkcja celu minimalizowana przez heurystyki
static inline double graph_cost(int n, adj adj){
    if (!is_connected(n,adj)){
        return COST_INF;
    }
    double lam[MAXN];
    jacobi(n,adj,lam);
    return integrality_cost(n,lam);
}

//obl calkowitoliczbowwe (dokladne)
//wspolczynniki wielomianu charakterystycznego
static inline void char_poly(int n, adj adj, i128 *c){
    i128 m[MAXN_EXACT][MAXN_EXACT], am[MAXN_EXACT][MAXN_EXACT];
    memset(m, 0, sizeof(m));
    c[n]=1;
    for (int k=1; k<=n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                i128 sum=0;
                for (int l=0; l<n; l++){
                    if (adj[i][l]){
                        sum+=m[l][j];
                    }
                }
                am[i][j]=sum;
            }
        }
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                m[i][j]=am[i][j]+(i==j?c[n-k+1]:0);
            }
        }
        i128 tr=0;
        for (int i=0; i<n; i++){
            for (int l=0; l<n; l++){
                if (adj[i][l]){
                    tr+=m[l][i];
                }
            }
        }
        c[n-k]=-tr/k;
    }
}

//sprawdza czy wielomain ma n pierwiastkow calk 
//dzielenie hornera
static inline int integer_spectrum(int n, const i128 *c, int maxdeg, int *spec){
    i128 p[MAXN_EXACT+1];
    for (int i=0; i<=n; i++){
        p[i]=c[i];
    }
    int deg=n;
    int found=0;
    for (int r=-maxdeg; r<=maxdeg; r++){
        while(deg>0){
            i128 q[MAXN_EXACT+1], rem=p[deg];
            for (int i=deg-1; i>=0; i--){
                q[i]=rem;
                rem=p[i]+rem*r;
            }
            if (rem!=0){
                break;
            }
            for (int i=0; i<deg; i++){
                p[i]=q[i];
            }
            deg--;
            spec[found++]=r;
        }
    }
    return found==n;
}

//weryfikacja wyniku
typedef struct {
    int n, k;
    int connected;
    int sweeps;
    double cost;
    int float_integral;
    int exact_integral;
    int integral_graph;
    int spec[MAXN_EXACT];
} Verdict;

//pelna weryfikacja jacobi(float)+wielomian charakterystyczny(dokladna)
static inline void verify_graph(int n, adj adj, double tol, Verdict *v){
    v->n=n;
    v->k=count_edges(n,adj);
    v->connected=is_connected(n,adj);
    double lam[MAXN];
    v->sweeps=jacobi(n,adj,lam);
    v->cost=integrality_cost(n,lam);
    v->float_integral=v->cost<tol;
    v->exact_integral=-1;
    if (n<=MAXN_EXACT){
        i128 c[MAXN_EXACT+1];
        int maxdeg=max_degree(n,adj);
        char_poly(n,adj,c);
        v->exact_integral=integer_spectrum(n,c,maxdeg,v->spec);
    }
    v->integral_graph=v->connected&&(v->exact_integral==1 || (v->exact_integral==-1 && v->float_integral));
}

//wypisuje widmo jako tab json lub null
static inline void print_spectrum_json(FILE *f, const Verdict *v){
    if (v->exact_integral != 1){
        fputs("null", f);
        return;
    }
    fputc('[', f);
    for (int i=0; i<v->n; i++){
        fprintf(f, "%s%d", i?",":"", v->spec[i]);
    }
    fputc(']',f);
}

//statystyki przeszukiwania
typedef struct {
    long long evals;
    long long steps;
    int timed_out;
} SearchStats;

//usuwa kraw rem i dodaje add
static inline void swap_edges(adj adj, Pair rem, Pair add){
    adj[rem.u][rem.v]=adj[rem.v][rem.u]=0;
    adj[add.u][add.v]=adj[add.v][add.u]=1;
}

static inline double local_search(int n, adj adj, Rng *rng, double deadline, double tol, SearchStats *st){
    Pair edges[MAX_PAIRS];
    Pair non_edges[MAX_PAIRS];
    double cur=graph_cost(n,adj);
    st->evals++;
    while (cur>=tol){
        int ne=0;
        int nn=0;
        for (int i=0; i<n; i++){
            for (int j=i+1; j<n; j++){
                if (adj[i][j]){
                    edges[ne++]= (Pair){i,j};
                } else {
                    non_edges[nn++]=(Pair){i,j};
                }
            }
        }
        if (rng){
            shuffle_pairs(rng, edges, ne);
            shuffle_pairs(rng, non_edges,nn);
        }
        double best=cur;
        int bi=-1;
        int bj=-1;
        int improved=0;
        int timeout=0;
        for (int a=0; a<ne && !improved && !timeout; a++){
            for (int b=0; b<nn; b++){
                if (now_s()>=deadline){
                    timeout=1;
                    break;
                }
                swap_edges(adj, edges[a], non_edges[b]);
                double c=graph_cost(n,adj);
                st->evals++;
                if (c<best-COST_EPS){
                    best=c;
                    bi=a;
                    bj=b;
                    if (rng){
                        improved=1;
                        break;
                    }
                }
                swap_edges(adj, non_edges[b], edges[a]);
            }
        }
        if (rng){
            if (improved){
                cur=best;
                st->steps++;
            }
        } else if (bi>=0){
            swap_edges(adj, edges[bi], non_edges[bj]);
            cur=best;
            improved=1;
            st->steps++;
        }
        if (timeout){
            st->timed_out=1;
            break;
        }
        if (!improved){
            break;
        }
    }
    return cur;
}

static inline void print_solution_json(const char *alg, int n, int k, long long seed, int restarts, const SearchStats *st, double cost, double time_s, const char *g6, const Verdict *v) {
    printf("{\n");
    printf("  \"algorithm\": \"%s\",\n", alg);
    printf("  \"n\": %d,\n  \"k\": %d,\n  \"seed\": %lld,\n", n, k, seed);
    printf("  \"restarts\": %d,\n", restarts);
    printf("  \"evaluations\": %lld,\n  \"steps\": %lld,\n", st->evals, st->steps);
    printf("  \"timed_out\": %d,\n", st->timed_out);
    printf("  \"time_s\": %.3f,\n", time_s);
    printf("  \"cost\": %.6e,\n", cost);
    printf("  \"graph6\": \"%s\",\n", g6);
    printf("  \"connected\": %d,\n", v->connected);
    printf("  \"exact_integral\": %d,\n", v->exact_integral);
    printf("  \"integral_graph\": %d,\n", v->integral_graph);
    printf("  \"spectrum\": ");
    print_spectrum_json(stdout, v);
    printf("\n}\n");
}
 
#endif