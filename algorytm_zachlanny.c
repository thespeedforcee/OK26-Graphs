#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAXN 62
//z weryfikator.c
static int jacobi(int n, int adj[MAXN][MAXN], double *lam){
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
//koszt całkowitości?
static double integrality_cost(int n, const double *lam){
    double c=0.0;
    for (int i=0; i<n; i++){
        c+=fabs(lam[i]-round(lam[i]));
    }
    return c;
}
//z generaotr.c
typedef struct{
    uint64_t s;
} Rng;

static uint64_t rng_next(Rng *r) {
    uint64_t z=(r->s+=0x9E3779B97F4A7C15ULL);
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;
    z=(z^(z>>27))*0x94D049BB133111EBULL;
    return z^(z>>31);
}

static int rng_int(Rng *r, int m){
    return (int)(rng_next(r)%(uint64_t)m);
}

static void generate(Rng *r, int n, int k, int adj[MAXN][MAXN]){
    memset(adj, 0, sizeof(int)*MAXN*MAXN);
    //losowa permutacja wierzchołków
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
    //kazdy kolejny wierzcholek laczony z losowym wczesniejszym
    for (int i=1; i<n; i++){
        int p=perm[rng_int(r,i)];
        adj[p][perm[i]]=adj[perm[i]][p]=1;
    }
    //lista brakujacych par
    int(*miss)[2]=malloc(sizeof(int[2])*(size_t)(n*(n-1)/2));
    int m=0;
    for (int u=0; u<n; u++){
        for (int v=u+1; v<n; v++){
            if (!adj[u][v]){
                miss[m][0]=u;
                miss[m][1]=v;
                m++;
            }
        }
    }
    //losowanie bez zwracania k-(n-1) par 
    for (int i=0; i<k-(n-1); i++){
        int j=i+rng_int(r, m-i);
        int a=miss[j][0];
        int b=miss[j][1];
        miss[j][0]=miss[i][0];
        miss[j][1]=miss[i][1];
        adj[a][b]=adj[b][a]=1;
    }
    free(miss);
}

double evaluate_graph(int n, int adj[MAXN][MAXN]){
    double lam[MAXN];
    int sweeps=jacobi(n,adj,lam);
    double cost=integrality_cost(n,lam);
    return cost;
}

void Algorytm_Zachlanny(int n, int k, int adj[MAXN][MAXN]){
    double current_cost=evaluate_graph(n, adj);
    printf("Początkowy koszt całkowitości: %.3e\n", current_cost);
    int iteration=0;
    int improved=1;
    while (improved && current_cost>1e-6){
        improved=0;
        int edges[MAXN*MAXN/2][2];
        int num_edges=0;
        int non_edges[MAXN*MAXN/2][2];
        int num_non_edges=0;
        for (int i=0; i<n; i++){
            for (int j=i+1; j<n; j++){
                if (adj[i][j]){
                    edges[num_edges][0]=i;
                    edges[num_edges][1]=j;
                    num_edges++;
                } else {
                    non_edges[num_non_edges][0]=i;
                    non_edges[num_non_edges][1]=j;
                    num_non_edges++;
                }
            }
        }
        double best_cost=current_cost;
        int best_remove=-1;
        int best_add=-1;
        for (int e=0; e<num_edges; e++){
            int u1=edges[e][0];
            int v1=edges[e][1];
            for (int ne=0; ne<num_non_edges; ne++){
                int u2=non_edges[ne][0];
                int v2=non_edges[ne][1];
                //usuń krawędź (u1,v1) i dodaj krawędź (u2,v2)
                adj[u1][v1]=adj[v1][u1]=0;
                adj[u2][v2]=adj[v2][u2]=1;
                double new_cost=evaluate_graph(n, adj);
                if (new_cost<best_cost){
                    best_cost=new_cost;
                    best_remove=e;
                    best_add=ne;
                }
                //przywróć graf
                adj[u1][v1]=adj[v1][u1]=1;
                adj[u2][v2]=adj[v2][u2]=0;
            }
        }
        if (best_cost<current_cost){
            //wykonaj najlepszą zmianę
            int u1=edges[best_remove][0];
            int v1=edges[best_remove][1];
            int u2=non_edges[best_add][0];
            int v2=non_edges[best_add][1];
            adj[u1][v1]=adj[v1][u1]=0;
            adj[u2][v2]=adj[v2][u2]=1;
            current_cost=best_cost;
            improved=1;
            iteration++;
            printf("Iteracja %d: koszt calkowitosci: %.3e\n", iteration, current_cost);
        }
    }
    printf("\n");
    if (current_cost<=1e-6){
        printf("Znaleziono graf calkowity po %d iteracjach.\n", iteration);
    } else {
        printf("Nie znaleziono grafu calkowitego po %d iteracjach. Ostateczny koszt calkowitolci: %.3e\n", iteration, current_cost);
    }
}

static void to_grapf6(int n, int adj[MAXN][MAXN], char *buf){
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

int main(int argc, char *argv[]){
    if (argc<3){
        fprintf(stderr, "Uzycie: %s n k [seed=-1]\n", argv[0]);
        return 1;
    }
    int n=atoi(argv[1]);
    int k=atoi(argv[2]);
    int adj[MAXN][MAXN];
    long long seed=(argc>3)?atoll(argv[3]):-1;
    if (n<2 || n>MAXN || k<n-1 || k>n*(n-1)/2){
        fprintf(stderr, "Blad:niepoprawne parametry (n=%d, k=%d)\n",n,k);
        return 1;
    }
    if (seed<0){
        seed=(long long)time(NULL); //ziarno ustawione raz tutaj
    }
    Rng rng={ (uint64_t)seed };
    fprintf(stderr,"Generuje graf: n=%d, k=%d, seed=%lld\n",n,k,seed);
    adj[MAXN][MAXN];
    generate(&rng,n,k,adj);
    char g6_start[512];
    to_grapf6(n,adj,g6_start);
    fprintf(stderr, "Poczatkowy graf: %s\n", g6_start);
    Algorytm_Zachlanny(n, k, adj);
    char g6_final[512];
    to_grapf6(n,adj,g6_final);
    printf("{\n");
    printf("  \"n\": %d,\n", n);
    printf("  \"k\": %d,\n", k);
    printf("  \"seed\": %lld,\n", seed);
    printf("  \"graph6_start\": \"%s\",\n", g6_start);
    printf("  \"graph6_end\": \"%s\"\n", g6_final);
    printf("}\n");
    return 0;
}