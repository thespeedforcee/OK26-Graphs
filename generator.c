#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#define MAXN 62

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
//koduje graf, n<=62; zwraca lancuch w buf.
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
        fprintf(stderr, "Uzycie: %s n k [seed=-1][rep=1]\n", argv[0]);
        return 1;
    }
    int n=atoi(argv[1]);
    int k=atoi(argv[2]);
    long long seed=(argc>3)?atoll(argv[3]):-1;
    int rep=(argc>4)?atoi(argv[4]):1;

    if (n<2 || n>MAXN || k<n-1 || k>n*(n-1)/2 || rep<1){
        fprintf(stderr, "Blad:niepoprawne parametry (n=%d, k=%d, rep=%d)\n",n,k,rep);
        return 1;
    }
    if (seed<0){
        seed=(long long)time(NULL); //ziarno ustawione raz tutaj
    }
    Rng rng={ (uint64_t)seed };
    fprintf(stderr,"Generuje %d grafow: n=%d, k=%d, seed=%lld\n",rep,n,k,seed);

    static int adj[MAXN][MAXN];
    char g6[512];
    printf("{\"n\":%d,\"k\":%d,\"seed\":%lld,\"graphs\":[",n,k,seed);
    for (int t=0; t<rep; t++){
        generate(&rng,n,k,adj);
        to_grapf6(n,adj,g6);
        printf("%s{\"graph6\":\"%s\",\"edges\":[",t? ",":"", g6);
        int first=1;
        for (int u=0; u<n; u++){
            for (int v=u+1; v<n; v++){
                if (adj[u][v]){
                    printf("%s[%d,%d]",first? "" : ",", u,v);
                    first=0;
                }
            }
        }
        printf("]}");
    }
    printf("]}\n");
    return 0;
}