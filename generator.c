#include "integral.h"

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
    seed=resolve_seed(seed);
    Rng rng={ (uint64_t)seed };
    fprintf(stderr,"Generuje %d grafow: n=%d, k=%d, seed=%lld\n",rep,n,k,seed);

    static adj adj;
    char g6[G6_BUF];
    printf("{\"n\":%d,\"k\":%d,\"seed\":%lld,\"graphs\":[",n,k,seed);
    for (int t=0; t<rep; t++){
        generate_connected(&rng,n,k,adj);
        g6_encode(n,adj,g6);
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