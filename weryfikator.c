#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MAXN 62
#define MAXN_EXACT 25

typedef __int128 i128;

//wczytuje graf do macierzy sasiedztwa, zwraca n lub -1
static int from_grapf6(const char *s, int adj[MAXN][MAXN]){
    int n=(unsigned char)s[0]-63;
    if (n<1 || n>MAXN){
        return -1;
    }
    memset(adj, 0, sizeof(int)*MAXN*MAXN);
    int need=(n*(n-1)/2+5)/6;
    if ((int)strlen(s)<1+need){
        return -1;
    }
    int pos=1;
    int bit=5;
    for (int j=1; j<n; j++){
        for (int i=0; i<j; i++){
            int v=((unsigned char)s[pos]-63) >> bit&1;
            adj[i][j]=adj[j][i]=v;
            if (--bit<0){
                bit=5;
                pos++;
            }
        }
    }
    return n;
}

//test spójności bfs
static int is_connected(int n, int adj[MAXN][MAXN]){
    int q[MAXN];
    int vis[MAXN] = {0};
    int h=0;
    int t=0;
    int cnt=1;
    q[t++]=0;
    vis[0]=1;
    while (h<t){
        int u=q[h++];
        for (int v=0; v<n; v++){
            if (adj[u][v] && !vis[v]){
                vis[v]=1;
                q[t++]=v;
                cnt++;
            }
        }
    }
    return cnt==n;
}

//metoda jacobiega cykliczna
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
//wspol. wielomianu charakterystycznego, dziel. przez k jest dokł
static void char_poly(int n, int adj[MAXN][MAXN], i128 *c){
    static i128 M[MAXN][MAXN], AM[MAXN][MAXN];
    memset(M, 0, sizeof(M));
    c[n]=1;
    for (int k=1; k<=n; k++){
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                i128 s=0;
                for (int l=0; l<n; l++){
                    if (adj[i][l]){
                        s+=M[l][j];
                    }
                }
                AM[i][j]=s;
            }
        }
        for (int i=0; i<n; i++){
            for (int j=0; j<n; j++){
                M[i][j]=AM[i][j] + (i==j?c[n-k+1]:0);
            }
        }
        i128 tr=0;
        for (int i=0; i<n; i++){
            for (int l=0; l<n; l++){
                if (adj[i][l]){
                    tr+=M[l][i];
                }
            }
        }
        c[n-k]=-tr/k;
    }
}
//sprawdza czy wielomian ma n pierwiastkow całk. dzielenie hornera
static int all_integer_roots(int n, i128 *c, int maxdeg, int *mult){
    i128 p[MAXN+1];
    for (int i=0; i<=n; i++){
        p[i]=c[i];
    }
    int deg=n;
    int found=0;
    for (int r=-maxdeg; r<=maxdeg; r++){
        mult[r+maxdeg]=0;
        while(deg>0){
            i128 q[MAXN+1], rem=p[deg];
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
            found++;
            mult[r+maxdeg]++;
        }
    }
    return found==n;
}

//przetwarza jeden graf i wypisuje json
static void process(const char *g6, double tol){
    static int adj[MAXN][MAXN];
    int n=from_grapf6(g6, adj);
    if (n<0){
        fprintf(stderr, "niepoprawny graf: %s\n", g6);
        return;
    }
    int k=0;
    int maxdeg=0;
    for (int i=0; i<n; i++){
        int d=0;
        for (int j=0; j<n; j++){
            d+=adj[i][j];
        }
        k+=d;
        if (d>maxdeg){
            maxdeg=d;
        }
    }
    k/=2;
    int conn=is_connected(n,adj);
    double lam[MAXN];
    int sweeps=jacobi(n,adj,lam);
    double cost=integrality_cost(n,lam);
    int float_ok=cost<tol;
    int exact_ok=-1;
    if(n<=MAXN_EXACT){
        i128 c[MAXN+1];
        int mult[2*MAXN+1];
        char_poly(n,adj,c);
        exact_ok=all_integer_roots(n,c,maxdeg,mult);
    }
    printf("{\"graph6\":\"%s\",\"n\":%d,\"k\":%d,\"connected\":%d,""\"cost\":%.3e,\"jacobi_sweeps\":%d,\"float_integral\":%d,""\"exact_integral\":%d,\"integral_graph\":%d}\n", g6,n,k,conn,cost,sweeps, float_ok,exact_ok,conn&&(exact_ok==1 || (exact_ok==-1 && float_ok)));
    if (exact_ok>=0 && exact_ok != float_ok){
        fprintf(stderr, "uwaga: rozbieznosc float vs exact: %s\n (float=%d, exact=%d)", g6, float_ok, exact_ok);

    }
}  

static long getline(char **line, size_t *cap, FILE *f){
    size_t len=0;
    int ch;
    if (!*line || *cap==0){
        *cap=4096;
        *line=malloc(*cap);
        if (!*line){
            return -1;
        }
    }
    while ((ch=fgetc(f))!=EOF){
        if (len+2>*cap){
            *cap*=2;
            char *tmp=realloc(*line,*cap);
            if (!tmp){
                return -1;
            }
            *line=tmp;
        }
        (*line)[len++]=(char)ch;
        if(ch=='\n'){
            break;
        }
    }
    if(len==0){
        return -1;
    }
    (*line)[len]='\0';
    return (long)len;
}

int main(int argc, char *argv[]){
    double tol=(argc>1)?atof(argv[1]):1e-6;
    char *line=NULL;
    size_t cap=0;
    while(getline(&line, &cap, stdin)>0){
        char *p=line;
        if (strstr(line, "\"graph6\":\"")){
            while ((p=strstr(p, "\"graph6\":\""))!=NULL){
                p+=10;
                char *e=strchr(p, '"');
                if (!e){
                    break;
                }
                *e='\0';
                process(p,tol);
                p=e+1;
            }
        } else {
            line[strcspn(line, "\r\n")]='\0';
            if (line[0]){
                process(line,tol);
            }
        }
    }
    free(line);
    return 0;
}