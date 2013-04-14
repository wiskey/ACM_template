#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;
#define N 220
int g[N][N];
class Dinic{
public:
    static const int INF = 1000000007, SIZE = 205;
    int c[SIZE][SIZE];
    int n,s,t,l[SIZE],e[SIZE];
    int flow(int maxf = INF){
        int left=maxf;
        while(build()) left-=push(s,left);
        return maxf-left;
    }
    int push(int x, int f){
        if(x==t) return f;
        int& y=e[x],sum=f;
        for(;y<n;y++) if(c[x][y]>0 && l[x]+1==l[y]){
            int cnt=push(y,min(sum,c[x][y]));
            c[x][y]-=cnt;
            c[y][x]+=cnt;
            sum-=cnt;
            if(!sum) return f;
        }
        return f-sum;
    }
    bool build(){
        int m=0;
        memset(l,255,sizeof(l));
        l[e[m++]=s]=0;
        for(int i=0;i<m;i++) for(int y=0;y<n;y++)
            if(c[e[i]][y]>0 && l[y]<0) l[e[m++]=y]=l[e[i]]+1;
        memset(e,0,sizeof(e));
        return l[t]>=0;
    }
}net;
void Init(int n) {
    memset(net.c, 0, sizeof(net.c));
    for (int i=0; i<n; i++) for (int j=0; j<n; j++)
        net.c[i][j] = g[i][j];
    net.s = 0;
    net.t = n-1;
    net.n = n;
}
int main() {
    freopen("in.txt", "r", stdin);

    int T;
    scanf("%d", &T);
    for (int cas=1; cas<=T; cas++) {
        printf("Case #%d:\n", cas);
        int n, m;
        scanf("%d %d", &n, &m);
        memset(g, 0, sizeof(g));
        int ans = 0, a, b, w;
        for (int i=1; i<=m; i++) {
            scanf("%d%d%d", &a, &b, &w);
            a--, b--;
            g[a][b] = g[b][a] = w;
            Init(n);
            int flow = net.flow();
            if (flow > ans) {
                printf("%d %d\n", i, flow-ans);
                ans = flow;
            }
        }
    }
    return 0;
}
