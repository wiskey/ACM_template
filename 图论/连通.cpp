#include <cstdio>
#include <cstring>
using namespace std;
const int MAXN = 1000;
const int MAXM = 100000;
int dfn[MAXN], low[MAXN];
struct node {
    int v, w, pre;
}edge[MAXM<<1];
int nEdge, pos[MAXN];
int n, m;

void connect(int u, int v, int w) {
    nEdge++;
    edge[nEdge].pre = pos[u];
    edge[nEdge].v = v;
    edge[nEdge].w = w;
    pos[u] = nEdge;
}
int Tarjanbfs(int cur, int &sig, int &tcc, int from) {
    sta[++head] = cur;
    dfn[cur] = ++sig;
    low[cur] = sig;
    for (int j=pos[cur]; j; j=edge[j].pre) {
        if (i^1 == from) continue;
        if (!dfn[edge[i].v) {
            Tarjanbfs(edge[i].v, sig, tcc, cur);
            low[cur] = min(low[cur], low[edge[i].v]);
            if (dfn[cur] <= low[edge[i].v]) {
                ans[tcc][0] = cur;
                int cow = 1;
                do {
                    ans[tcc][cow++] = sta[head];
                }while (sta[head--] != edge[i].v);
                if (cow > 2) ++tcc;
            }
        }
        else low[cur] = min(low[cur], low[edge[i].v]);
    }
    return 0;
}
int main() {
    int T;
    scanf("%d", &T);
    while (T--) {
        scanf("%d%d", &n, &m);
        int u, v, w;
        nEdge = 0;
        for (int i=0; i<m; i++) {
            scanf("%d%d%d", &u, &v, &w);
            connect(u, v, w);
            connect(v, u, w);
        }
    }
    return 0;
}








