#include <cstdio>
#include <cstring>
#include <queue>
#include <algorithm>
using namespace std;
const int MAXN = 100;
const int MAXM = 10000;
struct node {
    int v, w, pre;
} edge[MAXN];
int nEdge, pos[MAXN];
int n, m;

void connect(int u, int v, int w) {
    nEdge++;
    edge[nEdge].pre = pos[u];
    edge[nEdge].v = v;
    edge[nEdge].w = w;
    pos[u] = nEdge;
}
int spfa() {
    bool inq[MAXN];
    int d[MAXN];
    memset(d, 0x7f, sizeof(d));
    queue<int> Q;
    memset(inq, false, sizeof(inq));
    Q.push(1);
    d[1] = 0;
    inq[1] = true;
    while (!Q.empty()) {
        int v, now = Q.front();
        Q.pop();
        for (int j=pos[now]; j; j=edge[j].pre) {
            v = edge[j].v;
            if (d[now] + edge[j].w < d[v]) {
                d[v] = d[now] + edge[j].w;
                if (!inq[v]) {
                    inq[v] = true;
                    Q.push(v);
                }
            }
        }
        inq[now] = false;
    }
    return d[n];
}
int main() {
    while (scanf("%d%d", &n, &m) == 2) {
        int u, v, w;
        nEdge = 0;
        for (int i=0; i<m; i++) {
            scanf("%d%d%d", &u, &v, &w);
            connect(u, v, w);
            connect(v, u, w);
        }
        printf("%d\n", spfa());
    }
    return 0;
}
