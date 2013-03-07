#include <algorithm>
#include <queue>
#include <cstring>
#include <cstdio>
using namespace std;

const int maxn = 10000;
const int maxm = 1000000;
const int INF = 1000000000;

struct record {
    int v, f, c, next;
} edge[maxm];

int cas, ans, cl, n, m, s, t, aug;
int dist[maxn], pre[maxn], pos[maxn];
bool vis[maxn];

void connect(int a, int b, int f, int c) {
    cl++;
    edge[cl].next = pos[a];
    edge[cl].v = b;
    edge[cl].f = f;
    edge[cl].c = c;
    pos[a] = cl;
    cl++;
    edge[cl].next = pos[b];
    edge[cl].v = a;
    edge[cl].f = 0;
    edge[cl].c = -c;
    pos[b] = cl;
}

bool spfa() {
    queue<int> q;

    memset(vis, 0, sizeof (vis));
    for (int i = 0; i < n; i++) dist[i] = INF;

    q.push(s);
    dist[s] = 0;
    pre[s] = 0;
    while (!q.empty()) {
        int k = q.front();
        q.pop();
        vis[k] = false;
        for (int p=pos[k]; p; p=edge[p].next)
            if (edge[p].f>0 && edge[p].c+dist[k]<dist[edge[p].v]) {
                dist[edge[p].v] = edge[p].c + dist[k];
                pre[edge[p].v] = p;
                if (!vis[edge[p].v]) {
                    q.push(edge[p].v);
                    vis[edge[p].v] = true;
                }
            }
    }
    if (dist[t] == INF) return false;
    aug = INF;
    for (int p=pre[t]; p; p=pre[edge[p^1].v])
        aug = min(aug, edge[p].f);

    for (int p=pre[t]; p; p=pre[edge[p^1].v]) {
        edge[p].f -= aug;
        edge[p^1].f += aug;
    }
    ans += dist[t] * aug;
    return true;
}

int main() {
    scanf("%d", &cas);
    while (cas--) {
        cl = 1;
        scanf("%d%d%d%d", &n, &m, &s, &t);
        memset(pos, 0, sizeof(pos));
        int u, v, f, c;
        for (int i = 0; i < m; i++) {
            scanf("%d%d%d%d", &u, &v, &f, &c);
            connect(u, v, f, c);
        }
        ans = 0;
        while (spfa()) ;
        printf("%d\n", ans);
    }
    return 0;
}

