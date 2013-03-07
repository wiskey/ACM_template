/*
有向图的spfa，可以修改成无向图的
*/

#include <cstdio>
#include <cstring>
#include <vector>
#include <queue>
using namespace std;
#define MAXN 300005
#define inf 0x7f7f7f7f
struct node {
    int w, v;
    node(int _w, int _v): w(_w), v(_v) {}
};
vector<node> adj[MAXN];
queue<int> Q;
int d[MAXN];
int n, m, s, e;
bool inqueue[MAXN];

int main() {
    while (scanf("%d%d", &n, &m)==2) {
        for (int i=1; i<=n; i++) {
            adj[i].clear();
            d[i] = inf;
            inqueue[i] = false;
        }
        int u, v, w;
        for (int i=0; i<m; i++) {
            scanf("%d%d%d", &u, &v, &w);
            adj[u].push_back(node(w, v));
        }

        scanf("%d%d", &s, &e);   //输入的询问的点
        Q.push(s);
        d[s] = 0;
        inqueue[s] = true;
        while (!Q.empty()) {
            int now = Q.front();
            Q.pop();
            inqueue[now] = false;
            for (unsigned int i=0; i<adj[now].size(); i++)
                if (d[now]+adj[now][i].w < d[adj[now][i].v]) {
                    d[adj[now][i].v] = d[now]+adj[now][i].w;
                    if (!inqueue[adj[now][i].v]) {
                        Q.push(adj[now][i].v);
                        inqueue[adj[now][i].v] = true;
                    }
                }
        }
        if (d[e] == inf) printf("-1\n");//s->e的最短路，如果没有，输出-1
        else printf("%d\n", d[e]);
    }
    return 0;
}
