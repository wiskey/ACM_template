/*
 *无向图的桥及边的双连通分量，Tarjan算法O(E)
 */
#include <cstdio>
#include <cstring>
using namespace std;
#define MAXN 10000
#define MAXM 1000000

struct node {
    int v, w, pre;
} edge[MAXM];
int pos[MAXN], nEdge; //图的存储：链式前向星（池子法）

struct Bridge {
    int u, v;
} bridge[MAXM];  //用来记录桥
int tot; //桥的个数

int fa[MAXN], cc; //fa：各个点所属的缩点（连通块），cc连通块的个数
int dfn[MAXN], low[MAXN], time; //时间戳
int stack[MAXN], top;   //用于维护连通块的
int n, m;   //点的个数和边的条数

void connect(int u, int v, int w) {
    nEdge++;
    edge[nEdge].pre = pos[u];
    edge[nEdge].v = v;
    edge[nEdge].w = w;
    pos[u] = nEdge;
}

void tarjan(int cur, int from) {
    low[cur] = dfn[cur] = time++;
    stack[++top] = cur;
    for (int p=pos[cur]; p; p=edge[p].pre) {
        int v = edge[p].v;
        if (v == from) continue;  //注意一下这里
        if (!dfn[v]) {
            tarjan(v, cur);
            if (low[v] < low[cur]) low[cur] = low[v];
            if (low[v] > dfn[cur]) {
                bridge[tot].u = cur;
                bridge[tot++].v = v;
                cc++;
                do {
                    fa[stack[top]] = cc;
                } while (stack[top--] != v);
            }
        } else if (low[cur] > dfn[v]) low[cur] = dfn[v];
    }
}
int main() {
    scanf("%d%d", &n, &m);

    memset(pos, 0, sizeof(pos));
    nEdge = 0;
    int u, v, w;
    for (int i=0; i<m; i++) {
        scanf("%d%d%d", &u, &v, &w);
        connect(u, v, w);
        connect(v, u, w);
    }

    memset(dfn, 0, sizeof(dfn));
    memset(fa, -1, sizeof(fa));

    cc = tot = 0;
    for (int i=1; i<=n; i++)   //可以处理不连通的无向图，如果连通只需要一次即可
        if (!dfn[i]) {
            top = time = 1;
            tarjan(i, -1);
            ++cc;
            for (int j=1; j<=n; j++)   //特殊处理顶点的连通块
                if (dfn[j] && fa[j]==-1) fa[j] = cc;
        }

    for (int i=1; i<=n; i++)
        printf("%d ", fa[i]);  //输出每个节点所属于的连通块（缩点标号）
    printf("\n");

    for (int i=0; i<tot; i++)
        printf("%d %d\n", bridge[i].u, bridge[i].v); //输出所有的桥

    return 0;
}
