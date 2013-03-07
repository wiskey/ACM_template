#include <cstring>
#include <algorithm>
#include <cstdio>
using namespace std;

const int MAXM = 1000000;
const int MAXN = 10000;
const int INF = 1000000000;

struct record {
    int v, f, next;
} edge[MAXM];

int pos[MAXN], dis[MAXN], vh[MAXN], cl;
int his[MAXN], di[MAXN], pre[MAXN];

void connect(int a, int b, int f) {
    cl++;
    edge[cl].next = pos[a];
    edge[cl].v = b;
    edge[cl].f = f;
    pos[a] = cl;
    cl++;
    edge[cl].next = pos[b];
    edge[cl].v = a;
    edge[cl].f = 0; //若为无向边，则f = f
    pos[b] = cl;
}
int maxflow(int s, int t, int n) {
    vh[0] = n; //初始化GAP数组（默认所有点的距离标号均为0，则距离标号为0的点数量为n）
    for (int i = 0; i < n; i++) di[i] = pos[i]; //初始化当前弧
    int i = s, aug = INF, flow = 0; //初始化一些变量，flow为全局流量，aug为当前增广路的流量
    bool flag = false; //标记变量，记录是否找到了一条增广路（若没有找到则修正距离标号）
    while (dis[s] < n) {
        his[i] = aug; //保存当前流量
        flag = false;
        for (int p=di[i]; p; p=edge[p].next)
            if ((edge[p].f > 0) && (dis[edge[p].v] + 1 == dis[i])) {//利用距离标号判定可行弧
                flag = true; //发现可行弧
                di[i] = p; //更新当前弧
                aug = min(aug, edge[p].f); //更新当前流量
                pre[edge[p].v] = p; //记录前驱结点
                i = edge[p].v; //在弧上向前滑动
                if (i == t) {//遇到汇点，发现可增广路
                    flow += aug; //更新全局流量
                    while (i != s) {//减少增广路上相应弧的容量，并增加其反向边容量
                        edge[pre[i]].f -= aug;
                        edge[pre[i]^1].f += aug;
                        i = edge[pre[i]^1].v;
                    }
                    aug = INF;
                }
                break;
            }
        if (flag) continue; //若发现可行弧则继续，否则更新标号
        int min = n - 1;
        for (int p=pos[i]; p; p=edge[p].next)
            if ((edge[p].f > 0) && (dis[edge[p].v] < min)) {
                di[i] = p; //不要忘了重置当前弧
                min = dis[edge[p].v];
            }

        --vh[dis[i]];
        if (vh[dis[i]] == 0) break; //更新vh数组，若发现距离断层，则算法结束（GAP优化）
        dis[i] = min + 1;
        ++vh[dis[i]];
        if (i != s) {//退栈过程
            i = edge[pre[i]^1].v;
            aug = his[i];
        }
    }
    return flow;
}
int main() {
    int n, m, s, t, cas;
    scanf("%d", &cas);
    while (cas--) {
        scanf("%d%d%d%d", &n, &m, &s, &t);
        //初始化
        cl = 1;
        memset(dis, 0, sizeof(dis));
        memset(vh, 0, sizeof(vh));
        memset(pos, 0, sizeof(pos));
        //建图
        int u, v, f;
        for (int i = 0; i < m; i++) {
            scanf("%d%d%d", &u, &v, &f);
            connect(u, v, f);
        }
        int ans = maxflow(s, t, n);
        printf("%d\n", ans);
    }
    return 0;
}
