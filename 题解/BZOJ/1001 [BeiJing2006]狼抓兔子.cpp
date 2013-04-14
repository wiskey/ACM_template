/*
平面图的最大流（最小割）可以转化为对应偶图的最短路算法
集训队论文：2008day2  周冬《两极相通——浅析最大—最小定理在信息学竞赛中的应用》
*/

/**************************************************************
    Problem: 1001
    User: icpc046
    Language: C++
    Result: Accepted
    Time:4648 ms
    Memory:104260 kb
****************************************************************/
 
#include <cstdio>
#include <cstring>
#include <vector>
#include <queue>
using namespace std;
#define N 2000004
#define inf 0x3f3f3f3f
struct node {
    int v, w;
};
int s, t;
vector<node> g[N];
bool inq[N];
int d[N];
void Init(int n, int m) {
    int nn = n-1, mm = m-1;
    s = 0, t = 2*nn*mm + 1;
    int add = nn * mm;
    for (int i=s; i<=t; i++) g[i].clear();
    int w, a, b;
    for (int i=1; i<=n; i++)
        for (int j=1; j<m; j++) {
            scanf("%d", &w);
            if (i == 1) { a = j, b = s; }
            else if (i == n) { a = t, b = add + (i-2)*mm+j; }
            else { a = (i-1)*mm + j, b = add + (i-2)*mm+j; }
            g[a].push_back((node){b, w});
            g[b].push_back((node){a, w});
        }
    for (int i=1; i<n; i++)
        for (int j=1; j<=m; j++) {
            scanf("%d", &w);
            if (j == 1) { a = add + (i-1)*mm + j, b = t; }
            else if (j == m) { a = 0, b = i*mm; }
            else { b = (i-1)*mm + j -1, a = add + b + 1; }
            g[a].push_back((node){b, w});
            g[b].push_back((node){a, w});
        }
    for (int i=1; i<n; i++)
        for (int j=1; j<m; j++) {
            scanf("%d", &w);
            a = (i-1)*mm + j;
            b = add + a;
            g[a].push_back((node){b, w});
            g[b].push_back((node){a, w});
        }
}
int spfa() {
    queue<int> Q;
    memset(inq, false, sizeof(inq));
    memset(d, 0x3f, sizeof(d));
    Q.push(s);
    d[s] = 0;
    int now, u, v, w;
    while (!Q.empty()) {
        u = Q.front();
        Q.pop();
        inq[u] = false;
        for (size_t i=0; i<g[u].size(); i++) {
            v = g[u][i].v, w = g[u][i].w;
            if (d[u] + w < d[v]) {
                d[v] = d[u] + w;
                if (!inq[v]) { Q.push(v), inq[v] = true; }
            }
        }
    }
    return d[t];
}
int main() {
    int n, m;
    scanf("%d%d", &n, &m);
    if (n == 1) {
        int c, ret = inf;
        for (int i=0; i<m; i++) {
            scanf("%d", &c);
            if (ret > c) ret = c;
        }
        printf("%d\n", ret);
        return 0;
    }
    if (m == 1) {
        int c, ret = inf;
        for (int i=1; i<n; i++) {
            scanf("%d", &c);
            if (c < ret) ret = c;
        }
        printf("%d\n", ret);
        return 0;
    }
    Init(n, m);
    printf("%d\n", spfa());
 
    return 0;
}