#include <cstdio>
#include <cstring>
using namespace std;
const int MAXN = 500;
bool g[MAXN][MAXN], flag;
int dfn[MAXN], time;
int n, m;

void init() {
    scanf("%d%d", &n, &m);
    memset(g, false, sizeof(g));
    memset(dfn, 0, sizeof(dfn));
    flag = true;
    int u, v;
    for (int i=0; i<m; i++) {
        scanf("%d%d", &u, &v);
        g[u][v] = g[v][u] = true;;
    }
}
void dfs(int now) {
    dfn[now] = ++time;
    for (int i=1; i<=n; i++) {
        if (g[now][i]) {
            if (!dfn[i]) dfs(i);
            else if ((dfn[now] - dfn[i])%2 == 0)
                flag = false;
        }
    }
}
int main() {

    init();
    time = 0;
    for (int i=1; i<=n; i++)
        if (!dfn[i]) dfs(i);

    if (flag) printf("YES\n");
    else printf("NO\n");

    return 0;
}












