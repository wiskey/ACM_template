#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;

#define MAXN 305
#define INF 99999999

int w[MAXN][MAXN];
int lx[MAXN], ly[MAXN];
bool sx[MAXN], sy[MAXN];
int pre[MAXN], slack;
int n, m;

bool CrossPath(int u) {
    sx[u] = true;
    for(int i = 0; i < m; i++) {
        if(sy[i]) continue;
        int tmp = lx[u] + ly[i] - w[u][i];
        if(!tmp) {
            sy[i] = true;
            if(pre[i] == -1 || CrossPath(pre[i])) {
                pre[i] = u;
                return true;
            }
        } else if(slack > tmp) slack = tmp;
    }
    return false;
}
int KM() {
    int i, j, k, res;
    memset(pre, -1, sizeof(pre));
    memset(ly, 0, sizeof(ly));
    for(i = 0; i < n; i++)
        for(lx[i]=-INF, j = 0; j < m; j++)
            lx[i] = max(lx[i], w[i][j]);
    for(k = 0; k < n; k++) {
        while(1) {
            slack = INF;
            memset(sx, false, sizeof(sx));
            memset(sy, false, sizeof(sy));
            if(CrossPath(k)) break;
            for(i = 0; i < n; i++)
                if(sx[i]) lx[i] -= slack;
            for(i = 0; i < m; i++)
                if(sy[i]) ly[i] += slack;
        }
    }
    res = 0;
    for(i = 0; i < m; i++)
        res += w[pre[i]][i];
    return res;
}

int main() {
    int i, j, ans;
    while(scanf("%d", &n) != EOF) {
        m = n;
        for(i = 0; i < n; i++)
            for(j = 0; j < m; j++)
                scanf("%d", &w[i][j]);
        ans = KM();
        printf("%d\n", ans);
    }
    return 0;
}


