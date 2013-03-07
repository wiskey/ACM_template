#include <cstdio>
#include <cstring>
#include <vector>
using namespace std;
#define MAXN 102

int n;
vector<int> adj[MAXN];
int pre[MAXN];
bool used[MAXN];

bool CrossPath(int k) {
    for (unsigned int i=0; i<adj[k].size(); i++) {
        int j=adj[k][i];
        if (!used[j]) {
            used[j]=true;
            if (pre[j]==0 || CrossPath(pre[j])) {
                pre[j]=k;
                return true;
            }
        }
    }
    return false;
}
int hungary() {
    int match = 0;
    memset(pre, 0, sizeof(pre));
    for (int i=1; i<=n; i++) {
        memset(used, false, sizeof(used));
        if (CrossPath(i))
            match++;
    }
    return match;
}
int main() {
    int a, b, m;
    scanf("%d %d", &n, &m);
    for (int i=0; i<m; i++) {
        scanf("%d%d", &a, &b);
        adj[a].push_back(b);
        //adj[b].push_back(a);
    }
    int match = hungary();
    printf("%d\n", match);
    return 0;
}
