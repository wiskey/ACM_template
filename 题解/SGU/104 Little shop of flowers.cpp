/*
	DP: 带“路径”输出的
	权值存在负数，因此预先加上100之后变成正数，可以得到简化
*/
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>
using namespace std;
#define N 102
int f[N][N], a[N][N], p[N], F, V;

int main() {
    scanf("%d%d", &F, &V);
    for (int i=1; i<=F; i++)
        for (int j=1; j<=V; j++) {
            scanf("%d", &a[i][j]);
            a[i][j] += 100;
        }
    memset(f, 0, sizeof(f));
    memset(p, 0, sizeof(p));

    for (int i=1; i<=F; i++)
        for (int j=1; j<=V; j++)
            f[i][j] = max(f[i][j-1], f[i-1][j-1] + a[i][j]);
    printf("%d\n", f[F][V]-100*F);
    int j = V;
    for (int i=F; i>0; i--) {
        while (f[i][j] == f[i][j-1]) j--;
        p[i] = j--;
    }
    for (int i=1; i<=F; i++) printf("%d%c", p[i], (i==F)?'\n':' ');
    return 0;
}
