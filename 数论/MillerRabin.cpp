#include <cstdio>
#include <cstdlib>
const int MAX = 4;

long long PowerMod(long long a, long long b, long long k) {
    long long ret = 1, f = a;
    while (b) {
        if (b & 1) ret = ret * f % k;
        f = f * f % k;
        b >>= 1;
    }
    return ret;
}
bool MillerRabin(long long n) {
    int i;
    long long tmp;
    srand(100);
    for (i = 0; i < MAX; i++) {
        tmp = rand() % (n - 1) + 1;
        if (PowerMod(tmp, n - 1, n) != 1)
            break;
    }
    return (i == MAX);
}
int main() {
    int l, r;
    while (scanf("%d%d", &l, &r) == 2) {
        int sum = 0;
        if (l <= 2) l = 3, sum = 1;
        if (!(l&1)) l++;
        for (int i=l; i<=r; i+=2)
            if(MillerRabin(i)) sum++;
        printf("%d\n", sum);
    }
    return 0;
}
