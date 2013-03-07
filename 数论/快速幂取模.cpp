#include <cstdio>
long long PowerMod(long long a, long long b, long long k) {
    long long tmp = a, ret = 1;
    while (b) {
        if (b & 1) ret = ret * tmp % k;
        tmp = tmp * tmp % k;
        b >>= 1;
    }
    return ret;
}
int main() {
    long long a, b, k;

    while (scanf("%lld%lld%lld", &a, &b, &k) == 3)
        printf("%lld\n", PowerMod(a, b, k));
    return 0;
}
